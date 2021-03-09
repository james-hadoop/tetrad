///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.Vector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.*;

/**
 * Implements the continuous BIC score for FGES.
 *
 * @author Joseph Ramsey
 */
public class ZhangShenBoundScore implements Score {
    private double trueErrorVariance = 1;

    private int depth = -1;
    // The covariance matrix.
    private ICovarianceMatrix covariances;

    // The variables of the covariance matrix.
    private final List<Node> variables;

    // The sample size of the covariance matrix.
    private final int sampleSize;

    // True if verbose output should be sent to out.
    private boolean verbose = false;

    // Sample size or equivalent sample size.
    private double N;

    // A recpord of lambdas for each m0.
    private List<Double> lambdas;

    // The minimim probability lower bound for risk.
    private double riskBound = 0.05;

    // The data, if it is set.
    private Matrix data;

    // True if sume of squares should be calculated, false if estimated.
    private boolean calculateSquaredEuclideanNorms = false;

    // True if row subsets should be calculated.
    private boolean calculateRowSubsets = false;

    // A Map from nodes to the BICs for their minimal models.
    Map<Node, Double> bics = new HashMap<>();

    // A map from nodes to the varey's of their minimal models.
    Map<Node, Double> vareys = new HashMap<>();

    // A ,ap from nodes to the lambdas for their minimal models.
    Map<Node, Double> lambdass = new HashMap<>();

    // Depth for finding initial minimal models.
    private double epsilon = 1e-6;

    /**
     * Constructs the score using a covariance matrix.
     */
    public ZhangShenBoundScore(ICovarianceMatrix covariances, double epsilon, int depth) {
        this.epsilon = epsilon;
        this.depth = depth;

        if (covariances == null) {
            throw new NullPointerException();
        }

        setCovariances(covariances);
        this.variables = covariances.getVariables();
        this.sampleSize = covariances.getSampleSize();
    }

    /**
     * Constructs the score using a covariance matrix.
     */
    public ZhangShenBoundScore(DataSet dataSet, boolean calculateSquaredEuclideanNorms,
                               double epsilon, int depth, double trueErrorVariance) {
        this.epsilon = epsilon;
        this.depth = depth;

        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.calculateSquaredEuclideanNorms = calculateSquaredEuclideanNorms;
        this.variables = dataSet.getVariables();
        this.sampleSize = dataSet.getNumRows();

        if (!dataSet.existsMissingValue()) {
            DataSet _dataSet = DataUtils.center(dataSet);
            this.data = _dataSet.getDoubleData();

            setCovariances(new CovarianceMatrix(dataSet));
            calculateRowSubsets = false;
            setTrueStats();

            return;
        }

        calculateRowSubsets = true;
        DataSet _dataSet = DataUtils.center(dataSet);
        this.data = _dataSet.getDoubleData();
        setTrueStats();
    }

    private void setTrueStats() {
        for (Node y : variables) {

            List<Node> _adj = new ArrayList<>(variables);
            _adj.remove(y);

            List<Node> adj = new ArrayList<>();

            for (Node node : _adj) {
                Vector column = data.getColumn(variables.indexOf(y));
                Vector column1 = data.getColumn(variables.indexOf(node));
                if ((abs(StatUtils.correlation(column.toArray(), column1.toArray()))) >= epsilon) {
                    adj.add(node);
                }
            }

            DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), depth == -1 ? adj.size() : depth);
            int[] choice;

            double _score = POSITIVE_INFINITY;
            double _varey = NEGATIVE_INFINITY;
            double _lambda = NEGATIVE_INFINITY;
            List<Node> _pi = new ArrayList<>();

            while ((choice = gen.next()) != null) {
                List<Node> pi = GraphUtils.asList(choice, adj);

                int[] indices = indices(pi);
                double varey = getVarey(variables.indexOf(y), indices);
                double lambda = getLambda(pi.size());

                int pn = variables.size() - 1;
                double sigma2 = 1;
                double score = sampleSize * varey + pi.size() * 2 * (log(pn) + log(log(pn))) * trueErrorVariance;

                if (score < _score) {
                    _score = score;
                    _varey = varey;
                    _lambda = lambda;
                    _pi = pi;
                }
            }

            System.out.println(" var = " + y + " _varey = " + _varey + " _lambda = " + _lambda + " _score = " + _score + " _pi = " + _pi);

            bics.put(y, _score);
            vareys.put(y, 3.);//_varey);//covariances.getValue(variables.indexOf(y), variables.indexOf(y)));
            lambdass.put(y, getLambda(0));
        }
    }

    private int[] indices(List<Node> __adj) {
        int[] indices = new int[__adj.size()];
        for (int t = 0; t < __adj.size(); t++) indices[t] = variables.indexOf(__adj.get(t));
        return indices;
    }

    @Override
    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScoreDiff(x, y, new int[0]);
    }

    public double localScore(int i, int... parents) {
        final int p = parents.length;
        double sum;

        double varey = getVarey(i, parents);

        if (calculateSquaredEuclideanNorms) {
            sum = getSquaredEucleanNorm(i, parents);
        } else {
            sum = N * varey;
        }

        double varey2 = vareys.get(variables.get(i));
        double lambda = lambdass.get(variables.get(i));

        lambda = zhangShenLambda(variables.size() - 1, parents.length, 0);

        return -sum - lambda * p * varey2;
    }

    private double getVarey(int i, int[] parents) {
        int[] all = concat(i, parents);
        Matrix cov = getCov(getRows(i, parents), all, all);
        return getVarey(cov, parents);
    }

    private double getVarey(Matrix cov, int[] parents) {
        int[] pp = indexedParents(parents);

        Matrix covxx = cov.getSelection(pp, pp);
        Matrix covxy = cov.getSelection(pp, new int[]{0});

        Matrix b = adjustedCoefs(covxx.inverse().times(covxy));
        Matrix times = b.transpose().times(cov).times(b);
        return sqrt(times.get(0, 0));
    }

    private double getLambda(int m) {
        if (lambdas == null) {
            lambdas = new ArrayList<>();
        }

        if (lambdas.size() - 1 < m) {
            for (int t = lambdas.size(); t <= m; t++) {
                double lambda = zhangShenLambda(variables.size() - 1, t, riskBound);
                lambdas.add(lambda);
            }
        }

        return lambdas.get(m);
    }

    private double getSquaredEucleanNorm(int i, int[] parents) {
        int[] rows = new int[data.rows()];
        for (int t = 0; t < rows.length; t++) rows[t] = t;

        Matrix y = data.getSelection(rows, new int[]{i});
        Matrix x = data.getSelection(rows, parents);

        Matrix x2 = new Matrix(x.rows(), x.columns() + 1);
        for (int q = 0; q < x.rows(); q++) {
            for (int r = 0; r < x.columns(); r++) {
                x2.set(q, r, x.get(q, r));
            }
        }

        for (int q = 0; q < x.rows(); q++) {
            x2.set(q, x.columns(), 1);
        }

        x = x2;

        Matrix xT = x.transpose();
        Matrix xTx = xT.times(x);
        Matrix xTxInv = xTx.inverse();
        Matrix xTy = xT.times(y);
        Matrix b = xTxInv.times(xTy);

        Matrix yhat = x.times(b);

        double sum = 0.0;

        for (int q = 0; q < data.rows(); q++) {
            double diff = data.get(q, i) - yhat.get(q, 0);
            sum += diff * diff;
        }

        return sum;
    }


    @NotNull
    public Matrix adjustedCoefs(Matrix b) {
        Matrix byx = new Matrix(b.rows() + 1, 1);
        byx.set(0, 0, 1);
        for (int j = 0; j < b.rows(); j++) byx.set(j + 1, 0, -b.get(j, 0));
        return byx;
    }

    /**
     * Specialized scoring method for a single parent. Used to speed up the effect edges search.
     */
    public double localScore(int i, int parent) {
        return localScore(i, new int[]{parent});
    }

    /**
     * Specialized scoring method for no parents. Used to speed up the effect edges search.
     */
    public double localScore(int i) {
        return localScore(i, new int[0]);
    }

    public ICovarianceMatrix getCovariances() {
        return covariances;
    }

    public int getSampleSize() {
        return sampleSize;
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    @Override
    public List<Node> getVariables() {
        return variables;
    }

    @Override
    public Node getVariable(String targetName) {
        for (Node node : variables) {
            if (node.getName().equals(targetName)) {
                return node;
            }
        }

        return null;
    }

    @Override
    public int getMaxDegree() {
        return (int) ceil(log(sampleSize));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        int i = variables.indexOf(y);

        int[] k = indices(z);

        double v = localScore(i, k);

        return Double.isNaN(v);
    }

    @Override
    public Score defaultScore() {
        return new ZhangShenBoundScore(covariances, 0, -1);
    }

    private void setCovariances(ICovarianceMatrix covariances) {
//        this.covariances = new CorrelationMatrix(covariances);
        this.covariances = covariances;
        this.N = covariances.getSampleSize();
    }

    private static int[] append(int[] z, int x) {
        int[] _z = Arrays.copyOf(z, z.length + 1);
        _z[z.length] = x;
        return _z;
    }

    private static int[] indexedParents(int[] parents) {
        int[] pp = new int[parents.length];
        for (int j = 0; j < pp.length; j++) pp[j] = j + 1;
        return pp;
    }

    private static int[] concat(int i, int[] parents) {
        int[] all = new int[parents.length + 1];
        all[0] = i;
        System.arraycopy(parents, 0, all, 1, parents.length);
        return all;
    }

    private Matrix getCov(List<Integer> rows, int[] _rows, int[] cols) {
        if (rows == null) {
            return getCovariances().getSelection(_rows, cols);
        }

        Matrix cov = new Matrix(_rows.length, cols.length);

        for (int i = 0; i < _rows.length; i++) {
            for (int j = 0; j < cols.length; j++) {
                double mui = 0.0;
                double muj = 0.0;

                for (int k : rows) {
                    mui += data.get(k, _rows[i]);
                    muj += data.get(k, cols[j]);
                }

                mui /= rows.size() - 1;
                muj /= rows.size() - 1;

                double _cov = 0.0;

                for (int k : rows) {
                    _cov += (data.get(k, _rows[i]) - mui) * (data.get(k, cols[j]) - muj);
                }

                double mean = _cov / (rows.size());
                cov.set(i, j, mean);
            }
        }

        return cov;
    }

    private List<Integer> getRows(int i, int[] parents) {
        if (!calculateRowSubsets) {
            return null;
        }

        List<Integer> rows = new ArrayList<>();

        K:
        for (int k = 0; k < data.rows(); k++) {
            if (Double.isNaN(data.get(k, i))) continue;

            for (int p : parents) {
                if (Double.isNaN(data.get(k, p))) continue K;
            }

            rows.add(k);
        }

        return rows;
    }

    public void setRiskBound(double riskBound) {
        if (riskBound < 0 || riskBound > 1) throw new IllegalStateException(
                "Risk probability should be in [0, 1]: " + this.riskBound);
        this.riskBound = riskBound;
    }

    public double zhangShenLambda(int pn, int m0, double riskBound) {
        double high = 100000.0;
        double low = 0.0;

        while (high - low > 1e-10) {
            double lambda = (high + low) / 2.0;

            double p = getP(pn, m0, lambda);

            if (p < 1.0 - riskBound) {
                low = lambda;
            } else {
                high = lambda;
            }
        }

        return (high + low) / 2.0;
    }

    private double getP(int pn, int m0, double lambda) {
        return 2 - pow(1 + exp(-(lambda - 1) / 2.) * sqrt(lambda), pn - m0);
    }

    public void setTrueErrorVariance(double trueErrorVariance) {
        this.trueErrorVariance = trueErrorVariance;
    }
}


