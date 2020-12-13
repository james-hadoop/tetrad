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

import edu.cmu.tetrad.data.CorrelationMatrix;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.MatrixUtils;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.Vector;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Double.NaN;
import static java.lang.Math.*;

/**
 * Implements the continuous BIC score for FGES.
 *
 * @author Joseph Ramsey
 */
public class SemBicScore implements Score {

    private double ess;
    // The dataset.
    private DataSet dataSet;

    // The correlation matrix.
    private ICovarianceMatrix covariances;

    // The variables of the covariance matrix.
    private List<Node> variables;

    // The sample size of the covariance matrix.
    private final int sampleSize;

    // True if verbose output should be sent to out.
    private boolean verbose = false;

    // A  map from variable names to their indices.
    private final Map<Node, Integer> indexMap;

    // The penalty penaltyDiscount, 1 for standard BIC.
    private double penaltyDiscount = 1.0;

    // The structure prior, 0 for standard BIC.
    private double structurePrior = 0.0;

    // Equivalent sample size
    private Matrix matrix;

    // The rule type to use.
    private RuleType ruleType = RuleType.HIGH_DIMENSIONAL;
    private double[] mu;

    /**
     * Constructs the score using a covariance matrix.
     */
    public SemBicScore(ICovarianceMatrix covariances) {
        if (covariances == null) {
            throw new NullPointerException();
        }

        setCovariances(covariances);
        this.variables = covariances.getVariables();
        this.sampleSize = covariances.getSampleSize();
        this.indexMap = indexMap(this.variables);
    }

    /**
     * Constructs the score using a covariance matrix.
     */
    public SemBicScore(DataSet dataSet) {
        if (dataSet == null) {
            throw new NullPointerException();
        }

        double[][] columns = dataSet.getDoubleData().transpose().toArray();

        mu = new double[columns.length];

        for (int t = 0; t < columns.length; t++) {
            mu[t] = mean(columns[t]);
        }

        if (!dataSet.existsMissingValue()) {
//            setCovariances(new CorrelationMatrix(new CovarianceMatrix(dataSet, false)));
            setCovariances(new CovarianceMatrix(dataSet, false));
            this.variables = covariances.getVariables();
            this.sampleSize = covariances.getSampleSize();
            this.indexMap = indexMap(this.variables);


            return;
        }

        this.dataSet = dataSet;

        this.variables = dataSet.getVariables();
        this.sampleSize = dataSet.getNumRows();
        this.indexMap = indexMap(this.variables);
    }

    @Override
    public double localScoreDiff(int x, int y, int[] z) {

        if (ruleType == RuleType.NANDY) {
            return nandyBic(x, y, z);
        } else {
            return localScore(y, append(z, x)) - localScore(y, z);
        }
    }

    public double nandyBic(int x, int y, int[] z) {
        double sp1 = getStructurePrior(z.length + 1);
        double sp2 = getStructurePrior(z.length);

        Node _x = variables.get(x);
        Node _y = variables.get(y);
        List<Node> _z = getVariableList(z);

        List<Integer> rows = getRows(x, z);

        if (rows != null) {
            rows.retainAll(Objects.requireNonNull(getRows(y, z)));
        }

        double r = partialCorrelation(_x, _y, _z, rows);

        double c = getPenaltyDiscount();

        return -sampleSize * log(1.0 - r * r) - c * log(sampleSize)
                - 2.0 * (sp1 - sp2);
    }

    private int[] append(int[] z, int x) {
        int[] _z = Arrays.copyOf(z, z.length + 1);
        _z[z.length] = x;
        return _z;
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScoreDiff(x, y, new int[0]);
    }

    public double localScore(int i, int... parents) {
        List<Integer> rows = getRows(i, parents);

        final int p = parents.length;
        int k = 2 * p + 1;

        int[] all = concat(i, parents);

        // Only do this once.
        Matrix cov = getCov(rows, all, all);

        int[] pp = orderedParents(parents);

        Matrix covxx = cov.getSelection(pp, pp);
        Matrix covxy = cov.getSelection(pp, new int[]{0});

        // My calculation.
        Matrix b = covxx.inverse().times(covxy);
//        Matrix b2 = adjustedParents(p, b);
//        Matrix times = b2.transpose().times(cov).times(b2);
//        double s2 = times.get(0, 0);

        // Ricardo's calculation.
        double s2 = cov.get(0, 0);
        Vector _cxy = covxy.getColumn(0);
        Vector _b = b.getColumn(0);
        s2 -= _cxy.dotProduct(_b);

        double V = variables.size();
        double N = sampleSize / 2.;

        if (ruleType == RuleType.CHICKERING || ruleType == RuleType.NANDY) {

            // Standard BIC, with penalty discount and structure prior.
            double c = getPenaltyDiscount();

            return -N * log(s2) - k * c * log(N);// - 2 * getStructurePrior(p);

        } else if (ruleType == RuleType.HIGH_DIMENSIONAL) {

            // Pseudo-BIC formula, set Wikipedia page for BIC. With penalty discount and structure prior.

            // We will just take 6 * omega * (1 + gamma) to be a number >= 6. To be compatible with other scores,
            // we will use c + 5 for this value, where c is the penalty discount. So a penalty discount of 1 (the usual)
            // will correspond to 6 * omega * (1 + gamma) of 6, the minimum.

            double c = getPenaltyDiscount();

            return -N * log(s2) - k * c * 6 * log(V);// - 2 * getStructurePrior(p);

        } else {
            throw new IllegalStateException("That rule type is not implemented: " + ruleType);
        }
    }

    public int[] orderedParents(int[] parents) {
        int[] pp = new int[parents.length];
        for (int j = 0; j < pp.length; j++) pp[j] = j + 1;
        return pp;
    }

    public int[] concat(int i, int[] parents) {
        int[] all = new int[parents.length + 1];
        all[0] = i;
        System.arraycopy(parents, 0, all, 1, parents.length);
        return all;
    }

    @NotNull
    public Matrix adjustedParents(int p, Matrix b) {
        Matrix byx = new Matrix(p + 1, 1);
        byx.set(0, 0, 1);
        for (int j = 0; j < p; j++) byx.set(j + 1, 0, -b.get(j, 0));
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

    public double getPenaltyDiscount() {
        return penaltyDiscount;
    }

    public double getStructurePrior() {
        return structurePrior;
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

    public DataSet getDataSet() {
        return dataSet;
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }

    public void setStructurePrior(double structurePrior) {
        this.structurePrior = structurePrior;
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

    public void setVariables(List<Node> variables) {
        if (covariances != null) {
            covariances.setVariables(variables);
        }

        this.variables = variables;
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
        return (int) Math.ceil(log(sampleSize));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        int i = variables.indexOf(y);

        int[] k = new int[z.size()];

        for (int t = 0; t < z.size(); t++) {
            k[t] = variables.indexOf(z.get(t));
        }

        double v = localScore(i, k);

        return Double.isNaN(v);
    }

    private void setCovariances(ICovarianceMatrix covariances) {
        this.covariances = covariances;
        this.matrix = this.covariances.getMatrix();

        Matrix cor = new CorrelationMatrix(covariances).getMatrix();

//        for (int i = 0; i < cor.rows(); i++) {
//            for (int j = 0; j < cor.columns(); j++) {
//                cor.set(i, j, abs(cor.get(i, j)));
//            }
//        }
//
//        System.out.println("cor = " + cor);
//
//        System.out.println("cor.trace() = "+ cor.trace());


        cor = cor.minus(Matrix.identity(cor.rows()));

        double sum = 0.0;

        int count = 0;

        for (int i = 0; i < cor.rows(); i++) {
            for (int j = 0; j < cor.columns(); j++) {
                if (i == j) continue;
                double r = cor.get(i, j);
                r = abs(r);

                if (r >= 0) {
                    sum += r;
                }

                count++;
            }
        }

        double rho = sum / count;

//        double V = covariances.getSize();
        double N = covariances.getSampleSize();
//
//        double rho = sqrt((((cor.times(cor))).trace() / 2.) * (15. / (N * N)));

//        double rho = sqrt((cor.times(cor).trace())) * (1. / (V * V));

        ess = N / (1 + (N - 1) * rho);


//        double sum = 0.0;
//        int count = 0;
//
//        for (int i = 0; i < V; i++) {
//            for (int j = 0; j < V; j++) {
//                double rho = cor.get(i, j);
//
//                double f = N / (1 + (N - 1) * rho);
//
//                if (f > 0) {
//                    sum += f;
//                    count++;
//                }
//            }
//        }
//
//        ess = sum / (count);

        ess = N;

        System.out.println("N = " + N + " ess = " + ess + " rho = " + rho);


    }

    private double getStructurePrior(int parents) {
        if (abs(getStructurePrior()) <= 0) {
            return 0;
        } else {
            double p = (getStructurePrior()) / (variables.size());
            return -((parents) * Math.log(p) + (variables.size() - (parents)) * Math.log(1.0 - p));
        }
    }

    private List<Node> getVariableList(int[] indices) {
        List<Node> variables = new ArrayList<>();
        for (int i : indices) {
            variables.add(this.variables.get(i));
        }
        return variables;
    }

    private Map<Node, Integer> indexMap(List<Node> variables) {
        Map<Node, Integer> indexMap = new HashMap<>();

        for (int i = 0; variables.size() > i; i++) {
            indexMap.put(variables.get(i), i);
        }

        return indexMap;
    }

    /**
     * @return a string representation of this score.
     */
    public String toString() {
        NumberFormat nf = new DecimalFormat("0.00");
        return "SEM BIC Score penalty " + nf.format(penaltyDiscount);
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
                    mui += dataSet.getDouble(k, _rows[i]);
                    muj += dataSet.getDouble(k, cols[j]);
                }

                mui /= rows.size() - 1;
                muj /= rows.size() - 1;

                double _cov = 0.0;

                for (int k : rows) {
                    _cov += (dataSet.getDouble(k, _rows[i]) - mui) * (dataSet.getDouble(k, cols[j]) - muj);
                }

                double mean = _cov / (rows.size());
                cov.set(i, j, mean);
            }
        }

        return cov;
    }

    private List<Integer> getRows(int i, int[] parents) {
        if (dataSet == null) {
            return null;
        }

        List<Integer> rows = new ArrayList<>();

        K:
        for (int k = 0; k < dataSet.getNumRows(); k++) {
            if (Double.isNaN(dataSet.getDouble(k, i))) continue;

            for (int p : parents) {
                if (Double.isNaN(dataSet.getDouble(k, p))) continue K;
            }

            rows.add(k);
        }

        return rows;
    }

    private double partialCorrelation(Node x, Node y, List<Node> z, List<Integer> rows) {
        try {
            return StatUtils.partialCorrelation(MatrixUtils.convertCovToCorr(getCov(rows, indices(x, y, z))));
        } catch (Exception e) {
            return NaN;
        }
    }

    private int[] indices(Node x, Node y, List<Node> z) {
        int[] indices = new int[z.size() + 2];
        indices[0] = indexMap.get(x);
        indices[1] = indexMap.get(y);
        for (int i = 0; i < z.size(); i++) indices[i + 2] = indexMap.get(z.get(i));
        return indices;
    }

    private Matrix getCov(List<Integer> rows, int[] cols) {
        if (dataSet == null) {
            return matrix.getSelection(cols, cols);
        }

        Matrix cov = new Matrix(cols.length, cols.length);

        for (int i = 0; i < cols.length; i++) {
            for (int j = i + 1; j < cols.length; j++) {
                double mui = 0.0;
                double muj = 0.0;

                for (int k : rows) {
                    mui += dataSet.getDouble(k, cols[i]);
                    muj += dataSet.getDouble(k, cols[j]);
                }

                mui /= rows.size() - 1;
                muj /= rows.size() - 1;

                double _cov = 0.0;

                for (int k : rows) {
                    _cov += (dataSet.getDouble(k, cols[i]) - mui) * (dataSet.getDouble(k, cols[j]) - muj);
                }

                double mean = _cov / (rows.size());
                cov.set(i, j, mean);
                cov.set(j, i, mean);
            }
        }

        for (int i = 0; i < cols.length; i++) {
            double mui = 0.0;

            for (int k : rows) {
                mui += dataSet.getDouble(k, cols[i]);
            }

            mui /= rows.size();

            double _cov = 0.0;

            for (int k : rows) {
                _cov += (dataSet.getDouble(k, cols[i]) - mui) * (dataSet.getDouble(k, cols[i]) - mui);
            }

            double mean = _cov / (rows.size());
            cov.set(i, i, mean);
        }

        return cov;
    }

    public void setRuleType(RuleType ruleType) {
        this.ruleType = ruleType;
    }

    public enum RuleType {CHICKERING, NANDY, HIGH_DIMENSIONAL}
}


