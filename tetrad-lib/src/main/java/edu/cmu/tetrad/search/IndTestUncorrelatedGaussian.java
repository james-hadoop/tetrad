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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.NumberFormatUtil;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.text.NumberFormat;
import java.util.*;

import static java.lang.Math.*;

/**
 * Checks conditional independence of variable in a continuous data set using a conditional correlation test
 * for the nonlinear nonGaussian case.
 *
 * @author Joseph Ramsey
 */
public final class IndTestUncorrelatedGaussian implements IndependenceTest, ScoreForFact {

    /**
     * The variables of the covariance data, in order. (Unmodifiable list.)
     */
    private final List<Node> variables;

    /**
     * The significance level of the independence tests.
     */
    private double alpha;

    /**
     * Formats as 0.0000.
     */
    private static NumberFormat nf = NumberFormatUtil.getInstance().getNumberFormat();

    /**
     * Stores a reference to the data set passed in through the constructor.
     */
    private final DataSet dataSet;

    /**
     * True iff the fast FDR adjustment should be made.
     */
    private boolean fastFDR = false;

    /**
     * True if verbose output should be printed.
     */
    private boolean verbose = false;
    private double score = Double.NaN;

    IndTestFisherZ fisher;

    /**
     * Looks up the index of a record in the the sorted order for each variable z.
     */
    private final List<Map<Integer, Integer>> reverseLookup;

    /**
     * The matrix of data, N x M, where N is the number of samples, M the number
     * of variables, gotten from dataSet.
     */
    private final double[][] data;

    private int kernelRegressionSampleSize = 100;

    /**
     * The ith array gives indices into the ith variables in sorted order.
     */
    private final ArrayList<List<Integer>> sortedIndices;

    /**
     * Map from nodes to the indices.
     */
    private final Map<Node, Integer> indices;


    //==========================CONSTRUCTORS=============================//

    /**
     * Constructs a new Independence test which checks independence facts based on the correlation data implied by the
     * given data set (must be continuous). The given significance level is used.
     *
     * @param dataSet A data set containing only continuous columns.
     * @param alpha   The q level of the test.
     */
    public IndTestUncorrelatedGaussian(DataSet dataSet, double alpha) {
        if (!(dataSet.isContinuous())) {
            throw new IllegalArgumentException("Data set must be continuous.");
        }

        if (!(alpha >= 0 && alpha <= 1)) {
            throw new IllegalArgumentException("Q mut be in [0, 1]");
        }

        this.dataSet = DataUtils.getNonparanormalTransformed(dataSet);

        this.data = dataSet.getDoubleData().transpose().toArray();

        List<Node> nodes = dataSet.getVariables();

        this.variables = Collections.unmodifiableList(nodes);

        indices = new HashMap<>();

        for (int i = 0; i < variables.size(); i++) {
            indices.put(variables.get(i), i);
        }

        sortedIndices = new ArrayList<>();

        for (double[] x : data) {
            List<Integer> sorted = new ArrayList<>();
            for (int t = 0; t < x.length; t++) sorted.add(t);

            sorted.sort(Comparator.comparingDouble(o -> x[o]));
            sortedIndices.add(sorted);
        }

        reverseLookup = new ArrayList<>();

        for (int z2 = 0; z2 < data.length; z2++) {
            Map<Integer, Integer> m = new HashMap<>();

            for (int j = 0; j < data[z2].length; j++) {
                m.put(sortedIndices.get(z2).get(j), j);
            }

            reverseLookup.add(m);
        }

        fisher = new IndTestFisherZ(this.dataSet, alpha);
    }

    //==========================PUBLIC METHODS=============================//

    /**
     * Creates a new IndTestCramerT instance for a subset of the variables.
     */
    public IndependenceTest indTestSubset(List<Node> vars) {
        throw new UnsupportedOperationException();
    }

    public boolean isIndependent(Node x, Node y, List<Node> z) {
        int n = dataSet.getNumRows();

        if (fisher.isDependent(x, y, z)) {
            return false;
        } else {

            if (z.isEmpty()) {
                return bivariateGaussian(x, y, allRows(n));
            } else {

                int[] _z = new int[z.size()];

                for (int m = 0; m < z.size(); m++) {
                    _z[m] = indices.get(z.get(m));
                }

                for (int i = 0; i < n; i+=10) {
                    Set<Integer> js = getCloseZs(data, _z, i, getKernelRegressionSampleSize());

                    if (!bivariateGaussian(x, y, js)) {
                        return false;
                    }
                }

                return true;
            }
        }
    }

    private boolean zeroCorr(Node x, Node y, Set<Integer> rows) {
        List<Integer> _rows = new ArrayList<>(rows);

        double[] x1 = subset(x, _rows);
        double[] x2 = subset(y, _rows);

        double r = StatUtils.correlation(x1, x2);

        double z1 = 0.5 * sqrt(_rows.size() - 3) * (log(1 + r) - log(1 - r));
        double p2 = 2.0 * (1 - new NormalDistribution(0, 1).cumulativeProbability(abs(z1)));

        return p2 > alpha;
    }

    private boolean bivariateGaussian(Node x, Node y, Set<Integer> rows) {
        List<Integer> _rows = new ArrayList<>(rows);

        double[] x1 = subset(x, _rows);
        double[] x2 = subset(y, _rows);

        x1 = DataUtils.standardizeData(x1);
        x2 = DataUtils.standardizeData(x2);

        final int n = _rows.size();


        double s1 = StatUtils.sd(x1);
        double s2 = StatUtils.sd(x2);

        double r = StatUtils.correlation(x1, x2);

        double[] x1c = DataUtils.center(x1);
        double[] x2c = DataUtils.center(x2);

        TetradMatrix X = new TetradMatrix(2, n);

        for (int i = 0; i < n; i++) {
            X.set(0, i, x1c[i]);
            X.set(1, i, x2c[i]);
        }

        TetradMatrix A = new TetradMatrix(2, 2);
        A.set(0, 0, s1 * s1);
        A.set(0, 1, r * s1 * s2);
        A.set(1, 0, r * s1 * s2);
        A.set(1, 1, s2 * s2);

        TetradMatrix Y = A.inverse().times(X);

        double[] y1 = Y.getRow(0).toArray();
        double[] y2 = Y.getRow(1).toArray();

        final double m21 = moment(y1, y2, 2, 1);
        final double m12 = moment(y1, y2, 1, 2);
        final double m30 = moment(y1, y2, 3, 0);
        final double m03 = moment(y1, y2, 0, 3);

        double u3 = n * ((pow(m21, 2) + pow(m12, 2)) / 2.0 + ((pow(m30, 2) + pow(m03, 2)) / 6.0));

        final double m22 = moment(y1, y2, 2, 2);
        final double m31 = moment(y1, y2, 3, 1);
        final double m13 = moment(y1, y2, 1, 3);
        final double m04 = moment(y1, y2, 0, 4);
        final double m40 = moment(y1, y2, 4, 0);

        double u4 = n * (pow(m22 - 1, 2) / 4.0 + (pow(m31, 2) + pow(m13, 2)) / 6.0 + (pow(m04 - 1, 2) + pow(m40 - 1, 2)) / 24.0);

        double p = 1.0 - new ChiSquaredDistribution(4 * n).cumulativeProbability(u3);

//        System.out.println("u3 = " + u3);

        return p > alpha;

    }

    private double[] subset(Node x, List<Integer> rows) {
        double[] d = data[indices.get(x)];

        double[] d2 = new double[rows.size()];

        for (int i = 0; i < rows.size(); i++) {
            d2[i] = d[rows.get(i)];
        }

        return d2;
    }

    private double moment(double[] x, double[] y, int m, int n) {
        int N = x.length;
        double sum = 0.0;

        for (int j = 0; j < x.length; j++) {
            sum += Math.pow(x[j], m) * Math.pow(y[j], n);
        }

        return sum / N;
    }

    private Set<Integer> allRows(int n) {
        Set<Integer> rows = new HashSet<>();

        for (int i = 0; i < n; i++) {
            rows.add(i);
        }

        return rows;
    }

    private Set<Integer> getCloseZs(double[][] data, int[] _z, int i, int sampleSize) {
        Set<Integer> js = new HashSet<>();

        if (sampleSize > data[0].length) sampleSize = (int) ceil(0.8 * data.length);
        if (_z.length == 0) return new HashSet<>();

        int radius = 0;

        while (true) {
            for (int z1 : _z) {
                int q = reverseLookup.get(z1).get(i);

                if (q - radius >= 0 && q - radius < data[z1].length) {
                    final int r2 = sortedIndices.get(z1).get(q - radius);
                    js.add(r2);
                }

                if (q + radius >= 0 && q + radius < data[z1].length) {
                    final int r2 = sortedIndices.get(z1).get(q + radius);
                    js.add(r2);
                }

            }

            if (js.size() >= sampleSize) return js;

            radius++;
        }
    }

    public boolean isIndependent(Node x, Node y, Node... z) {
        return isIndependent(x, y, Arrays.asList(z));
    }

    public boolean isDependent(Node x, Node y, List<Node> z) {
        return !isIndependent(x, y, z);
    }

    public boolean isDependent(Node x, Node y, Node... z) {
        List<Node> zList = Arrays.asList(z);
        return isDependent(x, y, zList);
    }

    public double getPValue() {
        return 0.0;
    }

    /**
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     */
    public void setAlpha(double alpha) {
        if (alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("Significance out of range.");
        }

        this.alpha = alpha;
    }

    /**
     * Gets the getModel significance level.
     */
    public double getAlpha() {
        return this.alpha;
    }

    /**
     * @return the list of variables over which this independence checker is capable of determinine independence
     * relations-- that is, all the variables in the given graph or the given data set.
     */
    public List<Node> getVariables() {
        return this.variables;
    }

    /**
     * @return the variable with the given name.
     */
    public Node getVariable(String name) {
        for (Node node : variables) {
            if (node.getName().equals(name)) return node;
        }

        throw new IllegalArgumentException();
    }

    /**
     * @return the list of variable varNames.
     */
    public List<String> getVariableNames() {
        List<Node> variables = getVariables();
        List<String> variableNames = new ArrayList<>();
        for (Node variable1 : variables) {
            variableNames.add(variable1.getName());
        }
        return variableNames;
    }

    /**
     * If <code>isDeterminismAllowed()</code>, deters to IndTestFisherZD; otherwise throws
     * UnsupportedOperationException.
     */
    public boolean determines(List<Node> z, Node x) throws UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * @return the data set being analyzed.
     */
    public DataSet getData() {
        return dataSet;
    }

    @Override
    public ICovarianceMatrix getCov() {
        return null;
    }

    @Override
    public List<DataSet> getDataSets() {
        return null;
    }

    @Override
    public int getSampleSize() {
        return 0;
    }

    @Override
    public List<TetradMatrix> getCovMatrices() {
        return null;
    }

    @Override
    public double getScore() {
        return score;
    }

    @Override
    public double getScoreForFact(IndependenceFact fact) {
        throw new UnsupportedOperationException();
    }

    /**
     * @return a string representation of this test.
     */
    public String toString() {
        return "Conditional Correlation, q = " + nf.format(getAlpha());
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * The minimum sample size to use for the kernel regression.
     */
    public int getKernelRegressionSampleSize() {
        return kernelRegressionSampleSize;
    }

    public void setKernelRegressionSampleSize(int kernelRegressionSampleSize) {
        this.kernelRegressionSampleSize = kernelRegressionSampleSize;
    }
}



