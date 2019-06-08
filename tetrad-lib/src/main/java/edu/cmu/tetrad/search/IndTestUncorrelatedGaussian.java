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
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.covariance;
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

        this.alpha = alpha;

//        this.dataSet = dataSet;
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
//        if (fisher.isDependent(x, y, z)) return false;

        int n = dataSet.getNumRows();

        double[] _x = Arrays.copyOf(data[indices.get(x)], data[0].length);
        double[] _y = Arrays.copyOf(data[indices.get(y)], data[0].length);

        if (z.isEmpty()) {
//            _x = nonparanormal(_x);
//            _y = nonparanormal(_y);
            return bivariateIndependent(_x, _y, alpha);
        } else {

            int[] _z = new int[z.size()];

            for (int m = 0; m < z.size(); m++) {
                _z[m] = indices.get(z.get(m));
            }

            for (int i = 0; i < n; i += 20) {
                int j = i;//RandomUtil.getInstance().nextInt(n);
                List<Integer> js = getCloseZs(data, _z, j, getKernelRegressionSampleSize());
                List<Double> distances = new ArrayList<>();

                for (int k = 0; k < js.size(); k++) {
                    double d = distance(data, _z, j, js.get(k));
                    distances.add(d);
                }

                if (!(bivariateIndependent(subset(_x, js, distances), subset(_y, js, distances), alpha))) {
                    return false;
                }
            }

            return true;
        }

    }

    private boolean bivariateIndependent(double[] x, double[] y, double alpha) {
        double[] _x = nonparanormal(x);
        double[] _y = nonparanormal(y);

        return zeroCorr(_x, _y, alpha) && bivariateGaussian(_x, _y, alpha);
    }

    final static NormalDistribution normalDistribution = new NormalDistribution(0, 1);

    private boolean zeroCorr(double[] x, double[] y, double alpha) {
        try {
            final int n = x.length;
            double r = StatUtils.correlation(x, y);

            r += 0.1;

            if (r >= 1.0) r = 0.999999999;
            if (r <= -1.0) r = -0.999999999;

            double z = 0.5 * sqrt(n) * (log(1 + r) - log(1 - r));
            double p = 2.0 * (1 - normalDistribution.cumulativeProbability(abs(z)));

//            System.out.println("p = " + p + " z = " + z);
            if (abs(z) > 8) return true;
            return p > alpha;
        } catch (Exception e) {
            e.printStackTrace();
            throw e;
        }
    }

    private static boolean bivariateGaussian(double[] x, double[] y, double alpha) {

        final int n = x.length;

        double[] x1 = DataUtils.center(x);
        double[] x2 = DataUtils.center(y);

        double s1 = StatUtils.sd(x1);
        double s2 = StatUtils.sd(x2);

        if (s1 == 0) s1 = 1e-20;
        if (s2 == 0) s2 = 1e-20;

        double r = StatUtils.correlation(x1, x2);

        if (r >= 1.0) r = 0.999999999;
        if (r <= -1.0) r = -0.999999999;

        TetradMatrix X = new TetradMatrix(2, n);

        for (int i = 0; i < n; i++) {
            X.set(0, i, x1[i]);
            X.set(1, i, x2[i]);
        }

        TetradMatrix A = new TetradMatrix(2, 2);
        A.set(0, 0, s1 * s1);
        A.set(0, 1, r * s1 * s2);
        A.set(1, 0, r * s1 * s2);
        A.set(1, 1, s2 * s2);

        TetradMatrix Y;

        try {
            Y = A.inverse().times(X);
        } catch (SingularMatrixException e) {
            return true;
        }

        double[] y1 = Y.getRow(0).toArray();
        double[] y2 = Y.getRow(1).toArray();

        final double m21 = moment(y1, y2, 2, 1);
        final double m12 = moment(y1, y2, 1, 2);
        final double m30 = moment(y1, y2, 3, 0);
        final double m03 = moment(y1, y2, 0, 3);

        double m = 500;

        double u3 = m * ((pow(m21, 2) + pow(m12, 2)) / 2.0 + ((pow(m30, 2) + pow(m03, 2)) / 6.0));

//        final double m22 = moment(y1, y2, 2, 2);
//        final double m31 = moment(y1, y2, 3, 1);
//        final double m13 = moment(y1, y2, 1, 3);
//        final double m04 = moment(y1, y2, 0, 4);
//        final double m40 = moment(y1, y2, 4, 0);
//
//        double u4 = m * (pow(m22 - 1, 2) / 4.0 + (pow(m31, 2) + pow(m13, 2)) / 6.0 + (pow(m04 - 1, 2) + pow(m40 - 1, 2)) / 24.0);

//        System.out.println("u3 = " + u3);

        double p = 1.0 - new ChiSquaredDistribution(6).cumulativeProbability(u3);

        return p > alpha;
    }

    private static double[] subset(double[] x, List<Integer> rows, List<Double> distances) {
        double[] d2 = new double[rows.size()];

        for (int i = 0; i < rows.size(); i++) {
            d2[i] = x[rows.get(i)] * kernelEpinechnikov(distances.get(i), 1.0);
        }

        return d2;
    }

    private static double moment(double[] x, double[] y, int m, int n) {
        int N = x.length;
        double sum = 0.0;

        for (int j = 0; j < x.length; j++) {
            sum += Math.pow(x[j], m) * Math.pow(y[j], n);
        }

        return sum / N;
    }

    private List<Integer> allRows() {
        List<Integer> rows = new ArrayList<>();

        for (int i = 0; i < dataSet.getNumRows(); i++) {
            rows.add(i);
        }

        return rows;
    }

    private List<Integer> getCloseZs(double[][] data, int[] _z, int i, int sampleSize) {
        List<Integer> js = new ArrayList<>();

        if (sampleSize > data[0].length) sampleSize = data.length;
        if (_z.length == 0) return allRows();

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

    public static double[] nonparanormal(double[] x1) {
//        double std1 = StatUtils.sd(x1);
//        double mu1 = StatUtils.mean(x1);
        x1 = DataUtils.standardizeData(x1);

        double[] x = ranks(x1);
        int n = x1.length;
        double delta = 1.0 / (4.0 * Math.pow(n, 0.25) * Math.sqrt(Math.PI * Math.log(n)));


        delta = .999;

        for (int i = 0; i < x.length; i++) {
            x[i] /= n;
            if (x[i] > delta) x[i] = delta;
            if (x[i] < (1. - delta)) x[i] = 1. - delta;
            x[i] = normalDistribution.inverseCumulativeProbability(x[i]);
        }

//        double std = StatUtils.sd(x);
//
//        for (int i = 0; i < x.length; i++) {
//            x[i] /= std;
//            x[i] *= std1;
//            x[i] += mu1;
//        }

        return x;
    }

    private static double[] ranks(double[] x) {
        double[] ranks = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            double d = x[i];
            int count = 0;

            for (int k = 0; k < x.length; k++) {
                if (x[k] <= d) {
                    count++;
                }
            }

            ranks[i] = count;
        }

        return ranks;
    }

    private static double width = 0.5;

    private static double kernelEpinechnikov(double z, double h) {
        z /= width * h;
        if (abs(z) > 1) return 0.0;
        else return (/*0.75 **/ (1.0 - z * z));
    }

    // Euclidean distance.
    private static double distance(double[][] data, int[] z, int i, int j) {
        double sum = 0.0;

        for (int _z : z) {
            double d = (data[_z][i] - data[_z][j]) / 2.0;

            if (!Double.isNaN(d)) {
                sum += d * d;
            }
        }

        return sqrt(sum);
    }

    private static double nonparametricFisherZ(double[] _x, double[] _y) {

        // Testing the hypothesis that _x and _y are uncorrelated and assuming that 4th moments of _x and _y
        // are finite and that the sample is large.
        double[] __x = standardize(_x);
        double[] __y = standardize(_y);

        double r = covariance(__x, __y); // correlation
        int N = __x.length;

        // Non-parametric Fisher Z test.
        double z = 0.5 * sqrt(N) * (log(1.0 + r) - log(1.0 - r));

        return z / (sqrt((moment22(__x, __y))));
    }

    private static double moment22(double[] x, double[] y) {
        int N = x.length;
        double sum = 0.0;

        for (int j = 0; j < x.length; j++) {
            sum += x[j] * x[j] * y[j] * y[j];
        }

        return sum / N;
    }

    // Standardizes the given data array. No need to make a copy here.
    private static double[] standardize(double[] data) {
        double sum = 0.0;

        for (double d : data) {
            sum += d;
        }

        double mean = sum / data.length;

        for (int i = 0; i < data.length; i++) {
            data[i] = data[i] - mean;
        }

        double var = 0.0;

        for (double d : data) {
            var += d * d;
        }

        var /= (data.length);
        double sd = sqrt(var);

        for (int i = 0; i < data.length; i++) {
            data[i] /= sd;
        }

        return data;
    }

}



