///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (c) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
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
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.TetradLogger;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.*;

/**
 * Runs the FASK (Fast Adjacency Skewnmess) algorithm.
 *
 * @author Joseph Ramsey
 */
public final class Fask implements GraphSearch {

    // The score to be used for the FAS adjacency search.
    private final IndependenceTest test;

    // An initial graph to orient, skipping the adjacency step.
    private Graph initialGraph = null;

    // Elapsed time of the search, in milliseconds.
    private long elapsed = 0;

    // The data sets being analyzed. They must all have the same variables and the same
    // number of records.
    private final DataSet dataSet;

    // For the Fast Adjacency Search.
    private int depth = -1;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // A threshold for including extra adjacencies due to skewness. Default is 0 (no skew edges).
    private double skewEdgeThreshold = 0;

    // A theshold for making 2-cycles. Default is 0 (no 2-cycles.)
    private double twoCycleThreshold = 0;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if the nonlinear trend between X and Y should be removed.
    private boolean removeNonlinearTrend = false;

    // Conditioned correlations are checked to make sure they are different from zero (since if they
    // are zero, the FASK theory doesn't apply).
    private double zeroAlpha = 0.01;
    private boolean assumptionsSatisfied = false;
    private boolean twoCycle = false;
    private double lr;
    private double lrP;

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask(DataSet dataSet, IndependenceTest test) {
        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException("For FASK, the dataset must be entirely continuous");
        }

        this.dataSet = dataSet;
        this.test = test;
    }

    public Fask(DataSet dataSet, Graph initialGraph) {
        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException("For FASK, the dataset must be entirely continuous");
        }

        this.dataSet = dataSet;
        this.initialGraph = initialGraph;
        this.test = null;
    }

    //======================================== PUBLIC METHODS ====================================//

    /**
     * Runs the search on the concatenated data, returning a graph, possibly cyclic, possibly with
     * two-cycles. Runs the fast adjacency search (FAS, Spirtes et al., 2000) follows by a modification
     * of the robust skew rule (Pairwise Likelihood Ratios for Estimation of Non-Gaussian Structural
     * Equation Models, Smith and Hyvarinen), together with some heuristics for orienting two-cycles.
     *
     * @return the graph. Some of the edges may be undirected (though it shouldn't be many in most cases)
     * and some of the adjacencies may be two-cycles.
     */
    public Graph search() {
        long start = System.currentTimeMillis();
        NumberFormat nf = new DecimalFormat("0.000");

        DataSet dataSet = DataUtils.center(this.dataSet);

        List<Node> variables = dataSet.getVariables();
        double[][] colData = dataSet.getDoubleData().transpose().toArray();

        TetradLogger.getInstance().forceLogMessage("FASK v. 2.0");
        TetradLogger.getInstance().forceLogMessage("");
        TetradLogger.getInstance().forceLogMessage("# variables = " + dataSet.getNumColumns());
        TetradLogger.getInstance().forceLogMessage("N = " + dataSet.getNumRows());
        TetradLogger.getInstance().forceLogMessage("Skewness edge threshold = " + skewEdgeThreshold);
        TetradLogger.getInstance().forceLogMessage("2-cycle threshold = " + twoCycleThreshold);
        if (isRemoveNonlinearTrend()) {
            TetradLogger.getInstance().forceLogMessage("Removing nonlinear trend");
        }
        TetradLogger.getInstance().forceLogMessage("");

        Graph G0;

        if (getInitialGraph() != null) {
            TetradLogger.getInstance().forceLogMessage("Using initial graph.");

            Graph g1 = new EdgeListGraph(getInitialGraph().getNodes());

            for (Edge edge : getInitialGraph().getEdges()) {
                Node x = edge.getNode1();
                Node y = edge.getNode2();

                if (!g1.isAdjacentTo(x, y)) g1.addUndirectedEdge(x, y);
            }

            g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());

            G0 = g1;
        } else {
            TetradLogger.getInstance().forceLogMessage("Running FAS-Stable, alpha = " + test.getAlpha());

            FasStable fas = new FasStable(test);
            fas.setDepth(getDepth());
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            G0 = fas.search();
        }

        TetradLogger.getInstance().forceLogMessage("");

        SearchGraphUtils.pcOrientbk(knowledge, G0, G0.getNodes());

        Graph graph = new EdgeListGraph(variables);

        TetradLogger.getInstance().forceLogMessage("X\tY\tMethod\tLR\tEdge");

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                // Centered
                double[] x = colData[i];
                double[] y = colData[j];

                double[] corxyx = cov(x, y, x);
                double[] corxyy = cov(x, y, y);

                double c1 = corxyx[1];
                double c2 = corxyy[1];

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y)) || (skewEdgeThreshold > 0 && abs(c1 - c2) > getSkewEdgeThreshold())) {
                    double lrxy;

                    if (isRemoveNonlinearTrend()) {

                        // Will work either way, picking the one that's better for the causal pairs data.
                        // The reason is that reversing gives the opposite direction, but taking the
                        // residuals in the opposite direction reverses it again.
                        lrxy = leftRight(y, x, X, Y);
                    } else {
                        lrxy = leftRight(x, y, X, Y);
                    }

                    if (edgeForbiddenByKnowledge(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge_forbidden"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "<->" + Y
                        );
                        continue;
                    }

                    if (!isAssumptionsSatisfied()) {
                        System.out.println(X + "---" + Y + " Assumptions not satisfied");
//                        continue;
                    }

                    if (knowledgeOrients(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "-->" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "<--" + Y
                        );
                        graph.addDirectedEdge(Y, X);
                    } else if (abs(lrxy) == 0) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t0-coef"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + " " + Y
                        );
                    } else if (abs(lrxy) < twoCycleThreshold) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t2-cycle"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "<=>" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                        graph.addDirectedEdge(Y, X);
                    } else {
                        if (lrxy > 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lrxy)
                                    + "\t" + X + "-->" + Y
                            );
                            graph.addDirectedEdge(X, Y);
                        } else if (lrxy < 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lrxy)
                                    + "\t" + X + "<--" + Y
                            );
                            graph.addDirectedEdge(Y, X);
                        } else {
                            //
                        }
                    }
                }
            }
        }

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    private double leftRight(double[] x, double[] yPlusRx, Node X, Node Y) {
        x = Arrays.copyOf(x, x.length);
        yPlusRx = Arrays.copyOf(yPlusRx, yPlusRx.length);

        x = DataUtils.center(x);
        yPlusRx = DataUtils.center(yPlusRx);

        double b = 1;

        double[] sums = getSums(x, yPlusRx, b);

        double[] covx = cov(x, yPlusRx, x);
        double[] covy = cov(x, yPlusRx, yPlusRx);

        double exxx = covx[10];
        double exxy = covy[10];

        double diff = sums[0] / exxx - sums[1] / exxy;

        boolean assumptionsSatisfied = true;

        double n1 = covx[4];
        double n2 = covy[4];

        double c1 = covx[2];
        double c2 = covy[2];

        // Need to do this first.
        if (isNonzeroDiff(n1, n2, c1, c2, zeroAlpha, X, Y)) {
            assumptionsSatisfied = false;
        }

        if (isNonzeroCoef(correlation(x, yPlusRx), n1, zeroAlpha)) {
            assumptionsSatisfied = false;
        }

        if (isNonzeroSkewness(skewness(x), x.length, zeroAlpha)
                && isNonzeroSkewness(skewness(yPlusRx), yPlusRx.length, zeroAlpha)) {
            assumptionsSatisfied = false;
        }

        this.assumptionsSatisfied = assumptionsSatisfied;

        this.lr = diff;
        return diff;
    }

    private boolean isNonzeroDiff(double n1, double n2, double c1, double c2, double alpha, Node X, Node Y) {
        double z1 = 0.5 * (log(1.0 + c1) - log(1.0 - c1));
        double z2 = 0.5 * (log(1.0 + c2) - log(1.0 - c2));

        double zdiff = (z1 - z2) / (1.0 / n1 + 1.0 / n2);

        // One sided.
        double p = 1.0 - new NormalDistribution(0, 1)
                .cumulativeProbability(abs(zdiff));

        this.setLrP(p);

        return p > alpha;
    }

    private boolean isNonzeroCoef(double r, double n, double alpha) {
        double z = 0.5 * sqrt(n) * (log(1 + r) - log(1 - r));
        double p = 2 * (1 - new NormalDistribution(0, 1).cumulativeProbability(abs(z)));
        return p > alpha;
    }

    private boolean isNonzeroSkewness(double g1, double n, double alpha) {
        double G1 = (sqrt((n * (n - 1)) / (n - 2))) * g1;
        double se = sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)));
        double ratio = G1 / se;
        double p = 2 * (1 - new NormalDistribution(0, 1).cumulativeProbability(abs(ratio)));
        return p > alpha;
    }

    /**
     * @return The depth of search for the Fast Adjacency Search (FAS).
     */
    public int getDepth() {
        return depth;
    }

    /**
     * @param depth The depth of search for the Fast Adjacency Search (S). The default is -1.
     *              unlimited. Making this too high may results in statistical errors.
     */
    public void setDepth(int depth) {
        this.depth = depth;
    }

    /**
     * @return The elapsed time in milliseconds.
     */
    public long getElapsedTime() {
        return elapsed;
    }

    /**
     * @return the current knowledge.
     */
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    /**
     * @param knowledge Knowledge of forbidden and required edges.
     */
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public boolean isRemoveNonlinearTrend() {
        return removeNonlinearTrend;
    }

    public void setRemoveNonlinearTrend(boolean removeNonlinearTrend) {
        this.removeNonlinearTrend = removeNonlinearTrend;
    }

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public double getSkewEdgeThreshold() {
        return skewEdgeThreshold;
    }

    public void setSkewEdgeThreshold(double skewEdgeThreshold) {
        this.skewEdgeThreshold = skewEdgeThreshold;
    }

    public boolean isUseFasAdjacencies() {
        return useFasAdjacencies;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public void setTwoCycleThreshold(double twoCycleThreshold) {
        this.twoCycleThreshold = twoCycleThreshold;
    }

    public void setZeroAlpha(double zeroAlpha) {
        this.zeroAlpha = zeroAlpha;
    }

    public enum RegressionType {LINEAR, NONLINEAR}

    /**
     * Calculates the residuals of y regressed nonparametrically onto y. Left public
     * so it can be accessed separately.
     * <p>
     * Here we want residuals of x regressed onto y. I'll tailor the method to that.
     *
     * @return the nonlinear residuals of y regressed onto x.
     */
    public static double[] residuals(final double[] y, final double[] x, RegressionType regressionType) {
        double[] residuals;

        if (regressionType == RegressionType.LINEAR) {
            RegressionResult result = RegressionDataset.regress(y, new double[][]{x});
            residuals = result.getResiduals().toArray();
        } else {
            int N = y.length;
            residuals = new double[N];
            double[] sum = new double[N];
            double[] totalWeight = new double[N];
            double h = h1(x);

            for (int j = 0; j < N; j++) {
                double yj = y[j];

                for (int i = 0; i < N; i++) {
                    double d = distance(x, i, j);
                    double k = kernelGaussian(d, h);
                    sum[i] += k * yj;
                    totalWeight[i] += k;
                }
            }

            for (int i = 0; i < N; i++) {
                residuals[i] = y[i] - sum[i] / totalWeight[i];
            }
        }

        return residuals;
    }

    //======================================== PRIVATE METHODS ====================================//

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }

    private static double h1(double[] xCol) {
        int N = xCol.length;
        double w;

        if (N < 200) {
            w = 0.8;
        } else if (N < 1200) {
            w = 0.5;
        } else {
            w = 0.3;
        }

        return w;
    }

    private static double distance(double[] data, int i, int j) {
        double sum = 0.0;

        double d = (data[i] - data[j]) / 2.0;

        if (!Double.isNaN(d)) {
            sum += d * d;
        }

        return sqrt(sum);
    }

    private static double kernelGaussian(double z, double h) {
        z /= 3 * h;
        return exp(-z * z);
    }

    private static double[] cov(double[] x, double[] y, double[] condition) {
        double exy = 0.0;
        double exx = 0.0;
        double eyy = 0.0;

        double ex = 0.0;
        double ey = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] > 0) {
                exy += x[k] * y[k];
                exx += x[k] * x[k];
                eyy += y[k] * y[k];
                ex += x[k];
                ey += y[k];
                n++;
            }
        }

        exy /= n;
        exx /= n;
        eyy /= n;
        ex /= n;
        ey /= n;

        double sxy = exy - ex * ey;
        double sx = exx - ex * ex;
        double sy = eyy - ey * ey;

        return new double[]{sxy, sxy / sqrt(sx * sy), sx, sy, (double) n, ex, ey, sxy / sx, exy / sqrt(exx * eyy), exx, eyy};
    }

    private static double[] getSums(double[] yPlusRx, double[] x, double b) {
//        yPlusRx = Arrays.copyOf(yPlusRx, yPlusRx.length);
//        yPlusRx = residuals(yPlusRx, x, RegressionType.NONLINEAR);
//        x = residuals(x, yPlusRx, RegressionType.NONLINEAR);

//        RegressionResult result = RegressionDataset.regress(x, new double[][]{yPlusRx});
        double[] r_y_plux_rx = residuals(x, yPlusRx, RegressionType.LINEAR);// result.getResiduals().toArray();

        double[] y = new double[x.length];

        for (int k = 0; k < y.length; k++) {
            y[k] = yPlusRx[k] - b * r_y_plux_rx[k];
        }

//        RegressionResult result3 =  RegressionDataset.regress(x, new double[][]{y});
        double[] rxy = residuals(x, y, RegressionType.LINEAR);// result3.getResiduals().toArray();

        double sum1 = 0.0;
        double sum2 = 0.0;

        for (int i = 0; i < y.length; i++) {
            if (y[i] > 0) sum1 += y[i] * rxy[i];
            if (x[i] > 0) sum2 += y[i] * rxy[i];
        }

        return new double[]{sum1, sum2};
    }

    public boolean isAssumptionsSatisfied() {
        return assumptionsSatisfied;
    }

    public boolean isTwoCycle() {
        return twoCycle;
    }

    private double nonparametricFisherZ(double[] _x, double[] _y) {

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

    private double moment22(double[] x, double[] y) {
        int N = x.length;
        double sum = 0.0;

        for (int j = 0; j < x.length; j++) {
            sum += x[j] * x[j] * y[j] * y[j];
        }

        return sum / N;
    }

    // Standardizes the given data array. No need to make a copy here.
    private double[] standardize(double[] data) {
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

    public double getLr() {
        return lr;
    }

    public void setLr(double lr) {
        this.lr = lr;
    }

    public double getLrP() {
        return lrP;
    }

    public void setLrP(double lrP) {
        this.lrP = lrP;
    }
}






