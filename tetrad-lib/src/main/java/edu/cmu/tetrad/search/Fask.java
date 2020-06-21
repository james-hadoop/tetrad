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
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;

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

    private double[][] data;

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
    private boolean removeResiduals = false;

    // Conditioned correlations are checked to make sure they are different from zero (since if they
    // are zero, the FASK theory doesn't apply).
    private double zeroAlpha = 1.0;
    private boolean twoCycle = false;
    private double lr;

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

        DataSet dataSet = DataUtils.standardizeData(this.dataSet);
        this.data = dataSet.getDoubleData().transpose().toArray();

        List<Node> variables = dataSet.getVariables();
        double[][] colData = dataSet.getDoubleData().transpose().toArray();

        TetradLogger.getInstance().forceLogMessage("FASK v. 2.0");
        TetradLogger.getInstance().forceLogMessage("");
        TetradLogger.getInstance().forceLogMessage("# variables = " + dataSet.getNumColumns());
        TetradLogger.getInstance().forceLogMessage("N = " + dataSet.getNumRows());
        TetradLogger.getInstance().forceLogMessage("Skewness edge threshold = " + skewEdgeThreshold);
        TetradLogger.getInstance().forceLogMessage("2-cycle threshold = " + twoCycleThreshold);
        if (isRemoveResiduals()) {
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

                double lrOrig = leftRightOrig(x, y);//leftRightOrig(x, y) - leftRightOrig(y, x);

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y)) || (skewEdgeThreshold > 0 && abs(lrOrig) > getSkewEdgeThreshold())) {
                    double lrxy1 = leftRightMinnesota(x, y);
                    double lrxy2 = leftRightMinnesota(y, x);
                    this.lr = -lrxy1;

                    if (edgeForbiddenByKnowledge(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge_forbidden"
                                + "\t" + nf.format(lrxy1)
                                + "\t" + X + "<->" + Y
                        );
                        continue;
                    }

                    if (knowledgeOrients(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lrxy1)
                                + "\t" + X + "-->" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lrxy1)
                                + "\t" + X + "<--" + Y
                        );
                        graph.addDirectedEdge(Y, X);
                    } else if (abs(lrxy1) == 0) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t0-coef"
                                + "\t" + nf.format(lrxy1)
                                + "\t" + X + " " + Y
                        );
                    }
//                    else if (abs(lrOrig) < getSkewEdgeThreshold()) {
//                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tNon-edge"
//                                + "\t" + nf.format(abs(lrOrig))
//                                + "\t" + X + "<=>" + Y
//                        );
////                        graph.addDirectedEdge(X, Y);
////                        graph.addDirectedEdge(Y, X);
//                    }
                    else if (abs(lrOrig) < twoCycleThreshold) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t2-cycle"
                                + "\t" + nf.format(abs(lrOrig))
                                + "\t" + X + "<=>" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                        graph.addDirectedEdge(Y, X);
                    } else {
//                        if (lrxy1 > 0 == lrxy2 > 0) {
//                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
//                                    + "\t" + nf.format(lrxy1)
//                                    + "\t" + X + "-->" + Y
//                            );
//                            graph.addUndirectedEdge(X, Y);
//                        }
//                        else
                            if (lrxy1 > 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lrxy1)
                                    + "\t" + X + "-->" + Y
                            );
                            graph.addDirectedEdge(X, Y);
                        } else if (lrxy1 < 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lrxy1)
                                    + "\t" + X + "<--" + Y
                            );
                            graph.addDirectedEdge(Y, X);
                        }
                    }
                }
            }
        }

//        Graph G1 = new EdgeListGraph(graph);
//
//        for (int i = 0; i < variables.size(); i++) {
//            for (int j = i + 1; j < variables.size(); j++) {
//                Node X = variables.get(i);
//                Node Y = variables.get(j);
//
//                // Centered
//                double[] x = colData[i];
//                double[] y = colData[j];
//
//                if (bidirected(x, y, G1, X, Y)) {
//                    graph.removeEdges(X, Y);
//                    graph.addDirectedEdge(X, Y);
//                    graph.addDirectedEdge(Y, X);
//                }
//
//            }
//        }

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    private double leftRight(double[] x, double[] y) {
//        if (true) {
////            return robustSkew(x, y);
//            return skew(x, y);
//        }

        double diff;

        if (isRemoveResiduals()) {
            diff = getDiff(x, y) + getDiff(y, x);
        } else {
            diff = getDiff(y, x) - getDiff(x, y);
        }

        return diff;
    }

    private double leftRightOrig(double[] x, double[] y) {
        final double a = correlation(x, y);

//        if (a < 0) {// && sk_ey > 0) {
//            for (int i = 0; i < x.length; i++) x[i] *= -1;
//        }

//        x = correctSkewness(x);
//        y = correctSkewness(y);

        double[] cx = cov(x, y, x);
        double[] cy = cov(x, y, y);
        double lr = cx[8] - cy[8];

//        if (correlation(x, y) < 0) {
//            lr *= -1;
//        }

        return lr;
    }

    private double leftRightMinnesota(double[] x, double[] y) {
        x = correctSkewness(x);
        y = correctSkewness(y);

        final double cxyx = cov2(x, y, x);
        final double cxyy = cov2(x, y, y);
        final double cxxx = cov2(x, x, x);
        final double cyyx = cov2(y, y, x);
        final double cxxy = cov2(x, x, y);
        final double cyyy = cov2(y, y, y);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        if (sk_ey < 0) {
            lr *= -1;
        }

        final double a = correlation(x, y);

//        if (a < 0) {// && sk_ey > 0) {
//            lr *= -1;
//        }

        return -lr;
    }

    private boolean bidirected(double[] x, double[] y, Graph G0, Node X, Node Y) {
        x = DataUtils.standardizeData(x);
        y = DataUtils.standardizeData(y);

        if (G0.isAdjacentTo(X, Y)) return false;

        Set<Node> adjSet = new HashSet<>(G0.getAdjacentNodes(X));
        adjSet.addAll(G0.getAdjacentNodes(Y));
        List<Node> adj = new ArrayList<>(adjSet);
        adj.remove(X);
        adj.remove(Y);

        depth = 3;

        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), Math.min(depth, adj.size()));
        int[] choice;

        while ((choice = gen.next()) != null) {
            List<Node> _adj = GraphUtils.asList(choice, adj);
            double[][] _Z = new double[_adj.size()][];

            for (int f = 0; f < _adj.size(); f++) {
                Node _z = _adj.get(f);
                int column = dataSet.getColumn(_z);
                _Z[f] = data[column];//dataSet.getDoubleData().getColumn(column).toArray();
            }

            double pc = 0;
            double pc1 = 0;
            double pc2 = 0;

            try {
                pc = partialCorrelation(x, y, _Z, x, Double.NEGATIVE_INFINITY, +1);
                pc1 = partialCorrelation(x, y, _Z, x, 0, +1);
                pc2 = partialCorrelation(x, y, _Z, y, 0, +1);
            } catch (SingularMatrixException e) {
                System.out.println("Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                TetradLogger.getInstance().log("info", "Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                continue;
            }

            int nc = StatUtils.getRows(x, Double.NEGATIVE_INFINITY, +1).size();
            int nc1 = StatUtils.getRows(x, 0, +1).size();
            int nc2 = StatUtils.getRows(y, 0, +1).size();

            double z = 0.5 * (log(1.0 + pc) - log(1.0 - pc));
            double z1 = 0.5 * (log(1.0 + pc1) - log(1.0 - pc1));
            double z2 = 0.5 * (log(1.0 + pc2) - log(1.0 - pc2));

            double zv1 = (z - z1) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc1 - 3)));
            double zv2 = (z - z2) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc2 - 3)));

            boolean rejected1 = abs(zv1) > twoCycleThreshold;
            boolean rejected2 = abs(zv2) > twoCycleThreshold;

            boolean possibleTwoCycle = false;

            if (zv1 < 0 && zv2 > 0 && rejected1) {
                possibleTwoCycle = true;
            } else if (zv1 > 0 && zv2 < 0 && rejected2) {
                possibleTwoCycle = true;
            } else if (rejected1 && rejected2) {
                possibleTwoCycle = true;
            }

            if (!possibleTwoCycle) {
                return false;
            }
        }

        return true;
    }

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition, double threshold, double direction) throws SingularMatrixException {
        double[][] cv = StatUtils.covMatrix(x, y, z, condition, threshold, direction);
        TetradMatrix m = new TetradMatrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }


    private double robustSkew(double[] xData, double[] yData) {
//        if (true) {
//            xData = correctSkewnesses(xData);
//            yData = correctSkewnesses(yData);
//        }

        xData = Arrays.copyOf(xData, xData.length);
        yData = Arrays.copyOf(yData, yData.length);

        double[] xx = new double[xData.length];
        double[] yy = new double[yData.length];

        for (int i = 0; i < xData.length; i++) {
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            double xi = xData[i];
            double yi = yData[i];

            double s1 = g(xi) * yi;
            double s2 = xi * g(yi);

            xx[i] = s1;
            yy[i] = s2;
        }

        double mxx = mean(xx);
        double myy = mean(yy);

        return mxx - myy;
    }

    private double skew(double[] xData, double[] yData) {
//        if (true) {
//            xData = correctSkewnesses(xData);
//            yData = correctSkewnesses(yData);
//        }

        xData = Arrays.copyOf(xData, xData.length);
        yData = Arrays.copyOf(yData, yData.length);

        double[] xx = new double[xData.length];
        double[] yy = new double[yData.length];

        for (int i = 0; i < xData.length; i++) {
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            double xi = xData[i];
            double yi = yData[i];

            double s1 = xi * xi * yi;
            double s2 = xi * yi * yi;

            xx[i] = s1;
            yy[i] = s2;
        }

        double mxx = mean(xx);
        double myy = mean(yy);

        return mxx - myy;
    }

    private double g(double x) {
        return Math.log(Math.cosh(Math.max(x, 0)));
    }

    private double[] residuals(double[] _y, double[][] _x) {
        TetradMatrix y = new TetradMatrix(new double[][]{_y}).transpose();
        TetradMatrix x = new TetradMatrix(_x).transpose();

        TetradMatrix xT = x.transpose();
        TetradMatrix xTx = xT.times(x);
        TetradMatrix xTxInv = xTx.inverse();
        TetradMatrix xTy = xT.times(y);
        TetradMatrix b = xTxInv.times(xTy);

        TetradMatrix yHat = x.times(b);
        if (yHat.columns() == 0) yHat = y.copy();

        return y.minus(yHat).getColumn(0).toArray();
    }

    private double[] correctSkewnesses(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }

    private boolean isNonzeroDiff(double n1, double n2, double c1, double c2, double alpha, Node X, Node Y) {
        double z1 = 0.5 * (log(1.0 + c1) - log(1.0 - c1));
        double z2 = 0.5 * (log(1.0 + c2) - log(1.0 - c2));

        double zdiff = (z1 - z2) / (1.0 / n1 + 1.0 / n2);

        // One sided.
        double p = 1.0 - new NormalDistribution(0, 1)
                .cumulativeProbability(abs(zdiff));

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

    public boolean isRemoveResiduals() {
        return removeResiduals;
    }

    public void setRemoveResiduals(boolean removeResiduals) {
        this.removeResiduals = removeResiduals;
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
        z /= 2.5 * h;
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

    private double getDiff(double[] x, double[] y) {
        double[] yPrime = Arrays.copyOf(y, y.length);

        if (isRemoveResiduals()) {
            if (correlation(x, y) < 0) {
                for (int i = 0; i < x.length; i++) x[i] *= -1;
            }

            double[] r = residuals(x, y, RegressionType.LINEAR);

            for (int k = 0; k < yPrime.length; k++) {
                yPrime[k] -= r[k];
            }

//            yPrime = DataUtils.standardizeData(yPrime);

        }

        x = DataUtils.standardizeData(x);
        yPrime = DataUtils.standardizeData(yPrime);

        double[] rPrime = residuals(x, yPrime, RegressionType.LINEAR);

        rPrime = DataUtils.standardizeData(rPrime);

        double eyrxy = 0.0;
        double eyrxx = 0.0;
        double eyyx = 0.0;
        double eyyy = 0.0;
        int n1 = 0;
        int n2 = 0;

        for (int i = 0; i < yPrime.length; i++) {
            if (x[i] > 0) {
                eyrxx += yPrime[i] * rPrime[i];
                eyyx += yPrime[i] * yPrime[i];
                n1++;
            }

            if (yPrime[i] > 0) {
                eyrxy += yPrime[i] * rPrime[i];
                eyyy += yPrime[i] * yPrime[i];
                n2++;
            }
        }

        eyrxx /= n1;
        eyrxy /= n2;
        eyyx /= n1;
        eyyy /= n2;

//        double[] sums = new double[]{eyrxy, eyrxx, eyyy, eyyx, n2, n1};
        return eyrxy / eyyy - eyrxx / eyyx;
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

    private static double cov2(double[] x, double[] y, double[] condition) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }


    private double[] correctSkewness(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }

    public double getLr() {
        return lr;
    }


}






