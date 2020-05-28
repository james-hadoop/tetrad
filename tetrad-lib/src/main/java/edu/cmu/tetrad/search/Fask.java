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

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.correlation;
import static java.lang.Math.*;

/**
 * Fast adjacency search followed by robust skew orientation. Checks are done for adding
 * two-cycles. The two-cycle checks do not require non-Gaussianity. The robust skew
 * orientation of edges left or right does.
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
    private DataSet dataSet = null;

    // For the Fast Adjacency Search.
    private int depth = -1;

    // For the SEM BIC score, for the Fast Adjacency Search.
    private double penaltyDiscount = 1;

    // Alpha for orienting 2-cycles. Usually needs to be low.
    private double alpha = 1e-6;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // Data as a double[][].
    private final double[][] data;

    // Cutoff for T tests for 2-cycle tests.
    private double cutoff;

    // A threshold for including extra adjacencies due to skewness.
    private double extraEdgeThreshold = 0.3;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if skew adjacencies should be included in the output.
    private boolean useSkewAdjacencies = true;

    private double delta = -0.95;
    private int smoothSkewIntervals = 15;
    private int smoothSkewMinCount = 10;

    private Map<NodePair, Double> confidence = new HashMap<>();
    private boolean removeNonlinearTrend = false;
    private boolean omit = false;

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask(DataSet dataSet, IndependenceTest test) {
        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException("For FASK, the dataset must be entirely continuous");
        }

        this.dataSet = dataSet;
        this.test = test;
        data = dataSet.getDoubleData().transpose().toArray();
    }

    public Fask(DataSet dataSet, Graph initialGraph) {
        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException("For FASK, the dataset must be entirely continuous");
        }

        this.dataSet = dataSet;
        this.initialGraph = initialGraph;
        this.test = null;

        data = dataSet.getDoubleData().transpose().toArray();
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

        setCutoff(alpha);

        DataSet dataSet = DataUtils.standardizeData(this.dataSet);

        List<Node> variables = dataSet.getVariables();
        double[][] colData = dataSet.getDoubleData().transpose().toArray();
        Graph G0;

        if (getInitialGraph() != null) {
            Graph g1 = new EdgeListGraph(getInitialGraph().getNodes());

            for (Edge edge : getInitialGraph().getEdges()) {
                Node x = edge.getNode1();
                Node y = edge.getNode2();

                if (!g1.isAdjacentTo(x, y)) g1.addUndirectedEdge(x, y);
            }

            g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());

            G0 = g1;
        } else {
            FasStable fas = new FasStable(test);
            fas.setDepth(getDepth());
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            G0 = fas.search();

//            Fges fges = new Fges(new SemBicScore(new CovarianceMatrix(dataSet)));
//            fges.setFaithfulnessAssumed(false);
//            fges.setKnowledge(knowledge);
//            G0 = fges.search();
        }

        SearchGraphUtils.pcOrientbk(knowledge, G0, G0.getNodes());

//        System.out.println("Orientation");

        Graph graph = new EdgeListGraph(variables);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                // Centered
                final double[] x = colData[i];
                final double[] y = colData[j];

                double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
                double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y)) || (isUseSkewAdjacencies() && Math.abs(c1 - c2) > getExtraEdgeThreshold())) {
                    if (edgeForbiddenByKnowledge(X, Y)) {
                        // Don't add an edge.
                    } else if (knowledgeOrients(X, Y)) {
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        graph.addDirectedEdge(Y, X);
                    } else if (alpha > 0.0 && bidirected(x, y, G0, X, Y) == 1) {
                        Edge edge1 = Edges.directedEdge(X, Y);
                        Edge edge2 = Edges.directedEdge(Y, X);
                        graph.addEdge(edge1);
                        graph.addEdge(edge2);
                    } else {
                        boolean lrxy;

                        try {
                            lrxy = leftRight2(x, y, X, Y);
                        } catch (Exception e) {
                            return graph;
                        }

                        if (lrxy) {
                            graph.addDirectedEdge(X, Y);
                        } else {
                            graph.addDirectedEdge(Y, X);
                        }
                    }
                }
            }
        }

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

//    public Graph search2() {
//        DataSet dataSet = DataUtils.standardizeData(this.dataSet);
//        double[][] colData = dataSet.getDoubleData().transpose().toArray();
//
//        List<Node> variables = dataSet.getVariables();
//
////        Collections.sort(variables, new Comparator<Node>() {
////            @Override
////            public int compare(Node o1, Node o2) {
////                if (o1 == o2) return 0;
////                int i = variables.indexOf(o1);
////                int j = variables.indexOf(o2);
////                final double[] x = colData[i];
////                final double[] y = colData[j];
////                return leftRightMinnesota(x, y) ? +1 : -1;
////            }
////        });
//
//        Graph graph = new EdgeListGraph(variables);
//
//        for (int i = 0; i < variables.size(); i++) {
//            for  (int j = i + 1; j < variables.size(); j++) {
//                final double[] x = colData[i];
//                final double[] y = colData[j];
//
//                if (leftRightMinnesota(x, y)) {
//                    graph.addDirectedEdge(variables.get(i), variables.get(j));
//                } else  {
//                    graph.addDirectedEdge(variables.get(j), variables.get(i));
//                }
//            }
//        }
//
////        Knowledge2 knowledge = new Knowledge2();
////
////        for (int i = 0; i < variables.size(); i++) {
////            knowledge.addToTier(i + 1, variables.get(i).getName());
////        }
//
//        int numOfNodes = variables.size();
//        for (int i = 0; i < numOfNodes; i++) {
//            for (int j = i + 1; j < numOfNodes; j++) {
//                Node n1 = variables.get(i);
//                Node n2 = variables.get(j);
//
//                if (n1.getName().startsWith("E_") || n2.getName().startsWith("E_")) {
//                    continue;
//                }
//
//                Edge edge = graph.getEdge(n1, n2);
//                if (edge != null && edge.isDirected()) {
//                    knowledge.setForbidden(edge.getNode2().getName(), edge.getNode1().getName());
//                }
//            }
//        }
//
//        final SemBicScore score = new SemBicScore(new CovarianceMatrix(dataSet));
//        score.setPenaltyDiscount(penaltyDiscount);
//        Fges fges = new Fges(score);
//        fges.setKnowledge(knowledge);
//        return fges.search();
//    }

    private int bidirected(double[] x, double[] y, Graph G0, Node X, Node Y) {

        Set<Node> adjSet = new HashSet<>(G0.getAdjacentNodes(X));
        adjSet.addAll(G0.getAdjacentNodes(Y));
        List<Node> adj = new ArrayList<>(adjSet);
        adj.remove(X);
        adj.remove(Y);

        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), Math.min(depth, adj.size()));
        int[] choice;

        while ((choice = gen.next()) != null) {
            List<Node> _adj = GraphUtils.asList(choice, adj);
            double[][] _Z = new double[_adj.size()][];

            for (int f = 0; f < _adj.size(); f++) {
                Node _z = _adj.get(f);
                int column = dataSet.getColumn(_z);
                _Z[f] = data[column];
            }

            double pc = 0;
            double pc1 = 0;
            double pc2 = 0;

            try {
                pc = partialCorrelation(x, y, _Z, x, Double.NEGATIVE_INFINITY, +1);
                pc1 = partialCorrelation(x, y, _Z, x, 0, +1);
                pc2 = partialCorrelation(x, y, _Z, y, 0, +1);
            } catch (SingularMatrixException e) {
                System.out.println(" Singularity");
                return 0;
            } catch (org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException e) {
                System.out.println(" Not positive definite");
                return 0;
            }

            int nc = StatUtils.getRows(x, Double.NEGATIVE_INFINITY, +1).size();
            int nc1 = StatUtils.getRows(x, 0, +1).size();
            int nc2 = StatUtils.getRows(y, 0, +1).size();

            double z = 0.5 * (log(1.0 + pc) - log(1.0 - pc));
            double z1 = 0.5 * (log(1.0 + pc1) - log(1.0 - pc1));
            double z2 = 0.5 * (log(1.0 + pc2) - log(1.0 - pc2));

            double zv1 = (z - z1) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc1 - 3)));
            double zv2 = (z - z2) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc2 - 3)));

            boolean rejected1 = abs(zv1) > cutoff;
            boolean rejected2 = abs(zv2) > cutoff;

            boolean possibleTwoCycle = false;

            if (zv1 < 0 && zv2 > 0 && rejected1) {
                possibleTwoCycle = true;
            } else if (zv1 > 0 && zv2 < 0 && rejected2) {
                possibleTwoCycle = true;
            } else if (rejected1 && rejected2) {
                possibleTwoCycle = true;
            }

            if (!possibleTwoCycle) {
                return -1;
            }
        }

        return 1;
    }

    private boolean leftRight(double[] x, double[] y) {
        double left = cov(x, y, x) / (sqrt(cov(x, x, x) * cov(y, y, x)));
        double right = cov(x, y, y) / (sqrt(cov(x, x, y) * cov(y, y, y)));
        double lr = left - right;

        double r = StatUtils.correlation(x, y);
        double sx = StatUtils.skewness(x);
        double sy = StatUtils.skewness(y);

        r *= signum(sx) * signum(sy);
        lr *= signum(r);
        if (r < delta) lr *= -1;

        return lr > 0;
    }

    private boolean leftRight2(double[] x, double[] y, Node X, Node Y) {
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);

        double a = StatUtils.correlation(x, y);

        if (a < 0) {
            for (int i = 0; i < x.length; i++) x[i] *= -1;
        }

        if (isRemoveNonlinearTrend()) {
            double[] res = residuals(y, x);

            for (int i = 0; i < x.length; i++) {
                x[i] = x[i] - res[i];
            }
        }

        final double cxyx = cov(x, y, x);
        final double cxyy = cov(x, y, y);
        final double cxxx = cov(x, x, x);
        final double cyyx = cov(y, y, x);
        final double cxxy = cov(x, x, y);
        final double cyyy = cov(y, y, y);

        double lr = ((cxyx / sqrt(cxxx * cyyx)) - (cxyy / sqrt(cxxy * cyyy)));

        double n1 = cov2(x, y, x)[1];
        double n2 = cov2(x, y, y)[1];

        double c1 = cxyx / sqrt(cxxx * cyyx);
        double c2 = cxyy / sqrt(cxxy * cyyy);

        double z1 = 0.5 * sqrt(n1) * (log(1 + c1) - log(1 - c1));
        double z2 = 0.5 * sqrt(n2) * (log(1 + c2) - log(1 - c2));

        double zdiff = (z1 - z2) / sqrt((1. / (n1 - 3) + 1. / (n2 - 3)));

        double p = 2.0 * (1 - new TDistribution(x.length - 1)
                .cumulativeProbability(abs(zdiff)));

        confidence.put(new NodePair(X, Y), p);

        return lr > 0;
    }

    private boolean leftRightMinnesota(double[] x, double[] y) {
        x = correctSkewness(x);
        y = correctSkewness(y);

        final double cxyx = cov(x, y, x);
        final double cxyy = cov(x, y, y);
        final double cxxx = cov(x, x, x);
        final double cyyx = cov(y, y, x);
        final double cxxy = cov(x, x, y);
        final double cyyy = cov(y, y, y);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        final double a = correlation(x, y);

        if (a < 0 && sk_ey > 0) {
            lr *= -1;
        }

        return lr > delta;
    }

    private double[] correctSkewness(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }

    private static double cov(double[] x, double[] y, double[] condition) {
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

    private static double[] cov2(double[] x, double[] y, double[] condition) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return new double[]{exy / n, n};
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

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition, double threshold, double direction) throws SingularMatrixException {
        double[][] cv = StatUtils.covMatrix(x, y, z, condition, threshold, direction);
        TetradMatrix m = new TetradMatrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }

    /**
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     */
    private void setCutoff(double alpha) {
        if (alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("Significance out of range: " + alpha);
        }

        this.cutoff = StatUtils.getZForAlpha(alpha);
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
     * @return Returns the penalty discount used for the adjacency search. The default is 1,
     * though a higher value is recommended, say, 2, 3, or 4.
     */
    public double getPenaltyDiscount() {
        return penaltyDiscount;
    }

    /**
     * @param penaltyDiscount Sets the penalty discount used for the adjacency search.
     *                        The default is 1, though a higher value is recommended, say,
     *                        2, 3, or 4.
     */
    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }

    /**
     * @param alpha Alpha for orienting 2-cycles. Needs to be on the low side usually. Default 1e-6.
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
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

    //======================================== PRIVATE METHODS ====================================//

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public double getExtraEdgeThreshold() {
        return extraEdgeThreshold;
    }

    public void setExtraEdgeThreshold(double extraEdgeThreshold) {
        this.extraEdgeThreshold = extraEdgeThreshold;
    }

    public boolean isUseFasAdjacencies() {
        return useFasAdjacencies;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public boolean isUseSkewAdjacencies() {
        return useSkewAdjacencies;
    }

    public void setUseSkewAdjacencies(boolean useSkewAdjacencies) {
        this.useSkewAdjacencies = useSkewAdjacencies;
    }

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

    private static int smoothlySkewed(double[] x, double[] y, int numIntervals, int minCount) {
        double minP = 0;
        double maxX = 1;
        int right = 0;
        int left = 0;

        double interval = (maxX - minP) / numIntervals;

        for (int j = 0; j < numIntervals; j++) {
            double p1 = minP + j * interval;
            double p2 = minP + (j + 1) * interval;

            double top = StatUtils.quantile(y, p2);
            double bottom = StatUtils.quantile(y, p1);

            double max = StatUtils.max(x);
            double min = StatUtils.min(x);

            for (int i = 0; i < numIntervals; i++) {
                double b = (i + 1) * (max(abs(min), abs(max)) / numIntervals);
                int count1 = 0;
                int count2 = 0;

                for (int k = 0; k < x.length; k++) {
                    if (y[k] >= bottom && y[k] <= top && x[k] >= -b && x[k] < 0) {
                        count1++;
                    } else if (y[k] >= bottom && y[k] <= top && x[k] > 0 && x[k] <= b) {
                        count2++;
                    }
                }

                if (count1 > count2 && count1 > minCount) {
                    left++;
                } else if (count2 > count1 && count2 > minCount) {
                    right++;
                }
            }
        }

        return Integer.compare(left, right);
    }

    public static double skewness(double[] x, double[] y, double left, double right) {
        double secondMoment = 0.0;
        double thirdMoment = 0.0;

        double meany = StatUtils.mean(y);

        int count = 0;

        for (int i = 0; i < y.length; i++) {
            if (x[i] < left || x[i] > right) continue;
            if (Double.isNaN(y[i])) continue;
            double s = y[i] - meany;
            if (s == 0) continue;
            count++;
            secondMoment += s * s;
            thirdMoment += s * s * s;
        }

        if (secondMoment == 0) {
            secondMoment = 1e-5;
        }

        double ess = secondMoment / count;
        double esss = thirdMoment / count;

        return esss / Math.pow(ess, 1.5);
    }

    public int getSmoothSkewIntervals() {
        return smoothSkewIntervals;
    }

    public void setSmoothSkewIntervals(int smoothSkewIntervals) {
        this.smoothSkewIntervals = smoothSkewIntervals;
    }

    public int getSmoothSkewMinCount() {
        return smoothSkewMinCount;
    }

    public void setSmoothSkewMinCount(int smoothSkewMinCount) {
        this.smoothSkewMinCount = smoothSkewMinCount;
    }

    private double function(int index, double x) {
        double g = 1.0;

        for (int i = 1; i <= index; i++) {
            g *= x;
        }

        if (abs(g) == Double.POSITIVE_INFINITY) g = Double.NaN;

        return g;
    }

    public double getConfidence(Node X, Node Y) {
        return confidence.get(new NodePair(X, Y));
    }

    /**
     * Calculates the residuals of y regressed nonparametrically onto y. Left public
     * so it can be accessed separately.
     * <p>
     * Here we want residuals of x regressed onto y. I'll tailor the method to that.
     *
     * @return a double[2][] array. The first double[] array contains the residuals for y
     * and the second double[] array contains the resituls for x.
     */
    public static double[] residuals(final double[] y, final double[] x) {

        int N = y.length;

        double[] residuals = new double[N];

        double[] sum = new double[N];

        double[] totalWeight = new double[N];

        double h = h1(x);

        for (int j = 0; j < N; j++) {
            double yj = y[j];

            for (int i = 0; i < N; i++) {
                double d = distance(x, i, j);
                double k = kernelGaussian(d, 5, h);
                sum[i] += k * yj;
                totalWeight[i] += k;
            }
        }

        for (int i = 0; i < N; i++) {
            residuals[i] = y[i] - sum[i] / totalWeight[i];
        }

        return residuals;
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

    private static double kernelGaussian(double z, double width, double h) {
        z /= width * h;
        return Math.exp(-z * z);
    }

    public boolean isRemoveNonlinearTrend() {
        return removeNonlinearTrend;
    }

    public void setRemoveNonlinearTrend(boolean removeNonlinearTrend) {
        this.removeNonlinearTrend = removeNonlinearTrend;
    }

    public boolean isOmit() {
        return omit;
    }

    public void setOmit(boolean omit) {
        this.omit = omit;
    }
}






