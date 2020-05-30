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
import edu.cmu.tetrad.util.StatUtils;

import java.util.Arrays;
import java.util.List;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

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

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // A threshold for including extra adjacencies due to skewness. Defauls is 0 (no extra edges).
    private double skewEdgeThreshold = 0;

    // A theshold for making 2-cycles. Default is 0 (no 2-cycles.)
    private double twoCycleThreshold = 0;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

//    private Map<NodePair, Double> confidence = new HashMap<>();
    private boolean removeNonlinearTrend = false;

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

        Graph graph = new EdgeListGraph(variables);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                // Centered
                double[] x = colData[i];
                double[] y = colData[j];

                x = Arrays.copyOf(x, x.length);
                y = Arrays.copyOf(y, y.length);

                if (isRemoveNonlinearTrend()) {
                    double[] res = residuals(y, x);

                    for (int k = 0; k < x.length; k++) {
                        x[k] = x[k] - res[k];
                    }
                }

                double[] corxyx = StatUtils.cov(x, y, x, 0, +1);
                double[] corxyy = StatUtils.cov(x, y, y, 0, +1);

                double c1 = corxyx[1];
                double c2 = corxyy[1];

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y)) || (skewEdgeThreshold > 0 && abs(c1 - c2) > getSkewEdgeThreshold())) {
                    double lrxy = leftRight(x, y);

                    if (edgeForbiddenByKnowledge(X, Y)) continue;

                    if (knowledgeOrients(X, Y)) {
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        graph.addDirectedEdge(Y, X);
                    } else if (abs(lrxy) < twoCycleThreshold) {
                        graph.addDirectedEdge(X, Y);
                        graph.addDirectedEdge(Y, X);
                    } else {
                        if (lrxy > 0) {
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

    private double leftRight(double[] x, double[] y) {
        double[] covx = StatUtils.cov(x, y, x, 0, +1);
        double[] covy = StatUtils.cov(x, y, y, 0, +1);

        double a = StatUtils.correlation(x, y);

        if (a < 0) {
            for (int i = 0; i < x.length; i++) x[i] *= -1;
        }

        double lr = covx[8] - covy[8];

//        double n1 = covx[4];
//        double n2 = covy[4];
//
//        double s1 = covx[2];
//        double s2 = covy[2];
//
//        double c1 = covx[8];
//        double c2 = covy[8];
//
//        double z1 = 0.5 * (log(1 + c1) - log(1 - c1));
//        double z2 = 0.5 * (log(1 + c2) - log(1 - c2));
//
//        double zdiff = (z1 - z2) / sqrt(((s1 * s1) / (n1 - 3)) + ((s2 * s2) / (n2 - 3)));
//
//        double p = 2.0 * (1 - new NormalDistribution(0, 1)
//                .cumulativeProbability(abs(zdiff)));
//
//        confidence.put(new NodePair(X, Y), 1 - p);

        return lr;
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
//
//    public double getConfidence(Node X, Node Y) {
//        return confidence.get(new NodePair(X, Y));
//    }

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
                double k = kernelGaussian(d, 50, h);
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
}






