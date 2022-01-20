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

import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.GRASP;
import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradLogger;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.*;

/**
 * Runs the FASK (Fast Adjacency Skewness) algorithm. The reference is Sanchez-Romero, R., Ramsey, J. D.,
 * Zhang, K., Glymour, M. R., Huang, B., & Glymour, C. (2019). Estimating feedforward and feedback
 * effective connections from fMRI time series: Assessments of statistical methods. Network Neuroscience,
 * 3(2), 274-306, though it has been improved in some ways from that version, and some pairwise methods from
 * Hyvärinen, A., & Smith, S. M. (2013). Pairwise likelihood ratios for estimation of non-Gaussian structural
 * equation models. Journal of Machine Learning Research, 14(Jan), 111-152 have been included for
 * comparison (and potential use!--they are quite good!).
 * <p>
 * This method (and the Hyvarinen and Smith methods) make the assumption that the data are generated by
 * a linear, non-Gaussian causal process and attempts to recover the causal graph for that process. They
 * do not attempt to recover the parametrization of this graph; for this a separate estimation algorithm
 * would be needed, such as linear regression regressing each node onto its parents. A further assumption
 * is made, that there are no latent common causes of the algorithm. This is not a constraint on the pairwise
 * orientation methods, since they orient with respect only to the two variables at the endpoints of an edge
 * and so are happy with all other variables being considered latent with respect to that single edge. However,
 * if the built-in adjacency search is used (FAS-Stable), the existence of latents will throw this method
 * off.
 * <p>
 * As was shown in the Hyvarinen and Smith paper above, FASK works quite well even if the graph contains
 * feedback loops in most configurations, including 2-cycles. 2-cycles can be detected fairly well if the
 * FASK left-right rule is selected and the 2-cycle threshold set to about 0.1--more will be detected (or
 * hallucinated) if the threshold is set higher. As shown in the Sanchez-Romero reference above, 2-cycle
 * detection of the FASK algorithm using this rule is quite good.
 * <p>
 * Some edges may be undiscoverable by FAS-Stable; to recover more of these edges, a test related to the
 * FASK left-right rule is used, and there is a threshold for this test. A good default for this threshold
 * (the "skew edge threshold") is 0.3. For more of these edges, set this threshold to a lower number.
 * <p>
 * It is assumed that the data are arranged so the each variable forms a column and that there are no missing
 * values. The data matrix is assumed to be rectangular. To this end, the Tetrad DataSet class is used, which
 * enforces this.
 * <p>
 * Note that orienting a DAG for a linear, non-Gaussian model using the Hyvarinen and Smith pairwise rules
 * is alternatively known in the literature as Pairwise LiNGAM--see Hyvärinen, A., & Smith, S. M. (2013). Pairwise
 * likelihood ratios for estimation of non-Gaussian structural equation models. Journal of Machine Learning Research,
 * 14(Jan), 111-152. We include some of these methods here for comparison.
 *
 * @author Joseph Ramsey
 */
public final class Fask implements GraphSearch {

    // The method to use for finding the adjacencies.
    public enum AdjacencyMethod {FAS_STABLE, FAS_STABLE_CONCURRENT, FGES, EXTERNAL_GRAPH}

    // The left-right rule to use. Options include the FASK left-right rule and three left-right rules
    // from the Hyvarinen and Smith pairwise orientation paper: Robust Skew, Skew, and Tanh. In that
    // paper, "empirical" versions were given in which the variables are multiplied through by the
    // signs of the skewnesses; we follow this advice here (with good results). These others are provided
    // for comparison; in general they are quite good.
    public enum LeftRight {FASK1, FASK2, RSKEW, SKEW, TANH}

    // The score to be used for the FAS adjacency search.
    private final IndependenceTest test;

    // An initial graph to constrain the adjacency step.
    private Graph externalGraph = null;

    // Elapsed time of the search, in milliseconds.
    private long elapsed = 0;

    // The data sets being analyzed. They must all have the same variables and the same
    // number of records.
    private final DataSet dataSet;

    // For the Fast Adjacency Search, the maximum number of edges in a conditioning set.
    private int depth = -1;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // A threshold for including extra adjacencies due to skewness. Default is 0.3. For more edges, lower
    // this threshold.
    private double skewEdgeThreshold = 0;

    // A theshold for making 2-cycles. Default is 0 (no 2-cycles.) Note that the 2-cycle rule will only work
    // with the FASK left-right rule. Default is 0; a good value for finding a decent set of 2-cycles is 0.1.
    private double twoCycleScreeningCutoff = 0;

    // At the end of the procedure, two cycles marked in the graph (for having small LR differences) are then
    // tested statisstically to see if they are two-cycles, using this cutoff. To adjust this cutoff, set the
    // two cycle alpha to a number in [0, 1]. The default alpha  is 0.01.
    private double orientationCutoff;

    // The corresponding alpha.
    private double orientationAlpha;

    // Bias for orienting with negative coefficients.
    private double delta;

    // Whether X and Y should be adjusted for skewness. (Otherwise, they are assumed to have positive skewness.
    private boolean empirical = true;

    // True if FAS adjacencies should be included in the output, by default true.
    private boolean useFasAdjacencies = true;

    // By default, FAS Stable will be used for adjacencies, though this can be set.
    private AdjacencyMethod adjacencyMethod = AdjacencyMethod.FAS_STABLE;

    // The left right rule to use, default FASK.
    private LeftRight leftRight = LeftRight.RSKEW;

    // The graph resulting from search.
    private Graph graph;

    // Used for calculating coefficient values.
    private final RegressionDataset regressionDataset;

    double[][] D;

    /**
     * @param dataSet A continuous dataset over variables V.
     * @param test    An independence test over variables V. (Used for FAS.)
     */
    public Fask(DataSet dataSet, IndependenceTest test) {
        if (dataSet == null) {
            throw new NullPointerException("Data set not provided.");
        }

        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException("For FASK, the dataset must be entirely continuous");
        }

        this.dataSet = dataSet;
        this.test = test;

        regressionDataset = new RegressionDataset(dataSet);
        this.orientationCutoff = StatUtils.getZForAlpha(0.01);
        this.orientationAlpha = 0.01;
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

        List<Node> variables = dataSet.getVariables();
        double[][] lrs = getLrScores(); // Sets D.
//        D = dataSet.getDoubleData().transpose().toArray();

        for (int i = 0; i < variables.size(); i++) {
            System.out.println("Skewness of " + variables.get(i) + " = " + skewness(D[i]));
        }

        TetradLogger.getInstance().forceLogMessage("FASK v. 2.0");
        TetradLogger.getInstance().forceLogMessage("");
        TetradLogger.getInstance().forceLogMessage("# variables = " + dataSet.getNumColumns());
        TetradLogger.getInstance().forceLogMessage("N = " + dataSet.getNumRows());
        TetradLogger.getInstance().forceLogMessage("Skewness edge threshold = " + skewEdgeThreshold);
        TetradLogger.getInstance().forceLogMessage("Orientation Alpha = " + orientationAlpha);
        TetradLogger.getInstance().forceLogMessage("2-cycle threshold = " + twoCycleScreeningCutoff);
        TetradLogger.getInstance().forceLogMessage("");

        Graph G;

        if (adjacencyMethod == AdjacencyMethod.FAS_STABLE) {
            Fas fas = new Fas(test);
            fas.setStable(true);
            fas.setDepth(getDepth());
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            G = fas.search();
        } else if (adjacencyMethod == AdjacencyMethod.FAS_STABLE_CONCURRENT) {
            Grasp otherPermAlgs = new Grasp(new LinearGaussianBicScore(dataSet));
            otherPermAlgs.setKnowledge(knowledge);
            List<Node> order = otherPermAlgs.bestOrder(variables);
            G = otherPermAlgs.getGraph(false);

//            FasConcurrent fas = new FasConcurrent(test);
//            fas.setStable(true);
//            fas.setVerbose(false);
//            fas.setKnowledge(knowledge);
//            G = fas.search();
        } else if (adjacencyMethod == AdjacencyMethod.FGES) {
            Fges fas = new Fges(new ScoredIndTest(test));
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            G = fas.search();
        } else if (adjacencyMethod == AdjacencyMethod.EXTERNAL_GRAPH) {
            if (getExternalGraph() == null) throw new IllegalStateException("An external graph was not supplied.");

            Graph g1 = new EdgeListGraph(getExternalGraph().getNodes());

            for (Edge edge : getExternalGraph().getEdges()) {
                Node x = edge.getNode1();
                Node y = edge.getNode2();

                if (!g1.isAdjacentTo(x, y)) g1.addUndirectedEdge(x, y);
            }

            g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());

            G = g1;
        } else {
            throw new IllegalStateException("That method was not configured: " + adjacencyMethod);
        }

        G = GraphUtils.replaceNodes(G, dataSet.getVariables());

        TetradLogger.getInstance().forceLogMessage("");

        SearchGraphUtils.pcOrientbk(knowledge, G, G.getNodes());

        Graph graph = new EdgeListGraph(G.getNodes());

        TetradLogger.getInstance().forceLogMessage("X\tY\tMethod\tLR\tEdge");

        int V = variables.size();

        List<NodePair> twoCycles = new ArrayList<>();

        for (int i = 0; i < V; i++) {
            for (int j = i + 1; j < V; j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                // Centered
                double[] x = D[i];
                double[] y = D[j];

                double cx = correxp(x, y, x);
                double cy = correxp(x, y, y);

                if (G.isAdjacentTo(X, Y) || (abs(cx - cy) > skewEdgeThreshold)) {
                    double lr = lrs[i][j];// leftRight(x, y);

                    if (edgeForbiddenByKnowledge(X, Y) && edgeForbiddenByKnowledge(Y, X)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge_forbidden"
                                + "\t" + nf.format(lr)
                                + "\t" + X + "<->" + Y
                        );
                        continue;
                    }

                    if (knowledgeOrients(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lr)
                                + "\t" + X + "-->" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge"
                                + "\t" + nf.format(lr)
                                + "\t" + X + "<--" + Y
                        );
                        graph.addDirectedEdge(Y, X);
                    } else {
                        if (zeroDiff(i, j, D)) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t2-cycle Prescreen"
                                    + "\t" + nf.format(lr)
                                    + "\t" + X + "...TC?..." + Y
                            );

                            System.out.println(X + " " + Y + " lr = " + lr + " zero");
                            continue;
                        }

                        if (twoCycleScreeningCutoff > 0 && abs(faskLeftRightV2(x, y)) < twoCycleScreeningCutoff) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t2-cycle Prescreen"
                                    + "\t" + nf.format(lr)
                                    + "\t" + X + "...TC?..." + Y
                            );

                            twoCycles.add(new NodePair(X, Y));
                            System.out.println(X + " " + Y + " lr = " + lr + " zero");
                        }

                        if (lr > 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lr)
                                    + "\t" + X + "-->" + Y
                            );
                            graph.addDirectedEdge(X, Y);
                        } else if (lr < 0) {
                            TetradLogger.getInstance().forceLogMessage(Y + "\t" + X + "\tleft-right"
                                    + "\t" + nf.format(lr)
                                    + "\t" + Y + "-->" + X
                            );
                            graph.addDirectedEdge(Y, X);
                        }
                    }
                }
            }
        }

        if (twoCycleScreeningCutoff > 0 && orientationAlpha == 0) {
            for (NodePair edge : twoCycles) {
                Node X = edge.getFirst();
                Node Y = edge.getSecond();

                graph.removeEdges(X, Y);
                graph.addDirectedEdge(X, Y);
                graph.addDirectedEdge(Y, X);
                logTwoCycle(nf, variables, D, X, Y, "2-cycle Pre-screen");
            }
        } else if (twoCycleScreeningCutoff == 0 && orientationAlpha > 0) {
            for (Edge edge : graph.getEdges()) {
                Node X = edge.getNode1();
                Node Y = edge.getNode2();

                int i = variables.indexOf(X);
                int j = variables.indexOf(Y);

                if (twoCycleTest(i, j, D, graph, variables)) {
                    graph.removeEdges(X, Y);
                    graph.addDirectedEdge(X, Y);
                    graph.addDirectedEdge(Y, X);
                    logTwoCycle(nf, variables, D, X, Y, "2-cycle Tested");
                }
            }
        } else if (twoCycleScreeningCutoff > 0 && orientationAlpha > 0) {
            for (NodePair edge : twoCycles) {
                Node X = edge.getFirst();
                Node Y = edge.getSecond();

                int i = variables.indexOf(X);
                int j = variables.indexOf(Y);

                if (twoCycleTest(i, j, D, graph, variables)) {
                    graph.removeEdges(X, Y);
                    graph.addDirectedEdge(X, Y);
                    graph.addDirectedEdge(Y, X);
                    logTwoCycle(nf, variables, D, X, Y, "2-cycle Screened then Tested");
                }
            }
        }

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        this.graph = graph;

        double bic = new BicEst().getValue(null, graph, dataSet);
        graph.addAttribute("BIC", bic);

        return graph;
    }

    private void logTwoCycle(NumberFormat nf, List<Node> variables, double[][] d, Node X, Node Y, String type) {
        int i = variables.indexOf(X);
        int j = variables.indexOf(Y);

        double[] x = d[i];
        double[] y = d[j];

        double lr = leftRight(x, y);

        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t" + type
                + "\t" + nf.format(lr)
                + "\t" + X + "<=>" + Y
        );
    }

    /**
     * Returns the coefficient matrix for the search. If the search has not yet run, runs it,
     * then estimates coefficients of each node given its parents using linear regression and forms
     * the B matrix of coefficients from these estimates. B[i][j] != 0 means i->j with that coefficient.
     */
    public double[][] getB() {
        if (graph == null) search();

        List<Node> nodes = dataSet.getVariables();
        double[][] B = new double[nodes.size()][nodes.size()];

        for (int j = 0; j < nodes.size(); j++) {
            Node y = nodes.get(j);

            List<Node> pary = graph.getParents(y);
            RegressionResult result = regressionDataset.regress(y, pary);
            double[] coef = result.getCoef();

            for (int i = 0; i < pary.size(); i++) {
                B[nodes.indexOf(pary.get(i))][j] = coef[i + 1];
            }
        }

        return B;
    }

    /**
     * Returns a natrux matrix of left-right scores for the search. If lr = getLrScores(), then
     * lr[i][j] is the left right scores leftRight(data[i], data[j]);
     */
    public double[][] getLrScores() {
        List<Node> variables = dataSet.getVariables();
        double[][] D = DataUtils.standardizeData(dataSet).getDoubleData().transpose().toArray();

        double[][] lr = new double[variables.size()][variables.size()];

        for (int i = 0; i < variables.size(); i++) {
            for (int j = 0; j < variables.size(); j++) {
                lr[i][j] = leftRight(D[i], D[j]);
            }
        }

        this.D = D;

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

    public Graph getExternalGraph() {
        return externalGraph;
    }

    public void setExternalGraph(Graph externalGraph) {
        this.externalGraph = externalGraph;
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

    public void setTwoCycleScreeningCutoff(double twoCycleScreeningCutoff) {
        if (twoCycleScreeningCutoff < 0)
            throw new IllegalStateException("Two cycle screening threshold must be >= 0");
        this.twoCycleScreeningCutoff = twoCycleScreeningCutoff;
    }

    public void setOrientationAlpha(double orientationAlpha) {
        if (orientationAlpha < 0 || orientationAlpha > 1)
            throw new IllegalArgumentException("Two cycle testing alpha should be in [0, 1].");
        this.orientationCutoff = StatUtils.getZForAlpha(orientationAlpha);
        this.orientationAlpha = orientationAlpha;
    }

    public void setLeftRight(LeftRight leftRight) {
        this.leftRight = leftRight;
    }

    public void setAdjacencyMethod(AdjacencyMethod adjacencyMethod) {
        this.adjacencyMethod = adjacencyMethod;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

    public void setEmpirical(boolean empirical) {
        this.empirical = empirical;
    }

    public double leftRight(Node X, Node Y) {
        List<Node> variables = dataSet.getVariables();

        int i = -1;

        for (int k = 0; k < variables.size(); k++) {
            if (X.getName().equals(variables.get(k).getName())) i = k;
        }

        int j = -1;

        for (int k = 0; k < variables.size(); k++) {
            if (Y.getName().equals(variables.get(k).getName())) j = k;
        }

        double[] x = D[i];
        double[] y = D[j];

        return leftRight(x, y);

    }


    //======================================== PRIVATE METHODS ====================================//

    private double leftRight(double[] x, double[] y) {
        if (leftRight == LeftRight.FASK1) {
            return faskLeftRightV1(x, y);
        } else if (leftRight == LeftRight.FASK2) {
            return faskLeftRightV2(x, y);
        } else if (leftRight == LeftRight.RSKEW) {
            return robustSkew(x, y);
        } else if (leftRight == LeftRight.SKEW) {
            return skew(x, y);
        } else if (leftRight == LeftRight.TANH) {
            return tanh(x, y);
        }

        throw new IllegalStateException("Left right rule not configured: " + leftRight);
    }

    private double faskLeftRightV2(double[] x, double[] y) {
        double sx = skewness(x);
        double sy = skewness(y);
        double r = correlation(x, y);
        double lr = correxp(x, y, x) - correxp(x, y, y);

        if (empirical) {
            lr *= signum(sx) * signum(sy);
        }

//        lr *= signum(r);

        if (r < delta) {
            lr *= -1;
        }

        return lr;
    }

    private double faskLeftRightV1(double[] x, double[] y) {
        double left = cu(x, y, x) / (sqrt(cu(x, x, x) * cu(y, y, x)));
        double right = cu(x, y, y) / (sqrt(cu(x, x, y) * cu(y, y, y)));
        double lr = left - right;

        double r = StatUtils.correlation(x, y);
        double sx = StatUtils.skewness(x);
        double sy = StatUtils.skewness(y);

        if (empirical) {
            r *= signum(sx) * signum(sy);
        }

        lr *= signum(r);
        if (r < delta) lr *= -1;

        return lr;
    }

    private static double cu(double[] x, double[] y, double[] condition) {
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

    private double robustSkew(double[] x, double[] y) {

        if (empirical) {
            x = correctSkewness(x, skewness(x));
            y = correctSkewness(y, skewness(y));
        }

        double[] lr = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            lr[i] = g(x[i]) * y[i] - x[i] * g(y[i]);
        }

        return correlation(x, y) * mean(lr);
    }

    private double skew(double[] x, double[] y) {

        if (empirical) {
            x = correctSkewness(x, skewness(x));
            y = correctSkewness(y, skewness(y));
        }

        double[] lr = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            lr[i] = x[i] * x[i] * y[i] - x[i] * y[i] * y[i];
        }

        return correlation(x, y) * mean(lr);
    }

    private double tanh(double[] x, double[] y) {

        if (empirical) {
            x = correctSkewness(x, skewness(x));
            y = correctSkewness(y, skewness(y));
        }

        double[] lr = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            lr[i] = x[i] * Math.tanh(y[i]) - Math.tanh(x[i]) * y[i];
        }

        return correlation(x, y) * mean(lr);
    }

    private double g(double x) {
        return Math.log(Math.cosh(Math.max(x, 0)));
    }

    private boolean knowledgeOrients(Node X, Node Y) {
        return knowledge.isForbidden(Y.getName(), X.getName()) || knowledge.isRequired(X.getName(), Y.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node X, Node Y) {
        return knowledge.isForbidden(Y.getName(), X.getName()) && knowledge.isForbidden(X.getName(), Y.getName());
    }

    // Returns E(XY | Z > 0) / sqrt(E(XX | Z > 0) * E(YY | Z > 0)). Z is typically either X or Y.
    private static double correxp(double[] x, double[] y, double[] z) {
        return E(x, y, z) / sqrt(E(x, x, z) * E(y, y, z));
    }

    // Returns E(XY | Z > 0); Z is typically either X or Y.
    private static double E(double[] x, double[] y, double[] z) {
        double exy = 0.0;
        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (z[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }

    private double[] correctSkewness(double[] data, double sk) {
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * signum(sk);
        return data2;
    }

    private boolean twoCycleTest(int i, int j, double[][] D, Graph G0, List<Node> V) {
        Node X = V.get(i);
        Node Y = V.get(j);

        double[] x = D[i];
        double[] y = D[j];

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
                _Z[f] = D[column];
            }

            double pc, pc1, pc2;

            try {
                pc = partialCorrelation(x, y, _Z, x, Double.NEGATIVE_INFINITY);
                pc1 = partialCorrelation(x, y, _Z, x, 0);
                pc2 = partialCorrelation(x, y, _Z, y, 0);
            } catch (SingularMatrixException e) {
                System.out.println("Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                TetradLogger.getInstance().log("info", "Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                continue;
            }

            int nc = StatUtils.getRows(x, x, 0, Double.NEGATIVE_INFINITY).size();
            int nc1 = StatUtils.getRows(x, x, 0, +1).size();
            int nc2 = StatUtils.getRows(y, y, 0, +1).size();

            double z = 0.5 * (log(1.0 + pc) - log(1.0 - pc));
            double z1 = 0.5 * (log(1.0 + pc1) - log(1.0 - pc1));
            double z2 = 0.5 * (log(1.0 + pc2) - log(1.0 - pc2));

            double zv1 = (z - z1) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc1 - 3)));
            double zv2 = (z - z2) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc2 - 3)));

            boolean rejected1 = abs(zv1) > orientationCutoff;
            boolean rejected2 = abs(zv2) > orientationCutoff;

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

    private boolean zeroDiff(int i, int j, double[][] D) {
        double[] x = D[i];
        double[] y = D[j];

        double pc1, pc2;

        try {
            pc1 = partialCorrelation(x, y, new double[0][], x, 0);
            pc2 = partialCorrelation(x, y, new double[0][], y, 0);
        } catch (SingularMatrixException e) {
            throw new RuntimeException(e);
        }

        int nc1 = StatUtils.getRows(x, x, 0, +1).size();
        int nc2 = StatUtils.getRows(y, y, 0, +1).size();

        double z1 = 0.5 * (log(1.0 + pc1) - log(1.0 - pc1));
        double z2 = 0.5 * (log(1.0 + pc2) - log(1.0 - pc2));

        double zv = (z1 - z2) / sqrt((1.0 / ((double) nc1 - 3) + 1.0 / ((double) nc2 - 3)));

        return abs(zv) <= twoCycleScreeningCutoff;
    }

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition, double threshold) throws SingularMatrixException {
        double[][] cv = StatUtils.covMatrix(x, y, z, condition, threshold, 1);
        Matrix m = new Matrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }
}






