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
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradLogger;

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
    private boolean removeResiduals = false;

    // Conditioned correlations are checked to make sure they are different from zero (since if they
    // are zero, the FASK theory doesn't apply).
    private double lr;
    private double delta = 0;

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

        Graph graph = new EdgeListGraph(G0.getNodes());

        TetradLogger.getInstance().forceLogMessage("X\tY\tMethod\tLR\tEdge");

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                // Centered
                double[] x = colData[i];
                double[] y = colData[j];

                double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
                double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y)) || (skewEdgeThreshold > 0 && abs(c1 - c2) > getSkewEdgeThreshold())) {
                    double lrxy = leftRight(x, y);
                    double lryx = leftRight(y, x);

                    lr = 0;

                    if (edgeForbiddenByKnowledge(X, Y)) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tknowledge_forbidden"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "<->" + Y
                        );
                        continue;
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
                    }
                    else if (abs(lrxy) < twoCycleThreshold) {
                        TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\t2-cycle"
                                + "\t" + nf.format(lrxy)
                                + "\t" + X + "<=>" + Y
                        );
                        graph.addDirectedEdge(X, Y);
                        graph.addDirectedEdge(Y, X);
                    }
                    else {
                        if (lryx < 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lryx)
                                    + "\t" + X + "-->" + Y
                            );
                            this.lr = lrxy;
                            graph.addDirectedEdge(X, Y);
                        } else if (lrxy < 0) {
                            TetradLogger.getInstance().forceLogMessage(X + "\t" + Y + "\tleft-right"
                                    + "\t" + nf.format(lrxy)
                                    + "\t" + X + "<--" + Y
                            );
                            this.lr = lryx;
                            graph.addDirectedEdge(Y, X);
                        }
                        else {
                            graph.addUndirectedEdge(Y, X);
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
        double skx = skewness(x);
        double sky = skewness(y);
        double r = correlation(x, y);

        double lr = (E(x, y, x) - E(x, y, y));

        if (r * signum(skx) * signum(sky) < getDelta()) {
            lr *= -1;
        }

        return lr;
    }

    private static double E(double[] x, double[] y, double[] condition) {
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

    private double robustSkew(double[] xData, double[] yData) {

        if (true) {
            xData = correctSkewness(xData, skewness(xData));
            yData = correctSkewness(yData, skewness(yData));
        }

        double rho = correlation(xData, yData);

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

        return rho * (mxx - myy);
    }

    private double g(double x) {
        return Math.log(Math.cosh(Math.max(x, 0)));
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

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
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
        z /= 1 * h;
        return exp(-z * z);
    }

    public double getLr() {
        return lr;
    }

    private double[] correctSkewness(double[] data, double sk) {
        data = Arrays.copyOf(data, data.length);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(sk);
        return data2;
    }
}






