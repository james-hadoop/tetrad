package edu.cmu.tetrad.algcomparison.statistic.utils;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

/**
 * A confusion matrix for adjacencies--i.e. TP, FP, TN, FN for counts of adjacencies.
 *
 * @author jdramsey
 */
public class AdjacencyConfusion {
    private final int adjTp;
    private final int adjFp;
    private final int adjFn;
    private final int adjTn;

    public AdjacencyConfusion(Graph truth, Graph est) {
        est = GraphUtils.replaceNodes(est, truth.getNodes());

        int adjTp = 0;
        int adjFp = 0;
        int adjFn = 0;

        for (Edge edge : truth.getEdges()) {
            if (!est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                adjFn++;
            } else {
                adjTp++;
            }
        }

        for (Edge edge : est.getEdges()) {
            if (!truth.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                adjFp++;
            }
        }

        int allEdges = truth.getNumNodes() * (truth.getNumNodes() - 1) / 2;

        this.adjTp = adjTp;
        this.adjFp = adjFp;
        this.adjFn = adjFn;
        this.adjTn = allEdges - adjFn;
    }

    public int getAdjTp() {
        return adjTp;
    }

    public int getAdjFp() {
        return adjFp;
    }

    public int getAdjFn() {
        return adjFn;
    }

    public int getAdjTn() {
        return adjTn;
    }

}
