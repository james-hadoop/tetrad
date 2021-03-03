package edu.cmu.tetrad.algcomparison.statistic.utils;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.NodePair;

import java.util.HashSet;
import java.util.Set;

/**
 * A confusion matrix for adjacencies--i.e. TP, FP, TN, FN for counts of adjacencies.
 *
 * @author jdramsey
 */
public class AdjacencyConfusion {
    private int adjTp;
    private int adjFp;
    private int adjFn;
    private final int adjTn;

    public AdjacencyConfusion(Graph truth, Graph est) {
        adjTp = 0;
        adjFp = 0;
        adjFn = 0;

        est = GraphUtils.replaceNodes(est, truth.getNodes());

        Set<NodePair> trueEdges = new HashSet<>();
        Set<NodePair> estEdges = new HashSet<>();

        for (Edge edge : truth.getEdges()) {
            NodePair e = new NodePair(edge.getNode1(), edge.getNode2());
            trueEdges.add(e);
        }

        for (Edge edge : est.getEdges()) {
            NodePair e = new NodePair(edge.getNode1(), edge.getNode2());
            estEdges.add(e);
        }

        Set<NodePair> all = new HashSet<>(trueEdges);
        all.addAll(estEdges);

        for (NodePair edge : all) {
            if (estEdges.contains(edge) && !trueEdges.contains(edge)) {
                adjFp++;
            }

            if (!estEdges.contains(edge) && trueEdges.contains(edge)) {
                adjFn++;
            }

            if (estEdges.contains(edge) && trueEdges.contains(edge)) {
                adjTp++;
            }
        }

//        System.out.println("adjTp = " + adjTp + " adjFp = " + adjFp + " adjFn = " + adjFn
//             + " true edge = " + trueEdges.size() + " estEdge = " + estEdges.size() + " all edges = " + all.size());

        int allEdges = truth.getNumNodes() * (truth.getNumNodes() - 1) / 2;

        adjTn = allEdges - adjFn;
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
