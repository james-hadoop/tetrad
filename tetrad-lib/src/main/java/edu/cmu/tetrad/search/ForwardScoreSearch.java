package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.ArrayList;
import java.util.List;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FSS) algorithm.
 *
 * @author josephramsey
 */
public class ForwardScoreSearch {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;
    private boolean verbose = false;

    /**
     * Constructs a FSS search
     *
     * @param scorer the scorer used, by default FML BIC (for linear models). The score
     *               in general should be lower for better models.
     */
    public ForwardScoreSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        EdgeListGraph G0 = new EdgeListGraph(order);
        double s0 = scoreGraph(G0);

//        MOVES:
//        while (true) {
//            for (Edge edge : getRemoveMoves(graph)) {
//                graph.removeEdge(edge);
//                double score = scoreGraph(graph);
//
//                if (score < score0) {
//                    score0 = score;
////                    TetradLogger.getInstance().forceLogMessage("Decreases score (REMOVE): score = " + score);
//                    continue MOVES;
//                }
//
//                graph.addEdge(edge);
//            }
//
//            for (Edge edge : getAddMoves(order)) {
//                if (!graph.containsEdge(edge)) {
//                    graph.addEdge(edge);
//                    double _score = scoreGraph(graph);
//
//                    if (_score < score0) {
//                        score0 = _score;
//                        continue MOVES;
//                    }
//
//                    graph.removeEdge(edge);
//                }
//            }
//
//            break;
//        }

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Edge edge : getAddMoves(order)) {
                if (!G0.containsEdge(edge)) {
                    G0.addEdge(edge);
                    double s = scoreGraph(G0);

                    if (s < s0) {
                        s0 = s;
                        changed = true;
                    } else {
                        G0.removeEdge(edge);
                    }
                }
            }

            for (Edge edge : getRemoveMoves(G0)) {
                G0.removeEdge(edge);
                double s = scoreGraph(G0);

                if (s < s0) {
                    s0 = s;
                    changed = true;
                } else {
                    G0.addEdge(edge);
                }
            }
        }

        return G0;
    }

    private double scoreGraph(EdgeListGraph graph) {
        return scorer.score(graph);
    }

    private List<Edge> getRemoveMoves(Graph graph) {
        return new ArrayList<>(graph.getEdges());
    }

    private List<Edge> getAddMoves(List<Node> nodes) {
        List<Edge> edges = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
                edges.add(edge);
            }
        }

        return edges;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}
