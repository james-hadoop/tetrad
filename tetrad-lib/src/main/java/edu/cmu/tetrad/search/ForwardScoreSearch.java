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
        EdgeListGraph graph = new EdgeListGraph(order);
        double score0 = scoreGraph(graph);

        MOVES:
        while (true) {

            for (Edge edge : getAddMoves(graph)) {
                graph.addEdge(edge);
                double _score = scoreGraph(graph);

                if (_score < score0) {
                    score0 = _score;

                    if (verbose) {
                        TetradLogger.getInstance().forceLogMessage("Decreases score (ADD): score = " + _score);
                    }

                    continue MOVES;
                }

                graph.removeEdge(edge);
            }

            for (Edge edge : getRemoveMoves(graph)) {
                graph.removeEdge(edge);
                double score = scoreGraph(graph);

                if (score < score0) {
                    score0 = score;

                    if (verbose) {
                        TetradLogger.getInstance().forceLogMessage("Decreases score (REMOVE): score = " + score);
                    }

                    continue MOVES;
                } else {
                    graph.addEdge(edge);
                }
            }

            break;
        }

        return graph;
    }

    private double scoreGraph(EdgeListGraph graph) {
        return scorer.score(graph);
    }

    private List<Edge> getRemoveMoves(Graph graph) {
        return new ArrayList<>(graph.getEdges());
    }

    private List<Edge> getAddMoves(Graph graph) {
        List<Edge> edges = new ArrayList<>();
        List<Node> nodes = graph.getNodes();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (graph.isAdjacentTo(nodes.get(i), nodes.get(j))) {
                    continue;
                }

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
