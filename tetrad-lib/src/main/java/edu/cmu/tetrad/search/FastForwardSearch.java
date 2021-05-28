package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FFS) algorithm.
 *
 * @author josephramsey
 */
public class FastForwardSearch {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;
    private double score = Double.NaN;

    /**
     * Constructs a FFS search
     *
     * @param scorer the scorer used, by default FML BIC (for linear models). The score
     *               in general should be lower for better models.
     */
    public FastForwardSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        score = Double.NaN;

        EdgeListGraph G0 = new EdgeListGraph(order);
        double s0 = scorer.score(G0);

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Edge edge : getAddMoves(order)) {
                if (!G0.containsEdge(edge)) {
                    G0.addEdge(edge);
                    double s = scorer.score(edge);

                    if (s < s0) {
                        s0 = s;
                        changed = true;
                    } else {
                        G0.removeEdge(edge);
                        scorer.resetParameters(edge);
                    }
                }
            }

            for (Edge edge : getRemoveMoves(G0)) {
                G0.removeEdge(edge);
                double s = scorer.score(edge);

                if (s < s0) {
                    s0 = s;
                    changed = true;
                } else {
                    G0.addEdge(edge);
                    scorer.resetParameters(edge);
                }
            }
        }

        this.score = s0;

        return G0;
    }

    private List<Edge> getRemoveMoves(Graph graph) {
        return new ArrayList<>(graph.getEdges());
    }

    private List<Move> getMoves(List<Node> nodes, Graph G0) {
        List<Move> moves = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));

                if (G0.containsEdge(edge)) {
                    moves.add(new Move(edge, Type.REMOVE));
                } else {
                    moves.add(new Move(edge, Type.ADD));
                }
            }
        }

        return moves;
    }

    private List<Edge> getAddMoves(List<Node> nodes) {
        List<Edge> edges = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                edges.add(Edges.directedEdge(nodes.get(i), nodes.get(j)));
            }
        }

        return edges;
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }

    private enum Type {ADD, REMOVE}

    private class Move {
        Edge edge;
        Type type;

        public Move(Edge edge, Type type) {
            this.edge = edge;
            this.type = type;
        }
    }

    ;
}
