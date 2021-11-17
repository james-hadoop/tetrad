package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;

import java.util.ArrayList;
import java.util.LinkedList;
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
    private double score = Double.NaN;
    LinkedList<Move> lastMoves = new LinkedList<>();
    LinkedList<Move> currentMoves = new LinkedList<>();

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
        score = Double.NaN;

        lastMoves.clear();
        lastMoves.addAll(currentMoves);
        currentMoves.clear();

        EdgeListGraph G0 = new EdgeListGraph(order);
        double s0 = scorer.score(G0);

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Edge edge : getAddMoves(order, G0)) {
                if (!G0.containsEdge(edge)) {
                    Move move = lastMoves.isEmpty() ? null : lastMoves.removeFirst();

//                    if (move != null && move.getType() == Type.ADD) {
//                        double s = move.getScore();
//                        changed = true;
//                        currentMoves.add(new Move(edge, Type.ADD, s));
//                    } else {
                        G0.addEdge(edge);
                        double s = scorer.score(edge);

                        if (s < s0) {
                            s0 = s;
                            changed = true;
                            currentMoves.add(new Move(edge, Type.ADD, s));
                        } else {
                            G0.removeEdge(edge);
                            scorer.resetParameters(edge);
                        }

                        currentMoves.clear();
//                    }
                }
            }

            for (Edge edge : getRemoveMoves(order, G0)) {
                Move move = lastMoves.isEmpty() ? null : lastMoves.removeFirst();

//                if (move != null && move.getType() == Type.REMOVE) {
//                    double s = move.getScore();
//                    changed = true;
//                    currentMoves.add(new Move(edge, Type.REMOVE, s));
//                } else {
                    G0.removeEdge(edge);
                    double s = scorer.score(edge);

                    if (s < s0) {
                        s0 = s;
                        changed = true;
                        currentMoves.add(new Move(edge, Type.REMOVE, s));
                    } else {
                        G0.addEdge(edge);
                        scorer.resetParameters(edge);
                    }

                    currentMoves.clear();
//                }
            }
        }

        this.score = s0;

        return G0;
    }

    private List<Edge> getRemoveMoves(List<Node> nodes, Graph graph) {
        List<Edge> edges = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node n1 = nodes.get(i);
                Node n2 = nodes.get(j);

                Edge e = Edges.directedEdge(n1, n2);

                if (graph.getNumEdges() != 0 && !(graph.getAdjacentNodes(n1).size() > 1 || graph.getAdjacentNodes(n2).size() > 1)) {
                    continue;
                }

                if (graph.containsEdge(e)) {
                    edges.add(e);
                }
            }
        }

        return edges;
    }

    private List<Edge> getAddMoves(List<Node> nodes, Graph graph) {
        List<Edge> edges = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node n1 = nodes.get(i);
                Node n2 = nodes.get(j);

                if (graph.getNumEdges() != 0 && graph.getAdjacentNodes(n1).isEmpty() && graph.getAdjacentNodes(n2).isEmpty()) {
                    continue;
                }

                Edge e = Edges.directedEdge(n1, n2);

                if (!graph.containsEdge(e)) {
                    edges.add(e);
                }
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

    private static class Move {
        private final Edge edge;
        private final Type type;
        private final double score;

        public Move(Edge edge, Type type, double score) {
            this.edge = edge;
            this.type = type;
            this.score = score;
        }

        public Edge getEdge() {
            return edge;
        }

        public Type getType() {
            return type;
        }

        public double getScore() {
            return score;
        }
    }

    private enum Type{ADD, REMOVE}
}
