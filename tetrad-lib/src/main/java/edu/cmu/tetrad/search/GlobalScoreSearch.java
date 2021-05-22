package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.ArrayList;
import java.util.List;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class GlobalScoreSearch {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;

    /**
     * Constructs a GSS search
     *
     * @param scorer the scorer used, by default FML BIC (for linear models). The score
     *               in general should be lower for better models.
     */
    public GlobalScoreSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The order for the graph.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        Graph graph = new EdgeListGraph(order);
        double score0 = scorer.score(graph);

        boolean changed = true;

        MOVES:
        while (changed) {
            changed = false;

            for (Edge edge : getGraphEdges(graph)) {
                Edge edge2 = edge.reverse();
                Node x = edge2.getNode1();
                Node y = edge2.getNode2();

                List<Node> adj = graph.getAdjacentNodes(x);
                adj.retainAll(graph.getAdjacentNodes(y));

                for (Node z : adj) {

                    {
                        Graph graph1 = new EdgeListGraph(graph);


                        // z x y
                        graph1.removeEdge(x, y);
                        graph1.removeEdge(x, z);
                        graph1.removeEdge(z, y);
                        graph1.addDirectedEdge(z, x);
                        graph1.addDirectedEdge(z, y);

                        double score1 = scorer.score(graph);

                        if (!graph1.existsDirectedCycle()) {
                            System.out.println(edge2 + " Z-X-Y");

                            if (score1 < score0) {
                                graph = graph1;
                                score0 = score1;
                                changed = true;
                                TetradLogger.getInstance().forceLogMessage("Decreases score (FLIP-REMOVE): score = " + score1);
                            }
                        }
                    }

                    if (graph.isAdjacentTo(x, y)) {
                        Graph graph1 = new EdgeListGraph(graph);


                        // x y z
                        graph1.removeEdge(x, y);
                        graph1.removeEdge(x, z);
                        graph1.removeEdge(z, y);
                        graph1.addDirectedEdge(y, z);
                        graph1.addDirectedEdge(x, z);

                        double score1 = scorer.score(graph);

                        if (!graph1.existsDirectedCycle()) {
                            System.out.println(edge2 + " X Y Z");

                            if (score1 <= score0) {
                                graph = graph1;
                                score0 = score1;
                                changed = true;
                                TetradLogger.getInstance().forceLogMessage("Decreases score (FLIP-REMOVE): score = " + score1);
                            }
                        }
                    }

                    if (graph.isAdjacentTo(x, y)) {
                        Graph graph1 = new EdgeListGraph(graph);


                        // x z y
                        graph1.removeEdge(x, y);
                        graph1.removeEdge(x, z);
                        graph1.removeEdge(z, y);
                        graph1.addDirectedEdge(x, z);
                        graph1.addDirectedEdge(z, y);

                        double score1 = scorer.score(graph);

                        if (!graph1.existsDirectedCycle()) {
                            System.out.println(edge2 + " X Z Y");

                            if (score1 <= score0) {
                                graph = graph1;
                                score0 = score1;
                                changed = true;
                                TetradLogger.getInstance().forceLogMessage("Decreases score (FLIP-REMOVE): score = " + score1);
                            }
                        }
                    }

//                        {
//                            Graph graph1 = new EdgeListGraph(graph);
//
//                            // x z y
//                            graph1.removeEdge(x, y);
//                            graph1.removeEdge(x, z);
//                            graph1.removeEdge(z, y);
//                            graph1.addDirectedEdge(x, z);
//                            graph1.addDirectedEdge(z, y);
//
//                            double score1 = scorer.score(graph);
//
//                            if (score1 < score0) {
//                                graph = graph1;
//                                score0 = score1;
//                                changed = true;
//                            }
//                        }
                }

            }


            for (Edge edge : getGraphEdges(graph)) {
                graph.removeEdge(edge);
                double score = scorer.score(graph);

                if (score <= score0) {
                    score0 = score;
                    changed = true;
                    TetradLogger.getInstance().forceLogMessage("Decreases score (REMOVE): score = " + score);
                    continue MOVES;
                } else {
                    graph.addEdge(edge);
                }
            }


            for (Edge edge : getAddEdges(graph)) {
                Edge _edge = edge;

                double score1 = score(_edge, graph);
                double score2 = score(_edge.reverse(), graph);

                if (score2 < score1) {
                    _edge = _edge.reverse();
                    score1 = score2;
                }

                if (score1 < score0) {
                    graph.removeEdge(_edge.getNode1(), _edge.getNode2());
                    graph.addEdge(_edge);
                    score0 = score1;
                    changed = true;
                    TetradLogger.getInstance().forceLogMessage("Decreases score (ADD): score = " + score1);
                    continue  MOVES;
                } else {
                    graph.removeEdge(_edge.getNode1(), _edge.getNode2());
                }
            }
        }

        return graph;
    }

    private double score(Edge edge, Graph graph) {
        Edge edge0 = graph.getEdge(edge.getNode1(), edge.getNode2());

        if (edge0 != null) {
            graph.removeEdge(edge0);
        }

        if (graph.existsDirectedPathFromTo(edge.getNode2(), edge.getNode1())) {
            return Double.POSITIVE_INFINITY;
        }

        graph.addEdge(edge);
        double v = scorer.score(graph);
        graph.removeEdge(edge);

        if (edge0 != null) {
            graph.addEdge(edge0);
        }

        return v;
    }

    private Graph permute(Graph graph, List<Node> order) {
        Graph permuted = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            for (int j = i + 1; j < order.size(); j++) {
                if (graph.isAdjacentTo(order.get(i), order.get(j))) {
                    permuted.addDirectedEdge(order.get(i), order.get(j));
                }
            }
        }

        return permuted;
    }

    private List<Node> reverse(List<Node> order, Node n1, Node n2) {
        order = new ArrayList<>(order);

        int i1 = order.indexOf(n1);
        int i2 = order.indexOf(n2);

        order.remove(i1);
        order.add(i1, n2);

        order.remove(i2);
        order.add(i2, n1);

        return order;
    }

    private List<Edge> getGraphEdges(Graph graph) {
        return new ArrayList<>(graph.getEdges());
    }

    private List<Edge> getAddEdges(Graph graph) {
        List<Edge> edges = new ArrayList<>();
        List<Node> nodes = graph.getNodes();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (i == j) continue;

                if (graph.isAdjacentTo(nodes.get(i), nodes.get(j))) {
                    continue;
                }

                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
                edges.add(edge);
            }
        }

        return edges;
    }


}
