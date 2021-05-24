package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.ArrayList;
import java.util.List;

import static java.util.Collections.shuffle;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;

    /**
     * Constructs a GSS search
     *
     * @param scorer the scorer used, by default FML BIC (for linear models). The score
     *               in general should be lower for better models.
     */
    public BestOrderScoreSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The order for the graph.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        ForwardScoreSearch fss = new ForwardScoreSearch(scorer);

        List<Node> bestOrder = new ArrayList<>(order);
        Graph graph = fss.search(bestOrder);
//        double score = Double.POSITIVE_INFINITY;

//        for (int i = 0; i < 5; i++) {
//            List<Node> _order = new ArrayList<>(bestOrder);
//            shuffle(_order);
//            Graph _graph = fss.search(_order);
//            double _score = scorer.score(_graph);
//
//            if (_score < score) {
//                bestOrder = new ArrayList<>(_order);
//                graph = _graph;
//                score = _score;
//            }
//        }

        double score = scorer.score(graph);

        MOVES:
        while (true) {
            for (Edge edge : graph.getEdges()) {
//                if (graph.existsDirectedPathFromTo(edge.getNode2(), edge.getNode1())) {
//                    continue;
//                }

                List<Node> _order = reverse(bestOrder, edge.getNode1(), edge.getNode2());
                Graph _graph = fss.search(_order);
                double _score = scorer.score(_graph);

                if (_score < score) {
                    bestOrder = _order;
                    graph = _graph;
                    score = _score;
                    TetradLogger.getInstance().forceLogMessage("Decreases score: score = " + _score);
                    continue MOVES;
                }
            }

            break;
        }

        return graph;
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
}
