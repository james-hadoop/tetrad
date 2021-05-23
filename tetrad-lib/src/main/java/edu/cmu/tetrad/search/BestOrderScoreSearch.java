package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.ArrayList;
import java.util.List;

/**
 * Searches for a DAG by looking at select permutations of the variables and applying FSS. The DAG
 * returned is the one that optimizes the BIC score over these calls to FSS. Neighboring permutatinos
 * are considered for edges X->Y in the current graph, where reversing X and Y in the permutation and
 * making a call to FSS results in a better scoring graph.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;

    /**
     * Constructs a BOSS search
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
     * @param variables The variables for the graph.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> variables) {
        ForwardScoreSearch fss = new ForwardScoreSearch(scorer);

        List<Node> bestOrder = new ArrayList<>(variables);
        Graph graph = fss.search(bestOrder);
        double score = Double.POSITIVE_INFINITY;

        WHILE:
        while (true) {
            for (Edge edge : graph.getEdges()) {
                List<Node> _order = reverse(bestOrder, edge);
                Graph _graph = fss.search(_order);
                double _score = scorer.score(_graph);

                if (_score < score) {
                    bestOrder = _order;
                    graph = _graph;
                    score = _score;
                    TetradLogger.getInstance().forceLogMessage("Decreases score: score = " + _score);
                    continue WHILE;
                }
            }

            break;
        }

        return graph;
    }

    private List<Node> reverse(List<Node> order, Edge edge) {
        Node n1 = edge.getNode1();
        Node n2 = edge.getNode2();

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
