package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.List;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {

    // The score used (default FML BIC). Lower is better.
    private final FastForward fastForward;

    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;

    /**
     * Constructs a GSS search
     *
     * @param forward the fastForward algorithm used.
     */
    public BestOrderScoreSearch(FastForward forward) {
        this.fastForward = forward;
    }

    /**
     * Does the search.
     *
     * @param variables The variables for the graph.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> variables) {
        List<Node> b0 = new ArrayList<>(variables);
        Graph g0 = fastForward.search(b0);
        double s0 = fastForward.score();

        scoreOriginalOrder = s0;

        System.out.println("Original order = " + b0);

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node n : variables) {
                for (int j = 0; j < variables.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    Graph g1 = fastForward.search(b1);
                    double s1 = fastForward.score();

                    if (s1 < s0) {
                        s0 = s1;
                        b0 = b1;
                        g0 = g1;
                        changed = true;
                    }
                }
            }

            System.out.println("Updated order = " + b0);
        }

        scoreLearnedOrder = s0;

        return g0;
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }
}
