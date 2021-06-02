package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
        long start = System.currentTimeMillis();

        List<Node> b0 = new ArrayList<>(variables);
        double s0 = fastForward.score(b0);

        scoreOriginalOrder = s0;

        System.out.println("Original order = " + b0);

        boolean changed = true;

        int count = 0;

        while (changed) {
            changed = false;

            for (Node n : variables) {
                for (int j = 0; j < variables.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    count++;

                    double s1 = fastForward.score(b1);

                    if (s1 < s0) {
                        s0 = s1;
                        b0 = b1;
                        changed = true;
                    }
                }
            }

            System.out.println("Updated order = " + b0 + " count = " + count);
        }

        scoreLearnedOrder = s0;

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return fastForward.search(b0);
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }
}
