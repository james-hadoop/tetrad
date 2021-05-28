package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.Scorer;

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
    private final Scorer scorer;

    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;

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
     * @param variables The variables for the graph.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> variables) {
        FastForwardSearch ffs = new FastForwardSearch(scorer);

        List<Node> b0 = new ArrayList<>(variables);
        Graph g0 = ffs.search(b0);
        double s0 = scorer.score(g0);

        scoreOriginalOrder = s0;

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node n : variables) {
                for (int j = 0; j < b0.indexOf(n); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    Graph g1 = ffs.search(b1);
                    double s1 = ffs.score();

                    if (s1 < s0 - 0.0001) {
                        s0 = s1;
                        b0 = b1;
                        g0 = g1;
                        changed = true;
                    }
                }
            }
        }

        System.out.println("Best order = " + b0);

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
