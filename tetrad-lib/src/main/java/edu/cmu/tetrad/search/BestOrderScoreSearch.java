package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

import java.util.*;

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
    private final ForwardScore forwardScore;
    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;
    private Method method = Method.NONRECURSIVE;
    private int numRestarts = 1;

    /**
     * Constructs a GSS search
     *
     * @param forward the fastForward algorithm used.
     */
    public BestOrderScoreSearch(ForwardScore forward) {
        this.forwardScore = forward;
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
        double s0 = forwardScore.score(b0);

        scoreOriginalOrder = s0;

        System.out.println("Original order = " + b0);

        boolean changed = true;

        int count = 0;

        List<Node> br = new ArrayList<>(b0);

        Map<Node, Set<Node>> associates = getAssociates(variables);

        br = new ArrayList<>(br);
        double sr = forwardScore.score(br);

        for (int r = 0; r < numRestarts; r++) {
            if (method == Method.NONRECURSIVE) {

                // Until you can't do it anymore...
                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(br);
                    double s2 = sr;

                    // ...pick a variable v...
                    for (Node v : variables) {

                        // ...and move v to an optimal position.
                        for (Node w : associates.get(v)) {
                            int j = b2.indexOf(w);
                            List<Node> b1 = new ArrayList<>(b2);
                            b1.remove(v);
                            b1.add(j, v);

                            count++;

                            double s1 = forwardScore.score(b1);

                            if (s1 < s2) {
                                s2 = s1;
                                b2 = b1;
                            }
                        }
                    }

                    // Output the best order you find.
                    if (s2 < sr) {
                        br = b2;
                        sr = s2;
                        changed = true;
                    }
                }
            } else if (method == Method.RECURSIVE) {
                br = bossRecursive(br, sr, associates);
                sr = forwardScore.score(br);
            } else {
                throw new IllegalStateException("Unrecognized method: " + method);
            }

            if (sr < s0) {
                s0 = sr;
                b0 = new ArrayList<>(br);

                System.out.println("Updated order = " + b0 + " count = " + count);
            }

            shuffle(br);
            sr = forwardScore.score(br);
        }

        scoreLearnedOrder = s0;

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return forwardScore.search(b0);
    }

    private List<Node> bossRecursive(List<Node> b0, double s0, Map<Node, Set<Node>> associates) {

        List<Node> b2 = new ArrayList<>(b0);
        double s2 = s0;

        boolean found = false;

        for (Node v : b0) {
            for (Node w : associates.get(v)) {
                int j = b2.indexOf(w);
                List<Node> b1 = new ArrayList<>(b2);
                b1.remove(v);
                b1.add(j, v);

                double s1 = forwardScore.score(b1);

                if (s1 < s2) {
                    s2 = s1;
                    b2 = b1;
                    found = true;
                }
            }
        }

        if (found) {
            return bossRecursive(b2, s2, associates);
        } else {
            return b2;
        }

    }

    @NotNull
    private Map<Node, Set<Node>> getAssociates(List<Node> variables) {
        Map<Node, Set<Node>> associates = new HashMap<>();

        for (Node v : variables) {
            Set<Node> nodes = new HashSet<>();

            for (Node w : variables) {
                if (forwardScore.isAssociated(w, v)) {
                    nodes.add(w);
                }
            }

            associates.put(v,  nodes);
        }

        return associates;
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }

    public void setMethod(Method method) {
        this.method = method;
    }

    public void setNumRestarts(int numRestarts) {
        this.numRestarts = numRestarts;
    }

    public enum Method{NONRECURSIVE, RECURSIVE}
}
