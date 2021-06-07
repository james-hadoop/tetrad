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
public class GreedySparsestPermutation {

    // The score used (default FML BIC). Lower is better.
    private final ForwardScore forwardScore;
    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;
    private Method method = Method.NORECURSIVE;
    private int numRestarts = 1;

    /**
     * Constructs a GSS search
     *
     * @param forward the fastForward algorithm used.
     */
    public GreedySparsestPermutation(ForwardScore forward) {
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

        Map<Node, Set<Node>> associates = getAssociates(variables);

        List<Node> br = new ArrayList<>(b0);
        double sr = s0;

        for (int r = 0; r < numRestarts; r++) {
            if (method == Method.NORECURSIVE) {
                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(br);
                    double s2 = sr;

                    for (int i = 0; i < br.size() - 1; i++) {
                        if (!associates.get(br.get(i)).contains(br.get(i + 1))) continue;
                        List<Node> b1 = swap(br, i);

                        double s1 = forwardScore.score(b1);

                        if (s1 < s2) {
                            s2 = s1;
                            b2 = b1;
                        }
                    }

                    if (s2 < sr) {
                        sr = s2;
                        br = b2;
                        changed = true;
                    }
                }
            } else if (method == Method.RECURSIVE) {
                br = gspRecursive(br, sr);
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

    private List<Node> gspRecursive(List<Node> b0, double s0) {
        List<Node> b2 = new ArrayList<>(b0);
        double s2 = s0;
        boolean found = false;

        for (int i = 0; i < b2.size() - 1; i++) {
            List<Node> b1 = swap(b2, i);

            double s1 = forwardScore.score(b1);

            if (s1 < s2) {
                b2 = b1;
                s2 = s1;
                found = true;
            }
        }

        if (found) {
            return gspRecursive(b2, s2);
        } else {
            return b2;
        }
    }

    @NotNull
    private List<Node> swap(List<Node> b0, int i) {
        List<Node> b1 = new ArrayList<>(b0);

        Node n1 = b0.get(i);
        Node n2 = b0.get(i + 1);

        b1.remove(n1);
        b1.add(i + 1, n1);

        b1.remove(n2);
        b1.add(i, n2);
        return b1;
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

            associates.put(v, nodes);
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

    public enum Method{NORECURSIVE, RECURSIVE}
}
