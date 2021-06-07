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
    private int algorithm = 1;

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

        if (algorithm == 1) {
            Map<Node, Set<Node>> associates = getAssociates(variables);

            List<Node> br = new ArrayList<>(b0);
            double sr = s0;

            for (int R = 0; R < 10; R++) {
                List<Node> b3 = new ArrayList<>(variables);
                shuffle(b3);
                double s3 = forwardScore.score(b3);

                if (s3 < sr) {
                    sr = s3;
                    br = b3;
                }

                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(b0);
                    double s2 = sr;

                    for (int i = 0; i < br.size() - 1; i++) {
                        if (!associates.get(b0.get(i)).contains(b0.get(i + 1))) continue;
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

                if (sr < s0) {
                    s0 = sr;
                    b0 = br;

                    System.out.println("Updated order = " + b0 + " count = " + count);
                }
            }
        } else if (algorithm == 2) {
            List<Node> b1 = gspRecursive(b0, s0);
            double s1 = forwardScore.score(b0);

            if (s1 < s0) {
                s0 = s1;
                b0 = b1;

                System.out.println("Updated order = " + b0 + " count = " + count);
            }
        } else {
            throw new IllegalArgumentException("Expecting algorithm 1 (theirs) or 2 (mine), 3 (esp)");
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
//
//                for (Node r : nodes) {
//                    if (forwardScore.isAssociated(w, r)) {
//                        nodes.add(r);
//                    }
//
//                    for (Node s : nodes) {
//                        if (forwardScore.isAssociated(w, s)) {
//                            nodes.add(s);
//                        }
//                    }
//                }
            }

            associates.put(v, nodes);
        }

        return associates;
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

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }

    public void setAlgorithm(int anInt) {
        this.algorithm = anInt;
    }
}
