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
    private final FastForward fastForward;
    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;
    private int algorithm = 1;

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

        if (algorithm == 1) {
            Map<Node, Set<Node>> associates = getAssociates(variables);

            LinkedList<List<Node>> tabu = new LinkedList<>();

            int c1 = 0;
            boolean changed2 = true;

            while (true) {
                changed2 = false;

                List<Node> b3 = new ArrayList<>(variables);
                shuffle(b3);
                double s3 = fastForward.score(b3);

                if (s3 < s0) {
                    s0 = s3;
                    b0 = b3;
                }

                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(b0);
                    double s2 = s0;

                    for (int i = 0; i < b0.size() - 1; i++) {
//                        if (!associates.get(b0.get(i)).contains(b0.get(i + 1))) continue;
                        List<Node> b1 = swap(b0, i);

                        double s1 = fastForward.score(b1);

                        if (s1 < s2) {
                            s2 = s1;
                            b2 = b1;
                        }
                    }

                    if (s2 < s0) {
                        s0 = s2;
                        b0 = b2;
                        changed = true;
                    }
                }

                tabu.add(new ArrayList<>(b0));

                if (tabu.size() > 20) break;
            }
        } else if (algorithm == 2) {

            Map<Node, Set<Node>> associates = getAssociates(variables);

//            for (Node v :  variables) System.out.println("Node " + v + " # associates = " + associates.get(v).size());

            // Until you can't do it anymore...
            while (changed) {
                changed = false;

                List<Node> b2 = new ArrayList<>(b0);
                double s2 = s0;

                // ...pick a variable v...
                for (Node v : variables) {

                    // ...and move v to an optimal position.
                    for (Node w : associates.get(v)) {
                        int j = b2.indexOf(w);
                        List<Node> b1 = new ArrayList<>(b2);
                        b1.remove(v);
                        b1.add(j, v);

                        count++;

                        double s1 = fastForward.score(b1);

                        if (s1 < s2) {
                            s2 = s1;
                            b2 = b1;
                        }
                    }
                }

                // Output the best order you find.
                if (s2 < s0) {
                    b0 = b2;
                    s0 = s2;
                    changed = true;

                    System.out.println("Updated order = " + b0 + " count = " + count);
                }
            }
        } else if (algorithm == 3) {
            b0 = esp(b0, s0);
            s0 = fastForward.score(b0);
        } else if (algorithm == 4) {
            b0 = esp2(b0, s0, getAssociates(b0));
            s0 = fastForward.score(b0);


//            int r = 5;
//            Map<Node, Set<Node>> associates = getAssociates(b0);
//
//            for (int R = 0; R < r; R++) {
//                List<Node> b2 = new ArrayList<>(b0);
//                shuffle(b2);
//                List<Node> b1 = esp2(b2, s0, associates);
//                double s1 = fastForward.score(b1);
//
//                if (s1 < s0) {
//                    b0 = b1;
//                    s0 = s1;
//                }
//            }
        } else {
            throw new IllegalArgumentException("Expecting algorithm 1 (theirs) or 2 (mine), 3 (esp)");
        }

        scoreLearnedOrder = s0;

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return fastForward.search(b0);
    }

    private List<Node> esp2(List<Node> b0, double s0, Map<Node, Set<Node>> associates) {

        List<Node> b2 = new ArrayList<>(b0);
        double s2 = s0;

        boolean found = false;

        for (Node v : b0) {
            for (Node w : associates.get(v)) {
                int j = b2.indexOf(w);
                List<Node> b1 = new ArrayList<>(b2);
                b1.remove(v);
                b1.add(j, v);

                double s1 = fastForward.score(b1);

                if (s1 < s2) {
                    s2 = s1;
                    b2 = b1;
                    found = true;
                }
            }
        }

        if (found) {
            return esp(b2, s2);
        } else {
            return b2;
        }

    }

    private List<Node> esp(List<Node> b0, double s0) {
        List<Node> b2 = new ArrayList<>(b0);
        double s2 = s0;
        boolean found = false;

        for (int i = 0; i < b2.size() - 1; i++) {
            List<Node> b1 = swap(b2, i);

            double s1 = fastForward.score(b1);

            if (s1 < s2) {
                b2 = b1;
                s2 = s1;
                found = true;
            }
        }

        if (found) {
            return esp(b2, s2);
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
                if (fastForward.isAssociated(w, v)) {
                    nodes.add(w);
                }
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
