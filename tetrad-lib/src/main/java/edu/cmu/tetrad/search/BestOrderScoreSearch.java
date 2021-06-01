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
        List<Node> b0 = new ArrayList<>(variables);
        Graph g0 = fastForward.search(b0);
        double s0 = fastForward.score();

        scoreOriginalOrder = s0;

        System.out.println("Original order = " + b0);

        Set<List<Node>> tabu = new HashSet<>();

        boolean changed = true;

        int count = 0;

//        while (++count < 100000) {
//            if (RandomUtil.getInstance().nextDouble() < 0.01) {
//                b0 = new ArrayList<>(variables);
//                shuffle(b0);
//            }
//
//            changed = false;
//
//            for (int i = 0; i < b0.size() - 1; i++) {
//                int j = i + 1;
////                for (int j = 0; j < b0.size(); j++) {
//                List<Node> b1 = new ArrayList<>(b0);
//
//                Node n1 = b1.get(i);
//                Node n2 = b1.get(j);
//
//                b1.remove(n1);
//                b1.add(j, n1);
//
//                b1.remove(n2);
//                b1.add(i, n2);
//
//                if (tabu.contains(b1)) continue;
//                tabu.add(b1);
//
//                Graph g1 = fastForward.search(b1);
//                double s1 = fastForward.score();
//
//                if (s1 < s0) {
//                    s0 = s1;
//                    g0 = g1;
//                    changed = true;
//                }
//            }
////            }
//        }
//
//        int count = 0;
//
//        boolean changed = true

//        List<Node> b2 = b0;
//        Graph g2 = g0;
//        double s2 = s0;
//
//        WHILE:
        while (changed) {
            changed = false;

            for (Node n : variables) {
                for (int j = 0; j < variables.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    count++;

                    Graph g1 = fastForward.search(b1);
                    double s1 = fastForward.score();

                    if (s1 < s0) {

                        if (tabu.contains(b1)) continue;
                        tabu.add(b1);

                        s0 = s1;
                        b0 = b1;
                        g0 = g1;
                        changed = true;
                    }
                }
            }

//            b2 = b0;
//            g2 = g0;
//            s2 = s0;

            System.out.println("Updated order = " + b0 + " count = " + count);
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
