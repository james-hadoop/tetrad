package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.Scorer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
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
        ForwardScoreSearch fss = new ForwardScoreSearch(scorer);

        List<Node> b0 = new ArrayList<>(variables);
        Graph g0 = fss.search(b0);
        double s0 = scorer.score(g0);

        scoreOriginalOrder = s0;


        boolean changed = true;

        int[] pos = new int[variables.size()];
        List<Node> b2 = new ArrayList<>(variables);

        do {
            changed = false;

            for (Node n : variables) {
                double s5 = Double.POSITIVE_INFINITY;

                for (int j = 0; j < variables.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    Graph g1 = fss.search(b1);
                    double s1 = fss.score();

                    if (s1 < s5) {
                        pos[variables.indexOf(n)] = j;
                        s5 = s1;

//                        s0 = s1;
//                        b0 = b1;
//                        g0 = g1;
                        changed = true;
                    }
                }
            }

            System.out.println("POS " + Arrays.toString(pos));
            b2 = new ArrayList<>(b0);
            b0.sort(Comparator.comparingInt(o -> pos[variables.indexOf(o)]));
        } while (!b2.equals(b0));


//        for (int ii1 = 0; ii1 < variables.size(); ii1++) {
//            Node n1 = variables.get(ii1);
//
//            for (int i = 0; i < b0.indexOf(n1); i++) {
//                List<Node> b1 = new ArrayList<>(b0);
//                b1.remove(n1);
//                b1.add(i, n1);
//
//                for (int ii2 = ii1; ii2 < variables.size(); ii2++) {
//                    Node n2 = variables.get(ii2);
//
//                    for (int j = 0; j < b1.indexOf(n2); j++) {
//                        List<Node> b2 = new ArrayList<>(b1);
//                        b2.remove(n2);
//                        b2.add(j, n2);
//
//                        Graph g2 = fss.search(b2);
//                        double s2 = fss.score();
//
//                        if (s2 < s0) {
//                            s0 = s2;
//                            b0 = b2;
//                            g0 = g2;
////                            changed = true;
//                        }
//                    }
//                }
//            }
//        }

        System.out.println("Best order = " + b0);

        scoreLearnedOrder = s0;

        return fss.search(b0);

//        return g0;
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }
}
