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

        while (changed) {
            changed = false;

//            for (int i = 0; i < variables.size(); i++) {
//                for (int j = i + 1; j < variables.size(); j++) {
//                    List<Node> b1 = new ArrayList<>(b0);
//
//                    Node n1 = b0.get(i);
//                    Node n2 = b0.get(j);
//
//                    b1.remove(i);
//                    b1.add(i, n2);
//
//                    b1.remove(j);
//                    b1.add(j, n1);
//
//                    Graph g1 = fss.search(b1);
//                    double s1 = fss.score();
//
//                    if (s1 < s0) {
//                        s0 = s1;
//                        b0 = b1;
//                        g0 = g1;
//                        changed = true;
//                    }
//
//
//                }
//            }

            for (Node n : variables) {
                for (int j = 0; j < b0.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    Graph g1 = fss.search(b1);
                    double s1 = fss.score();

                    if (s1 < s0) {
                        s0 = s1;
                        b0 = b1;
                        g0 = g1;
                        changed = true;
                    }
                }
            }

            System.out.println("Interim order: " + b0);
        }

        System.out.println("Best order = " + b0);

        scoreLearnedOrder = s0;

        return g0;
    }

    /**
     * Does the search.
     *
     * @param variables The variables for the graph.
     * @return The estimated DAG.
     */
    public Graph search2(List<Node> variables) {
        ForwardScoreSearch fss = new ForwardScoreSearch(scorer);

        List<Node> b0 = new ArrayList<>(variables);
        Graph g0 = fss.search(b0);
        double s0 = scorer.score(g0);

        scoreOriginalOrder = s0;

        int[] pos = new int[variables.size()];
        List<Node> b2 = new ArrayList<>(variables);

        do {
            for (Node n : variables) {
                double s5 = Double.POSITIVE_INFINITY;

                for (int j = 0; j < variables.size(); j++) {
                    List<Node> b1 = new ArrayList<>(b0);
                    b1.remove(n);
                    b1.add(j, n);

                    fss.search(b1);
                    double s1 = fss.score();

                    if (s1 < s5) {
                        pos[variables.indexOf(n)] = j;
                        s5 = s1;
                    }
                }
            }

            System.out.println("POS " + Arrays.toString(pos));
            b2 = new ArrayList<>(b0);
            b0.sort(Comparator.comparingInt(o -> pos[variables.indexOf(o)]));
        } while (!b2.equals(b0));

        System.out.println("Best order = " + b0);

        scoreLearnedOrder = s0;

        Graph g2 = fss.search(b0);

        scoreLearnedOrder = fss.score();

        return g2;
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }

}
