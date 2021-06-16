package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.*;


/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {
    private final Score score;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public Graph search(List<Node> initialOrder, int method) {
        List<Node> order;

        switch (method) {
            case 1:
                order = getBestOrdering1(initialOrder);
                break;
            case 2:
                order = getBestOrdering2(initialOrder);
                break;
            case 3:
                order = getBestOrdering3(initialOrder);
                break;
            default:
                throw new IllegalArgumentException("Expecting 1 or 2.");
        }

        K3 k3 = new K3(score);

        List<Set<Node>> pis = k3.scoreResult(order, true).getPis();

        Graph G1 = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            for (Node z : pis.get(i)) {
                G1.addDirectedEdge(z, order.get(i));
            }
        }

        return SearchGraphUtils.patternForDag(G1);
    }

    public List<Node> getBestOrdering1(List<Node> variables) {
        TeyssierScorer scorer = new TeyssierScorer(score);

        long start = System.currentTimeMillis();

        List<Node> b0 = new ArrayList<>(variables);
        double s0 = scorer.score(b0);

        System.out.println("Original order = " + b0);

        boolean changed = true;
        List<Node> br = new ArrayList<>(b0);
        Map<Node, Set<Node>> associates = scorer.getAssociates(variables);
        double sr = scorer.score(br);

        // Take each variable in turn and move it to the best-scoring position. Repeat until
        // no more variables can be repositioned this way.
        while (true) {
            while (changed) {
                changed = false;

                List<Node> b2 = new ArrayList<>(br);
                double s2 = sr;

                for (Node v : variables) {
                    for (Node w : associates.get(v)) {
                        int j = b2.indexOf(w);
                        List<Node> b1 = new ArrayList<>(b2);
                        b1.remove(v);
                        b1.add(j, v);

                        double s1 = scorer.score(b1);

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

            if (sr < s0) {
                s0 = sr;
                b0 = new ArrayList<>(br);
                System.out.println("Updated order = " + b0);
            } else {
                break;
            }
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return b0;
    }

    // Using Teyssier and Kohler's neighbor swaps.
    public List<Node> getBestOrdering2(List<Node> initialOrder) {
        long start = System.currentTimeMillis();
        System.out.println("Original order = " + initialOrder);

        TeyssierScorer scorer = new TeyssierScorer(score);
        double best = scorer.score(initialOrder);

        // Take each variable in turn and move it to the best-scoring position. Repeat until
        // no more variables can be repositioned this way.
        while (true) {
            for (Node node : initialOrder) {
                scorer.moveTo(node, 0);
                double bestScore = scorer.score();
                int bestIndex = scorer.indexOf(node);
                boolean moved;

                do {
                    moved = scorer.moveRight(node);

                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        bestIndex = scorer.indexOf(node);
                    }
                } while (moved);

                scorer.moveTo(node, bestIndex);
            }

            if (scorer.score() < best) {
                best = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return scorer.getOrder();
    }

    public List<Node> getBestOrdering3(List<Node> initialOrder) {
        TeyssierScorer scorer = new TeyssierScorer(score);
        double s0 = scorer.score(initialOrder);
        greedySolun(scorer, s0);
        return scorer.getOrder();
    }

    private void greedySolun(TeyssierScorer scorer, double s0) {
        for (Node v : scorer.getOrder()) {
            scorer.moveTo(v, 0);

            for (int i = 0; i < scorer.getOrder().size(); i++) {
                scorer.moveTo(v, i);

                if (scorer.score() < s0) {
                    greedySolun(scorer, scorer.score());
                }
            }
        }
    }
}
