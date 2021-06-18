package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

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
    private final Score score;
    private boolean cachingScores = true;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public Graph search(List<Node> initialOrder) {
        List<Node> order = getBestOrdering(initialOrder);

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

    // Using Teyssier and Kohler's neighbor swaps.
    public List<Node> getBestOrdering(List<Node> initialOrder) {
        long start = System.currentTimeMillis();
        System.out.println("Original order = " + initialOrder);

        TeyssierScorer scorer = new TeyssierScorer(score);
        scorer.setCachingScores(cachingScores);

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = scorer.score(initialOrder);

        while (true) {
            for (Node node : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark();

                while (true) {
                    boolean moved = scorer.moveLeft(node);

                    if (!moved) break;

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark();
                    }
                }

                scorer.restoreBookmark();
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return scorer.getOrder();
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }
}
