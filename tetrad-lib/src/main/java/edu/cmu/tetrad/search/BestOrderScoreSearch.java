package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

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
    private final Score score;
    private boolean cachingScores = true;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public Graph search(List<Node> initialOrder) {
        Graph graph = bossSearch(initialOrder);
        return SearchGraphUtils.patternForDag(graph);
    }

    @NotNull
    private Graph getGraph(TeyssierScorer scorer) {
        List<Node> order = scorer.getOrder();
        Graph G1 = new EdgeListGraph(order);

        for (int p = 0; p < order.size(); p++) {
            for (Node z : scorer.getMb(p)) {
                G1.addDirectedEdge(z, order.get(p));
            }
        }

        return SearchGraphUtils.patternForDag(G1);
    }

    // Using Teyssier and Kohler's neighbor swaps.
    private Graph bossSearch(List<Node> initialOrder) {
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

        return getGraph(scorer);
    }

    private Graph bossSearch2(List<Node> initialOrder) {
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
                scorer.moveTo(node, 0);
                double bestScore = scorer.score();
                scorer.bookmark();

                while (true) {
                    boolean moved = scorer.moveRight(node);

                    if (!moved) break;

                    if (scorer.score() < bestScore) {
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

        return getGraph(scorer);
    }

    private Graph bossSearch3(List<Node> initialOrder) {
        long start = System.currentTimeMillis();
        System.out.println("Original order = " + initialOrder);

        TeyssierScorer scorer = new TeyssierScorer(score);
        scorer.setCachingScores(cachingScores);
        List<Node> running = new ArrayList<>();

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        for (Node m : initialOrder) {
            running.add(m);
            double bestScore = scorer.score(running);

//            for (Node node : running) {
            scorer.bookmark();

            for (int k = running.size() - 1; k >= 0; k--) {
                Node n = running.get(k);

                while (true) {
                    boolean moved = scorer.moveLeft(n);

                    if (!moved) break;

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark();
                    }
                }
            }

            scorer.restoreBookmark();
//            }

//            if (scorer.score() < overall) {
//                overall = scorer.score();
//                System.out.println("Updated order = " + scorer.getOrder());
//            }
//            else {
//                break;
//            }

            running = scorer.getOrder();
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return getGraph(scorer);
    }

    private void printReversible(TeyssierScorer scorer) {
        Graph dag = getGraph(scorer);

        for (Edge edge : dag.getEdges()) {
            if (!edge.isDirected()) continue;

            Node m = edge.getNode1();
            Node n = edge.getNode2();

            double score1 = scorer.score();
            scorer.bookmark();
            scorer.swap(m, n);
            double score2 = scorer.score();
            scorer.restoreBookmark();

            if (score2 - score1 <= 0) {
                System.out.println("Swapping " + n + " and " + m + " does not increase the score; diff = "
                        + (score2 - score1));
            }
        }
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }
}
