package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;

import static java.util.Collections.shuffle;


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
    private int numStarts = 0;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public Graph search(List<Node> initialOrder) {
        shuffle(initialOrder);


        long start = System.currentTimeMillis();
        System.out.println("Original order = " + initialOrder);

        TeyssierScorer scorer = new TeyssierScorer(score);
        scorer.setCachingScores(cachingScores);

        double best = scorer.score(initialOrder);
        List<Node> bestP = scorer.getOrder();

        for (int r = 0; r < numStarts; r++) {
            bossSearchPromotion(scorer);

            if (scorer.score() < best) {
                best = scorer.score();
                bestP = scorer.getOrder();
            }

            List<Node> shuffled = new ArrayList<>(scorer.getOrder());
            shuffle(shuffled);
            scorer.score(shuffled);

//            sparsestPermutation(scorer);


        }

        scorer.score(bestP);

        long stop = System.currentTimeMillis();

        System.out.println("Final " + scorer.getOrder());

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return SearchGraphUtils.patternForDag(getGraph(scorer));
    }

    // Using Teyssier and Kohler's neighbor swaps.
    public void bossSearchPromotion(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = Double.POSITIVE_INFINITY;

        while (true) {
            for (Node node : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark(1);

                while (true) {
                    boolean moved = scorer.moveLeft(node);

                    if (!moved) break;

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                }

                scorer.restoreBookmark(1);

//                if (scorer.weaklyBetterThanBookmark(2)) {
//                    System.out.println("WB node = " + node + " order = " + scorer.getOrder());
//                }
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }
    }


    private void bossSearchOrig(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = scorer.score();

        while (true) {
            for (Node node : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark(1);

                while (true) {
                    boolean moved = scorer.moveLeft(node);

                    if (!moved) break;

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                }

                scorer.restoreBookmark(1);
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }
    }

    public void bossSearchAllIndices(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = scorer.score();

        while (true) {
            for (Node node : scorer.getOrder()) {
                scorer.moveTo(node, 0);
                double bestScore = scorer.score();
                scorer.bookmark(1);

                while (true) {
                    boolean moved = scorer.moveRight(node);

                    if (!moved) break;

                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                }

                scorer.restoreBookmark(1);
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }
    }

    private void bossSearchPromotionPairs(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = scorer.score();

        while (true) {
            double bestScore = scorer.score();

            for (Node v : scorer.getOrder()) {
                scorer.bookmark(1);

                while (true) {
                    boolean moved = scorer.moveLeft(v);

                    if (!moved) break;

                    int indexv = scorer.indexOf(v);

                    for (Node w : scorer.getOrder()) {
//                        int indexw = scorer.indexOf(w);
//                        if (indexw <= indexv) continue;

//                        if (visited.contains(w)) continue;
//                        if (v == w) continue;

                        scorer.bookmark(2);

                        while (true) {
                            boolean movedw = scorer.moveLeft(w);

                            if (!movedw) break;

                            if (scorer.score() < bestScore) {
//                                System.out.println("v = " + v + " w = " + w);
                                bestScore = scorer.score();
                                scorer.bookmark(2);
                                break;
                            }
                        }

                        scorer.restoreBookmark(2);
                    }

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                }

                scorer.restoreBookmark(1);
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
//                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }
    }

    private void bossSearchPromotionTriples(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = scorer.score();

        while (true) {
            double bestScore = scorer.score();

            for (Node v : scorer.getOrder()) {
                scorer.bookmark(1);

                while (true) {
                    boolean moved = scorer.moveLeft(v);
                    if (!moved) break;

                    for (Node w : scorer.getOrder()) {
                        scorer.bookmark(2);

                        while (true) {
                            boolean movedw = scorer.moveLeft(w);
                            if (!movedw) break;

                            for (Node r : scorer.getOrder()) {
                                scorer.bookmark(3);

                                while (true) {
                                    boolean movedr = scorer.moveLeft(r);
                                    if (!movedr) break;

                                    if (scorer.score() < bestScore) {
                                        bestScore = scorer.score();
                                        scorer.bookmark(3);
                                    }
                                }

                                scorer.restoreBookmark(3);
                            }

                            if (scorer.score() < bestScore) {
                                bestScore = scorer.score();
                                scorer.bookmark(2);
                            }
                        }

                        scorer.restoreBookmark(2);
                    }

                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                }

                scorer.restoreBookmark(1);
            }

            if (scorer.score() < overall) {
                overall = scorer.score();
                System.out.println("Updated order = " + scorer.getOrder());
            } else {
                break;
            }
        }
    }

    private void printReversible(TeyssierScorer scorer) {
        Graph dag = getGraph(scorer);

        for (Edge edge : dag.getEdges()) {
            if (!edge.isDirected()) continue;

            Node m = edge.getNode1();
            Node n = edge.getNode2();

            double score1 = scorer.score();
            scorer.bookmark(1);
            scorer.swap(m, n);
            double score2 = scorer.score();
            scorer.restoreBookmark(1);

            if (score2 - score1 <= 0) {
                System.out.println("Swapping " + n + " and " + m + " does not increase the score; diff = "
                        + (score2 - score1));
            }
        }
    }

    private void sparsestPermutation(TeyssierScorer scorer) {
        double minScore = Double.POSITIVE_INFINITY;
        List<Node> minP = null;

        List<Node> variables = scorer.getOrder();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);

            double score = scorer.score(p);
            if (score < minScore) {
                minScore = score;
                minP = p;
            }
        }

        scorer.score(minP);
    }


    @NotNull
    private Graph getGraph(TeyssierScorer scorer) {
        List<Node> order = scorer.getOrder();
        Graph G1 = new EdgeListGraph(order);

        for (int p = 0; p < order.size(); p++) {
//            System.out.println(order.get(p) + " mb = " + scorer.getMb(p));
            for (Node z : scorer.getMb(p)) {
                G1.addDirectedEdge(z, order.get(p));
            }
        }

        return SearchGraphUtils.patternForDag(G1);
    }


    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public void setNumStarts(int numStarts) {
        this.numStarts = numStarts;
    }
}
