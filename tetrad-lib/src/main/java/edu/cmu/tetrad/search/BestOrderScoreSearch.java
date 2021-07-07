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
    private Score score;
    private IndependenceTest test;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private Method method = Method.PROMOTION;
    private boolean verbose = false;
    private boolean returnCpdag = false;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public BestOrderScoreSearch(IndependenceTest test) {
        this.test = test;
    }

    public Graph search(List<Node> variables) {
        long start = System.currentTimeMillis();

        TeyssierScorer scorer;

        if (score != null) {
            scorer = new TeyssierScorer(score);
        } else {
            scorer = new TeyssierScorer(test);
        }

        scorer.setCachingScores(cachingScores);

//        double best = scorer.score(variables);
//        List<Node> bestP =scorer.getOrder();

        double best = Double.NEGATIVE_INFINITY;
        List<Node> bestP = null;

        for (int r = 0; r < numStarts; r++) {
            List<Node> shuffled = new ArrayList<>(variables);
            shuffle(shuffled);
            scorer.score(shuffled);

            if (method == Method.PROMOTION) {
                bossSearchPromotion(scorer);
//            } else if (method == Method.PROMOTION_PAIRS) {
//                bossSearchPromotionPairs(scorer);
            } else if (method == Method.ALL_INDICES) {
                bossSearchAllIndices(scorer);
            } else if (method == Method.SP) {
                sp(scorer);
            } else if (method == Method.ESP) {
                esp(scorer);
            } else if (method == Method.TSP) {
                tsp(scorer);
            }

            if (bestP == null || scorer.score() < best) {
                best = scorer.score();
                bestP = scorer.getOrder();
            }
        }

        scorer.score(bestP);

        long stop = System.currentTimeMillis();

        if (verbose) {
            System.out.println("Final " + scorer.getOrder());
            System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
        }

        if (returnCpdag) {
            return SearchGraphUtils.patternForDag(getGraph(scorer));
        } else {
            return getGraph(scorer);
        }
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
            for (Node v : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark(1);

                do {
                    if (scorer.score() <= bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                } while (scorer.moveLeft(v));

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

    public void bossSearchAllIndices(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        double overall = Double.POSITIVE_INFINITY;

        while (true) {
            for (Node v : scorer.getOrder()) {
                scorer.moveTo(v, 0);
                double bestScore = scorer.score();
                scorer.bookmark(1);

                do {
                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        scorer.bookmark(1);
                    }
                } while (scorer.moveRight(v));

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

//    private void bossSearchPromotionPairs(TeyssierScorer scorer) {
//
//        // Take each variable in turn and try moving it to each position to the left (i.e.,
//        // try promoting it in the causal order). Score each causal order by building
//        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
//        // Markov blanket algorithm, to find putative parents for each node.) Finally,
//        // place the node in whichever position yielded the highest score. Do this
//        // for each variable in turn. Once you're done, do it all again, until no more
//        // variables can be relocated.
//        double overall = Double.POSITIVE_INFINITY;
//
//        while (true) {
//            for (Node v : scorer.getOrder()) {
////                scorer.moveTo(v, 0);
//                double bestScore = scorer.score();
//                scorer.bookmark(1);
//
//                do {
//                    for (Node w : scorer.getOrder()) {
//                        if (v == w) continue;
////                        scorer.moveTo(w, 0);
//
//                        do {
//                            if (scorer.score() <= bestScore) {
//                                bestScore = scorer.score();
//                                scorer.bookmark(1);
//                            }
//                        } while (scorer.moveLeft(w));
//                    }
//
////                    scorer.restoreBookmark(1);
//
//                    if (scorer.score() >= bestScore) break;
//                } while (scorer.moveLeft(v));
//
//                scorer.restoreBookmark(1);
//            }
//
//            if (scorer.score() < overall) {
//                overall = scorer.score();
////                System.out.println("Updated order = " + scorer.getOrder());
//            } else {
//                break;
//            }
//        }
//    }

//    private void bossSearchPromotionTriples(TeyssierScorer scorer) {
//
//        // Take each variable in turn and try moving it to each position to the left (i.e.,
//        // try promoting it in the causal order). Score each causal order by building
//        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
//        // Markov blanket algorithm, to find putative parents for each node.) Finally,
//        // place the node in whichever position yielded the highest score. Do this
//        // for each variable in turn. Once you're done, do it all again, until no more
//        // variables can be relocated.
//        double bestScore = scorer.score();
//        scorer.bookmark(1);
//
//        for (Node v : scorer.getOrder()) {
//            double bestScoreV = bestScore;
//
//            do {
//                for (Node w : scorer.getOrder()) {
//                    do {
//                        for (Node r : scorer.getOrder()) {
//                            do {
//                                if (scorer.score() <= bestScore) {
//                                    bestScore = scorer.score();
//                                    scorer.bookmark(1);
//                                }
//                            } while (scorer.moveLeft(r));
//                        }
//
//                        if (bestScoreV >= bestScore) break;
//                    } while (scorer.moveLeft(w));
//                }
//
//                if (bestScoreV >= bestScore) break;
//            } while (scorer.moveLeft(v));
//        }
//
//        scorer.restoreBookmark(1);
//    }

    public void esp(TeyssierScorer scorer) {
        EspVisit visit = new EspVisit(scorer.getOrder(), scorer.score());
        EspVisit visit0;

        do {
            visit0 = visit;
            visit = espVisit(scorer, visit);
        } while (visit.getScore() != visit0.getScore());

        scorer.score(visit.getOrder());
    }

    public EspVisit espVisit(TeyssierScorer scorer, EspVisit visit) {
        for (Node v : visit.getOrder()) {
            scorer.score(visit.getOrder());
            if (scorer.moveRight(v)) {
                if (scorer.score() < visit.getScore()) {
                    return espVisit(scorer, new EspVisit(scorer.getOrder(), scorer.score()));
                }
            }
        }

        return visit;
    }

    public void tsp(TeyssierScorer scorer) {
        Graph graph = getGraph(scorer);
        double score = graph.getNumEdges();
        TspVisit visit = tspVisit(scorer, new TspVisit(graph, scorer.getOrder(), score));
        scorer.score(visit.getOrder());
    }

    public TspVisit tspVisit(TeyssierScorer scorer, TspVisit visit) {

        WHILE:
        while (true) {
            for (Edge edge : visit.getGraph().getEdges()) {
                Node v = edge.getNode1();
                Node w = edge.getNode2();

                if (visit.getGraph().getParents(v).equals(visit.getGraph().getParents(w))) {
                    scorer.swap(v, w);
                    List<Node> newOrdering = scorer.getOrder();
                    Graph newGraph = getGraph(scorer);

                    if (newGraph.getNumEdges() <= visit.getNumEdges()) {
                        visit = new TspVisit(newGraph, newOrdering, newGraph.getNumEdges());
                        continue WHILE;
                    }
                }
            }

            break;
        }

        return visit;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public boolean isReturnCpdag() {
        return returnCpdag;
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

    public void sp(TeyssierScorer scorer) {
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

        return G1;
//        return SearchGraphUtils.patternForDag(G1);
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public void setNumStarts(int numStarts) {
        this.numStarts = numStarts;
    }

    public Method getMethod() {
        return method;
    }

    public void setMethod(Method method) {
        this.method = method;
    }

    public enum Method {PROMOTION, /*PROMOTION_PAIRS,*/ ALL_INDICES, SP, ESP, TSP}

    private static class EspVisit {
        private final List<Node> order;
        private final double score;

        public EspVisit(List<Node> order, double score) {
            this.order = order;
            this.score = score;
        }

        public List<Node> getOrder() {
            return order;
        }

        public double getScore() {
            return score;
        }
    }

    private static class TspVisit {
        private final Graph graph;
        private final List<Node> order;
        private final double numEdges;

        public TspVisit(Graph graph, List<Node> order, double numEdges) {
            this.graph = graph;
            this.order = order;
            this.numEdges = numEdges;
        }

        public Graph getGraph() {
            return graph;
        }

        public double getNumEdges() {
            return numEdges;
        }

        public List<Node> getOrder() {
            return order;
        }
    }
}
