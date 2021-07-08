package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    private int gspDepth = -1;

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

        double best = Double.NEGATIVE_INFINITY;
        List<Node> bestP = null;

        for (int r = 0; r < numStarts; r++) {
            List<Node> shuffled = new ArrayList<>(variables);
            shuffle(shuffled);
            scorer.score(shuffled);

            if (method == Method.PROMOTION) {
                bossSearchPromotion(scorer);
            } else if (method == Method.ALL_INDICES) {
                bossSearchAllIndices(scorer);
            } else if (method == Method.SP) {
                sp(scorer);
            } else if (method == Method.ESP) {
                esp(scorer);
            } else if (method == Method.GSP) {
                gsp(scorer);
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
            } else {
                break;
            }
        }
    }

    public void esp(TeyssierScorer scorer) {
        EspVisit visit = new EspVisit(scorer.getOrder(), scorer.score());
        EspVisit visit0;

        do {
            visit0 = visit;
            visit = espVisit(scorer, visit);
        } while (visit.getScore() != visit0.getScore());

        scorer.score(visit.getOrder());
    }

    public EspVisit espVisit(TeyssierScorer scorer, EspVisit G) {
        EspVisit Gtau = G;

        for (int i = 0; i < G.getOrder().size() - 1; i++) {
            Node v = G.getOrder().get(i);
            Node w = G.getOrder().get(i + 1);

            scorer.score(G.getOrder());

            scorer.swap(v, w);

            if (scorer.score() < Gtau.getScore()) {
                Gtau = espVisit(scorer, new EspVisit(scorer.getOrder(), scorer.score()));
            }
        }

        return Gtau;
    }

    public void gsp(TeyssierScorer scorer) {
        GspVisit visit0;
        GspVisit visit1 = new GspVisit(scorer);

        do {
            visit0 = visit1;
            visit1 = gspVisit(scorer, visit0, this.gspDepth == -1 ? Integer.MAX_VALUE : this.gspDepth, 0);
        } while (visit1 != visit0);

        scorer.score(visit1.getOrder());
    }

    public GspVisit gspVisit(TeyssierScorer scorer, GspVisit visit0, int maxDepth, int depth) {
        for (Edge edge : visit0.getGraph().getEdges()) {
            Node v = Edges.getDirectedEdgeTail(edge);
            Node w = Edges.getDirectedEdgeHead(edge);

            List<Node> parentsv = visit0.getGraph().getParents(v);
            List<Node> parentsw = visit0.getGraph().getParents(w);

            parentsw.remove(v);

            if (parentsv.equals(parentsw)) {
                scorer.swap(v, w);
                GspVisit visit2 = new GspVisit(scorer);

                if (depth < maxDepth && visit2.getNumEdges() <= visit0.getNumEdges()) {
                    if (visit2.getNumEdges() == visit0.getNumEdges()) depth++;
                    else depth = 0;

                    GspVisit visit1 = gspVisit(scorer, visit2, maxDepth, depth);

                    if (visit1.getNumEdges() < visit0.getNumEdges()) {
                        scorer.swap(v, w);
                        return visit1;
                    }
                }

                scorer.swap(v, w);
            }
        }

        return visit0;
    }

    public boolean isVerbose() {
        return verbose;
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
            for (Node z : scorer.getMb(p)) {
                G1.addDirectedEdge(z, order.get(p));
            }
        }

        return G1;
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

    public void setGspDepth(int gspDepth) {
        this.gspDepth = gspDepth;

    }

    public enum Method {PROMOTION, ALL_INDICES, SP, ESP, GSP}

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

    private static class GspVisit {
        private final Graph graph;
        private final List<Node> order;

        public GspVisit(TeyssierScorer scorer) {
            this.graph = getGraph(scorer);
            this.order = scorer.getOrder();
        }

        public Graph getGraph() {
            return graph;
        }

        public double getNumEdges() {
            return this.graph.getNumEdges();
        }

        public List<Node> getOrder() {
            return order;
        }

        private Graph getGraph(TeyssierScorer scorer) {
            List<Node> order = scorer.getOrder();
            Graph G1 = new EdgeListGraph(order);

            for (int p = 0; p < order.size(); p++) {
                for (Node z : scorer.getMb(p)) {
                    G1.addDirectedEdge(z, order.get(p));
                }
            }

            return G1;
        }
    }
}
