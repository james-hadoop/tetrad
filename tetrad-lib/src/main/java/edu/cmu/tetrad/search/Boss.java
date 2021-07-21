package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.OrderedPair;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

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
public class Boss {
    private final List<Node> variables;
    private Score score;
    private IndependenceTest test;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private Method method = Method.BOSS_PROMOTION;
    private boolean verbose = false;
    private boolean doFinalOrientation = false;

    public Boss(Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
    }

    public Boss(IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(List<Node> order) {
        long start = System.currentTimeMillis();
        order = new ArrayList<>(order);

        TeyssierScorer scorer;

        if (score != null) {
            scorer = new TeyssierScorer(score);
        } else {
            scorer = new TeyssierScorer(test);
        }

        scorer.setCachingScores(cachingScores);

        double best = Double.POSITIVE_INFINITY;
        List<Node> bestPerm = null;

        for (int r = 0; r < numStarts; r++) {
            if (r > 0) scorer.shuffleVariables();
            scorer.score(order);
            List<Node> perm;

            if (method == Method.BOSS_PROMOTION) {
                perm = bossSearchPromotion(scorer);
            } else if (method == Method.BOSS_ALL_INDICES) {
                perm = bossSearchAllIndices(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else {
                throw new IllegalArgumentException("Unrecognized method: " + method);
            }

            if (scorer.score() < best) {
                best = scorer.score();
                bestPerm = perm;
            }
        }

        long stop = System.currentTimeMillis();

        if (verbose) {
            System.out.println("Final order = " + scorer.getOrder());
            System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
        }

        return bestPerm;
    }

    public List<Node> bossSearchPromotion(TeyssierScorer scorer) {

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.
        if (doFinalOrientation) {
            while (true) {
                Set<EqualsEntry> equals = promoteLoop(scorer);
                double score = scorer.score();
                finalOrientation(scorer, equals);
                if (scorer.score() == score) break;
            }

            return scorer.getOrder();
        } else {
            promoteLoop(scorer);
            return scorer.getOrder();
        }
    }

    public List<Node> bossSearchAllIndices(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        if (doFinalOrientation) {
            while (true) {
                Set<EqualsEntry> equals = allIndicesLoop(scorer);
                double score = scorer.score();
                finalOrientation(scorer, equals);
                if (scorer.score() == score) break;
            }

            return scorer.getOrder();
        } else {
            allIndicesLoop(scorer);
            return scorer.getOrder();
        }
    }

    private Set<EqualsEntry> promoteLoop(TeyssierScorer scorer) {
        boolean reduced;
        Set<EqualsEntry> equals = new HashSet<>();

        do {
            reduced = false;

            for (Node v : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark();

                while (scorer.promote(v)) {
                    if (scorer.score() < bestScore) {
                        reduced = true;
                        scorer.bookmark();
                        equals.clear();
                    } else {
                        equals.add(new EqualsEntry(v, scorer.indexOf(v), scorer.score()));
                    }
                }

                scorer.goToBookmark();
            }

            if (verbose) {
                System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Promotion)");
            }
        } while (reduced);

        return equals;
    }

    private Set<EqualsEntry> allIndicesLoop(TeyssierScorer scorer) {
        boolean reduced;
        Set<EqualsEntry> equals2 = new HashSet<>();

        do {
            reduced = false;
            List<Node> order = scorer.getOrder();
            scorer.bookmark();

            for (Node v : order) {
                double bestScore = scorer.score();
                scorer.moveToFirst(v);

                if (scorer.score() < bestScore) {
                    bestScore = scorer.score();
                    reduced = true;
                    equals2.clear();
                    scorer.bookmark();
                } else {
                    equals2.add(new EqualsEntry(v, scorer.indexOf(v), scorer.score()));
                }

                while (scorer.demote(v)) {
                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        reduced = true;
                        equals2.clear();
                        scorer.bookmark();
                    } else {
                        equals2.add(new EqualsEntry(v, scorer.indexOf(v), scorer.score()));
                    }
                }

                scorer.goToBookmark();
            }

            if (verbose) {
                System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (All indices)");
            }
        } while (reduced);

        return equals2;
    }

    private void finalOrientation(TeyssierScorer scorer, Set<EqualsEntry> equals) {
        scorer.bookmark();

        Graph graph = scorer.getGraph(false);
        double bestScore = scorer.score();

        Set<OrderedPair<Node>> pairs = new HashSet<>();

        for (EqualsEntry e : equals) {
            Node v = e.getV();

            for (Node r : graph.getAdjacentNodes(v)) {
                for (Node r2 : graph.getAdjacentNodes(v)) {
                    if (graph.isDefCollider(r, v, r2)) {
                        pairs.add(new OrderedPair<>(v, r));
                        pairs.add(new OrderedPair<>(v, r2));
                    } else if (graph.isDirectedFromTo(v, r) && graph.isDirectedFromTo(v, r2)) {
                        pairs.add(new OrderedPair<>(v, r));
                        pairs.add(new OrderedPair<>(v, r2));
                    }
                }
            }
        }

        for (OrderedPair<Node> pair : pairs) {
            Node r = pair.getSecond();
            scorer.bookmark();
            scorer.moveTo(r, 0);

            while (scorer.demote(r)) {
                if (scorer.score() < bestScore) {
                    scorer.bookmark();
                    bestScore = scorer.score();
                    scorer.bookmark();

                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Final Orientation)");
                    }
                }
            }

            scorer.goToBookmark();
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Final Orientation)");
        }
    }

    public List<Node> esp(TeyssierScorer scorer) {
        EspVisit visit = new EspVisit(scorer.getOrder(), scorer.score());
        EspVisit visit0;

        do {
            visit0 = visit;
            visit = espVisit(scorer, visit);
        } while (visit.getScore() != visit0.getScore());

        return visit.getOrder();
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

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public List<Node> sp(TeyssierScorer scorer) {
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


        return minP;
    }

    @NotNull
    public Graph getGraph(List<Node> order, boolean pattern) {
        TeyssierScorer scorer;

        if (score != null) {
            scorer = new TeyssierScorer(score);
        } else {
            scorer = new TeyssierScorer(test);
        }

        scorer.score(order);
        return scorer.getGraph(pattern);
    }

    public void setCacheScores(boolean cachingScores) {
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

    public List<Node> getVariables() {
        return this.variables;
    }

    public void setDoFinalOrientation(boolean doFinalOrientation) {
        this.doFinalOrientation = doFinalOrientation;
    }

    public enum Method {BOSS_PROMOTION, BOSS_ALL_INDICES, SP}

    public static class EqualsEntry {
        private final Node v;
        private final int index;
        private final double score;

        public EqualsEntry(Node v, int index, double score) {
            this.v = v;
            this.index = index;
            this.score = score;
        }

        public Node getV() {
            return v;
        }

        public int getIndex() {
            return index;
        }

        public double getScore() {
            return score;
        }

        public String toString() {
            return "(" + v + "/" + index /*+ "," + score */ + ")";
        }

        public int hashCode() {
            return 17 * v.hashCode() + index;
        }

        public boolean equals(Object o) {
            if (o == null) return false;
            if (!(o instanceof EqualsEntry)) return false;
            EqualsEntry e = (EqualsEntry) o;

            return v == e.v && index == e.index && score == e.score;
        }
    }

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
}
