package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

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
public class Gsp {
    private Score score;
    private IndependenceTest test;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private boolean verbose = false;
    private boolean returnCpdag = false;
    private int gspDepth = -1;

    public Gsp(Score score) {
        this.score = score;
    }

    public Gsp(IndependenceTest test) {
        this.test = test;
    }

    public Graph search(List<Node> order) {
        long start = System.currentTimeMillis();

        TeyssierScorer scorer;

        if (score != null) {
            scorer = new TeyssierScorer(score);
        } else {
            scorer = new TeyssierScorer(test);
        }

        scorer.setCachingScores(cachingScores);

        double best = Double.POSITIVE_INFINITY;
        List<Node> bestP = null;

        scorer.score(order);

        for (int r = 0; r < numStarts; r++) {
//            if (r > 0) scorer.shuffleVariables();

            gsp(scorer);

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

    private void gsp(TeyssierScorer scorer) {
        int maxDepth = this.gspDepth == -1 ? Integer.MAX_VALUE : this.gspDepth;

        List<Node> order;
        List<Node> best = scorer.getOrder();

        do {
            order = gspVisit(scorer, scorer.getNumEdges(), maxDepth,
                    1, new HashSet<>(), null, null);
            if (order != null) {
                best = order;
            }
        } while (order != null);

        scorer.score(best);
    }

    private List<Node> gspVisit(TeyssierScorer scorer, int e, int maxDepth, int depth, Set<Node> path,
                                Node vv, Node ww) {
        if (depth > maxDepth) return null;
        if (scorer.getNumEdges() < e) return scorer.getOrder();
        if (scorer.getNumEdges() > e) return null;

//        path.add(vv);
        path.add(ww);
        Graph graph0 = getGraph(scorer);
        scorer.bookmark();

        for (Edge edge : graph0.getEdges()) {
            scorer.goToBookmark();

            Node v = Edges.getDirectedEdgeTail(edge);
            Node w = Edges.getDirectedEdgeHead(edge);

//            if (path.contains(v)) continue;
            if (path.contains(w)) continue;

            Set<Node> parentsv = new HashSet<>(graph0.getParents(v));
            Set<Node> parentsw = new HashSet<>(graph0.getParents(w));
            parentsw.remove(v);

            if (parentsv.equals(parentsw)) {
                scorer.swap(v, w);

                List<Node> order = gspVisit(scorer, e, maxDepth, depth + 1, path, v, w);

                if (order != null) {
                    return order;
                }
            }
        }

        path.remove(ww);
        return null;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
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

    public void setGspDepth(int gspDepth) {
        this.gspDepth = gspDepth;
    }

    public void setReturnCpdag(boolean returnCpdag) {
        this.returnCpdag = returnCpdag;
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
