package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
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
    private boolean breakTies = false;
    private int gspDepth = -1;

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
            } else if (method == Method.GSP) {
                perm = gsp(scorer);
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

        /*
         Take each variable in turn and try moving it to each position to the left (i.e.,
         try promoting it in the causal order). Score each causal order by building
         a DAG using something like the K2 method. (We actually use Grow-Shrink, a
         Markov blanket algorithm, to find putative parents for each node.) Finally,
         place the node in whichever position yielded the highest score. Do this
         for each variable in turn. Once you're done, do it all again, until no more
         variables can be relocated.
        */
        if (breakTies) {
            while (true) {
                promoteLoop(scorer);
                double score = scorer.score();
                breakTie(scorer);
                if (scorer.score() == score) break;
            }
        } else {
            promoteLoop(scorer);
        }

        return scorer.getOrder();
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

        if (breakTies) {
//            while (true) {
            allIndicesLoop(scorer);
//                double score = scorer.score();
            breakTie(scorer);
////                if (scorer.score() == score) break;
//            }
        } else {
            allIndicesLoop(scorer);
        }

        return scorer.getOrder();
    }

    private void promoteLoop(TeyssierScorer scorer) {
        boolean reduced;

        do {
            reduced = false;

            for (Node v : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.bookmark();

                while (scorer.promote(v)) {
                    if (scorer.score() < bestScore) {
                        reduced = true;
                        scorer.bookmark();
                    }
                }

                scorer.goToBookmark();
            }

            if (verbose) {
                System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Promotion)");
            }
        } while (reduced);

    }

    private void allIndicesLoop(TeyssierScorer scorer) {
        boolean reduced;

        do {
            reduced = false;
            scorer.bookmark();

            for (Node v : scorer.getOrder()) {
                double bestScore = scorer.score();
                scorer.moveToFirst(v);

                if (scorer.score() < bestScore) {
                    bestScore = scorer.score();
                    reduced = true;
                    scorer.bookmark();
                }

                while (scorer.demote(v)) {
                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        reduced = true;
                        scorer.bookmark();
                    }
                }

                scorer.goToBookmark();
            }

            if (verbose) {
                System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (All indices)");
            }
        } while (reduced);

    }

    private void breakTie(TeyssierScorer scorer) {
        String label = "Break Tie";

//        System.out.println("Permutation: " + scorer.getOrder());

        List<Node> order = scorer.getOrder();
        double score = scorer.score();

        // Pick a v.

        I:
        for (int i = 0; i < order.size(); i++) {
            Node v = order.get(i);

            // ...and promote it to another position that yields the same score.
            for (int j = 0; j < order.indexOf(v); j++) {
                if (i == j) continue;
                scorer.moveTo(v, j);

                if (scorer.score() > score) {
                    continue;
                }

                if (scorer.score() < score) {
                    return;
                }

//                System.out.println("Moving " + v + " from " + i + " to " + j
//                        + " same score = " + scorer.score() + " " + scorer.getOrder());

                scorer.bookmark();

                // Then pick an r != v
                for (int k = 0; k < scorer.getOrder().size(); k++) {
                    Node r = scorer.getOrder().get(k);
                    if (r == v) continue;

                    Set<Node> vMb = scorer.getMb(scorer.indexOf(v));
                    Set<Node> rMb = scorer.getMb(scorer.indexOf(r));

                    if (!vMb.contains(r) && !rMb.contains(v)) {
                        continue;
                    }

                    scorer.swap(v, r);

                    if (scorer.score() < score) {
                        return;
                    }

                    scorer.goToBookmark();
                }
            }
        }

        scorer.goToBookmark();

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (" + label + ")");
        }
    }

    private List<Node> gsp(TeyssierScorer scorer) {
        int maxDepth = this.gspDepth == -1 ? Integer.MAX_VALUE : this.gspDepth;

        List<Node> order;
        List<Node> best = scorer.getOrder();
        double s = scorer.score(best);



        while (true) {
            order = gspVisit(scorer, maxDepth,10, new HashSet<>(), null);
            double score = scorer.score(order);

            breakTie(scorer);

            if (score == s) {
                break;
            } else {
                s = score;
            }
        }

        scorer.score(best);
        return best;
    }

    private List<Node> gspVisit(TeyssierScorer scorer, int maxDepth, int depth, Set<Node> path,
                                Node ww) {
        maxDepth = maxDepth == -1 ? 10000 : maxDepth;
        if (depth > maxDepth) return scorer.getOrder();

        List<Node> order = scorer.getOrder();
        double s = scorer.score();

        path.add(ww);
        Graph graph0 = getGraph(scorer.getOrder(), false);
        scorer.bookmark();

        for (Edge edge : graph0.getEdges()) {
            if (Edges.isUndirectedEdge(edge)) continue;

            Node v = Edges.getDirectedEdgeTail(edge);
            Node w = Edges.getDirectedEdgeHead(edge);

            if (path.contains(w)) continue;

            Set<Node> mbv = scorer.getMb(scorer.indexOf(v));
            Set<Node> mbw = scorer.getMb(scorer.indexOf(w));
            mbw.remove(v);

            if (mbv.equals(mbw)) {
                System.out.println("Swapping " + v + " and " + w + " " + graph0.getEdge(v, w));

                scorer.swap(v, w);

                if (scorer.score() < s) {
                    return gspVisit(scorer, maxDepth, depth + 1, path, w);
                }
            }

            scorer.goToBookmark();
        }

        path.remove(ww);
        return order;
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

    public void setBreakTies(boolean breakTies) {
        this.breakTies = breakTies;
    }

    public void setGspDepth(int gspDepth) {
        this.gspDepth = gspDepth;
    }

    public enum Method {BOSS_PROMOTION, BOSS_ALL_INDICES, SP, GSP}
}
