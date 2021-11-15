package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static java.lang.Double.NEGATIVE_INFINITY;


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
    private Method method = Method.BOSS;
    private boolean verbose = false;
    private int gspDepth = -1;
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();
    private boolean firstRunUseDataOrder = false;
    private int depth = -1;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.GrowShrinkMb;
    private TeyssierScorer scorer;

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

        if (score != null && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(score);
//            scorer.setParentCalculation(parentCalculation);
        } else if (test != null) {
            scorer = new TeyssierScorer(test);
//            scorer.setParentCalculation(parentCalculation);
        } else {
            throw new IllegalArgumentException("Need a score (not GraphScore) or a test.");
        }

        scorer.setKnowledge(knowledge);
        scorer.setScoreType(scoreType);

        scorer.setCachingScores(cachingScores);

        List<Node> bestPerm = new ArrayList<>(order);
        scorer.evaluate(bestPerm);
        double best = scorer.score();

        for (int r = 0; r < numStarts; r++) {
            if (firstRunUseDataOrder) {
                if (r > 0)
                    scorer.shuffleVariables();
            } else {
                scorer.shuffleVariables();
            }

            List<Node> _order = scorer.getOrder();
            makeValidKnowledgeOrder(_order);
            scorer.evaluate(_order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else if (method == Method.GSP) {
                perm = gsp(scorer);
            } else {
                throw new IllegalArgumentException("Unrecognized method: " + method);
            }

            scorer.evaluate(perm);

            if (scorer.score() > best) {
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

    public int getNumEdges() {
        return scorer.getNumEdges();
    }

    private void makeValidKnowledgeOrder(List<Node> order) {
        if (!knowledge.isEmpty()) {
            order.sort((o1, o2) -> {
                if (o1.getName().equals(o2.getName())) {
                    return 0;
                } else if (knowledge.isRequired(o1.getName(), o2.getName())) {
                    return 1;
                } else if (knowledge.isRequired(o2.getName(), o1.getName())) {
                    return -1;
                } else if (knowledge.isForbidden(o2.getName(), o1.getName())) {
                    return -1;
                } else if (knowledge.isForbidden(o1.getName(), o2.getName())) {
                    return 1;
                } else {
                    return 1;
                }
            });
        }
    }

    private boolean satisfiesKnowledge(List<Node> order) {
        if (!knowledge.isEmpty()) {
            for (int i = 0; i < order.size(); i++) {
                for (int j = i + 1; j < order.size(); j++) {
                    if (knowledge.isForbidden(order.get(i).getName(), order.get(j).getName())) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public List<Node> boss(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal from). Score each causal from by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated. Then try to remove edges X--Y by reorienting edge
        // in triangles X--V--Y, repeatedly, until no more such edges can be removed.
        // Apparently, it suffices in the oracle to consider only on such triangle at
        // a time, though to hedge bets we allow the user the user to specify a maximum
        // number of triangles to consider per such edge.

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

//        do {
        bossLoop(scorer);
//        } while (triangleLoop1(scorer));

        return scorer.getOrder();
    }


    private void bossLoop0(TeyssierScorer scorer) {
        scorer.bookmark();

        double s0;

        do {
            s0 = scorer.score();

            for (Node v : scorer.getOrder()) {
                scorer.moveToLast(v);

                double s = NEGATIVE_INFINITY;

                if (scorer.score() > s0) {
                    if (satisfiesKnowledge(scorer.getOrder())) {
                        scorer.bookmark();
                        s = scorer.score();
                    }
                }

                while (scorer.promote(v)) {
                    if (scorer.score() > s && scorer.score() > s0) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            scorer.bookmark();
                            s = scorer.score();
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (scorer.score() > s0);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

    private void bossLoop(TeyssierScorer scorer) {
        scorer.bookmark();

        double s0;

        do {
            s0 = scorer.score();
//            List<Node> _best = null;

            double s = NEGATIVE_INFINITY;

//            V:
            for (Node v : scorer.getOrder()) {
//                BestMut bestMut = getBestMut(scorer, v, s);
//
////                System.out.println(bestMut.bestMut.size() + " score = " + bestMut.score);
//
//                for (int w = bestMut.bestMut.size() - 1; w >= 0; w--) {
//                    List<Node> best = bestMut.bestMut.get(w);
//                    scorer.evaluate(best);
//
////                    if (scorer.score() > s) {
//                    if (satisfiesKnowledge(scorer.getOrder())) {
//                        scorer.bookmark();
//                        s = bestMut.score;
//                        _best = best;
////                        break V;
//                    }
////                    }
//                }

                for (int i = scorer.size() - 1; i >= 0; i--) {
                    scorer.moveTo(v, i);

                    if (scorer.score() > s && scorer.score() > s0) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            scorer.bookmark();
                            s = scorer.score();
                        }
                    }
                }

//                scorer.evaluate(_best);
                scorer.goToBookmark();
            }
        } while (scorer.score() > s0);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

    private BestMut getBestMut(TeyssierScorer scorer, Node v, double score) {
        List<List<Node>> best = new ArrayList<>();
        double s = NEGATIVE_INFINITY;

        for (int i = scorer.size() - 1; i >= 0; i--) {
            scorer.moveTo(v, i);

            if (scorer.score() > s && scorer.score() >= score) {
                best.clear();
                s = scorer.score();
            }

            if (scorer.score() >= score) {
                best.add(scorer.getOrder());
            }
        }

        BestMut bestMut = new BestMut();
        bestMut.bestMut = best;
        bestMut.score = s;

//        for (List<Node> o : best) {
//            scorer.evaluate(o);
//            System.out.println("scoreee = " + scorer.score());
//        }

        return bestMut;
    }

    private void bossLoop2(TeyssierScorer scorer) {
        scorer.bookmark();
        double score;

        do {
            score = scorer.score();

            for (Node v : scorer.getOrder()) {
                double score2 = scorer.score();

                while (scorer.promote(v)) {
                    if (scorer.score() > score2 && scorer.score() > score) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            scorer.bookmark();
                            score2 = scorer.score();
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (scorer.score() > score);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

//    private void bossLoop(TeyssierScorer scorer) {
//        scorer.bookmark();
//        double s0 = NEGATIVE_INFINITY;
//        List<Node> best = null;
//
//        while (true) {
//            for (Node v : scorer.getOrder()) {
//                double s = scorer.score();
//                scorer.bookmark();
//
//                for (int i = 0; i < scorer.size(); i++) {
//                    scorer.moveTo(v, i);
//
//                    if (scorer.score() >= s) {
//                        if (satisfiesKnowledge(scorer.getOrder())) {
//                            s = scorer.score();
//                            scorer.bookmark();
//                        }
//                    }
//                }
//
//                scorer.goToBookmark();
//            }
//
//            if (scorer.score() > s0) {
//                s0 = scorer.score();
//                best = scorer.getOrder();
//            } else {
//                break;
//            }
//
////            if (scorer.score() <= s0) break;
//        }
//
//        if (verbose) {
//            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
//        }
//
//        scorer.evaluate(best);
//    }

//    public List<Node> bossSearchAllIndices(TeyssierScorer scorer) {
//
//        // Take each variable in turn and try moving it to each position to the left (i.e.,
//        // try promoting it in the causal order). Score each causal order by building
//        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
//        // Markov blanket algorithm, to find putative parents for each node.) Finally,
//        // place the node in whichever position yielded the highest score. Do this
//        // for each variable in turn. Once you're done, do it all again, until no more
//        // variables can be relocated.
//        double overall = Double.POSITIVE_INFINITY;
//        List<Node> best = null;
//
//        while (true) {
//            for (Node v : scorer.getOrder()) {
//                scorer.moveTo(v, 0);
//                double bestScore = scorer.score();
//                scorer.bookmark(1);
//
//                while (scorer.moveRight(v)) {
//                    if (scorer.score() <= bestScore) {
//                        bestScore = scorer.score();
//                        scorer.bookmark(1);
//                    }
//                }
//
//                scorer.goToBookmark(1);
//            }
//
//            if (scorer.score() < overall) {
//                overall = scorer.score();
//                best = scorer.getOrder();
//                if (verbose) {
//                    System.out.println("Updated order = " + scorer.getOrder());
//                }
//            } else {
//                break;
//            }
//        }
//
//        return best;
//    }

    private boolean triangleLoop1(TeyssierScorer scorer) {
        List<Node> b = scorer.getOrder();

        for (Node v : b) {
            for (int _r1 = 0; _r1 < b.size(); _r1++) {
                for (int _r2 = _r1 + 1; _r2 < b.size(); _r2++) {
                    Node r1 = scorer.get(_r1);
                    Node r2 = scorer.get(_r2);

                    if (!adjacent(r1, r2, scorer)) continue;
//                    if (!covered(r2, r2, scorer)) continue;
//
                    // r1->v<-r2, try r1->v->r2, r1<-v->r1, r1<-v<-r2
                    if (scorer.getParents(v).contains(r1) && scorer.getParents(v).contains(r2)) {
                        double score = scorer.score();
                        scorer.bookmark();

                        flip(r1, v, scorer);
                        flip(r2, v, scorer);

                        if (scorer.score() > score && satisfiesKnowledge(scorer.getOrder())) {
                            return true;
                        }

                        scorer.goToBookmark();

//                        flip(v, r1, scorer);
//
//                        if (scorer.score() > score && satisfiesKnowledge(scorer.getOrder())) {
//                            return true;
//                        }
//
//                        scorer.goToBookmark();
//
//                        flip(v, r2, scorer);
//
//                        if (scorer.score() > score && satisfiesKnowledge(scorer.getOrder())) {
//                            return true;
//                        }
//
//                        scorer.goToBookmark();
                    }

                    // r1<-v->r2, try r1->v<-r1
                    if ((scorer.getParents(r1).contains(v) && scorer.getParents(r2).contains(v))) {
                        double score = scorer.score();
                        scorer.bookmark();

                        flip(v, r1, scorer);
                        flip(v, r2, scorer);

                        if (scorer.score() > score && satisfiesKnowledge(scorer.getOrder())) {
                            return true;
                        }

                        scorer.goToBookmark();
                    }
                }
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (" + "Simultaneous Moves" + ")");
        }

        return false;
    }

    private void flip(Node x, Node y, TeyssierScorer scorer) {
//        scorer.swap(x, y);

        if (scorer.getParents(x).contains(y)) {
            scorer.moveTo(x, scorer.indexOf(y));
        } else if (scorer.getParents(y).contains(x)) {
            scorer.moveTo(y, scorer.indexOf(x));
        }
    }

    private boolean triangle(Node v, Node r1, Node r2, TeyssierScorer scorer) {
        return adjacent(v, r1, scorer) && adjacent(v, r2, scorer) && adjacent(r1, r2, scorer);
    }

    private boolean adjacent(Node v, Node r1, TeyssierScorer scorer) {
        return scorer.getParents(v).contains(r1) || scorer.getParents(r1).contains(v);
    }

    private List<Node> gsp(TeyssierScorer scorer) {
        double score = scorer.score();

        while (true) {
            List<Node> order = gspLoop(scorer, 1, gspDepth == -1 ? Integer.MAX_VALUE : gspDepth);
            scorer.evaluate(order);
            if (scorer.score() <= score) break;
            score = scorer.score();
        }

        return scorer.getOrder();
    }

    private List<Node> gspLoop(final TeyssierScorer scorer, final int depth, final int gspDepth) {
        List<Node> _order = scorer.getOrder();

        if (depth > gspDepth) {
            return _order;
        }

        double score = scorer.score();

        for (Node w : scorer.getOrder()) {
            for (Node v : scorer.getParents(w)) {
                if (covered(v, w, scorer)) {
                    scorer.swap(v, w);

                    if (scorer.score() > score) {
                        return gspLoop(scorer, depth + 1, gspDepth);
                    }

                    scorer.evaluate(_order);
                }
            }
        }

        return _order;
    }

    // Checks if v-->w is covered.
    private boolean covered(Node v, Node w, TeyssierScorer scorer) {
        Set<Node> pw = new HashSet<>(scorer.getParents(w));
        pw.remove(v);
        Set<Node> pv = new HashSet<>(scorer.getParents(v));
        return pv.equals(pw);
    }

    public List<Node> sp(TeyssierScorer scorer) {
        double minScore = NEGATIVE_INFINITY;
        List<Node> minP = null;

        List<Node> variables = scorer.getOrder();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);
            scorer.evaluate(p);

            if (scorer.score() > minScore) {
                minScore = scorer.score();
                minP = p;
            }
        }

        return minP;
    }

    @NotNull
    public Graph getGraph(boolean cpdag) {
        DataModel dataModel;

        if (score != null) {
            dataModel = score.getData();
        } else {
            dataModel = test.getData();
        }

        Graph graph1 = scorer.getGraph(cpdag);

        if (dataModel != null) {
            double bic = new BicEst(score).getValue(null, graph1, dataModel);
            graph1.addAttribute("BIC", bic);
        }

        return graph1;
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

    public void setScoreType(TeyssierScorer.ScoreType scoreType) {
        this.scoreType = scoreType;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public void setGspDepth(int gspDepth) {
        this.gspDepth = gspDepth;
    }

    public void setFirstRunUseDataOrder(boolean firstRunUseDataOrder) {
        this.firstRunUseDataOrder = firstRunUseDataOrder;
    }

    public void setMaxNumTriangles(int maxNumTriangles) {
        if (maxNumTriangles < -1) throw new IllegalArgumentException("Depth must be >= -1.");
    }

    public int getDepth() {
        return depth;
    }

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");

        this.depth = depth;
    }

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public enum Method {BOSS, SP, GSP}

    private class BestMut {
        List<List<Node>> bestMut;
        double score;
    }
}
