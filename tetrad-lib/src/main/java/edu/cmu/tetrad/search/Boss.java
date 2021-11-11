package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodePair;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.*;

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

        if (score != null && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(score);
        } else if (test != null) {
            scorer = new TeyssierScorer(test);
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

        do {
            bossLoop(scorer);
        } while (extraSwapsLoop(scorer, depth == -1 ? Integer.MAX_VALUE : depth));

        return scorer.getOrder();
    }

    private void bossLoop(TeyssierScorer scorer) {
        scorer.bookmark();

        double score;

        do {
            score = scorer.score();

            for (Node v : scorer.getOrder()) {
                scorer.moveToLast(v);

                double score2 = NEGATIVE_INFINITY;

                if (scorer.score() > score) {
                    if (satisfiesKnowledge(scorer.getOrder())) {
                        scorer.bookmark();
                        score2 = scorer.score();
                    }
                }

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

    private boolean extraSwapsLoop(TeyssierScorer scorer, int depth) {
        List<Node> b = scorer.getOrder();

        for (Node r1 : b) {
            for (Node r2 : scorer.getParents(r1)) {
                List<NodePair> pairs = new ArrayList<>();

                for (Node v : b) {
                    if (triangle(v, r1, r2, scorer)) {
                        pairs.add(new NodePair(v, r1));
                        pairs.add(new NodePair(v, r2));
                    }
                }

                Collections.shuffle(pairs);

                if (doPairs(scorer, pairs, depth)) return true;

                List<NodePair> _pairs = new ArrayList<>(pairs);
                Collections.shuffle(_pairs);

                if (doPairs(scorer, _pairs, depth)) return true;
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (" + "Simultaneous Moves" + ")");
        }

        return false;
    }

    private boolean doPairs(TeyssierScorer scorer, List<NodePair> pairs, int depth) {
//        System.out.println("doing # pairs = " + pairs);
        if (depth > pairs.size()) depth = pairs.size();

        ChoiceGenerator iterator = new ChoiceGenerator(pairs.size(), depth);
        int[] choice;

        while ((choice = iterator.next()) != null) {
            if (choice.length == 0) continue;
            double score = scorer.score();
            scorer.bookmark();

            for (int i : choice) {
                flip(scorer, pairs.get(i).getFirst(), pairs.get(i).getSecond());
            }

            if (scorer.score() > score && satisfiesKnowledge(scorer.getOrder())) {
                return true;
            }

            scorer.goToBookmark();
        }

        return false;
    }

    private void removeEdgesInTrianglesLoop(TeyssierScorer scorer) {
//        while (doNumTriangles(scorer)) ;
//
//        doNumTriangles(scorer);
//        doNumTriangles(scorer);
//        doNumTriangles(scorer);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (" + "Simultaneous Moves" + ")");
        }
    }

    private boolean doNumTriangles(TeyssierScorer scorer) {
        boolean changed = false;

        for (Node x : variables) {
            for (Node y : scorer.getAdjacentNodes(x)) {
                if (!scorer.adjacent(x, y)) continue;

                List<Node> xy = new ArrayList<>();
                xy.add(x);
                xy.add(y);

                List<Node> t = new ArrayList<>();

                for (Node v : variables) {
                    if (triangle(x, y, v, scorer)) {
                        t.add(v);
                    }
                }

                if (t.isEmpty()) continue;

                DepthChoiceGenerator gen = new DepthChoiceGenerator(t.size(), t.size());
                int[] choice;

                W:
                while ((choice = gen.next()) != null) {
                    List<Node> c = GraphUtils.asList(choice, t);

                    for (Node t1 : c) {
                        ChoiceGenerator gen2 = new ChoiceGenerator(2, 2);
                        int[] choice2;

                        while ((choice2 = gen2.next()) != null) {
                            scorer.bookmark();

                            List<Node> c2 = GraphUtils.asList(choice2, xy);
                            flip(scorer, c2.get(choice2[0]), t1);

                            if (!scorer.getParents(y).contains(x)) {
                                changed = true;
                                break W;
                            }
                        }

                        scorer.goToBookmark();
                    }
                }
            }
        }

        return changed;
    }

    private void flip(TeyssierScorer scorer, Node x, Node y) {
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
        Set<List<Node>> path = new HashSet<>();
        gspLoop(scorer, path, 1, gspDepth == -1 ? Integer.MAX_VALUE : gspDepth);
        return scorer.getOrder();
    }

    private double gspLoop(final TeyssierScorer scorer, final Set<List<Node>> path, final int depth,
                           final int gspDepth) {
        if (depth > gspDepth) return scorer.score();
        List<Node> order = scorer.getOrder();
        double score = scorer.score();

        for (Node w : scorer.getOrder()) {
            for (Node v : scorer.getParents(w)) {
                if (covered(scorer, v, w)) {
                    scorer.swap(v, w);

                    if (scorer.score() >= score && !path.contains(order)) {
                        path.add(order);

                        double _score = gspLoop(scorer, path, score == scorer.score() ? depth + 1 : 1, gspDepth);

                        if (_score > score) {
                            return _score;
                        }

//                        path.remove(order);
                    }
                }
            }
        }

        scorer.evaluate(order);
        return score;
    }

    // Checks if v-->w is covered.
    private boolean covered(TeyssierScorer scorer, Node v, Node w) {
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
    public Graph getGraph(List<Node> order, boolean cpdag) {
        TeyssierScorer scorer;
        DataModel dataModel;

        if (score != null) {
            scorer = new TeyssierScorer(score);
            scorer.setKnowledge(knowledge);
            dataModel = score.getData();
        } else {
            scorer = new TeyssierScorer(test);
            scorer.setKnowledge(knowledge);
            dataModel = test.getData();
        }

        scorer.setScoreType(scoreType);

        scorer.evaluate(order);

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

    public enum Method {BOSS, SP, GSP}
}
