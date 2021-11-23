package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
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
    private int triangleDepth = -1;
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
            scorer.setParentCalculation(parentCalculation);
        } else if (test != null) {
            scorer = new TeyssierScorer(test);
            scorer.setParentCalculation(parentCalculation);
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
            } else if (method == Method.AGSP) {
                perm = agsp(scorer);
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

        while (true) {
            bossLoop(scorer);
            double s = scorer.score();
            triangleLoop(scorer);
            if (scorer.score() <= s) break;
        }

        return scorer.getOrder();
    }

    private void bossLoop(TeyssierScorer scorer) {
        double s0;
        scorer.bookmark();

        do {
            s0 = scorer.score();

            double s = NEGATIVE_INFINITY;
            for (Node v : scorer.getOrder()) {
                for (int i = scorer.size() - 1; i >= 0; i--) {
                    scorer.moveTo(v, i);

                    if (scorer.score() > s) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            s = scorer.score();
                            scorer.bookmark();
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

    private void bossLoop2(TeyssierScorer scorer) {

        double score0 = NEGATIVE_INFINITY;

        while (scorer.score() > score0) {
            score0 = scorer.score();

            List<NodePair> allPairs = new ArrayList<>();

            for (int i = scorer.size() - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                    Node w = scorer.get(i);
                    Node v = scorer.get(j);

                    NodePair pair = new NodePair(w, v);
                    allPairs.add(pair);
                }
            }

            for (NodePair pair : allPairs) {
                Node w = pair.getFirst();
                Node v = pair.getSecond();

                double score = scorer.score();
                scorer.swap(w, v);

                if (!(scorer.score() > score)) {
                    scorer.swap(w, v);
                }
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

    private void triangleLoop(TeyssierScorer scorer) {
        if (triangleDepth == 0) return;

        Set<NodePair> path = new HashSet<>();
        int _depth = triangleDepth == -1 ? Integer.MAX_VALUE : triangleDepth;
        double s = scorer.score();

        for (int i = scorer.size() - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                Node x = scorer.get(i);
                Node y = scorer.get(j);
                if (!scorer.adjacent(x, y)) continue;

                List<NodePair> ZZ = new ArrayList<>();

                for (Node z : scorer.getOrder()) {
                    if (scorer.triangle(x, y, z)) {
                        ZZ.add(new NodePair(x, z));
                        ZZ.add(new NodePair(y, z));
                    }
                }

                List<Node> order = scorer.getOrder();

                triangleVisit(scorer, ZZ, path, _depth);

                if (scorer.score() > s) {
                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges() + " Score = "
                                + scorer.score() + " (Triangle)");
                    }

                    return;
                }

                scorer.evaluate(order);
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Triangle)");
        }
    }

    private void triangleVisit(TeyssierScorer scorer, List<NodePair> ZZ, Set<NodePair> path, int depth) {
        if (path.size() > depth) return;

        double score = scorer.score();
        List<Node> order = scorer.getOrder();

        for (NodePair z : ZZ) {
            if (path.contains(z)) continue;

            scorer.swap(z.getFirst(), z.getSecond());
            path.add(z);

            triangleVisit(scorer, ZZ, path, depth);

            path.remove(z);

            if (scorer.score() > score) {
                return;
            }

            scorer.evaluate(order);
        }
    }

    private List<Node> gsp(TeyssierScorer scorer) {
        Set<NodePair> path = new HashSet<>();

        List<NodePair> edges = new ArrayList<>();

        for (int i = 0; i < scorer.size(); i++) {
            for (int j = i + 1; j < scorer.size(); j++) {
                edges.add(new NodePair(scorer.get(i), scorer.get(j)));
            }
        }

        double score = scorer.score();
        double _score = score;

        do {
            score = _score;
            List<Node> order = gspLoop(scorer, 1, gspDepth == -1 ? Integer.MAX_VALUE : gspDepth, path, edges, null);

            path.clear();

            scorer.evaluate(order);
            List<Node> order2 = gspLoopOneEqualStep(scorer, 1, gspDepth == -1 ? Integer.MAX_VALUE : gspDepth, path, edges, null);

            if (order2 == null) {
                scorer.evaluate(order);
                break;
            }

            scorer.evaluate(order2);
            _score = scorer.score();
        } while (_score > score);

        return scorer.getOrder();
    }

    private List<Node> gspLoop(TeyssierScorer scorer, int depth, int gspDepth, Set<NodePair> path, List<NodePair> edges, NodePair edge) {
        if (edge != null) {
            if (path.contains(edge)) return scorer.getOrder();
            path.add(edge);
        }

        List<Node> order = scorer.getOrder();

        scorer.evaluate(order);
        List<Node> newOrder = scorer.getOrder();
        double newScore = scorer.score();

        if (depth > gspDepth) {
            return scorer.getOrder();
        }

        double score = scorer.score();
        List<Node> best = scorer.getOrder();

        for (NodePair _edge : edges) {
            Node x = _edge.getFirst();
            Node y = _edge.getSecond();

            scorer.evaluate(order);

            if (covered(x, y, scorer)) {
                flip(scorer, x, y);

                if (scorer.score() > newScore) {
                    List<Node> _order = gspLoop(scorer, depth + 1, gspDepth, path, edges, edge);

                    scorer.evaluate(_order);

                    if (scorer.score() > score) {
                        newScore = scorer.score();
                        newOrder = scorer.getOrder();
                    }
                }
            }
        }

        path.remove(edge);
        return newOrder;
    }

    private List<Node> gspLoopOneEqualStep(TeyssierScorer scorer, int depth, int gspDepth, Set<NodePair> path, List<NodePair> edges, NodePair edge) {
        if (edge != null) {
            if (path.contains(edge)) return scorer.getOrder();
            path.add(edge);
        }

        List<Node> order = scorer.getOrder();

        scorer.evaluate(order);
//        List<Node> newOrder = scorer.getOrder();
        double newScore = scorer.score();

        for (NodePair _edge : edges) {
            Node x = _edge.getFirst();
            Node y = _edge.getSecond();

            scorer.evaluate(order);

            if (covered(x, y, scorer)) {
                flip(scorer, x, y);

                if (scorer.score() == newScore) {
                    return scorer.getOrder();
                }
            }
        }

//        path.remove(edge);
        return null;
    }

    private List<Node> gsp2(TeyssierScorer scorer) {
        Graph complete = GraphUtils.completeGraph(scorer.getGraph(false));
        List<Edge> edges = new ArrayList<>(complete.getEdges());
        edges.addAll(new ArrayList<>(edges));

        double score0 = NEGATIVE_INFINITY;

        while (scorer.score() > score0) {
            score0 = scorer.score();

            for (Edge edge : edges) {
                Node w = edge.getNode1();
                Node v = edge.getNode2();

                if (covered(w, v, scorer)) {
                    System.out.println(w + "---" + v + " covered");

                    double score = scorer.score();
                    scorer.swap(w, v);

                    if (scorer.score() >= score) {
                        System.out.println(scorer.score() + " = score");
                        continue;
                    }

                    scorer.swap(w, v);
                }
            }
        }

        return scorer.getOrder();
    }

    private List<Node> gsp3(TeyssierScorer scorer) {
        List<NodePair> edges = new ArrayList<>();

        for (int i = 0; i < scorer.size(); i++) {
            for (int j = i + 1; j < scorer.size(); j++) {
                edges.add(new NodePair(scorer.get(i), scorer.get(j)));
            }
        }

        return gsp3Visit(scorer, null, scorer.getOrder(), edges, new ArrayList<>(), 1);
    }

    private List<Node> gsp3Visit(TeyssierScorer scorer, NodePair edge, List<Node> order,
                                 List<NodePair> edges, List<NodePair> path, int depth) {
        System.out.println("depth = " + depth);
        if (depth > (gspDepth == -1 ? Integer.MAX_VALUE : gspDepth)) return scorer.getOrder();

        if (edge != null) {
            if (path.contains(edge)) return order;
            path.add(edge);
        }

        scorer.evaluate(order);
        List<Node> newOrder = scorer.getOrder();
        double newScore = scorer.score();

        for (NodePair _edge : edges) {
            Node x = _edge.getFirst();
            Node y = _edge.getSecond();

            scorer.evaluate(order);

            if (covered(x, y, scorer)) {
                double score = scorer.score();
                flip(scorer, x, y);

                if (scorer.score() >= newScore) {
                    List<Node> _order = gsp3Visit(scorer, _edge, scorer.getOrder(), edges, path, depth + 1);
                    scorer.evaluate(_order);

                    if (scorer.score() > newScore) {
                        newScore = scorer.score();
                        newOrder = scorer.getOrder();
                        return scorer.getOrder();
                    }
                }
            }
        }

        path.remove(edge);
        return newOrder;
    }

    private void flip(TeyssierScorer scorer, Node x, Node y) {
        if (!scorer.adjacent(x, y)) throw new IllegalArgumentException("Expecting adjacent X, Y");

        if (scorer.getParents(x).contains(y)) {
            scorer.moveTo(x, scorer.indexOf(y));
//            System.out.println("Moving " + x + " to before " + y + ": " + scorer.getOrder() + " count = " + -scorer.score());
        } else if (scorer.getParents(y).contains(x)) {
//            System.out.println("Moving " + y + " to before " + x + ": " + scorer.getOrder() + " count = " + -scorer.score());
            scorer.moveTo(y, scorer.indexOf(x));
        }
    }

    private void flip2(TeyssierScorer scorer, Node x, Node y) {
        if (!scorer.adjacent(x, y)) throw new IllegalArgumentException("Expecting adjacent X, Y");
        if (!scorer.getParents(y).contains(x)) {
            throw new IllegalArgumentException("Expecting " + x + "-->" + y);
        }

//        if (scorer.getParents(x).contains(y)) {
//        scorer.moveTo(x, scorer.indexOf(y));
//            System.out.println("Moving " + x + " to before " + y + ": " + scorer.getOrder() + " count = " + -scorer.score());
//        } else if (scorer.getParents(y).contains(x)) {
//            System.out.println("Moving " + y + " to before " + x + ": " + scorer.getOrder() + " count = " + -scorer.score());
        scorer.moveTo(y, scorer.indexOf(x));
//        }
    }

    private List<Node> agsp(@NotNull TeyssierScorer scorer) {
        double oldScore;
        double newScore;
        int round = 0;

        do {
            if (verbose) {
                System.out.println("AGSP Round = " + ++round);
            }

            oldScore = scorer.score();

            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < scorer.size(); i++) {
                    for (int j = i + 1; j < scorer.size(); j++) {
                        Node x = scorer.get(i);
                        Node y = scorer.get(j);

                        if (!scorer.adjacent(y, x)) continue;

                        if (scorer.getParents(x).contains(y)) {
                            Node r = x;
                            x = y;
                            y = r;
                        }

                        scorer.bookmark();

                        double s1 = scorer.score();

                        if (scorer.indexOf(x) < scorer.indexOf(y)) {
                            scorer.moveTo(y, scorer.indexOf(x));
                        } else {
                            scorer.moveTo(y, scorer.indexOf(x) + 1);
                        }

                        double s2 = scorer.score();

                        if (s2 >= s1) {
                            continue;
                        }

                        scorer.goToBookmark();
                    }
                }
            }

            newScore = scorer.score();
        } while (newScore > oldScore);

        return scorer.getOrder();
    }

    // Checks if v-->w is covered.
    private boolean covered(Node v, Node w, TeyssierScorer scorer) {
        if (!scorer.adjacent(v, w)) throw new IllegalArgumentException(
                "When checking if <X, Y> is covered, X and Y should be adjacent.");
        Set<Node> pw = new HashSet<>(scorer.getParents(w));
        pw.remove(v);
        Set<Node> pv = new HashSet<>(scorer.getParents(v));
        pv.remove(w);
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

    public void setTriangleDepth(int triangleDepth) {
        if (triangleDepth < -1) throw new IllegalArgumentException("Depth should be >= -1.");

        this.triangleDepth = triangleDepth;
    }

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public enum Method {BOSS, SP, AGSP}
}
