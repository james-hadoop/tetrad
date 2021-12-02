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
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Double.*;
import static java.util.Collections.reverse;


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
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();
    private boolean useDataOrder = false;
    private int triangleDepth = -1;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.GrowShrinkMb;
    private TeyssierScorer scorer;

    public Boss(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
    }

    public Boss(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
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

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (!useDataOrder) {
                scorer.shuffleVariables();
            }

            List<Node> _order = scorer.getOrder();
            makeValidKnowledgeOrder(_order);
            scorer.evaluate(_order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
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

    public List<Node> boss(@NotNull TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal from). Score each causal from by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated. Then try to remove edges X--Y by reorienting edge
        // in triangles X--V--Y, repeatedly, until no more such edges can be removed.
        // Apparently, it suffices in the oracle to consider only on such triangle at
        // a time, though to hedge bets we allow the user to specify a maximum
        // number of triangles to consider per such edge.

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        double s0;
        double s;

        do {
            betterMutation(scorer);
            s0 = scorer.score();
            betterTriMutation2(scorer);
            s = scorer.score();
        } while (s > s0);

        return scorer.getOrder();
    }

    private void betterMutation1(@NotNull TeyssierScorer scorer) {
        double s = NEGATIVE_INFINITY;
        scorer.bookmark();

        while (scorer.score() > s) {
            for (Node v : scorer.getOrder()) {
                for (int i = 0; i < scorer.size(); i++) {
//                for (int i = scorer.size() - 1; i >= 0; i--) {
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
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

    private void betterMutation(@NotNull TeyssierScorer scorer) {
        double s0;
        double s = scorer.score();
        scorer.bookmark();

        do {
            s0 = s;

            for (Node x : scorer.getOrder()) {
                s = NEGATIVE_INFINITY;

                for (int i = scorer.size() - 1; i >= 0; i--) {
                    scorer.moveTo(x, i);

                    if (scorer.score() > s) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            s = scorer.score();
                            scorer.bookmark();
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (s > s0);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }
    }

    private void betterTriMutation(@NotNull TeyssierScorer scorer) {
        if (triangleDepth == 0) return;

        Set<NodePair> path = new HashSet<>();
        int _depth = triangleDepth == -1 ? Integer.MAX_VALUE : triangleDepth;
        MultiSwapRet ret0 = new MultiSwapRet(scorer.getOrder(), scorer.score());

        List<Node> order = scorer.getOrder();
        reverse(order);

        for (int i = 0; i < scorer.size(); i++) {
            Node x = order.get(i);

            for (int j = 1; j < i; j++) {
                Node y = order.get(j);

                if (!scorer.adjacent(x, y)) continue;

                List<NodePair> Z = new ArrayList<>();

                for (Node z : order) {
                    if (scorer.triangle(x, y, z)) {
                        Z.add(new NodePair(x, z));
                        Z.add(new NodePair(y, z));
                    }
                }

                MultiSwapRet best = ret0;

                for (NodePair w : Z) {
                    MultiSwapRet ret = multiSwap(scorer, w, Z, path, _depth);
                    if (ret.score > best.score) {
                        path.remove(w);
                        best = ret;
                        scorer.evaluate(best.order);
                    }
                }

                if (best.score > ret0.score) {
                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges() + " Score = "
                                + scorer.score() + " (Triangle)");
                    }

                    scorer.evaluate(best.order);
                    return;
                }

                scorer.evaluate(ret0.order);
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Triangle)");
        }
    }

    private MultiSwapRet multiSwap(TeyssierScorer scorer, @NotNull NodePair z, List<NodePair> Z,
                                   Set<NodePair> path, int depth) {
        if (path.size() > depth) return new MultiSwapRet(scorer.getOrder(), scorer.score());
        if (path.contains(z)) return new MultiSwapRet(scorer.getOrder(), scorer.score());
        MultiSwapRet ret0 = new MultiSwapRet(scorer.getOrder(), scorer.score());

        path.add(z);
        scorer.swap(z.getFirst(), z.getSecond());

        for (NodePair w : Z) {
            MultiSwapRet ret1 = multiSwap(scorer, w, Z, path, depth);
            if (ret1.score > ret0.score) {
                path.remove(z);
                return ret1;
            }
        }

        scorer.swap(z.getFirst(), z.getSecond());
        path.remove(z);
        return ret0;
    }

    private void betterTriMutation2(@NotNull TeyssierScorer scorer) {
        if (triangleDepth == 0) return;

        int _depth = triangleDepth == -1 ? Integer.MAX_VALUE : triangleDepth;
        MultiSwapRet ret0 = new MultiSwapRet(scorer.getOrder(), scorer.score());

        for (int i = 0; i < scorer.size(); i++) {
            Node x = scorer.get(i);

            for (int j = 1; j < i; j++) {
                Node y = scorer.get(j);

                if (!scorer.adjacent(x, y)) continue;

                Set<Node> W = new HashSet<>();

                for (Node z : scorer.getOrder()) {
                    if (scorer.triangle(x, y, z)) {
                        W.add(x);
                        W.add(y);
                        W.add(z);
                    }
                }

                MultiSwapRet ret = multiSwap2(scorer, W, _depth);

                if (ret.score > ret0.score) {
                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges() + " Score = "
                                + scorer.score() + " (Triangle)");
                    }

                    scorer.evaluate(ret.order);
                    return;
                }

                scorer.evaluate(ret0.order);
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Triangle)");
        }
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    private MultiSwapRet multiSwap2(TeyssierScorer scorer, Set<Node> W, int depth) {
        List<Node> pi = scorer.getOrder();
        depth = Math.min(depth, W.size());
        double s0 = scorer.score();
        List<Node> WW = new ArrayList<>(W);

        int[] indices = new int[WW.size()];
        for (int i = 0; i < WW.size(); i++) indices[i] = scorer.index(WW.get(i));

        MultiSwapRet best = new MultiSwapRet(scorer.getOrder(), scorer.score());

        ChoiceGenerator gen = new ChoiceGenerator(WW.size(), depth);
        int[] choice;

        while ((choice = gen.next()) != null) {
            PermutationGenerator permGen = new PermutationGenerator(depth);
            int[] perm;

            while ((perm = permGen.next()) != null) {
                for (int i = 0; i < depth; i++) {
                    scorer.moveTo(WW.get(choice[i]), indices[choice[perm[i]]]);
                }

                if (scorer.score() > s0) {
                    best = new MultiSwapRet(scorer.getOrder(), scorer.score());
                }
            }
        }

        return best;

//        scorer.evaluate(pi);
//        return new MultiSwapRet(scorer.getOrder(), scorer.score());
    }

    private static class MultiSwapRet {
        public MultiSwapRet(List<Node> pi, double score) {
            this.order = pi;
            this.score = score;
        }

        public List<Node> order;
        public double score;
    }

    private List<Node> gasp(@NotNull TeyssierScorer scorer) {
        if (scorer.size() < 2) return scorer.getOrder();

        List<NodePair> list = new ArrayList<>();

        for (int i = 0; i < scorer.size(); i++) {
            for (int j = i + 1; j < scorer.size(); j++) {
                list.add(new NodePair(scorer.get(i), scorer.get(j)));
            }
        }

        reverse(list);
        Deque<NodePair> deque = new LinkedList<>(list);

        Map<NodePair, Integer> counts = new HashMap<>();
        for (NodePair pair : deque) counts.put(pair, 0);


        do {
            NodePair pair = deque.removeFirst();
            deque.addLast(pair);
            counts.put(pair, counts.get(pair) + 1);

            if (counts.get(deque.peek()) > 2) {
                scorer.goToBookmark();
                continue;
            }

            if (scorer.adjacent(pair.getFirst(), pair.getSecond())) {
                scorer.bookmark();

                double sOld = scorer.score();
                scorer.reverse(pair.getFirst(), pair.getSecond());
                double sNew = scorer.score();

                if (sNew >= sOld) {
                    continue;
                }

                scorer.goToBookmark();
            }
        } while (counts.get(deque.peek()) <= 5);

        return scorer.getOrder();
    }

    public List<Node> sp(@NotNull TeyssierScorer scorer) {
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

    public void setUseDataOrder(boolean useDataOrder) {
        this.useDataOrder = useDataOrder;
    }

    public void setTriangleDepth(int triangleDepth) {
        if (triangleDepth < -1) throw new IllegalArgumentException("Depth should be >= -1.");

        this.triangleDepth = triangleDepth;
    }

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public enum Method {BOSS, SP, GASP}

}
