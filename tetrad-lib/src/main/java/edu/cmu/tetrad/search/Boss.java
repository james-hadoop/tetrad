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

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.util.Collections.shuffle;


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
    private int depth = -1;
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
//        scorer.evaluate(bestPerm);
        double best = NEGATIVE_INFINITY;// scorer.score();

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (!useDataOrder) {
                shuffle(order);
//                scorer.shuffleVariables();
            }

            List<Node> _order = new ArrayList<>(order);//scorer.getOrder();
            makeValidKnowledgeOrder(_order);
            scorer.evaluate(_order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer, start);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
            } else if (method == Method.QUICK_GASP) {
                perm = quickGasp(scorer);
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
            System.out.println("Elapsed time = " + (stop - start) / 1000.0 + " s");
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

    public List<Node> boss(@NotNull TeyssierScorer scorer, long start) {

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        double s0;
        double s;

        do {
            betterMutation(scorer, start);
            s0 = scorer.score();
            betterTriMutation(scorer, start);
            s = scorer.score();
        } while (s > s0);

        return scorer.getOrder();
    }

    private void betterMutation(@NotNull TeyssierScorer scorer, long start) {
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
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (Single Moves)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }
    }

    private void betterTriMutation(@NotNull TeyssierScorer scorer, long start) {
        if (depth == 0) return;

        int _depth = depth == -1 ? Integer.MAX_VALUE : depth;
        MultiSwapRet ret0 = new MultiSwapRet(scorer.getOrder(), scorer.score());

        for (int i = 0; i < scorer.size(); i++) {
            Node x = scorer.get(i);

            for (int j = 1; j < i; j++) {
                Node y = scorer.get(j);

                if (!scorer.adjacent(x, y)) continue;

                Set<Node> W = new HashSet<>();

                for (Node z : scorer.getOrder()) {
                    if (scorer.adjacent(x, z) && scorer.adjacent(y, z)) {
//                        W.add(x);
//                        W.add(y);
                        W.add(z);
                    }
                }

                MultiSwapRet ret = multiSwap(scorer, x, y, W, _depth);

                if (ret.score > ret0.score) {
                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges()
                                + " Score = "
                                + scorer.score()
                                + " (Triangle)"
                                + " Elapsed = " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
                    }

                    scorer.evaluate(ret.order);
                    return;
                }

                scorer.evaluate(ret0.order);
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = "
                    + scorer.score()
                    + " (Triangle)"
                    + " Elapsed = " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    private MultiSwapRet multiSwap(TeyssierScorer scorer, Node x, Node y, Set<Node> W, int depth) {
        depth = Math.min(depth, W.size());
        List<Node> WW = new ArrayList<>(W);

        int[] indices = new int[WW.size()];
        for (int i = 0; i < WW.size(); i++) indices[i] = scorer.index(WW.get(i));

        MultiSwapRet best = new MultiSwapRet(scorer.getOrder(), scorer.score());
        List<Node> pi = scorer.getOrder();

        ChoiceGenerator gen = new ChoiceGenerator(WW.size(), depth);
        int[] choice;

        while ((choice = gen.next()) != null) {
            int[] choice2 = Arrays.copyOf(choice, choice.length + 2);
            choice2[choice.length] = scorer.index(x);
            choice2[choice.length + 1] = scorer.index(y);

            List<Node> myNodes = new ArrayList<>();
            List<Integer> myIndices = new ArrayList<>();

            int min = Integer.MAX_VALUE;
            int max = Integer.MIN_VALUE;

            for (int i : choice2) {
                Node node = scorer.get(i);
                myNodes.add(node);
                int index = scorer.index(node);
                myIndices.add(index);

                if (index > max) max = index;
                if (index < min) min = index;
            }

            PermutationGenerator permGen = new PermutationGenerator(depth);
            int[] perm;

            while ((perm = permGen.next()) != null) {
                for (int w = 0; w < depth; w++) {
                    pi.set(myIndices.get(perm[w]), myNodes.get(w));
                }

                scorer.evaluate(pi, 0, scorer.size() - 1);

                if (scorer.score() > best.score) {
                    best = new MultiSwapRet(scorer.getOrder(), scorer.score());
                }
            }
        }

        return best;
    }

    private List<Node> gasp(@NotNull TeyssierScorer scorer) {
        Graph newGraph = scorer.getGraph(true);
        Graph oldGraph;

        int _depth = depth == -1 ? 100 : depth;

        for (int k = 0; k < _depth; k++) {
            if (verbose) {
                System.out.println("Round " + (k + 1));
            }

            do {
                oldGraph = newGraph;
                List<NodePair> pairs = scorer.getAdjacencies();
                shuffle(pairs);

                for (NodePair pair : pairs) {
                    scorer.bookmark();

                    double sOld = scorer.score();
                    scorer.tuck(pair.getFirst(), pair.getSecond());
                    double sNew = scorer.score();

                    if (sNew < sOld) {
                        scorer.goToBookmark();
                    }
                }

                newGraph = scorer.getGraph(true);
            } while (!newGraph.equals(oldGraph));
        }

        return scorer.getOrder();
    }


    private List<Node> quickGasp(@NotNull TeyssierScorer scorer) {
        int _depth = depth == -1 ? 100 : depth;

        for (int k = 0; k < _depth; k++) {
            System.out.println("Round " + (k + 1));

            List<NodePair> pairs = scorer.getAdjacencies();
            shuffle(pairs);

            for (NodePair pair : pairs) {
                scorer.bookmark();

                double sOld = scorer.score();
                scorer.tuck(pair.getFirst(), pair.getSecond());
                double sNew = scorer.score();

                if (sNew < sOld) {
                    scorer.goToBookmark();
                }
            }
        }

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

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");

        this.depth = depth;
    }

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }


    public enum Method {BOSS, SP, GASP, QUICK_GASP}

    private static class MultiSwapRet {
        public List<Node> order;
        public double score;

        public MultiSwapRet(List<Node> pi, double score) {
            this.order = pi;
            this.score = score;
        }
    }
}
