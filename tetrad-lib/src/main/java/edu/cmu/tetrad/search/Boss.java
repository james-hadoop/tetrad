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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Math.min;
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
        double best = NEGATIVE_INFINITY;// scorer.score();

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (!useDataOrder) {
                shuffle(order);
            }

            List<Node> _order = new ArrayList<>(order);
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
        int depth = this.depth == -1 ? 100 : this.depth;

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        List<Node> pi, pip = scorer.getOrder();

//        do {
        pi = pip;
        pi = betterMutation(scorer, pi, start);
        pip = betterTriMutation(scorer, pi, depth, start);
        pip = betterMutation(scorer, pip, start);
//        } while (scorer.evaluate(pip) > scorer.evaluate(pi));

        return pip;
    }

    private List<Node> betterMutation(@NotNull TeyssierScorer scorer, List<Node> pi, long start) {
        double s;
        double sp = scorer.evaluate(pi);
        scorer.bookmark();

        do {
            s = sp;

            for (Node x : scorer.getOrder()) {
                sp = NEGATIVE_INFINITY;

                for (int i = scorer.size() - 1; i >= 0; i--) {
                    scorer.moveTo(x, i);

                    if (scorer.score() > sp) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            sp = scorer.score();
                            scorer.bookmark();
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (sp > s);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (Single Moves)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
        }

        return scorer.getOrder();
    }

    private List<Node> betterTriMutation(@NotNull TeyssierScorer scorer, List<Node> pi, int depth, long start) {
        List<NodePair> pairs = scorer.getAdjacencies();
        List<Node> pip = pi;

        for (NodePair pair : pairs) {
            Node x = pair.getFirst();
            Node y = pair.getSecond();

            double score = scorer.evaluate(pip);

            Set<Node> W = new HashSet<>();

            for (Node z : pip) {
                if (scorer.triangle(x, y, z)) {
                    W.add(z);
                }
            }

            List<Node> pipp = multiSwap(scorer, pi, x, y, W, depth);

            if (scorer.evaluate(pipp) > score) {
                pip = pipp;
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = "
                    + scorer.score()
                    + " (Triangle)"
                    + " Elapsed = " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return pip;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    private List<Node> multiSwap(TeyssierScorer scorer, List<Node> pi, Node x, Node y, Set<Node> W, int depth) {
        List<Node> WW = new ArrayList<>(W);
        WW.add(x);
        WW.add(y);

        depth = min(depth, WW.size());

        List<Node> pipp = pi;

        ChoiceGenerator choiceGenerator = new ChoiceGenerator(WW.size(), depth);
        int[] choice;

        while ((choice = choiceGenerator.next()) != null) {
            PermutationGenerator permGen = new PermutationGenerator(depth);
            int[] perm;

            while ((perm = permGen.next()) != null) {
                List<Node> delta = new ArrayList<>();
                for (int j : perm) delta.add(WW.get(choice[perm[j]]));

                List<Node> pip = subMutation(pi, delta);

                if (scorer.evaluate(pip) > scorer.evaluate(pipp)) {
                    pipp = pip;
                }
            }
        }

        return pipp;
    }

    private List<Node> subMutation(List<Node> pi, List<Node> delta) {
        List<Node> pi2 = new ArrayList<>(pi);
        int i = -1;
        for (int j = 0; j < pi.size(); j++) {

            if (delta.contains(pi.get(j))) {
                ++i;
                pi2.set(j, delta.get(i));
            }

            if (i == pi.size() - 1) {
                return pi2;
            }
        }

        return pi;
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
}
