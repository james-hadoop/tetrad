package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
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
            scorer.score(_order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = esp(scorer);//, start);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
            } else if (method == Method.QUICK_GASP) {
                perm = quickGasp(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else {
                throw new IllegalArgumentException("Unrecognized method: " + method);
            }

            scorer.score(perm);

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

        List<Node> pi = scorer.getOrder();
        pi = betterMutation(scorer, pi, start);
        pi = betterTriMutation(scorer, pi, depth, start);
        pi = betterMutation(scorer, pi, start);

        return pi;
    }

    private List<Node> betterMutation(@NotNull TeyssierScorer scorer, List<Node> pi, long start) {
        double s;
        double sp = scorer.score(pi);
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

            double score = scorer.score(pip);

            Set<Node> W = new HashSet<>();

            for (Node z : pip) {
                if (scorer.triangle(x, y, z)) {
                    W.add(z);
                }
            }

            List<Node> pipp = multiSwap(scorer, pi, x, y, W, depth);

            if (scorer.score(pipp) > score) {
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

                if (scorer.score(pip) > scorer.score(pipp)) {
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

    public List<Node> tuck(List<Node> pi, Node x, Node y) {
        List<Node> pip = new ArrayList<>(pi);
        pip.remove(y);
        pip.add(pi.indexOf(x), y);
        return pip;
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
        double maxScore = NEGATIVE_INFINITY;
        List<Node> maxP = null;

        List<Node> variables = scorer.getOrder();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);
            if (scorer.score(p) > maxScore) {
                maxScore = scorer.score(p);
                maxP = p;
            }
        }

        return maxP;
    }

    public List<Node> swap(List<Node> pi, Node m, Node n) {
        List<Node> pip = new ArrayList<>(pi);

        int i = pip.indexOf(n);
        int j = pip.indexOf(m);

        pip.set(i, m);
        pip.set(j, n);

        return pip;
    }

    public List<Node> esp(@NotNull TeyssierScorer scorer) {
        List<Node> pi1, pi2 = scorer.getOrder(), pi3;

        Set<NodePair> path = new HashSet<>();

//        for (int k = 0; k < 5; k++) {
//            pi1 = pi2;
        pi1 = espVisit1(scorer, pi2, null, path);
        pi2 = espVisit2(scorer, pi1, path);
//        }
//        while (scorer.score(pi2) > scorer.score(pi1));

        return pi2;
    }

    public List<Node> espVisit1(@NotNull TeyssierScorer scorer, List<Node> pi, NodePair pair, Set<NodePair> path) {
        if (pair != null && path.contains(pair)) return pi;
        if (pair != null) path.add(pair);

        for (int i = 0; i < pi.size() - 1; i++) {
            List<Node> pip = swap(pi, pi.get(i), pi.get(i + 1));

            if (scorer.score(pip) > scorer.score(pi)) {
                List<Node> pipp = espVisit1(scorer, pip, new NodePair(pi.get(i), pi.get(i + 1)), path);

                if (scorer.score(pipp) > scorer.score(pip)) {
                    pi = pipp;
                }
            }
        }

        if (pair != null) path.remove(pair);
        return pi;
    }

    public List<Node> espVisit2(@NotNull TeyssierScorer scorer, List<Node> pi, Set<NodePair> path) {
        for (int i = 0; i < pi.size() - 1; i++) {
            Node X = pi.get(i);
            Node Y = pi.get(i + 1);
            List<Node> pip = swap(pi, X, Y);

            if (scorer.score(pip) == scorer.score(pi)) {
                NodePair pair = new NodePair(X, Y);
                if (path.contains(pair)) continue;
                path.add(pair);
                List<Node> pipp = espVisit2(scorer, pip, path);
                path.remove(pair);

                if (scorer.score(pipp) == scorer.score(pip)) {
                    pi = pipp;
                }
            }
        }

        return pi;
    }

    public List<Node> tsp(@NotNull TeyssierScorer scorer) {
        return tspVisit(scorer, scorer.getOrder(), null, new HashSet<>());
    }

    public List<Node> tspVisit(@NotNull TeyssierScorer scorer, List<Node> pi, NodePair pair, Set<NodePair> path) {
        if (pair != null && path.contains(pair)) return pi;
        if (pair != null) path.add(pair);

        List<Node> pip = new ArrayList<>(pi);

        for (Node Y : pi) {
            scorer.score(pip);
            for (Node X : scorer.getParents(Y)) {
                if (covered(scorer, X, Y)) {
                    List<Node> pipp = tuck(pip, X, Y);

                    if (scorer.score(pipp) >= scorer.score(pip)) {
                        pip = pipp;

                        List<Node> pippp = tspVisit(scorer, pipp, new NodePair(X, Y), path);

                        if (scorer.score(pippp) >= scorer.score(pip)) {
                            pip = pippp;
                        }
                    }
                }
            }
        }

        if (pair != null) path.remove(pair);
        return pip;
    }

    private boolean covered(TeyssierScorer scorer, Node x, Node y) {
        Set<Node> px = scorer.getParents(x);
        Set<Node> py = scorer.getParents(y);
        px.remove(y);
        return px.equals(py);
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
