package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
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
    private int maxPermSize = 4;
    private int gspDepth = -1;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.GrowShrinkMb;
    private TeyssierScorer scorer;
    private boolean doFgesFirst = false;
    private boolean useScore = true;

    public Boss(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
    }

    public Boss(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public Boss(@NotNull IndependenceTest test, Score score) {
        this.test = test;
        this.score = score;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();
        List<Node> _order = null;

        if (useScore && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(test, score);
            scorer.setUseScore(useScore);
            scorer.setParentCalculation(parentCalculation);

            if (doFgesFirst) {
                Fges fges = new Fges(score);
                fges.setKnowledge(knowledge);
                fges.setTurning(true);
                fges.setVerbose(verbose);
                Graph g = fges.search();
                _order = g.getCausalOrdering();
            }
        } else {
            scorer = new TeyssierScorer(test, score);
            scorer.setParentCalculation(parentCalculation);
            scorer.setUseScore(useScore);

            if (score != null && doFgesFirst) {
                Fges fges = new Fges(score);
                fges.setKnowledge(knowledge);
                fges.setTurning(true);
                fges.setVerbose(verbose);
                Graph g = fges.search();
                _order = g.getCausalOrdering();
            }
        }

        scorer.setKnowledge(knowledge);
        scorer.setScoreType(scoreType);

        scorer.setCachingScores(cachingScores);

        List<Node> bestPerm = new ArrayList<>(order);
        double best = NEGATIVE_INFINITY;

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (_order == null) {
                order = new ArrayList<>(order);

                if (!useDataOrder) {
                    shuffle(order);
                }

                makeValidKnowledgeOrder(order);
            }

            scorer.score(order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer, start);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
            } else if (method == Method.QUICK_GASP) {
                perm = quickGasp(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else if (method == Method.ESP) {
                perm = esp(scorer);
            } else if (method == Method.TSP) {
                perm = tsp(scorer);
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
        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        int _maxPermSize = maxPermSize < 0 ? Integer.MAX_VALUE : maxPermSize;

        List<Node> pi = scorer.getOrder();

        pi = betterMutation(scorer, pi, start);
        pi = correctTriMutation(scorer, pi, _maxPermSize, start);

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
                    + " (betterMutation)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
        }

        return scorer.getOrder();
    }

    private List<Node> correctTriMutation(@NotNull TeyssierScorer scorer, List<Node> pi, int depth, long start) {
        DepthChoiceGenerator choiceGenerator = new DepthChoiceGenerator(pi.size(), depth);
        int[] choice;

        while ((choice = choiceGenerator.next()) != null) {
            if (choice.length < 3) continue;

            if (triangleCondition(scorer, choice)) {
                PermutationGenerator permGen = new PermutationGenerator(choice.length);
                int[] perm;

                while ((perm = permGen.next()) != null) {
                    List<Node> pip = subMutation(pi, choice, perm);

                    if (scorer.score(pip) > scorer.score(pi)) {
                        pi = pip;
                    }
                }
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = "
                    + scorer.score()
                    + " (correctTriMutation)"
                    + " Elapsed = " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return pi;
    }

    private boolean triangleCondition(TeyssierScorer scorer, int[] choice) {
        for (int i = 0; i < choice.length; i++) {

            J:
            for (int j = i + 1; j < choice.length; j++) {
                Node x = scorer.get(choice[i]);
                Node y = scorer.get(choice[j]);
                if (!scorer.adjacent(x, y)) continue;

                for (int k = 0; k < choice.length; k++) {
                    if (k == i || k == j) continue;
                    Node z = scorer.get(choice[k]);

                    if (!(scorer.adjacent(z, x) && scorer.adjacent(z, y)))
                        continue J;
                }

                return true;
            }
        }

        return false;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    private List<Node> subMutation(List<Node> pi, int[] choice, int[] perm) {
        List<Node> pi2 = new ArrayList<>(pi);

        int[] choice2 = Arrays.copyOf(choice, choice.length);
        Arrays.sort(choice2);

        for (int i = 0; i < choice.length; i++) {
            pi2.set(choice2[i], pi.get(choice[perm[i]]));
        }

        return pi2;
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

        int _depth = gspDepth == -1 ? 100 : gspDepth;

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
        int _depth = gspDepth == -1 ? 100 : gspDepth;

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
        List<Node> pi = scorer.getOrder();

        Set<OrderedPair<Node>> path = new HashSet<>();

        for (int i = 0; i < 5; i++) {
            pi = espVisit1(scorer, pi, null, path);
//            pi = espVisit2(scorer, pi, null, path);
        }

        return pi;
    }

    public List<Node> espVisit1(@NotNull TeyssierScorer scorer,
                                List<Node> pi, OrderedPair<Node> pair,
                                Set<OrderedPair<Node>> path) {
        if (pair != null && path.contains(pair)) return pi;
        if (pair != null) path.add(pair);

        for (int i = 0; i < pi.size() - 1; i++) {
            if (pair != null) scorer.swap(pi.get(i), pi.get(i + 1));
            OrderedPair<Node> pathPair = new OrderedPair<>(pi.get(i), pi.get(i + 1));
            List<Node> pip = espVisit1(scorer, pi, pathPair, path);

            if (scorer.score(pip) >= scorer.score(pi)) {
                pi = pip;
            }

            if (pair != null) scorer.swap(pi.get(i), pi.get(i + 1));
        }

        if (pair != null) path.remove(pair);
        return pi;
    }

//    public List<Node> espVisit2(@NotNull TeyssierScorer scorer, List<Node> pi, NodePair
//            pair, Set<NodePair> path) {
//        if (pair != null && path.contains(pair)) return pi;
//        if (pair != null) path.add(pair);
//
//        for (int i = 0; i < pi.size() - 1; i++) {
//            List<Node> pip = swap(pi, pi.get(i), pi.get(i + 1));
//
//            if (scorer.score(pip) > scorer.score(pi)) {
//                pi = pip;
//                OrderedPair<Node> pathPair = new OrderedPair<>(pi.get(i), pi.get(i + 1));
//                List<Node> pipp = espVisit1(scorer, pip, pathPair, path);
//
//                if (scorer.score(pipp) == scorer.score(pip)) {
//                    pi = pipp;
//                } else if (scorer.score(pipp) > scorer.score(pip)) {
//                    return pi;
//                }
//            }
//        }
//
//        if (pair != null) path.remove(pair);
//        return pi;
//    }


    public List<Node> tsp(@NotNull TeyssierScorer scorer) {
        return tspVisit(scorer, scorer.getOrder(), null, new HashSet<>());
    }

    public List<Node> tspVisit(@NotNull TeyssierScorer scorer, List<Node> pi, NodePair
            pair, Set<NodePair> path) {
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

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public void setMaxPermSize(int maxPermSize) {
        if (maxPermSize < -1) throw new IllegalArgumentException(
                "Max perm size should be >= 2. This the maximum number of variables that will be " +
                        "permuted in place. If 2, the permutation step will be ignored.");
        this.maxPermSize = maxPermSize;
    }

    public void setGspDepth(int gspDepth) {
        if (gspDepth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.gspDepth = gspDepth;
    }

    public void setDoFgesFirst(boolean doFgesFirst) {
        this.doFgesFirst = doFgesFirst;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }

    public enum Method {BOSS, SP, GASP, QUICK_GASP, ESP, TSP}
}
