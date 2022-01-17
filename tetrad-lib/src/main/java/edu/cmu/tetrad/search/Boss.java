package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.OrderedPair;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Math.min;
import static java.util.Collections.shuffle;


/**
 * Implements various permutation algorithms, including BOSS and GASP.
 *
 * @author josephramsey
 */
public class Boss {
    private final List<Node> variables;
    private long start;
    private Score score;
    private IndependenceTest test;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private Method method = Method.BOSS1;
    private boolean verbose = false;
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();
    private boolean useDataOrder = false;
    private int depth = 4;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.Pearl;
    private TeyssierScorer scorer;
    private boolean useScore = true;
    private int maxPermSize = 4;
    private int numRounds = 50;
    private boolean useTuck = true;

    public Boss(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
        this.useScore = true;
    }

    public Boss(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
        this.useScore = false;
    }

    public Boss(@NotNull IndependenceTest test, Score score) {
        this.test = test;
        this.score = score;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();

        if (useScore && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(test, score);
            scorer.setUseScore(true);
            scorer.setParentCalculation(parentCalculation);
        } else {
            scorer = new TeyssierScorer(test, score);
            scorer.setParentCalculation(parentCalculation);
            scorer.setUseScore(false);
        }

        scorer.setKnowledge(knowledge);
        scorer.setScoreType(scoreType);
        scorer.clearBookmarks();

        scorer.setCachingScores(cachingScores);

        List<Node> bestPerm = new ArrayList<>(order);
        double best = NEGATIVE_INFINITY;

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (!useDataOrder) {
                shuffle(order);
            }

            this.start = System.currentTimeMillis();

            makeValidKnowledgeOrder(order);

            scorer.score(order);

            if (verbose) {
                System.out.println("Using " + method);
            }

            List<Node> perm;

            if (method == Method.BOSS1) {
                perm = boss1(scorer);
            } else if (method == Method.BOSS2) {
                perm = boss2(scorer);
            } else if (method == Method.GRASP) {
                perm = grasp(scorer);
            } else if (method == Method.RCG) {
                perm = rcg2(scorer);
            } else if (method == Method.ESP) {
                perm = esp(scorer);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else if (method == Method.SES) {
                perm = ses(scorer);
            } else if (method == Method.SHG) {
                perm = shuffledGrasp(scorer);
            } else if (method == Method.slowSES) {
                perm = ses(scorer);
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

    public List<Node> boss1(@NotNull TeyssierScorer scorer) {
        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        betterMutation(scorer);
        betterTriMutation(scorer, (maxPermSize < 0 ? Integer.MAX_VALUE : maxPermSize));

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (boss1)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    public List<Node> boss2(@NotNull TeyssierScorer scorer) {
        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        rcg(scorer);
        betterTriMutation(scorer, (maxPermSize < 0 ? Integer.MAX_VALUE : maxPermSize));

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (boss2)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    public List<Node> esp(@NotNull TeyssierScorer scorer) {
        if (depth <= 0) throw new IllegalArgumentException("Form ESP, max depth should be > 0");

        double sOld;
        double sNew = scorer.score();

        do {
            sOld = sNew;
            espDfs(scorer, sOld, (depth < 0 ? 100 : depth), 1);
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (ESP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    public List<Node> gasp(@NotNull TeyssierScorer scorer) {
        if (depth < 0) throw new IllegalArgumentException("Form GRASP, max depth should be >= 0");
        scorer.clearBookmarks();

        double sOld;
        double sNew = scorer.score();

        do {
            sOld = sNew;
            graspDfs(scorer, sOld, (depth < 0 ? Integer.MAX_VALUE : depth), 0, true);
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GASP))"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();

    }

    public List<Node> grasp(@NotNull TeyssierScorer scorer) {
        if (depth < 0) throw new IllegalArgumentException("Form GRASP, max depth should be >= 0");
        scorer.clearBookmarks();

        double sOld;
        double sNew = scorer.score();

        do {
            sOld = sNew;
            graspDfs(scorer, sOld, (depth < 0 ? Integer.MAX_VALUE : depth), 0, false);
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRASP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    public List<Node> shuffledGrasp(@NotNull TeyssierScorer scorer) {
        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        double sOld;
        double sNew = scorer.score();

        do {
            sOld = sNew;
            List<OrderedPair<Node>> edges = scorer.getEdges();
            shuffle(edges);
            shuffledGraspDfs(scorer, sOld, depth, 1, edges, 0);
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRASP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    public List<Node> ses(@NotNull TeyssierScorer scorer) {
        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        double sNew = scorer.score();
        double sOld;

        List<int[]> ops = new ArrayList<>();
        for (int i = 0; i < scorer.size(); i++) {
            for (int j = (useTuck ? i + 1 : 0); j < scorer.size(); j++) {
                if (i != j) {
                    ops.add(new int[]{i, j});
                }
            }
        }

        do {
            sOld = sNew;
            shuffle(ops);
            sesDfs(scorer, sOld, depth, 1, ops, new HashSet<>(), new HashSet<>());
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getOrder();
    }

    private void sesDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
                        List<int[]> ops, Set<Set<Node>> branchHistory, Set<Set<Set<Node>>> dfsHistory) {
        for (int[] op : ops) {
            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;

            if (currentDepth > 1 && !scorer.coveredEdge(x, y)) continue;

            Set<Set<Node>> current = new HashSet<>(branchHistory);
            Set<Node> adj = new HashSet<>();
            adj.add(x);
            adj.add(y);

            if (current.contains(adj)) continue;
            current.add(adj);

            if (op[0] < op[1] && scorer.coveredEdge(x, y)) {
                if (dfsHistory.contains(current)) continue;
                dfsHistory.add(current);
            }

            scorer.bookmark(currentDepth);
            scorer.moveTo(y, op[0]);

            if (violatesKnowledge(scorer.getOrder())) {
                scorer.goToBookmark(currentDepth);
            } else {
                double sNew = scorer.score();

                if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) {
//                if (sNew == sOld && currentDepth < depth && op[0] < op[1]) {
//                if (sNew == sOld && currentDepth < depth) {
                    sesDfs(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
                    sNew = scorer.score();
                }

                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
                    }
                    break;
                }
            }
        }
    }

    // "Random Carnival Game"
    public List<Node> rcg(@NotNull TeyssierScorer scorer) {
        if (numRounds <= 0) throw new IllegalArgumentException("For RCG, #rounds should be > 0");
        scorer.clearBookmarks();
        NumberFormat nf = new DecimalFormat("0.0");

        if (verbose) {
            System.out.println("\nInitial # edges = " + scorer.getNumEdges());
        }

        scorer.bookmark(1);

        int maxRounds = numRounds < 0 ? Integer.MAX_VALUE : numRounds;
        int unimproved = 0;

        for (int r = 1; r <= maxRounds; r++) {
            if (verbose) {
                System.out.println("### Round " + (r));
            }

            List<OrderedPair<Node>> pairs = new ArrayList<>();

            for (Node y : scorer.getOrder()) {
                for (Node x : scorer.getParents(y)) {
                    pairs.add(new OrderedPair<>(x, y));
                }
            }

            shuffle(pairs);

            int numImprovements = 0;
            int numEquals = 0;

            int visited = 0;
            int numPairs = pairs.size();

            for (OrderedPair<Node> pair : pairs) {
                visited++;

                Node x = pair.getFirst();
                Node y = pair.getSecond();
                if (!scorer.adjacent(x, y)) continue;

                double s0 = scorer.score();
                scorer.bookmark(0);

                // 'tuck' operation.
//                scorer.tuck(x, y);
                scorer.moveTo(y, scorer.index(x));

                if (violatesKnowledge(scorer.getOrderShallow())) {
                    scorer.goToBookmark(0);
                    continue;
                }

                double sNew = scorer.score();

                if (sNew < s0) {
                    scorer.goToBookmark(0);
                    continue;
                }

                scorer.bookmark(1);

                if (sNew > s0) {
                    numImprovements++;
                }

                if (verbose) {
                    if (sNew == s0) {
                        numEquals++;
                    }

                    if (sNew > s0) {
                        System.out.println("Round " + (r) + " # improvements = " + numImprovements
                                + " # unimproved = " + numEquals
                                + " # edges = " + scorer.getNumEdges() + " progress this round = "
                                + nf.format(100 * visited / (double) numPairs) + "%");
                    }
                }
            }

            if (numImprovements == 0) {
                unimproved++;
            }

            if (unimproved >= depth) {
                break;
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " #round = " + numRounds
                    + " (RCG)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        scorer.goToBookmark(1);
        return scorer.getOrder();
    }


    public List<Node> rcg2(@NotNull TeyssierScorer scorer) {
        if (numRounds <= 0) throw new IllegalArgumentException("For RCG, #rounds should be > 0");
        scorer.clearBookmarks();
        NumberFormat nf = new DecimalFormat("0.0");

        if (verbose) {
            System.out.println("\nInitial # edges = " + scorer.getNumEdges());
        }

        scorer.bookmark(1);

        int maxRounds = numRounds < 0 ? Integer.MAX_VALUE : numRounds;
        int unimproved = 0;

        for (int r = 1; r <= maxRounds; r++) {
            if (verbose) {
                System.out.println("### Round " + (r));
            }

            List<OrderedPair<Node>> pairs = new ArrayList<>();

            for (Node y : scorer.getOrder()) {
                for (Node x : scorer.getParents(y)) {
                    pairs.add(new OrderedPair<>(x, y));
                }
            }

            shuffle(pairs);

            int numImprovements = 0;
            int numEquals = 0;

            int visited = 0;
            int numPairs = pairs.size();

            for (OrderedPair<Node> pair : pairs) {
                visited++;

                Node x = pair.getFirst();
                Node y = pair.getSecond();
                if (!scorer.adjacent(x, y)) continue;

                double s0 = scorer.score();
                scorer.bookmark(0);

                // 'tuck' operation.
//                scorer.tuck(x, y);
                scorer.moveTo(y, scorer.index(x));

                if (violatesKnowledge(scorer.getOrderShallow())) {
                    scorer.goToBookmark(0);
                    continue;
                }

                double sNew = scorer.score();

                if (sNew < s0) {
                    scorer.goToBookmark(0);
                    continue;
                }

                scorer.bookmark(1);

                if (sNew > s0) {
                    numImprovements++;
                }

                if (verbose) {
                    if (sNew == s0) {
                        numEquals++;
                    }

                    if (sNew > s0) {
                        System.out.println("Round " + (r) + " # improvements = " + numImprovements
                                + " # unimproved = " + numEquals
                                + " # edges = " + scorer.getNumEdges() + " progress this round = "
                                + nf.format(100 * visited / (double) numPairs) + "%");
                    }
                }
            }


            for (OrderedPair<Node> pair : pairs) {
                visited++;

                Node x = pair.getFirst();
                Node y = pair.getSecond();
                if (!scorer.adjacent(x, y)) continue;

                double s0 = scorer.score();
                scorer.bookmark(0);

                // 'tuck' operation.
//                scorer.tuck(x, y);
                scorer.moveTo(x, scorer.index(y));

                if (violatesKnowledge(scorer.getOrderShallow())) {
                    scorer.goToBookmark(0);
                    continue;
                }

                double sNew = scorer.score();

                if (sNew < s0) {
                    scorer.goToBookmark(0);
                    continue;
                }

                scorer.bookmark(1);

                if (sNew > s0) {
                    numImprovements++;
                }

                if (verbose) {
                    if (sNew == s0) {
                        numEquals++;
                    }

                    if (sNew > s0) {
                        System.out.println("Round " + (r) + " # improvements = " + numImprovements
                                + " # unimproved = " + numEquals
                                + " # edges = " + scorer.getNumEdges() + " progress this round = "
                                + nf.format(100 * visited / (double) numPairs) + "%");
                    }
                }
            }

            if (numImprovements == 0) {
                unimproved++;
            }

            if (unimproved >= depth) {
                break;
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " #round = " + numRounds
                    + " (RCG)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        scorer.goToBookmark(1);
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
            scorer.score(p);

            if (scorer.score(p) > maxScore) {
                maxScore = scorer.score(p);
                maxP = p;
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (SP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
        }

        return maxP;
    }

    public void betterMutation(@NotNull TeyssierScorer scorer) {
        List<Node> pi = scorer.getOrder();
        double s;
        double sp = scorer.score(pi);
        scorer.bookmark();

        do {
            s = sp;

            for (Node x : scorer.getOrder()) {
                sp = NEGATIVE_INFINITY;
                int index = scorer.index(x);

                for (int i = index; i >= 0; i--) {
                    scorer.moveTo(x, i);

                    if (scorer.score() > sp) {
                        if (!violatesKnowledge(scorer.getOrder())) {
                            sp = scorer.score();
                            scorer.bookmark();
                        }
                    }
                }

                scorer.goToBookmark();
                scorer.bookmark();

                for (int i = index; i < scorer.size(); i++) {
                    scorer.moveTo(x, i);

                    if (scorer.score() > sp) {
                        if (!violatesKnowledge(scorer.getOrder())) {
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
    }

    private void betterTriMutation(@NotNull TeyssierScorer scorer, int maxPermSize) {
        if (maxPermSize < 0) throw new IllegalArgumentException("maxPermSize should be >= 0");

        List<Node> pi = scorer.getOrder();

        W:
        while (true) {
            for (int length = 2; length <= maxPermSize; length++) {
                boolean triangleConditionNeverTrue = true;

                ChoiceGenerator choiceGenerator = new ChoiceGenerator(pi.size(), min(length, pi.size()));
                int[] choice;

                while ((choice = choiceGenerator.next()) != null) {
                    if (triangleCondition(scorer, choice)) {
                        triangleConditionNeverTrue = false;
                        PermutationGenerator permGen = new PermutationGenerator(choice.length);
                        int[] perm;

                        while ((perm = permGen.next()) != null) {
                            List<Node> pip = subMutation(pi, choice, perm);

                            if (scorer.score(pip) > scorer.score(pi)) {
                                pi = pip;
                                continue W;
                            }
                        }
                    }
                }

                if (triangleConditionNeverTrue) {
                    if (verbose) System.out.println("Breaking at " + length);
                    break;
                }
            }

            break;
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " max perm size = " + maxPermSize
                    + " Score = "
                    + scorer.score()
                    + " (betterTriMutationA)"
                    + " Elapsed = " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }
    }

    @NotNull
    public Graph getGraph(boolean cpDag) {
        if (scorer == null) throw new IllegalArgumentException("Please run algorithm first.");
        Graph graph = scorer.getGraph(cpDag);
        graph.addAttribute("# edges", graph.getNumEdges());
        return graph;
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
                "Max perm size should be >= -1. This the maximum number of variables that will be " +
                        "permuted in place. If 2, the permutation step will be ignored. -1 = unlimited.");
        this.maxPermSize = maxPermSize;
    }

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.depth = depth;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }

    private boolean violatesKnowledge(List<Node> order) {
        if (!knowledge.isEmpty()) {
            for (int i = 0; i < order.size(); i++) {
                for (int j = i + 1; j < order.size(); j++) {
                    if (knowledge.isForbidden(order.get(i).getName(), order.get(j).getName())) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private void espDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth) {
        for (int i = 0; i < scorer.size() - 1; i++) {
            List<Node> pi = scorer.getOrder();
            scorer.swap(scorer.get(i), scorer.get(i + 1));

            if (violatesKnowledge(scorer.getOrder())) {
                scorer.score(pi);
                continue;
            }

            double sNew = scorer.score();

            if (sNew == sOld && currentDepth < depth) {
                espDfs(scorer, sNew, depth, currentDepth + 1);
                sNew = scorer.score();
            }

            if (sNew <= sOld) {
                scorer.score(pi);
            } else {
                break;
            }
        }
    }

    private void graspDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
                          boolean checkCovering) {
        for (OrderedPair<Node> adj : scorer.getEdges()) {
            Node x = adj.getFirst();
            Node y = adj.getSecond();
            if (checkCovering && !scorer.coveredEdge(x, y)) continue;
            scorer.bookmark(currentDepth);
            scorer.tuck(x, y);
//            System.out.println(scorer.getOrder() + " score = " + scorer.getNumEdges());

            if (violatesKnowledge(scorer.getOrder())) {
                scorer.goToBookmark(currentDepth);
                continue;
            }

            double sNew = scorer.score();

            if (sNew == sOld && currentDepth < depth) {
//                System.out.println("equals " + scorer.getOrder() + " sNew = " + sNew);
                graspDfs(scorer, sNew, depth, currentDepth + 1, checkCovering);
                sNew = scorer.score();
            }

            if (sNew <= sOld) {
                scorer.goToBookmark(currentDepth);
            } else {
                break;
            }
        }
    }

    private void shuffledGraspDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
                                  List<OrderedPair<Node>> edges, int i) {
        double sNew = sOld;

        while (i < edges.size() && sNew <= sOld) {
            OrderedPair<Node> adj = edges.get(i++);
            Node x = adj.getFirst();
            Node y = adj.getSecond();
            scorer.bookmark(currentDepth);
            scorer.tuck(x, y);

            if (violatesKnowledge(scorer.getOrder())) {
                scorer.goToBookmark(currentDepth);
            } else {
                sNew = scorer.score();

                if (sNew == sOld && currentDepth < depth) {
                    shuffledGraspDfs(scorer, sNew, depth, currentDepth + 1, edges, i);
                    sNew = scorer.score();
                }

                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                }
            }
        }
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


    private List<Node> subMutation(List<Node> pi, int[] choice, int[] perm) {
        List<Node> pi2 = new ArrayList<>(pi);

        int[] choice2 = Arrays.copyOf(choice, choice.length);
        Arrays.sort(choice2);

        for (int i = 0; i < choice.length; i++) {
            pi2.set(choice2[i], pi.get(choice[perm[i]]));
        }

        return pi2;
    }

    public void setNumRounds(int numRounds) {
        this.numRounds = numRounds;
    }



    public void setUseTuck(boolean useTuck) {
        this.useTuck = useTuck;
    }

    public enum Method {GASP, ESP, GRASP, RCG, BOSS1, BOSS2, SP, SES, SHG, slowSES}
}