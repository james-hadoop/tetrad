package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Math.min;
import static java.util.Collections.shuffle;


/**
 * Implements the GRASP algorithms, with various execution flags.
 *
 * @author bryanandrews
 * @author josephramsey
 */
public class Grasp3 {
    private final List<Node> variables;
    private Score score;
    private IndependenceTest test;
    private IKnowledge knowledge = new Knowledge2();
    private TeyssierScorer scorer;
    private long start;

    // flags
    private boolean checkCovering = false;
    private boolean useTuck = false;
    private boolean breakAfterImprovement = true;
    private boolean ordered = true;
    private boolean useScore = true;
    private boolean usePearl = true;
    private boolean verbose = false;
    private boolean cachingScores = true;
    private boolean useDataOrder = false;
    private boolean useEdgeRecursion;

    {
        useEdgeRecursion = false;
    }

    // other params
    private int depth = 4;
    private int numStarts = 1;

    public Grasp3(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
        this.useScore = true;
    }

    public Grasp3(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
        this.useScore = false;
    }

    public Grasp3(@NotNull IndependenceTest test, Score score) {
        this.test = test;
        this.score = score;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();

        scorer = new TeyssierScorer(test, score);
        scorer.setUsePearl(usePearl);

        if (usePearl) {
            scorer.setUseScore(false);
        } else {
            scorer.setUseScore(useScore && !(score instanceof GraphScore));
        }

        scorer.setKnowledge(knowledge);
        scorer.clearBookmarks();

        scorer.setCachingScores(cachingScores);

        List<Node> bestPerm = new ArrayList<>(order);
        float best = Float.NEGATIVE_INFINITY;

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (!useDataOrder) {
                shuffle(order);
            }

            this.start = System.currentTimeMillis();

            makeValidKnowledgeOrder(order);

            scorer.score(order);

            List<Node> perm;

            if (ordered) {
                boolean _checkCovering = checkCovering;
                boolean _useTuck = useTuck;
                setCheckCovering(true);
                setUseTuck(true);
                grasp2(scorer);
                setCheckCovering(false);
                grasp2(scorer);
                setUseTuck(false);
                perm = grasp2(scorer);
                setUseTuck(_useTuck);
                setCheckCovering(_checkCovering);
            } else {
                perm = grasp2(scorer);
            }

//            perm = betterMutation(scorer);

            scorer.score(perm);

//            betterTriMutation(scorer, 5);

            if (scorer.score() > best) {
                best = scorer.score();
                bestPerm = perm;
            }
        }

        long stop = System.currentTimeMillis();

        if (verbose) {
            System.out.println("Final order = " + scorer.getPi());
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

    double sNew = Double.NaN;
    double sOld = Double.NaN;

    public List<Node> grasp(@NotNull TeyssierScorer scorer) {

        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        if (Double.isNaN(sNew)) {
            sNew = scorer.score();
        }

        List<int[]> ops = new ArrayList<>();
        for (int i = 0; i < scorer.size(); i++) {
            for (int j = (/*useTuck ? i + 1 :*/ 0); j < scorer.size(); j++) {
                if (i != j) {
                    ops.add(new int[]{i, j});
                }
            }
        }

        do {
            sOld = sNew;
            shuffle(ops);
            graspDfs(scorer, sOld, depth, 1, ops, new HashSet<>(), new HashSet<>());
            sNew = scorer.score();
        } while (sNew > sOld);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

    private void graspDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
                          List<int[]> ops, Set<Set<Node>> branchHistory, Set<Set<Set<Node>>> dfsHistory) {
        for (int[] op : ops) {
            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;
            if (checkCovering && !scorer.coveredEdge(x, y)) continue;

//            if (currentDepth > 1 && !scorer.coveredEdge(x, y))
//                continue; // Uncomment to only tuck on covered edges within DFS

            Set<Set<Node>> current = new HashSet<>(branchHistory);
            Set<Node> adj = new HashSet<>();
            adj.add(x);
            adj.add(y);

            if (current.contains(adj))
                continue; // Uncomment to prevent tucks between variables that have already been tucked
            current.add(adj);

            if (op[0] < op[1] && scorer.coveredEdge(x, y)) {
                if (dfsHistory.contains(current))
                    continue; // Uncomment to prevent sets of tucks that have already been performed (possibly in a different order)
                dfsHistory.add(current);
            }

            scorer.bookmark(currentDepth);
            scorer.moveTo(y, op[0]);

            if (violatesKnowledge(scorer.getPi())) {
                scorer.goToBookmark(currentDepth);
            } else {
                sNew = scorer.score();

                if (sNew == sOld && currentDepth < depth) {
                    graspDfs(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
                    sNew = scorer.score();
                }

                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
                    }

                    if (breakAfterImprovement) {
                        break;
                    } else {
                        sOld = sNew;
                        dfsHistory.clear();
                        current.clear();
                    }
                }
            }
        }
    }

    public List<Node> grasp2(@NotNull TeyssierScorer scorer) {
        if (true) return grasp(scorer);

        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        if (Double.isNaN(sNew)) {
            sNew = scorer.score();
        }

        int eNew = scorer.getNumEdges();
        int eOld;

        List<int[]> ops = new ArrayList<>();
        for (int i = 0; i < scorer.size(); i++) {
            for (int j = (useTuck ? i + 1 : 0); j < scorer.size(); j++) {
                if (i != j) {
                    ops.add(new int[]{i, j});
                }
            }
        }

        int itr_depth = 1;

        do {
//            System.out.printf("\nIterative depth: %d\n\n", itr_depth);
            do {
                sOld = sNew;
                eOld = eNew;
                shuffle(ops);

//                if (ordered) {
//                    ops.sort((o1, o2) -> {
//                        Node x1 = scorer.get(o1[0]);
//                        Node y1 = scorer.get(o1[1]);
//                        Node x2 = scorer.get(o2[0]);
//                        Node y2 = scorer.get(o2[1]);
//
//                        boolean covered1 = scorer.coveredEdge(x1, y1);
//                        boolean covered2 = scorer.coveredEdge(x2, y2);
//
//                        if (covered1 && !covered2) return 1;
//                        else if (covered2 && !covered1) return -1;
//                        else return 0;
//                    });
//                }

                graspDfs2(scorer, sOld, eOld, itr_depth, 1, ops, new HashSet<>());
                sNew = scorer.score();
            } while (sNew > sOld);
        } while(itr_depth++ < depth);

//        List<Node> pi2 = scorer.getPi();
//        reverse(pi2);
//        scorer.score(pi2);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRASP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

    private void graspDfs2(@NotNull TeyssierScorer scorer, double sOld, int eOld, int depth, int currentDepth,
                           List<int[]> ops, Set<Set<Node>> dfsHistory) {
        for (int[] op : ops) {
            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;
            if (checkCovering && !scorer.coveredEdge(x, y)) continue;
            if (!useEdgeRecursion && currentDepth > 1 && !scorer.coveredEdge(x, y)) continue;

            Set<Node> tuck = new HashSet<>();
            tuck.add(x);
            tuck.add(y);

            if (dfsHistory.contains(tuck)) continue;

            scorer.bookmark(currentDepth);
            scorer.moveTo(y, op[0]);

            if (violatesKnowledge(scorer.getPi())) {
                scorer.goToBookmark(currentDepth);
            } else {
                sNew = scorer.score();
                int eNew = scorer.getNumEdges();

                if ((!useEdgeRecursion ? (sNew == sOld) : (sNew <= sOld && eNew <= eOld)) && currentDepth < depth) {
                    dfsHistory.add(tuck);
                    graspDfs2(scorer, sOld, eNew, depth, currentDepth + 1, ops, dfsHistory);
                    dfsHistory.remove(tuck);
                    sNew = scorer.score();
                    eNew = scorer.getNumEdges();
                }

                if (sNew < sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        dfsHistory.add(tuck);
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", eNew, sNew - sOld, dfsHistory);
                        dfsHistory.remove(tuck);
                    }

                    if (breakAfterImprovement) {
                        break;
                    } else {
                        sOld = sNew;
                        eOld = eNew;
                    }
                }
            }
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

    public List<Node> getVariables() {
        return this.variables;
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

    public void setUseTuck(boolean useTuck) {
        this.useTuck = useTuck;
    }

    public void setBreakAfterImprovement(boolean breakAfterImprovement) {
        this.breakAfterImprovement = breakAfterImprovement;
    }

    public void setOrdered(boolean ordered) {
        this.ordered = ordered;
    }

    public void setUsePearl(boolean usePearl) {
        this.usePearl = usePearl;
    }

    public void setCheckCovering(boolean checkCovering) {
        this.checkCovering = checkCovering;
    }

    public void setUseEdgeRecursion(boolean useEdgeRecursion) {
        this.useEdgeRecursion = useEdgeRecursion;
    }

    private void betterTriMutation(@NotNull TeyssierScorer scorer, int maxPermSize) {
        if (maxPermSize < 0) throw new IllegalArgumentException("maxPermSize should be >= 0");

        List<Node> pi = scorer.getPi();

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

    private List<Node> subMutation(List<Node> pi, int[] choice, int[] perm) {
        List<Node> pi2 = new ArrayList<>(pi);

        int[] choice2 = Arrays.copyOf(choice, choice.length);
        Arrays.sort(choice2);

        for (int i = 0; i < choice.length; i++) {
            pi2.set(choice2[i], pi.get(choice[perm[i]]));
        }

        return pi2;
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

    public List<Node> betterMutation(@NotNull TeyssierScorer scorer) {
        List<Node> pi = scorer.getPi();
        double s;
        double sp = scorer.score(pi);

        System.out.println("Original score = " + sp + " pi = " + pi + " graph = " + scorer.getGraph(false));
        scorer.bookmark(0);
        scorer.bookmark(1);

//        do {
            s = sp;

            for (Node x : scorer.getPi()) {
                scorer.bookmark(0);

                sp = NEGATIVE_INFINITY;
                int index = scorer.index(x);

                for (int i = index; i >= 0; i--) {
//                    scorer.moveTo(x, i);
                    scorer.tuck(x, scorer.get(i));

                    if (scorer.score() >= sp) {
                        if (!violatesKnowledge(scorer.getPi())) {
                            sp = scorer.score();
                            scorer.bookmark(1);

                            System.out.println("Tucking " + x + " to index " + i + " score = " + sp + " " + scorer.getPi());
                        }
                    }
                }

                scorer.goToBookmark(0);
                scorer.bookmark(0);

                for (int i = index; i < scorer.size(); i++) {
//                    scorer.moveTo(x, i);
                    scorer.tuck(x, scorer.get(i));

                    if (scorer.score() >= sp) {
                        if (!violatesKnowledge(scorer.getPi())) {
                            sp = scorer.score();
//                            scorer.bookmark(1);

                            System.out.println("Tucking " + x + " to index " + i + " score = " + sp + " " + scorer.getPi());
                        }
                    }
                }

                scorer.goToBookmark(0);
            }
//        } while (sp >= s);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (betterMutation)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
        }

        return scorer.getPi();
    }


}