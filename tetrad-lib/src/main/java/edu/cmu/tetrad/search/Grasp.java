package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static java.util.Collections.shuffle;


/**
 * Implements the GRASP algorithms, with various execution flags.
 *
 * @author josephramsey
 * @author bryanandrews
 */
public class Grasp {
    private final List<Node> variables;
    private Score score;
    private IndependenceTest test;
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();
    private TeyssierScorer scorer;
    private int depth = 4;
    private int numStarts = 1;
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

    public Grasp(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
        this.useScore = true;
    }

    public Grasp(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
        this.useScore = false;
    }

    public Grasp(@NotNull IndependenceTest test, Score score) {
        this.test = test;
        this.score = score;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();

        if (useScore && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(test, score);
            scorer.setUseScore(true);
            scorer.setUsePearl(usePearl);
        } else {
            scorer = new TeyssierScorer(test, score);
            scorer.setUsePearl(usePearl);

            if (usePearl) {
                scorer.setUseScore(false);
            } else {
                scorer.setUseScore(useScore);
            }
        }

        scorer.setKnowledge(knowledge);
        scorer.setScoreType(scoreType);
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

            List<Node> perm = ses(scorer);
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

    public List<Node> ses(@NotNull TeyssierScorer scorer) {
        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        float sNew = scorer.score();
        float sOld;

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

            if (ordered) {
                sesDfsOrdered(scorer, sOld, depth, 1, ops, new HashSet<>(), new HashSet<>());
            } else {
                sdsDfs(scorer, sOld, depth, 1, ops, new HashSet<>(), new HashSet<>());
            }
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

//    private void sesDfs(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
//                        List<int[]> ops, Set<Set<Node>> branchHistory, Set<Set<Set<Node>>> dfsHistory){
//        for (int[] op : ops) {
//            Node x = scorer.get(op[0]);
//            Node y = scorer.get(op[1]);
//
//            if (!scorer.adjacent(x, y)) continue;
//
//            if (currentDepth > 1 && !scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck on covered edges within DFS
////            if (!scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck covered edges
//
//            Set<Set<Node>> current = new HashSet<>(branchHistory);
//            Set<Node> adj = new HashSet<>();
//            adj.add(x);
//            adj.add(y);
//
//            if (current.contains(adj)) continue;
//            current.add(adj);
//
//            if (op[0] < op[1] && scorer.coveredEdge(x, y)) {
//                if (dfsHistory.contains(current)) continue;
//                dfsHistory.add(current);
//            }
//
//            scorer.bookmark(currentDepth);
//            scorer.moveTo(y, op[0]);
//
//            if (violatesKnowledge(scorer.getOrder())) {
//                scorer.goToBookmark(currentDepth);
//            } else {
//                double sNew = scorer.score();
//
//                if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) { // Uncomment to only DFS on covered edges
////                if (sNew == sOld && currentDepth < depth && op[0] < op[1]) { // Uncomment to only DFS on forward tucks
////                if (sNew == sOld && currentDepth < depth) { // Uncomment to DFS on anything
//                    sesDfs(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
//                    sNew = scorer.score();
//                }
//
//
//                if (sNew <= sOld) {
//                    scorer.goToBookmark(currentDepth);
//                } else {
//                    System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
//                    break; // Uncomment this to break after first improvement
////                    sOld = sNew; // Uncomment this to run rounds to completion
////                    dfsHistory.clear(); // Uncomment this to run rounds to completion
////                    current.clear(); // Uncomment this to run rounds to completion
//                }
////
////                if (sNew <= sOld) {
////                    scorer.goToBookmark(currentDepth);
////                } else {
////                    System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
////                    break; // Comment this out to run rounds to completion
////                }
//            }
//        }
//    }

    private void sdsDfs(@NotNull TeyssierScorer scorer, float sOld, int depth, int currentDepth, List<int[]> ops, Set<Set<Node>> branchHistory,
                        Set<Set<Set<Node>>> dfsHistory) {
        for (int[] op : ops) {
            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;

            if (checkCovering && !scorer.coveredEdge(x, y)) continue;

//            if (currentDepth > 1 && !scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck on covered edges within DFS
//            if (doCovered && !scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck covered edges
//            else if (!doCovered && scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck covered edges

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
                float sNew = scorer.score();

                if (!breakAfterImprovement) {
//                    if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) { // Uncomment to only DFS on covered edges
//                    if (sNew == sOld && currentDepth < depth && op[0] < op[1]) { // Uncomment to only DFS on forward tucks
                    if (sNew == sOld && currentDepth < depth) { // Uncomment to DFS on anything
                        sdsDfs(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
                        sNew = scorer.score();
                    }
                } else {

                    if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) { // Uncomment to only DFS on covered edges
//                    if (sNew == sOld && currentDepth < depth && op[0] < op[1]) { // Uncomment to only DFS on forward tucks
//                    if (sNew == sOld && currentDepth < depth) { // Uncomment to DFS on anything
                        sesDfsOrdered(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
                        sNew = scorer.score();
                    }
                }


                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
                    }

                    if (!breakAfterImprovement) {
                        sOld = sNew; // Uncomment this to run rounds to completion
                        dfsHistory.clear(); // Uncomment this to run rounds to completion
                        current.clear(); // Uncomment this to run rounds to completion
                    } else {
                        break; // Uncomment this to break after first improvement
                    }
                }
//
//                if (sNew <= sOld) {
//                    scorer.goToBookmark(currentDepth);
//                } else {
//                    System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
//                    break; // Comment this out to run rounds to completion
//                }
            }
        }
    }

    private void sesDfsOrdered(@NotNull TeyssierScorer scorer, float sOld, int depth, int currentDepth,
                               List<int[]> ops, Set<Set<Node>> branchHistory, Set<Set<Set<Node>>> dfsHistory) {
        sOld = sdsDfsOrderedVisit(scorer, sOld, depth, currentDepth, ops, branchHistory, dfsHistory,
                true, checkCovering);
        sdsDfsOrderedVisit(scorer, sOld, depth, currentDepth, ops, branchHistory, dfsHistory,
                false, checkCovering);
    }

    private float sdsDfsOrderedVisit(@NotNull TeyssierScorer scorer, float sOld, int depth, int currentDepth, List<int[]> ops, Set<Set<Node>> branchHistory,
                                     Set<Set<Set<Node>>> dfsHistory, boolean doCovered, boolean checkCovered) {
        for (int[] op : ops) {
            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;

            if (checkCovered && !scorer.coveredEdge(x, y)) continue;

//            if (currentDepth > 1 && !scorer.coveredEdge(x, y)) continue; // Uncomment to only tuck on covered edges within DFS
            if (doCovered == !scorer.coveredEdge(x, y)) continue;

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
                float sNew = scorer.score();

//                if (rcgConfig) {
////                    if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) { // Uncomment to only DFS on covered edges
////                    if (sNew == sOld && currentDepth < depth && op[0] < op[1]) { // Uncomment to only DFS on forward tucks
//                    if (sNew == sOld && currentDepth < depth) { // Uncomment to DFS on anything
//                        graspDfsOrdered(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory,
//                                checkCovered);
//                        sNew = scorer.score();
//                    }
//                } else {

                if (sNew == sOld && currentDepth < depth && dfsHistory.contains(current)) { // Uncomment to only DFS on covered edges
//                    if (sNew == sOld && currentDepth < depth && op[0] < op[1]) { // Uncomment to only DFS on forward tucks
//                    if (sNew == sOld && currentDepth < depth) { // Uncomment to DFS on anything
                    sesDfsOrdered(scorer, sNew, depth, currentDepth + 1, ops, current, dfsHistory);
                    sNew = scorer.score();
                    sOld = sNew;
                }
//                }

                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
                    }

                    if (!breakAfterImprovement) {
                        sOld = sNew; // Uncomment this to run rounds to completion
                        dfsHistory.clear(); // Uncomment this to run rounds to completion
                        current.clear(); // Uncomment this to run rounds to completion
                    } else {
                        break; // Uncomment this to break after first improvement
                    }
                }
//
//                if (sNew <= sOld) {
//                    scorer.goToBookmark(currentDepth);
//                } else {
//                    System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", scorer.getNumEdges(), sNew - sOld, current);
//                    break; // Comment this out to run rounds to completion
//                }
            }
        }
        return sOld;
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
}