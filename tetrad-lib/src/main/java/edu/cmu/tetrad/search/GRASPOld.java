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

import static java.lang.Float.NEGATIVE_INFINITY;
import static java.util.Collections.shuffle;


/**
 * Implements the GRASP permutation algorithms.
 *
 * @author josephramsey
 */
public class GRASPOld {
    private final List<Node> variables;
    private final Score score;
    private long start;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private boolean verbose = false;
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();
    private boolean useDataOrder = false;
    private int depth = 4;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.Pearl;
    private TeyssierScorer scorer;
    private boolean useTuck = false;

    public GRASPOld(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();

        scorer = new TeyssierScorer(null, score);
        scorer.setUseScore(true);
        scorer.setUsePearl(false);

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

            List<Node> perm;
            perm = ses(scorer);
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

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.depth = depth;
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
//        this.useTuck = useTuck;
    }
}