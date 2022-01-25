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
 * @author bryanandrews
 * @author josephramsey
 */
public class Grasp {
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

    // other params
    private int depth = 4;
    private int numStarts = 1;

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
                grasp(scorer);
                setCheckCovering(false);
                grasp(scorer);
                setUseTuck(false);
                perm = grasp(scorer);
                setUseTuck(_useTuck);
                setCheckCovering(_checkCovering);
            } else {
                perm = grasp(scorer);
            }

            scorer.score(perm);

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
//        if (true) {
//            return grasp2(scorer);
//        }

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
//                shuffle(ops);

                if (ordered) {
                    ops.sort((o1, o2) -> {
                        Node x1 = scorer.get(o1[0]);
                        Node y1 = scorer.get(o1[1]);
                        Node x2 = scorer.get(o2[0]);
                        Node y2 = scorer.get(o2[1]);

                        boolean covered1 = scorer.coveredEdge(x1, y1);
                        boolean covered2 = scorer.coveredEdge(x2, y2);

                        if (covered1 && !covered2) return 1;
                        else if (covered2 && !covered1) return -1;
                        else return 0;
                    });
                }

                graspDfs2(scorer, sOld, eOld, itr_depth, 1, ops, new HashSet<>());
                sNew = scorer.score();
            } while (sNew > sOld);
        } while(itr_depth++ < depth);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

    boolean useEdgeRecursion = true;

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

                if (sNew <= sOld) {
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
}