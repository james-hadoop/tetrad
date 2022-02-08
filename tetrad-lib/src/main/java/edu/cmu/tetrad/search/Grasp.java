package edu.cmu.tetrad.search;

import com.google.common.collect.Sets;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.util.Collections.shuffle;


/**
 * Implements the GRASP algorithms, with various execution flags.
 *
 * @author bryanandrews
 * @author josephramsey
 */
public class Grasp {
    private final List<Node> variables;
    private double sNew = Double.NaN;
    private Score score;
    private IndependenceTest test;
    private IKnowledge knowledge = new Knowledge2();
    private TeyssierScorer scorer;
    private long start;

    // flags
    private boolean checkCovering = false;
    private boolean useForwardTuckOnly = false;
    private boolean breakAfterImprovement = true;
    private boolean ordered = false;
    private boolean useScore = true;
    private boolean usePearl = false;
    private boolean verbose = false;
    private boolean cachingScores = true;
    private boolean useDataOrder = false;
    private boolean graspAlg = false;

    // other params
    private int depth = 4;
    private int uncoveredDepth = 0;
    private int numStarts = 1;
    private int timeout = -1;

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
        start = System.currentTimeMillis();

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

        scorer.score(order);

        for (int r = 0; r < (useDataOrder ? 1 : numStarts); r++) {
            if (timeOut()) break;

            if (!useDataOrder) {
                shuffle(order);
            }

            makeValidKnowledgeOrder(order);

            scorer.score(order);

            List<Node> perm;

            if (graspAlg) {
                if (ordered) {
                    boolean _checkCovering = checkCovering;
                    boolean _useTuck = useForwardTuckOnly;
                    setCheckCovering(true);
                    setUseForwardTuckOnly(true);
                    grasp1(scorer);
                    setCheckCovering(false);
                    grasp1(scorer);
                    setUseForwardTuckOnly(false);
                    perm = grasp1(scorer);
                    setUseForwardTuckOnly(_useTuck);
                    setCheckCovering(_checkCovering);
                } else {
                    perm = grasp1(scorer);
                }
            } else {
                perm = grasp4(scorer);
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


    public List<Node> grasp1(@NotNull TeyssierScorer scorer) {
        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        if (Double.isNaN(sNew)) {
            sNew = scorer.score();
        }

        List<int[]> ops = new ArrayList<>();
        for (int i = 0; i < scorer.size(); i++) {
            for (int j = (useForwardTuckOnly ? i + 1 : 0); j < scorer.size(); j++) {
                if (i != j) {
                    ops.add(new int[]{i, j});
                }
            }
        }

        double sOld;

        do {
            if (timeOut()) break;

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
            if (timeOut()) break;

            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;
            if (checkCovering && !scorer.coveredEdge(x, y)) continue;

            if (currentDepth > 1 && !scorer.coveredEdge(x, y))
                continue; // Uncomment to only tuck on covered edges within DFS

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
        int overallDepth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        int uncoveredDepth = this.uncoveredDepth < 0 ? Integer.MAX_VALUE : this.uncoveredDepth;
        scorer.clearBookmarks();

        double sNew = scorer.score();
        double sOld;

        int eNew = scorer.getNumEdges();
        int eOld;

        List<int[]> ops = new ArrayList<>();
        for (int i = 0; i < scorer.size(); i++) {
            for (int j = (useForwardTuckOnly ? i + 1 : 0); j < scorer.size(); j++) {
                if (i != j) {
                    ops.add(new int[]{i, j});
                }
            }
        }

        int itrUncoveredDepth = 0;

        do {
            if (timeOut()) break;

            do {
                if (timeOut()) break;

                sOld = sNew;
                eOld = eNew;
                shuffle(ops);
                graspDfs2(scorer, sOld, eOld, new int[]{overallDepth, itrUncoveredDepth}, 1,
                        ops, new HashSet<>(), scorer.getSkeleton(), new HashSet<>());
                sNew = scorer.score();
            } while (sNew > sOld);
        } while (itrUncoveredDepth++ < uncoveredDepth);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

    private void graspDfs2(@NotNull TeyssierScorer scorer, double sOld, int eOld, int[] depth, int currentDepth,
                           List<int[]> ops, Set<Set<Node>> tuckHistory, Set<Set<Node>> skeleton,
                           Set<Set<Set<Node>>> skeletonHistory) {

        boolean forbidUncoveredTucks = currentDepth > depth[1];
        boolean allowUncoveredTuckRecursion = currentDepth < depth[1];
        boolean allowRecursion = currentDepth < depth[0];

        for (int[] op : ops) {
            if (timeOut()) break;

            Node x = scorer.get(op[0]);
            Node y = scorer.get(op[1]);

            if (!scorer.adjacent(x, y)) continue;
            boolean coveredTuck = scorer.coveredEdge(x, y) && op[0] < op[1];
            if (!coveredTuck && forbidUncoveredTucks) continue;

            Set<Node> tuck = new HashSet<>();
            tuck.add(x);
            tuck.add(y);

            if (tuckHistory.contains(tuck)) continue;

            scorer.bookmark(currentDepth);
            scorer.moveTo(y, op[0]);

            Set<Set<Node>> symmDiff = Sets.symmetricDifference(skeleton, scorer.getSkeleton());
            if (!coveredTuck && skeletonHistory.contains(symmDiff)) {
                scorer.goToBookmark(currentDepth);
                continue;
            }
            skeletonHistory.add(symmDiff);

            if (violatesKnowledge(scorer.getPi())) {
                scorer.goToBookmark(currentDepth);
            } else {
                double sNew = scorer.score();
                int eNew = scorer.getNumEdges();

                if (allowRecursion && (allowUncoveredTuckRecursion ? eNew <= eOld && sNew <= sOld : coveredTuck && sNew == sOld)) {
                    tuckHistory.add(tuck);
                    graspDfs2(scorer, sOld, eNew, depth, currentDepth + 1, ops, tuckHistory, skeleton, skeletonHistory
                    );

                    tuckHistory.remove(tuck);
                    sNew = scorer.score();
                    eNew = scorer.getNumEdges();
                }

                if (sNew <= sOld) {
                    scorer.goToBookmark(currentDepth);
                } else {
                    if (verbose) {
                        tuckHistory.add(tuck);
                        System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s\n", eNew, sNew - sOld, tuckHistory);
                        tuckHistory.remove(tuck);
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

    public List<Node> grasp4(@NotNull TeyssierScorer scorer) {
        int depth = this.depth < 1 ? Integer.MAX_VALUE : this.depth;
        scorer.clearBookmarks();

        double sNew = scorer.score();
        double sOld;

        do {
            sOld = sNew;
            graspDfs4(scorer, sOld, depth, 1, new HashSet<>(), new HashSet<>());
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

    private void graspDfs4(@NotNull TeyssierScorer scorer, double sOld, int depth, int currentDepth,
                           Set<Set<Node>> tucks, Set<Set<Set<Node>>> dfsHistory) {
        for(Node y : scorer.getShuffledVariables()) {
            Set<Node> ancestors = scorer.getAncestors(y);
            List<Node> parents = new ArrayList<>(scorer.getParents(y));
            shuffle(parents);

            for(Node x : parents) {

                boolean covered = scorer.coveredEdge(x, y);
                Set<Node> tuck = new HashSet<>();
                tuck.add(x);
                tuck.add(y);

                if (covered && tucks.contains(tuck)) continue;

                int i = scorer.index(x);
                scorer.bookmark(currentDepth);

                boolean first = true;
                List<Node> Z = new ArrayList<>(scorer.getOrderShallow().subList(i + 1, scorer.index(y)));
                Iterator<Node> zItr = Z.iterator();
                do {
                    if(first) {
                        scorer.moveTo(y, i);
                        first = false;
                    } else {
                        Node z = zItr.next();
                        if (ancestors.contains(z)) {
                            scorer.moveTo(z, i++);
                        } else continue;
                    }

                    if (violatesKnowledge(scorer.getPi())) continue;

                    sNew = scorer.score();
                    if (sNew > sOld) {
                        if (verbose) {
                            System.out.printf("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s %s \n",
                                    scorer.getNumEdges(), sNew - sOld, tucks, tuck);
                        }
                        return;
                    }

                    if (covered && currentDepth < depth) {
                        tucks.add(tuck);
                        if (!dfsHistory.contains(tucks)) {
                            dfsHistory.add(new HashSet<>(tucks));
                            graspDfs4(scorer, sOld, depth, currentDepth + 1, tucks, dfsHistory);
                        }
                        tucks.remove(tuck);
                    }

                    if (scorer.score() > sOld) return;

                } while (zItr.hasNext());

                scorer.goToBookmark(currentDepth);
            }
        }
    }

    private boolean timeOut() {
        double time = (System.currentTimeMillis() - start) / 1000.0;
        return timeout > -1 && time > timeout;
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

    /**
     * Need to set this before setting uncovered depth.
     *
     * @param depth Overall recusion depth. Uncovered depth needs to be less than or equal to this.
     */
    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.depth = depth;
    }

    /**
     * Need to set this after setting depth. Uncovered depth must be <= depth.
     *
     * @param uncoveredDepth Maximum depth for uncovered edges once all covered edges have been reached in the
     *                       depth first recursion.
     */
    public void setUncoveredDepth(int uncoveredDepth) {
        if (uncoveredDepth > this.depth) {
            throw new IllegalArgumentException("Uncovered depth cannot exceed depth");
        }

        this.uncoveredDepth = uncoveredDepth;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }


    public void setUseForwardTuckOnly(boolean useForwardTuckOnly) {
        this.useForwardTuckOnly = useForwardTuckOnly;
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


    public void setGraspAlg(boolean graspAlg) {
        this.graspAlg = graspAlg;
    }

    public void setTimeout(int timeout) {
        this.timeout = timeout;
    }
}