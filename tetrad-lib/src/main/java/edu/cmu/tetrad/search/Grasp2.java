package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.OrderedPair;
import edu.cmu.tetrad.util.TetradLogger;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static java.util.Collections.shuffle;


/**
 * Implements the GRASP algorithms, with various execution flags.
 *
 * @author bryanandrews
 * @author josephramsey
 */
public class Grasp2 {
    private final List<Node> variables;
    double sNew = Double.NaN;
    int numRounds = 8;
    private Score score;
    private IndependenceTest test;
    private IKnowledge knowledge = new Knowledge2();
    private TeyssierScorer2 scorer;
    private long start;
    // flags
    private boolean useScore = true;
    private boolean useVermaPearl = false;
    private boolean ordered = false;
    private boolean verbose = false;
    private boolean cachingScores = true;
    private int uncoveredDepth = 1;
    private int nonSingularDepth = 1;
    // other params
    private int depth = 4;
    private int numStarts = 1;

    public Grasp2(@NotNull Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
        this.useScore = true;
    }

    public Grasp2(@NotNull IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
        this.useScore = false;
    }

    public Grasp2(@NotNull IndependenceTest test, Score score) {
        this.test = test;
        this.score = score;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(@NotNull List<Node> order) {
        long start = System.currentTimeMillis();
        order = new ArrayList<>(order);

        scorer = new TeyssierScorer2(test, score);
        scorer.setUseVermaPearl(useVermaPearl);

        if (useVermaPearl) {
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

        for (int r = 0; r < numStarts; r++) {
            if (r > 0) {
                shuffle(order);
            }

            this.start = System.currentTimeMillis();

            makeValidKnowledgeOrder(order);

            scorer.score(order);

            List<Node> perm = grasp(scorer);

            scorer.score(perm);

            if (scorer.score() > best) {
                best = scorer.score();
                bestPerm = perm;
            }
        }

        scorer.score(bestPerm);

        long stop = System.currentTimeMillis();

        if (verbose) {
            TetradLogger.getInstance().forceLogMessage("Final order = " + scorer.getPi());
            TetradLogger.getInstance().forceLogMessage("Elapsed time = " + (stop - start) / 1000.0 + " s");
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

    public List<Node> grasp(@NotNull TeyssierScorer2 scorer) {
//        if (true) grasp4(scorer);

        scorer.clearBookmarks();
//        scorer.initCache();
        List<int[]> depths = new ArrayList<>();

        // GRaSP-TSP
        if (ordered && uncoveredDepth != 0 && nonSingularDepth != 0) {
            depths.add(new int[]{depth < 1 ? Integer.MAX_VALUE : depth, 0, 0});
        }

        // GRaSP-ESP
        if (ordered && nonSingularDepth != 0) {
            depths.add(new int[]{depth < 1 ? Integer.MAX_VALUE : depth,
                    uncoveredDepth < 0 ? Integer.MAX_VALUE : uncoveredDepth, 0});
        }

        // GRaSP
        depths.add(new int[]{depth < 1 ? Integer.MAX_VALUE : depth,
                uncoveredDepth < 0 ? Integer.MAX_VALUE : uncoveredDepth,
                nonSingularDepth < 0 ? Integer.MAX_VALUE : nonSingularDepth});

        double sNew = scorer.score();
        double sOld;

        for (int[] depth : depths) {
            do {
//                scorer.updateCache();
                sOld = sNew;
                graspDfs(scorer, sOld, depth, 1, new HashSet<>(), new HashSet<>());
                sNew = scorer.score();
            } while (sNew > sOld);
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " (GRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

    private void graspDfs(@NotNull TeyssierScorer2 scorer, double sOld, int[] depth, int currentDepth,
                          Set<Set<Node>> tucks, Set<Set<Set<Node>>> dfsHistory) {

        for (Node y : scorer.getShuffledVariables()) {
//            scorer.resetCacheIfTooBig(1000000);

            Set<Node> ancestors = scorer.getAncestors(y);
            List<Node> parents = new ArrayList<>(scorer.getParents(y));
            shuffle(parents);
            for (Node x : parents) {

                boolean covered = scorer.coveredEdge(x, y);
                boolean singular = true;
                Set<Node> tuck = new HashSet<>();
                tuck.add(x);
                tuck.add(y);

                if (covered && tucks.contains(tuck)) continue;
                if (currentDepth > depth[1] && !covered) continue;

                int i = scorer.index(x);
                scorer.bookmark(currentDepth);

                boolean first = true;
                List<Node> Z = new ArrayList<>(scorer.getOrderShallow().subList(i + 1, scorer.index(y)));
                Iterator<Node> zItr = Z.iterator();
                do {
                    if (first) {
                        scorer.moveTo(y, i);
                        first = false;
                    } else {
                        Node z = zItr.next();
                        if (ancestors.contains(z)) {
                            if (scorer.getParents(z).contains(x)) {
                                singular = false;
                            }
                            scorer.moveTo(z, i++);
                        }
                    }
                } while (zItr.hasNext());

                if (currentDepth > depth[2] && !singular) {
                    scorer.goToBookmark(currentDepth);
                    continue;
                }

                if (violatesKnowledge(scorer.getPi())) continue;

                sNew = scorer.score();
                if (sNew > sOld) {
                    if (verbose) {
                        String msg = String.format("Edges: %d \t|\t Score Improvement: %f \t|\t Tucks Performed: %s %s",
                                scorer.getNumEdges(), sNew - sOld, tucks, tuck);
                        TetradLogger.getInstance().forceLogMessage(msg);
                    }
                    return;
                }

                if (sNew == sOld && currentDepth < depth[0]) {
                    tucks.add(tuck);
                    if (currentDepth > depth[1]) {
                        if (!dfsHistory.contains(tucks)) {
                            dfsHistory.add(new HashSet<>(tucks));
                            graspDfs(scorer, sOld, depth, currentDepth + 1, tucks, dfsHistory);
                        }
                    } else {
                        graspDfs(scorer, sOld, depth, currentDepth + 1, tucks, dfsHistory);
                    }
                    tucks.remove(tuck);
                }

                if (scorer.score() > sOld) return;

                scorer.goToBookmark(currentDepth);
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

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.depth = depth;
    }

    public void setUncoveredDepth(int uncoveredDepth) {
        if (depth < -1) throw new IllegalArgumentException("Uncovered depth should be >= -1.");
        this.uncoveredDepth = uncoveredDepth;
    }

    public void setNonSingularDepth(int nonSingularDepth) {
        if (depth < -1) throw new IllegalArgumentException("Non-singular depth should be >= -1.");
        this.nonSingularDepth = nonSingularDepth;
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

    public void setOrdered(boolean ordered) {
        this.ordered = ordered;
    }

    public void setUseVermaPearl(boolean useVermaPearl) {
        this.useVermaPearl = useVermaPearl;
    }

    public List<Node> quickGrasp(@NotNull TeyssierScorer2 scorer) {
        if (numRounds <= 0) throw new IllegalArgumentException("For quickGRaSP, #rounds should be > 0");
        scorer.clearBookmarks();

        NumberFormat nf = new DecimalFormat("0.00");

        if (verbose) {
            System.out.println("\nInitial # edges = " + scorer.getNumEdges());
        }

        scorer.bookmark(0);
        double s0 = scorer.score();

        int unimproved = 0;
        int rounds = 0;
        int maxRounds = numRounds < 0 ? Integer.MAX_VALUE : numRounds;

        while (unimproved < depth && ++rounds < maxRounds) {
            if (verbose) {
                System.out.println("### Round " + (rounds));
            }

            List<OrderedPair<Node>> edges = scorer.getEdges();

            int numImprovements = 0;
            int numEquals = 0;

            for (int w = 0; w < edges.size(); w++) {
                OrderedPair<Node> pair = edges.get(w);
                Node x = pair.getFirst();
                Node y = pair.getSecond();

                scorer.bookmark(1);

                if (!scorer.getParents(y).contains(x)) {
                    scorer.goToBookmark(1);
                    continue;
                }

                // 'tuck' operation.
                scorer.moveTo(y, scorer.index(x));

//                scorer.restartCacheIfTooBig(100000);

                if (violatesKnowledge(scorer.getPi())) {
                    scorer.goToBookmark(1);
                    continue;
                }

                double sNew = scorer.score();

                if (sNew < s0) {
                    scorer.goToBookmark(1);
                }

                if (verbose) {
                    if (sNew > s0) {
                        numImprovements++;
                    }

                    if (sNew == s0) {
                        numEquals++;
                    }

                    if (sNew > s0) {
                        System.out.println("Round " + (rounds) + " # improvements = " + numImprovements
                                + " # unimproved = " + numEquals
                                + " # edges = " + scorer.getNumEdges() + " progress this round = " + nf.format(100D * ((w + 1) / (double) edges.size())) + "%");
                    }
                }

                scorer.bookmark(0);
                s0 = scorer.score();
            }

            if (numImprovements == 0) {
                unimproved++;
            } else {
                unimproved = 0;
            }
        }

        scorer.goToBookmark(0);

//        if (quickGraphDoFinalGrasp) {
//            grasp(scorer);
//        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges()
                    + " Score = " + scorer.score()
                    + " #round = " + numRounds
                    + " (quickGRaSP)"
                    + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " s"));
        }

        return scorer.getPi();
    }

}