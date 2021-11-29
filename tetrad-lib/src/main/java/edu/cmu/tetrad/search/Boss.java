package edu.cmu.tetrad.search;

import edu.cmu.tetrad.algcomparison.statistic.BicEst;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodePair;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Double.NEGATIVE_INFINITY;


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
    private boolean firstRunUseDataOrder = false;
    private int triangleDepth = -1;
    private TeyssierScorer.ParentCalculation parentCalculation = TeyssierScorer.ParentCalculation.GrowShrinkMb;
    private TeyssierScorer scorer;
    private Set<TriangleProblem> triangleProblems = new HashSet<>();

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
        scorer.evaluate(bestPerm);
        double best = scorer.score();

        for (int r = 0; r < numStarts; r++) {
            if (firstRunUseDataOrder) {
                if (r > 0)
                    scorer.shuffleVariables();
            } else {
                scorer.shuffleVariables();
            }

            List<Node> _order = scorer.getOrder();
            makeValidKnowledgeOrder(_order);
            scorer.evaluate(_order);

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else if (method == Method.GASP) {
                perm = gasp(scorer);
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
            System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
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

    public List<Node> boss(@NotNull TeyssierScorer scorer) {
        triangleProblems = new HashSet<>();

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal from). Score each causal from by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated. Then try to remove edges X--Y by reorienting edge
        // in triangles X--V--Y, repeatedly, until no more such edges can be removed.
        // Apparently, it suffices in the oracle to consider only on such triangle at
        // a time, though to hedge bets we allow the user to specify a maximum
        // number of triangles to consider per such edge.

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        while (true) {

            bossLoop(scorer);

            double s = scorer.score();
            triangleLoop(scorer);
            if (scorer.score() == s) break;
        }

        return scorer.getOrder();
    }

    private void bossLoop(@NotNull TeyssierScorer scorer) {
        double s = scorer.score();
        scorer.bookmark();

        do {
            for (Node v : scorer.getOrder()) {
                for (int i = scorer.size() - 1; i >= 0; i--) {
                    scorer.moveTo(v, i);

                    if (scorer.score() > s) {
                        if (satisfiesKnowledge(scorer.getOrder())) {
                            s = scorer.score();
                            scorer.bookmark();
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (scorer.score() > s);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        scorer.score();
    }

    private void triangleLoop(@NotNull TeyssierScorer scorer) {
        if (triangleDepth == 0) return;

        Set<NodePair> path = new HashSet<>();
        int _triangleDepth = triangleDepth == -1 ? Integer.MAX_VALUE : triangleDepth;
        double s = scorer.score();
        List<Node> order = scorer.getOrder();

        for (int i = scorer.size() - 1; i >= 0; i--) {
            Node x = scorer.get(i);

            for (int j = i - 1; j >= 0; j--) {
                Node y = scorer.get(j);

                List<NodePair> ZZ = new ArrayList<>();

                for (Node z : scorer.getOrder()) {
                    if (scorer.triangle(x, y, z)) {
                        ZZ.add(new NodePair(x, z));
                        ZZ.add(new NodePair(y, z));
                    }
                }

                triangleVisit(scorer, ZZ, path, _triangleDepth);

                if (scorer.score() > s) {
                    if (verbose) {
                        System.out.println("# Edges = " + scorer.getNumEdges() + " Score = "
                                + scorer.score() + " (Triangle)");
                    }

                    return;
                }

                scorer.evaluate(order);
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Triangle)");
        }
    }

    private void triangleVisit(TeyssierScorer scorer, List<NodePair> ZZ, Set<NodePair> path, int depth) {
        if (path.size() > depth) return;

        double s = scorer.score();
        List<Node> pi = scorer.getOrder();

        for (NodePair z : ZZ) {
            if (path.contains(z)) continue;

            scorer.swap(z.getFirst(), z.getSecond());

            if (scorer.score() > s) {
                return;
            }

            path.add(z);

            triangleVisit(scorer, ZZ, path, depth);

            path.remove(z);

            if (scorer.score() > s) {
                return;
            }

            scorer.evaluate(pi);
        }
    }

    private List<Node> gasp1(@NotNull TeyssierScorer scorer) {
        if (scorer.size() < 2) return scorer.getOrder();

        List<NodePair> list = new ArrayList<>();

        for (int i = 0; i < scorer.size(); i++) {
            for (int j = i + 1; j < scorer.size(); j++) {
                list.add(new NodePair(scorer.get(i), scorer.get(j)));
            }
        }

        Collections.reverse(list);

        for (int k = 0; k < 3; k++) {
            for (NodePair pair : list) {
                if (scorer.adjacent(pair.getFirst(), pair.getSecond())) {
                    scorer.bookmark();
                    double sOld = scorer.score();

                    scorer.reverse(pair.getFirst(), pair.getSecond());
                    double sNew = scorer.score();

                    if (sNew >= sOld) {

                        // Accept change.
                        continue;
                    }

                    // Change is bad.
                    scorer.goToBookmark();
                }
            }
        }

        return scorer.getOrder();
    }

    private List<Node> gasp(@NotNull TeyssierScorer scorer) {
        if (scorer.size() < 2) return scorer.getOrder();

        List<NodePair> list = new ArrayList<>();

        for (int i = 0; i < scorer.size(); i++) {
            for (int j = i + 1; j < scorer.size(); j++) {
                list.add(new NodePair(scorer.get(i), scorer.get(j)));
            }
        }

        Collections.reverse(list);
        Deque<NodePair> deque = new LinkedList<>(list);

        Map<NodePair, Integer> counts = new HashMap<>();
        for (NodePair pair : deque) counts.put(pair, 0);


        do {
            NodePair pair = deque.removeFirst();
            deque.addLast(pair);
            counts.put(pair, counts.get(pair) + 1);

            if (counts.get(deque.peek()) > 2) {
                scorer.goToBookmark();
                continue;
            }

            if (scorer.adjacent(pair.getFirst(), pair.getSecond())) {
                scorer.bookmark();

                double sOld = scorer.score();
                scorer.reverse(pair.getFirst(), pair.getSecond());
                double sNew = scorer.score();

                if (sNew >= sOld) {
                    continue;
                }

                scorer.goToBookmark();
            }
        } while (counts.get(deque.peek()) <= 5);

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

    public void setFirstRunUseDataOrder(boolean firstRunUseDataOrder) {
        this.firstRunUseDataOrder = firstRunUseDataOrder;
    }

    public void setTriangleDepth(int triangleDepth) {
        if (triangleDepth < -1) throw new IllegalArgumentException("Depth should be >= -1.");

        this.triangleDepth = triangleDepth;
    }

    public void setParentCalculation(TeyssierScorer.ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public enum Method {BOSS, SP, GASP}

    private static class TriangleProblem {
        public Set<Node> XX;
        public Set<NodePair> ZZ;

        public TriangleProblem(Node x, Node y, List<NodePair> ZZ) {
            Set<Node> XX = new HashSet<>();
            XX.add(x);
            XX.add(y);
            this.XX = XX;
            this.ZZ = new HashSet<>(ZZ);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof TriangleProblem)) return false;
            TriangleProblem that = (TriangleProblem) o;
            return XX.equals(that.XX) && ZZ.equals(that.ZZ);
        }

        @Override
        public int hashCode() {
            return Objects.hash(XX, ZZ);
        }
    }
}
