package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.NumberFormatUtil;
import edu.cmu.tetrad.util.TetradLogger;
import org.jetbrains.annotations.NotNull;

import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.*;

import static edu.cmu.tetrad.graph.Edges.directedEdge;
import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Math.min;
import static java.util.Collections.shuffle;


/**
 * Implements the GRASP algorithms, with various execution flags.
 *
 * @author bryanandrews
 * @author josephramsey
 */
public class Boss {
    private final List<Node> variables;
    private Score score;
    private IndependenceTest test;
    private IKnowledge knowledge = new Knowledge2();
    private TeyssierScorer scorer;
    private long start;
    // flags
    private boolean useScore = true;
    private boolean usePearl;
    private boolean cachingScores = true;
    private boolean useDataOrder = true;

    private boolean verbose = true;

    // other params
    private int depth = 4;
    private int numStarts = 1;

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
        order = new ArrayList<>(order);

        this.scorer = new TeyssierScorer(this.test, this.score);
        this.scorer.setUseVermaPearl(this.usePearl);

        if (this.usePearl) {
            this.scorer.setUseScore(false);
        } else {
            this.scorer.setUseScore(this.useScore && !(this.score instanceof GraphScore));
        }

        this.scorer.setKnowledge(this.knowledge);
        this.scorer.clearBookmarks();

        this.scorer.setCachingScores(this.cachingScores);

        List<Node> bestPerm = null;
        double best = NEGATIVE_INFINITY;

        this.scorer.score(order);

        for (int r = 0; r < this.numStarts; r++) {
            if (Thread.interrupted()) break;

            if ((r == 0 && !this.useDataOrder) || r > 0) {
                shuffle(order);
            }

            this.start = System.currentTimeMillis();

            makeValidKnowledgeOrder(order);

            this.scorer.score(order);
            double s1, s2;

            do {
                betterMutation(scorer);
                s1 = scorer.score();
                this.graph = scorer.getGraph(true);
                bes();
                s2 = scorer.score(this.graph.getCausalOrdering());
            } while (s2 > s1);

            if (this.scorer.score() > best) {
                best = this.scorer.score();
                bestPerm = scorer.getPi();
            }
        }

        this.scorer.score(bestPerm);

        long stop = System.currentTimeMillis();

        if (this.verbose) {
            TetradLogger.getInstance().forceLogMessage("Final order = " + this.scorer.getPi());
            TetradLogger.getInstance().forceLogMessage("Elapsed time = " + (stop - start) / 1000.0 + " s");
        }

        return bestPerm;
    }

    public void betterMutation(@NotNull TeyssierScorer scorer) {
        List<Node> pi = scorer.getPi();
        double s;
        double sp = scorer.score(pi);

        do {
            s = sp;
            scorer.bookmark();

            for (Node k : scorer.getPi()) {
                sp = NEGATIVE_INFINITY;

                for (int j = 0; j < scorer.size(); j++) {
                    scorer.moveTo(k, j);

                    if (scorer.score() > sp) {
                        if (!violatesKnowledge(scorer.getPi())) {
                            sp = scorer.score();
                            scorer.bookmark();

                            if (verbose) {
                                System.out.print("\r# Edges = " + scorer.getNumEdges()
                                        + " Score = " + scorer.score()
                                        + " (betterMutation)"
                                        + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
                            }
                        }
                    }
                }

                scorer.goToBookmark();
            }
        } while (sp > s);
    }

    public void betterMutation2(@NotNull TeyssierScorer scorer) {
        double s;
        double sp;

        sp = scorer.score();
        do {
            s = sp;
            scorer.bookmark();

            for (int j = scorer.size() - 1; j >= 0; j--) {
                Node _j = scorer.get(j);
                for (int k = j - 1; k >= 0; k--) {

                    tuck(_j, k, scorer);

                    if (scorer.score() > sp) {
                        if (!violatesKnowledge(scorer.getPi())) {
                            sp = scorer.score();

                            if (verbose) {
                                System.out.print("\r# Edges = " + scorer.getNumEdges()
                                        + " Score = " + scorer.score()
                                        + " (betterMutation)"
                                        + " Elapsed " + ((System.currentTimeMillis() - start) / 1000.0 + " sp"));
                            }

                            scorer.bookmark();
                            scorer.moveTo(_j, scorer.index(_j) + 1);
                        }
                    }
                }

                scorer.goToBookmark();
                scorer.bookmark();
            }
        } while (sp > s);
    }

    private void tuck(Node k, int j, TeyssierScorer scorer) {
        if (!scorer.adjacent(k, scorer.get(j))) return;
        if (scorer.coveredEdge(k, scorer.get(j))) return;
        if (j >= scorer.index(k)) return;

        Set<Node> ancestors = scorer.getAncestors(k);
        for (int i = j+1; i < scorer.index(k); i++) {
            if (ancestors.contains(scorer.get(i))) {
                scorer.moveTo(scorer.get(i), j);
            }
        }
        scorer.moveTo(k, j);
    }

    public int getNumEdges() {
        return this.scorer.getNumEdges();
    }

    private void makeValidKnowledgeOrder(List<Node> order) {
        if (!this.knowledge.isEmpty()) {
            order.sort((o1, o2) -> {
                if (o1.getName().equals(o2.getName())) {
                    return 0;
                } else if (this.knowledge.isRequired(o1.getName(), o2.getName())) {
                    return 1;
                } else if (this.knowledge.isRequired(o2.getName(), o1.getName())) {
                    return -1;
                } else if (this.knowledge.isForbidden(o2.getName(), o1.getName())) {
                    return -1;
                } else if (this.knowledge.isForbidden(o1.getName(), o2.getName())) {
                    return 1;
                } else {
                    return 1;
                }
            });
        }
    }

    @NotNull
    public Graph getGraph(boolean cpDag) {
        if (this.scorer == null) throw new IllegalArgumentException("Please run algorithm first.");
        Graph graph = this.scorer.getGraph(cpDag);

        NumberFormat nf = NumberFormatUtil.getInstance().getNumberFormat();
        graph.addAttribute("score ", nf.format(this.scorer.score()));
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
        return this.verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
        this.test.setVerbose(verbose);
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public void setDepth(int depth) {
        if (depth < -1) throw new IllegalArgumentException("Depth should be >= -1.");
        this.depth = depth;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }

    private boolean violatesKnowledge(List<Node> order) {
        if (!this.knowledge.isEmpty()) {
            for (int i = 0; i < order.size(); i++) {
                for (int j = i + 1; j < order.size(); j++) {
                    if (this.knowledge.isForbidden(order.get(i).getName(), order.get(j).getName())) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    public void setUseRaskuttiUhler(boolean usePearl) {
        this.usePearl = usePearl;
    }

    public void setUseDataOrder(boolean useDataOrder) {
        this.useDataOrder = useDataOrder;
    }

    private Graph graph;
    private final SortedSet<Arrow> sortedArrowsBack = new ConcurrentSkipListSet<>();
    private boolean meekVerbose = false;
    private ConcurrentMap<Node, Integer> hashIndices;
    private final Map<Edge, ArrowConfigBackward> arrowsMapBackward = new ConcurrentHashMap<>();
    private int arrowIndex = 0;


    private void buildIndexing(List<Node> nodes) {
        this.hashIndices = new ConcurrentHashMap<>();

        int i = -1;

        for (Node n : nodes) {
            this.hashIndices.put(n, ++i);
        }
    }

    private void bes() {
        buildIndexing(variables);

        reevaluateBackward(new HashSet<>(variables));

        while (!sortedArrowsBack.isEmpty()) {
            Arrow arrow = sortedArrowsBack.first();
            sortedArrowsBack.remove(arrow);

            Node x = arrow.getA();
            Node y = arrow.getB();

            if (!graph.isAdjacentTo(x, y)) {
                continue;
            }

            Edge edge = graph.getEdge(x, y);

            if (edge.pointsTowards(x)) {
                continue;
            }

            if (!getNaYX(x, y).equals(arrow.getNaYX())) {
                continue;
            }

            if (!new HashSet<>(graph.getParents(y)).equals(new HashSet<>(arrow.getParents()))) {
                continue;
            }

            if (!validDelete(x, y, arrow.getHOrT(), arrow.getNaYX())) {
                continue;
            }

            Set<Node> complement = new HashSet<>(arrow.getNaYX());
            complement.removeAll(arrow.getHOrT());

            double _bump = deleteEval(x, y, complement,
                    arrow.parents, hashIndices);

            delete(x, y, arrow.getHOrT(), _bump, arrow.getNaYX());

            Set<Node> process = revertToCPDAG();
            process.add(x);
            process.add(y);
            process.addAll(graph.getAdjacentNodes(x));
            process.addAll(graph.getAdjacentNodes(y));

            reevaluateBackward(new HashSet<>(process));
        }
    }

    private void delete(Node x, Node y, Set<Node> H, double bump, Set<Node> naYX) {
        Edge oldxy = graph.getEdge(x, y);

        Set<Node> diff = new HashSet<>(naYX);
        diff.removeAll(H);

        graph.removeEdge(oldxy);

        int numEdges = graph.getNumEdges();
        if (numEdges % 1000 == 0) {
            System.out.println("Num edges (backwards) = " + numEdges);
        }

        if (verbose) {
            int cond = diff.size() + graph.getParents(y).size();

            String message = (graph.getNumEdges()) + ". DELETE " + x + " --> " + y
                    + " H = " + H + " NaYX = " + naYX
                    + " degree = " + GraphUtils.getDegree(graph)
                    + " indegree = " + GraphUtils.getIndegree(graph)
                    + " diff = " + diff + " (" + bump + ") "
                    + " cond = " + cond;
            TetradLogger.getInstance().forceLogMessage(message);
        }

        for (Node h : H) {
            if (graph.isParentOf(h, y) || graph.isParentOf(h, x)) {
                continue;
            }

            Edge oldyh = graph.getEdge(y, h);

            graph.removeEdge(oldyh);

            graph.addEdge(directedEdge(y, h));

            if (verbose) {
                TetradLogger.getInstance().forceLogMessage("--- Directing " + oldyh + " to "
                        + graph.getEdge(y, h));
            }

            Edge oldxh = graph.getEdge(x, h);

            if (Edges.isUndirectedEdge(oldxh)) {
                graph.removeEdge(oldxh);

                graph.addEdge(directedEdge(x, h));

                if (verbose) {
                    TetradLogger.getInstance().forceLogMessage("--- Directing " + oldxh + " to "
                            + graph.getEdge(x, h));
                }
            }
        }
    }


    private double deleteEval(Node x, Node y, Set<Node> complement, Set<Node> parents,
                              Map<Node, Integer> hashIndices) {
        Set<Node> set = new HashSet<>(complement);
        set.addAll(parents);
        set.remove(x);

        return -scoreGraphChange(x, y, set, hashIndices);
    }

    private double scoreGraphChange(Node x, Node y, Set<Node> parents,
                                    Map<Node, Integer> hashIndices) {
        int xIndex = hashIndices.get(x);
        int yIndex = hashIndices.get(y);

        if (x == y) {
            throw new IllegalArgumentException();
        }

        if (parents.contains(y)) {
            throw new IllegalArgumentException();
        }

        int[] parentIndices = new int[parents.size()];

        int count = 0;
        for (Node parent : parents) {
            parentIndices[count++] = hashIndices.get(parent);
        }

        return score.localScoreDiff(xIndex, yIndex, parentIndices);
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    private Set<Node> revertToCPDAG() {
        MeekRules rules = new MeekRules();
        rules.setKnowledge(getKnowledge());
        rules.setAggressivelyPreventCycles(true);
        rules.setVerbose(meekVerbose);
        return rules.orientImplied(graph);
    }

    private boolean validDelete(Node x, Node y, Set<Node> H, Set<Node> naYX) {
        boolean violatesKnowledge = false;

        if (existsKnowledge()) {
            for (Node h : H) {
                if (knowledge.isForbidden(x.getName(), h.getName())) {
                    violatesKnowledge = true;
                }

                if (knowledge.isForbidden(y.getName(), h.getName())) {
                    violatesKnowledge = true;
                }
            }
        }

        Set<Node> diff = new HashSet<>(naYX);
        diff.removeAll(H);
        return isClique(diff) && !violatesKnowledge;
    }

    private boolean existsKnowledge() {
        return !knowledge.isEmpty();
    }

    private boolean isClique(Set<Node> nodes) {
        List<Node> _nodes = new ArrayList<>(nodes);
        for (int i = 0; i < _nodes.size(); i++) {
            for (int j = i + 1; j < _nodes.size(); j++) {
                if (!graph.isAdjacentTo(_nodes.get(i), _nodes.get(j))) {
                    return false;
                }
            }
        }

        return true;
    }

    private Set<Node> getNaYX(Node x, Node y) {
        List<Node> adj = graph.getAdjacentNodes(y);
        Set<Node> nayx = new HashSet<>();

        for (Node z : adj) {
            if (z == x) {
                continue;
            }
            Edge yz = graph.getEdge(y, z);
            if (!Edges.isUndirectedEdge(yz)) {
                continue;
            }
            if (!graph.isAdjacentTo(z, x)) {
                continue;
            }
            nayx.add(z);
        }

        return nayx;
    }

    private void reevaluateBackward(Set<Node> toProcess) {
        class BackwardTask extends RecursiveTask<Boolean> {
            private final Node r;
            private final List<Node> adj;
            private final Map<Node, Integer> hashIndices;
            private final int chunk;
            private final int from;
            private final int to;

            private BackwardTask(Node r, List<Node> adj, int chunk, int from, int to,
                                 Map<Node, Integer> hashIndices) {
                this.adj = adj;
                this.hashIndices = hashIndices;
                this.chunk = chunk;
                this.from = from;
                this.to = to;
                this.r = r;
            }

            @Override
            protected Boolean compute() {
                if (to - from <= chunk) {
                    for (int _w = from; _w < to; _w++) {
                        final Node w = adj.get(_w);
                        Edge e = graph.getEdge(w, r);

                        if (e != null) {
                            if (e.pointsTowards(r)) {
                                calculateArrowsBackward(w, r);
                            } else if (e.pointsTowards(w)) {
                                calculateArrowsBackward(r, w);
                            } else {
                                calculateArrowsBackward(w, r);
                                calculateArrowsBackward(r, w);
                            }
                        }
                    }

                } else {
                    int mid = (to - from) / 2;

                    List<BackwardTask> tasks = new ArrayList<>();

                    tasks.add(new BackwardTask(r, adj, chunk, from, from + mid, hashIndices));
                    tasks.add(new BackwardTask(r, adj, chunk, from + mid, to, hashIndices));

                    invokeAll(tasks);
                }

                return true;
            }
        }

        for (Node r : toProcess) {
            List<Node> adjacentNodes = new ArrayList<>(toProcess);
            ForkJoinPool.commonPool().invoke(new BackwardTask(r, adjacentNodes, getChunkSize(adjacentNodes.size()), 0,
                    adjacentNodes.size(), hashIndices));
        }
    }

    private int getChunkSize(int n) {
        int chunk = n / Runtime.getRuntime().availableProcessors();
        if (chunk < 100) chunk = 100;
        return chunk;
    }

    private void calculateArrowsBackward(Node a, Node b) {
        if (existsKnowledge()) {
            if (!getKnowledge().noEdgeRequired(a.getName(), b.getName())) {
                return;
            }
        }

        Set<Node> naYX = getNaYX(a, b);
        Set<Node> parents = new HashSet<>(graph.getParents(b));

        List<Node> _naYX = new ArrayList<>(naYX);

        ArrowConfigBackward config = new ArrowConfigBackward(naYX, parents);
        ArrowConfigBackward storedConfig = arrowsMapBackward.get(directedEdge(a, b));
        if (storedConfig != null && storedConfig.equals(config)) return;
        arrowsMapBackward.put(directedEdge(a, b), new ArrowConfigBackward(naYX, parents));

        int _depth = min(depth, _naYX.size());

        final DepthChoiceGenerator gen = new DepthChoiceGenerator(_naYX.size(), _depth);//_naYX.size());
        int[] choice;
        Set<Node> maxComplement = null;
        double maxBump = Double.NEGATIVE_INFINITY;

        while ((choice = gen.next()) != null) {
            Set<Node> complement = GraphUtils.asSet(choice, _naYX);
            double _bump = deleteEval(a, b, complement, parents, hashIndices);

            if (_bump > maxBump) {
                maxBump = _bump;
                maxComplement = complement;
            }
        }

        if (maxBump > 0) {
            Set<Node> _H = new HashSet<>(naYX);
            _H.removeAll(maxComplement);
            addArrowBackward(a, b, _H, naYX, parents, maxBump);
        }
    }

    private void addArrowBackward(Node a, Node b, Set<Node> hOrT, Set<Node> naYX,
                                  Set<Node> parents, double bump) {
        Arrow arrow = new Arrow(bump, a, b, hOrT, null, naYX, parents, arrowIndex++);
        sortedArrowsBack.add(arrow);
    }

    private static class ArrowConfigBackward {
        private Set<Node> nayx;
        private Set<Node> parents;

        public ArrowConfigBackward(Set<Node> nayx, Set<Node> parents) {
            this.setNayx(nayx);
            this.setParents(parents);
        }

        public void setNayx(Set<Node> nayx) {
            this.nayx = nayx;
        }

        public Set<Node> getParents() {
            return parents;
        }

        public void setParents(Set<Node> parents) {
            this.parents = parents;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ArrowConfigBackward that = (ArrowConfigBackward) o;
            return nayx.equals(that.nayx) && parents.equals(that.parents);
        }

        @Override
        public int hashCode() {
            return Objects.hash(nayx, parents);
        }
    }


    private static class Arrow implements Comparable<Arrow> {

        private final double bump;
        private final Node a;
        private final Node b;
        private final Set<Node> hOrT;
        private final Set<Node> naYX;
        private final Set<Node> parents;
        private final int index;
        private Set<Node> TNeighbors;

        Arrow(double bump, Node a, Node b, Set<Node> hOrT, Set<Node> capTorH, Set<Node> naYX,
              Set<Node> parents, int index) {
            this.bump = bump;
            this.a = a;
            this.b = b;
            this.setTNeighbors(capTorH);
            this.hOrT = hOrT;
            this.naYX = naYX;
            this.index = index;
            this.parents = parents;
        }

        public double getBump() {
            return bump;
        }

        public Node getA() {
            return a;
        }

        public Node getB() {
            return b;
        }

        Set<Node> getHOrT() {
            return hOrT;
        }

        Set<Node> getNaYX() {
            return naYX;
        }

        // Sorting by bump, high to low. The problem is the SortedSet contains won't add a new element if it compares
        // to zero with an existing element, so for the cases where the comparison is to zero (i.e. have the same
        // bump, we need to determine as quickly as possible a determinate ordering (fixed) ordering for two variables.
        // The fastest way to do this is using a hash code, though it's still possible for two Arrows to have the
        // same hash code but not be equal. If we're paranoid, in this case we calculate a determinate comparison
        // not equal to zero by keeping a list. This last part is commened out by default.
        public int compareTo(@NotNull Arrow arrow) {

            final int compare = Double.compare(arrow.getBump(), getBump());

            if (compare == 0) {
                return Integer.compare(getIndex(), arrow.getIndex());
            }

            return compare;
        }

        public String toString() {
            return "Arrow<" + a + "->" + b + " bump = " + bump
                    + " t/h = " + hOrT
                    + " TNeighbors = " + getTNeighbors()
                    + " parents = " + parents
                    + " naYX = " + naYX + ">";
        }

        public int getIndex() {
            return index;
        }

        public Set<Node> getTNeighbors() {
            return TNeighbors;
        }

        public void setTNeighbors(Set<Node> TNeighbors) {
            this.TNeighbors = TNeighbors;
        }

        public Set<Node> getParents() {
            return parents;
        }
    }

}