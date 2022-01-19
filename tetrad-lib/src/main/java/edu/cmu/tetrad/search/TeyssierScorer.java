package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Math.floor;
import static java.util.Collections.shuffle;


/**
 * Implements a scorer as in Teyssier, M., & Koller, D. (2012). Ordering-based search: A simple and effective
 * algorithm for learning Bayesian networks. arXiv preprint arXiv:1207.1429. You give it a score function
 * and a variable ordering, and it computes the score. You can move any variable left or right, and it will
 * keep track of the score using the Teyssier and Kohler method. You can move a variable to a new position,
 * and you can bookmark a state and come back to it.
 *
 * @author josephramsey
 */
public class TeyssierScorer {
    private final List<Node> variables;
    private final Map<Node, Integer> variablesHash;
    private final Score score;
    private final IndependenceTest test;
    private final Map<Integer, LinkedList<Node>> bookmarkedOrders = new HashMap<>();
    private final Map<Integer, LinkedList<Pair>> bookmarkedScores = new HashMap<>();
    private final Map<Integer, Map<Node, Integer>> bookmarkedOrderHashes = new HashMap<>();
    private final Map<Integer, Float> bookmarkedRunningScores = new HashMap<>();
    private Map<Node, Map<Set<Node>, Float>> cache = new HashMap<>();
    private Map<Node, Integer> orderHash;
    private LinkedList<Node> order;
    private LinkedList<Pair> scores;
    private IKnowledge knowledge = new Knowledge2();
    private LinkedList<Set<Node>> prefixes;
    private ScoreType scoreType = ScoreType.SCORE;
    private float runningScore = 0F;

    private boolean useScore = true;
    private boolean usePearl = false;
    private boolean cachingScores = true;

    public TeyssierScorer(IndependenceTest test, Score score) {
        this.score = score;
        this.test = test;

        if (score != null) {
            this.variables = score.getVariables();
            this.order = new LinkedList<>(this.variables);
        } else if (test != null) {
            this.variables = test.getVariables();
            this.order = new LinkedList<>(this.variables);
        } else {
            throw new IllegalArgumentException("Need both a score and a test,");
        }

        this.orderHash = new HashMap<>();
        nodesHash(orderHash, order);

        this.variablesHash = new HashMap<>();
        nodesHash(variablesHash, variables);

        if (score instanceof GraphScore) {
            useScore = false;
        }
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public float score(List<Node> order) {
        runningScore = 0F;
        this.order = new LinkedList<>(order);
        this.scores = new LinkedList<>();

        for (int i1 = 0; i1 < order.size(); i1++) {
            this.scores.add(null);
        }

        this.prefixes = new LinkedList<>();
        for (int i1 = 0; i1 < order.size(); i1++) this.prefixes.add(null);
        initializeScores();
        return score();
    }

    public float score() {
        return runningScore;
//        return sum();
    }

    private float sum() {
        float score = 0;

        for (int i = 0; i < order.size(); i++) {
            float score1 = scores.get(i).getScore();

//            score += (Float.isNaN(score1) || Float.isInfinite(score1)) ? 0F : score1;
            score += score1;
        }

        return score;
    }

    public void tuck(Node x, Node y) {

        if (index(x) < index(y)) {
            moveTo(y, index(x));
        } else if (index(x) > index(y)) {
            moveTo(x, index(y));
        }

//        // x->y
//        if (getParents(y).contains(x)) {
//            moveTo(y, index(x));
//        }
//
//        // y -> x
//        else if (getParents(x).contains(y)) {
//            moveTo(x, index(y));
//        }
    }

    public void moveTo(Node v, int toIndex) {
        if (!order.contains(v)) return;

        int vIndex = index(v);

        if (vIndex == toIndex) return;

        if (lastMoveSame(vIndex, toIndex)) return;

        order.remove(v);
        order.add(toIndex, v);

        if (toIndex < vIndex) {
            updateScores(toIndex, vIndex);
        } else {
            updateScores(vIndex, toIndex);
        }
    }

    public boolean swap(Node m, Node n) {
        int i = orderHash.get(m);
        int j = orderHash.get(n);

        order.set(i, n);
        order.set(j, m);

        if (!validKnowledgeOrder(order)) {
            order.set(i, m);
            order.set(j, n);
            return false;
        }

        if (i < j) {
            updateScores(i, j);
        } else {
            updateScores(j, i);
        }

        return true;
    }

    public boolean coveredEdge(Node x, Node y) {
        if (!adjacent(x, y)) return false;
        Set<Node> px = getParents(x);
        Set<Node> py = getParents(y);
        px.remove(y);
        py.remove(x);
        return px.equals(py);
    }

    public List<Node> getOrder() {
        return new ArrayList<>(order);
    }

    public List<Node> getOrderShallow() {
        return order;
    }

    public int index(Node v) {
        if (!orderHash.containsKey(v)) {
            System.out.println();
        }

        Integer integer = orderHash.get(v);

        if (integer == null)
            throw new IllegalArgumentException("First 'evaluate' a permutation containing variable "
                    + v + ".");

        return integer;
    }

    public Set<Node> getParents(int p) {
        return new HashSet<>(scores.get(p).getParents());
    }

    public Set<Node> getParents(Node v) {
        return getParents(index(v));
    }

    public Set<Node> getAdjacentNodes(Node v) {
        Set<Node> adj = new HashSet<>();

        for (Node w : order) {
            if (getParents(v).contains(w) || getParents(w).contains(v)) {
                adj.add(w);
            }
        }

        return adj;
    }

    public Graph getGraph(boolean cpDag) {
        List<Node> order = getOrder();
        Graph G1 = new EdgeListGraph(variables);

        for (int p = 0; p < order.size(); p++) {
            for (Node z : getParents(p)) {
                G1.addDirectedEdge(z, order.get(p));
            }
        }

        GraphUtils.replaceNodes(G1, variables);

        if (cpDag) {
            return SearchGraphUtils.cpdagForDag(G1);
        } else {
            return G1;
        }
    }

    public List<NodePair> getAdjacencies() {
        List<Node> order = getOrder();
        Set<NodePair> pairs = new HashSet<>();

        for (int i = 0; i < order.size(); i++) {
            for (int j = 0; j < i; j++) {
                Node x = order.get(i);
                Node y = order.get(j);

                if (adjacent(x, y)) {
                    pairs.add(new NodePair(x, y));
                }
            }
        }

        return new ArrayList<>(pairs);
    }

    public List<OrderedPair<Node>> getEdges() {
        List<Node> order = getOrder();
        List<OrderedPair<Node>> edges = new ArrayList<>();

        for (Node y : order) {
            for (Node x : getParents(y)) {
                edges.add(new OrderedPair<>(x, y));
            }
        }

        return edges;
    }

    public Pair pearlParents(int p) {
        Node x = order.get(p);
        Set<Node> parents = new HashSet<>();
        Set<Node> prefix = getPrefix(p);

        for (Node y : prefix) {
            Set<Node> minus = new HashSet<>(prefix);
            minus.remove(y);

            if (test.isDependent(x, y, new ArrayList<>(minus))) {
                parents.add(y);
            }
        }

        return new Pair(parents, -parents.size());
    }

    private Pair pearlScore(int p) {
        Node x = order.get(p);
        Set<Node> parents = new HashSet<>();
        Set<Node> prefix = getPrefix(p);

        float s1 = score(x, new HashSet<>(prefix));

        for (Node y : prefix) {
            List<Node> minus = new ArrayList<>(prefix);
            minus.remove(y);

            float s2 = score(x, new HashSet<>(minus));

            if (s2 > s1) {
                parents.add(y);
            }
        }

        return new Pair(parents, parents.size());
    }

    public void bookmark(int key) {
        bookmarkedOrders.put(key, new LinkedList<>(order));
        bookmarkedScores.put(key, new LinkedList<>(scores));
        bookmarkedRunningScores.put(key, runningScore);
        bookmarkedOrderHashes.put(key, new HashMap<>(orderHash));
    }

    public void bookmark() {
        bookmark(Integer.MIN_VALUE);
    }

    public void goToBookmark(int key) {
        if (!bookmarkedOrders.containsKey(key)) {
            throw new IllegalArgumentException("That key was not bookmarked recently.");
        }

        order = bookmarkedOrders.remove(key);
        scores = bookmarkedScores.remove(key);
        runningScore = bookmarkedRunningScores.remove(key);
        orderHash = bookmarkedOrderHashes.remove(key);
    }

    public void clearBookmarks() {
        bookmarkedOrders.clear();
        bookmarkedScores.clear();
        bookmarkedRunningScores.clear();
        bookmarkedOrderHashes.clear();
    }

    public void goToBookmark() {
        goToBookmark(Integer.MIN_VALUE);
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public int size() {
        return order.size();
    }

    public int getNumEdges() {
        int numEdges = 0;

        for (int p = 0; p < order.size(); p++) {
            numEdges += getParents(p).size();
        }

        return numEdges;
    }

    public void shuffleVariables() {
        order = new LinkedList<>(order);
        shuffle(order);
        score(order);
    }

    public void setScoreType(ScoreType scoreType) {
        this.scoreType = scoreType;
    }

    public Node get(int j) {
        return order.get(j);
    }

    public boolean adjacent(Node a, Node c) {
        if (a == c) return false;
        return getParents(a).contains(c) || getParents(c).contains(a);
    }

    public boolean collider(Node a, Node b, Node c) {
        return getParents(b).contains(a) && getParents(b).contains(c);
    }

    public boolean triangle(Node x, Node y, Node z) {
        return adjacent(x, y) && adjacent(y, z) && adjacent(x, z);
    }

    public void setUseScore(boolean useScore) {
        if (!(this.score instanceof GraphScore)) {
            this.useScore = useScore;
        }
    }

    public boolean clique(List<Node> W) {
        for (int i = 0; i < W.size(); i++) {
            for (int j = i + 1; j < W.size(); j++) {
                if (!adjacent(W.get(i), W.get(j))) {
                    return false;
                }
            }
        }

        return true;
    }

    public void resetCacheIfTooBig(int maxSize) {
        if (cache.size() > maxSize) {
            cache = new HashMap<>();
            System.out.println("Clearing cacche...");
            System.gc();
        }
    }

    private boolean validKnowledgeOrder(List<Node> order) {
        for (int i = 0; i < order.size(); i++) {
            for (int j = i + 1; j < order.size(); j++) {
                if (knowledge.isForbidden(order.get(i).getName(), order.get(j).getName())) {
                    return false;
                }
            }
        }

        return true;
    }

    private void initializeScores() {
        for (int i1 = 0; i1 < order.size(); i1++) this.prefixes.set(i1, null);

        for (int i = 0; i < order.size(); i++) {
            recalculate(i);
            orderHash.put(order.get(i), i);
        }

        updateScores(0, order.size() - 1);
    }

    private void updateScores(int i1, int i2) {
        for (int i = i1; i <= i2; i++) {
            recalculate(i);
            orderHash.put(order.get(i), i);
        }
    }

    private float score(Node n, Set<Node> pi) {
        if (cachingScores) {
            cache.computeIfAbsent(n, w -> new HashMap<>());
            Float score = cache.get(n).get(pi);

            if (score != null) {
                return score;
            }
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variablesHash.get(p);
        }

        float v = (float) this.score.localScore(variablesHash.get(n), parentIndices);
//        v = Math.floor(1.e10 * v) * 1.e-10;
        v = Math.nextDown(v);

        if (cachingScores) {
            cache.computeIfAbsent(n, w -> new HashMap<>());
            cache.get(n).put(new HashSet<>(pi), v);
        }

        return v;
    }

    private Set<Node> getPrefix(int i) {
        Set<Node> prefix = new HashSet<>();

        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }

        return prefix;
    }

    private void recalculate(int p) {
        if (prefixes.get(p) == null || !prefixes.get(p).containsAll(getPrefix(p))) {
            Pair p1 = scores.get(p);
            Pair p2 = getParentsInternal(p);
            updateRunningScore(p1, p2);
            scores.set(p, p2);
        }
    }

    private void updateRunningScore(Pair p1, Pair p2) {
        float s1 = p1 == null ? 0F : p1.getScore();
        float s2 = p2 == null ? 0F : p2.getScore();
        runningScore -= s1;
        runningScore += s2;
    }

    private void nodesHash(Map<Node, Integer> nodesHash, List<Node> variables) {
        for (int i = 0; i < variables.size(); i++) {
            nodesHash.put(variables.get(i), i);
        }
    }

    private boolean lastMoveSame(int i1, int i2) {
        if (i1 <= i2) {
            for (int i = i1; i <= i2; i++) {
                if (!getPrefix(i).equals(prefixes.get(i))) return false;
            }
        } else {
            for (int i = i2; i <= i1; i++) {
                if (!getPrefix(i).equals(prefixes.get(i))) return false;
            }
        }

        return true;
    }

    @NotNull
    private Pair getGrowShrinkScore(int p) {
        Node n = order.get(p);

        Set<Node> parents = new HashSet<>();
        boolean changed = true;

        float sMax = score(n, new HashSet<>());
        List<Node> prefix = new ArrayList<>(getPrefix(p));

        // Grow-shrink
        while (changed) {
            changed = false;

            // Let z be the node that maximizes the score...
            Node z = null;

            for (Node z0 : prefix) {
                if (parents.contains(z0)) continue;

                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;
                parents.add(z0);

                float s2 = score(n, parents);

                if (s2 >= sMax) {
                    sMax = s2;
                    z = z0;
                }

                parents.remove(z0);
            }

            if (z != null) {
                parents.add(z);
                changed = true;
            }

        }

        boolean changed2 = true;

        while (changed2) {
            changed2 = false;

            Node w = null;

            for (Node z0 : new HashSet<>(parents)) {
                parents.remove(z0);

                float s2 = score(n, parents);

                if (s2 >= sMax) {
                    sMax = s2;
                    w = z0;
                }

                parents.add(z0);
            }

            if (w != null) {
                parents.remove(w);
                changed2 = true;
            }
        }

        if (scoreType == ScoreType.Edge) {
            return new Pair(parents, -parents.size());
        } else if (scoreType == ScoreType.SCORE) {
            return new Pair(parents, Float.isNaN(sMax) ? Float.POSITIVE_INFINITY : sMax);
        } else {
            throw new IllegalStateException("Unexpected score type: " + scoreType);
        }
    }

    private Pair getGrowShrinkIndependent(int p) {
        Node n = order.get(p);

        Set<Node> parents = new HashSet<>();

        Set<Node> prefix = getPrefix(p);

        boolean changed1 = true;

        while (changed1) {
            changed1 = false;

            for (Node z0 : prefix) {
                if (parents.contains(z0)) continue;
                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;

                if (test.isDependent(n, z0, new ArrayList<>(parents))) {
                    parents.add(z0);
                    changed1 = true;
                }
            }
        }

        boolean changed2 = true;

        while (changed2) {
            changed2 = false;

            for (Node z1 : new HashSet<>(parents)) {
                Set<Node> _p = new HashSet<>(parents);
                _p.remove(z1);

                if (test.isIndependent(n, z1, new ArrayList<>(_p))) {
                    parents.remove(z1);
                    changed2 = true;
                }
            }
        }

        return new Pair(parents, -parents.size());
    }

    private Pair getParentsInternal(int p) {
        if (usePearl) {
//            if (!useScore) {
            return pearlParents(p);
//            }
//            else {
//                return pearlScore(p);
//            }
        } else {
            if (useScore) {
                return getGrowShrinkScore(p);
            } else {
                return getGrowShrinkIndependent(p);
            }
        }
    }

    public void setUsePearl(boolean usePearl) {
        this.usePearl = usePearl;
    }

    public enum ScoreType {Edge, SCORE}

    public enum ParentCalculation {GrowShrinkMb, Pearl}

    private static class Pair {
        private final Set<Node> parents;
        private final float score;

        private Pair(Set<Node> parents, float score) {
            this.parents = parents;
            this.score = score;
        }

        public Set<Node> getParents() {
            return parents;
        }

        public float getScore() {
            return score;
        }

        public int hashCode() {
            return parents.hashCode() + (int) floor(10000D * score);
        }

        public boolean equals(Object o) {
            if (o == null) return false;
            if (!(o instanceof Pair)) return false;
            Pair thatPair = (Pair) o;
            return this.parents.equals(thatPair.parents) && this.score == thatPair.score;
        }
    }

//    public static class ScoreKey {
//        private final Set<Node> pi;
//
//        public ScoreKey(Set<Node> pi) {
//            this.pi = new HashSet<>(pi);
//        }
//
//        public int hashCode() {
//            return 11 * pi.hashCode();
//        }
//
//        public boolean equals(Object o) {
//            if (!(o instanceof ScoreKey)) {
//                return false;
//            }
//
//            ScoreKey spec = (ScoreKey) o;
//            return this.pi.equals(spec.pi);
//        }
//
//        public Set<Node> getPi() {
//            return pi;
//        }
//    }
}
