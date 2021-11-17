package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.util.Collections.shuffle;


/**
 * Implements a scorer as in Teyssier, M., & Koller, D. (2012). Ordering-based search: A simple and effective
 * algorithm for learning Bayesian networks. arXiv preprint arXiv:1207.1429. You give it a score function
 * and a variable ordering, and it computes the score. You can move any variable left or right and it will
 * keep track of the score using the Teyssier and Kohler method. You can move a vairable to a new position,
 * and you can bookmark a state and come baqck to it.
 *
 * @author josephramsey
 */
public class TeyssierScorer {
    private final Map<ScoreKey, Pair> cache = new HashMap<>();
    private final List<Node> variables;
    private final Map<Node, Integer> variablesHash;
    private ParentCalculation parentCalculation = ParentCalculation.GrowShrinkMb;
    private LinkedList<Node> bookmarkedOrder = new LinkedList<>();
    private LinkedList<Pair> bookmarkedScores = new LinkedList<>();
    private HashMap<Node, Integer> bookmarkedNodesHash = new HashMap<>();
    private Score score;
    private IndependenceTest test;
    private LinkedList<Node> order;
    private LinkedList<Pair> scores;
    private boolean cachingScores = true;
    private IKnowledge knowledge = new Knowledge2();
    private LinkedList<Set<Node>> prefixes;
    private ScoreType scoreType = ScoreType.Edge;
    private Map<Node, Integer> orderHash = new HashMap<>();

    public TeyssierScorer(Score score) {
        this.score = score;
        this.order = new LinkedList<>(score.getVariables());

        this.orderHash = new HashMap<>();
        nodesHash(orderHash, order);

        this.variables = score.getVariables();
        this.variablesHash = new HashMap<>();
        nodesHash(variablesHash, variables);
    }

    public TeyssierScorer(IndependenceTest test) {
        this.test = test;
        this.variables = test.getVariables();
        this.variablesHash = new HashMap<>();
        nodesHash(variablesHash, variables);
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public void evaluate(List<Node> order) {
        this.order = new LinkedList<>(order);
        initializeScores();
        this.bookmarkedOrder = new LinkedList<>(this.order);
        this.bookmarkedScores = new LinkedList<>(this.scores);
    }

    private void nodesHash(Map<Node, Integer> nodesHash, List<Node> variables) {
        for (int i = 0; i < variables.size(); i++) {
            nodesHash.put(variables.get(i), i);
        }
    }

    public double score() {
        return sum();
    }

    public boolean promote(Node v) {
        int index = orderHash.get(v);
        if (index == 0) return false;

        Node v1 = order.get(index - 1);
        Node v2 = order.get(index);

        order.set(index - 1, v2);
        order.set(index, v1);

        updateScores(index - 1, index);

        return true;
    }

    public boolean demote(Node v) {
        int index = orderHash.get(v);
        if (index >= size() - 1) return false;
        if (index == -1) return false;

        Node v1 = order.get(index);
        Node v2 = order.get(index + 1);

        order.set(index, v2);
        order.set(index + 1, v1);

        updateScores(index, index + 1);

        return true;
    }

    public void moveTo(Node v, int toIndex) {
        int vindex = indexOf(v);

        order.remove(v);
        order.add(toIndex, v);

        if (vindex < toIndex) {
            updateScores(vindex, toIndex);
        } else {
            updateScores(toIndex, vindex);
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

    public void moveToFirst(Node v) {
        int vindex = indexOf(v);

        order.remove(v);
        order.addFirst(v);

        updateScores(0, vindex);
    }

    public void moveToLast(Node v) {
        int vindex = indexOf(v);

        order.remove(v);
        order.addLast(v);

        updateScores(vindex, order.size() - 1);
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

    public List<Node> getOrder() {
        return new ArrayList<>(order);
    }

    public int indexOf(Node v) {
        if (!orderHash.containsKey(v)) {
            System.out.println();
        }

        return orderHash.get(v);
    }

    private double sum() {
        double score = 0;

        for (int i = 0; i < order.size(); i++) {
            double score1 = scores.get(i).getScore();

            score += (Double.isNaN(score1) || Double.isInfinite(score1)) ? 0 : score1;
//            score += score1;
        }

        return score;
    }

    private void initializeScores() {
        this.scores = new LinkedList<>();
        for (int i1 = 0; i1 < order.size(); i1++) this.scores.add(null);

        this.prefixes = new LinkedList<>();
        for (int i1 = 0; i1 < order.size(); i1++) this.prefixes.add(null);

        updateScores(0, order.size() - 1);
    }

    private void updateScores(int i1, int i2) {
        for (int i = i1; i <= i2; i++) {
            recalculate(i);
            orderHash.put(order.get(i), i);
        }
    }

    private double score(Node n, Set<Node> pi) {
        if (cachingScores) {
            ScoreKey key = new ScoreKey(n, pi);
            Pair pair = cache.get(key);

            if (pair != null) {
                return pair.getScore();
            }
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variablesHash.get(p);
        }

        double v = this.score.localScore(variablesHash.get(n), parentIndices);

        if (cachingScores) {
            ScoreKey key = new ScoreKey(n, pi);
            cache.put(key, new Pair(pi, v));
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
        if (!getPrefix(p).equals(prefixes.get(p))) {
            scores.set(p, getParentsInternal(p));
        }
    }

    public Set<Node> getParents(int p) {
        return new HashSet<>(scores.get(p).getParents());
    }

    public Set<Node> getParents(Node v) {
        return new HashSet<>(scores.get(indexOf(v)).getParents());
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

    public Graph getGraph(boolean cpdag) {
        List<Node> order = getOrder();
        Graph G1 = new EdgeListGraph(variables);

        for (int p = 0; p < order.size(); p++) {
            for (Node z : getParents(p)) {
                G1.addDirectedEdge(z, order.get(p));
            }
        }

        GraphUtils.replaceNodes(G1, variables);

        if (cpdag) {
            return SearchGraphUtils.cpdagForDag(G1);
        } else {
            return G1;
        }
    }

    private Pair getParentsInternal(int p) {
        if (parentCalculation == ParentCalculation.GrowShrinkMb) {
            if (test != null) {
                return getGrowShrinkIndep(p);
            } else {
                return getGrowShrinkScore(p);
            }
        } else if (parentCalculation == ParentCalculation.VermaPearl) {
            if (test != null) {
                return vermaPearl(p);
            } else {
                return vermaPearlScore(p);
            }
        } else {
            throw new IllegalStateException("Unrecognized parent calculation: " + parentCalculation);
        }
    }

    private Pair vermaPearl(int p) {
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

    private Pair vermaPearlScore(int p) {
        Node x = order.get(p);
        Set<Node> parents = new HashSet<>();
        Set<Node> prefix = getPrefix(p);

        double s1 = score(x, new HashSet<>(prefix));

        for (Node y : prefix) {
            List<Node> minus = new ArrayList<>(prefix);
            minus.remove(y);

            double s2 = score(x, new HashSet<>(minus));

            if (s2 < s1) {
                parents.add(y);
            }
        }

        return new Pair(parents, parents.size());
    }

    private Pair getGrowShrinkIndep(int p) {
        Node n = order.get(p);

        Set<Node> parents = new HashSet<>();

        Set<Node> prefix = getPrefix(p);

//        System.out.println("GS Indep target = " + n + " prefix = " + prefix + " order = " + order);

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
//
////        System.out.println("MB(" + n + ", " + prefix + ") = " + parents + " order = " + order);

        return new Pair(parents, -parents.size());
    }

    @NotNull
    private Pair getGrowShrinkScore(int p) {
        Node n = order.get(p);

        Set<Node> parents = new HashSet<>();
        boolean changed = true;

        double sMax = score(n, new HashSet<>());
        Set<Node> prefix = getPrefix(p);

        // Grow-shrink
        while (changed) {
            changed = false;

            // Let z be the node that maximizes the score...
            Node z = null;

            for (Node z0 : prefix) {
                if (parents.contains(z0)) continue;

                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;
                parents.add(z0);

                double s2 = score(n, parents);

                if (s2 > sMax) {
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

                double s2 = score(n, parents);

                if (s2 > sMax) {
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
            return new Pair(parents, Double.isNaN(sMax) ? Double.POSITIVE_INFINITY : sMax);
        } else {
            throw new IllegalStateException("Unexpected score type: " + scoreType);
        }
    }

    public void bookmark() {
        this.bookmarkedOrder = new LinkedList<>(order);
        this.bookmarkedScores = new LinkedList<>(scores);
        this.bookmarkedNodesHash = new HashMap<>(orderHash);
    }

    public void goToBookmark() {
        this.order = new LinkedList<>(bookmarkedOrder);
        this.scores = new LinkedList<>(bookmarkedScores);
        this.orderHash = new HashMap<>(bookmarkedNodesHash);
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
        evaluate(order);
    }

    public void setScoreType(ScoreType scoreType) {
        this.scoreType = scoreType;
    }

    public Node get(int j) {
        return order.get(j);
    }

    public boolean adjacent(Node a, Node c) {
        return getParents(indexOf(a)).contains(c) || getParents(indexOf(c)).contains(a);
    }

    public boolean collider(Node a, Node b, Node c) {
        return getParents(b).contains(a) && getParents(b).contains(c);
    }

    public void setParentCalculation(ParentCalculation parentCalculation) {
        this.parentCalculation = parentCalculation;
    }

    public double scoreAppended(Node v) {
        order.addLast(v);

        orderHash.put(v, order.indexOf(v));

        Pair pair = getParentsInternal(indexOf(v));
        order.removeLast();
        orderHash.remove(v);
        return pair.score;
    }

    public boolean triangle(Node x, Node y, Node z) {
        return adjacent(x, y) && adjacent(y, z) && adjacent(x, z);
    }

    public enum ScoreType {Edge, SCORE}

    public enum ParentCalculation {GrowShrinkMb, VermaPearl}

    private static class Pair {
        private final Set<Node> parents;
        private final double score;

        private Pair(Set<Node> parents, double score) {
            this.parents = parents;
            this.score = score;
        }

        public Set<Node> getParents() {
            return parents;
        }

        public double getScore() {
            return score;
        }
    }

    public static class ScoreKey {
        private final Node y;
        private final Set<Node> pi;

        public ScoreKey(Node y, Set<Node> pi) {
            this.y = y;
            this.pi = new HashSet<>(pi);
        }

        public int hashCode() {
            return 3 * y.hashCode() + 7 * pi.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof ScoreKey)) {
                return false;
            }

            ScoreKey spec = (ScoreKey) o;
            return y.equals(spec.y) && this.pi.equals(spec.pi);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPi() {
            return pi;
        }
    }
}
