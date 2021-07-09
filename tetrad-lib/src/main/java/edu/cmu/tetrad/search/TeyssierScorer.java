package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

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
    private final Map<Integer, Pair[]> bookmarkedMovingScores = new HashMap<>();
    private final Map<Integer, Node> bookmarkedNodeMoving = new HashMap<>();
    private Score score;
    private IndependenceTest test;
    private List<Node> variables;
    private List<Node> order;
    private Pair[] scores;
    private Node nodeMoving = null;
    private Pair[] movingScores = null;
    private Map<Integer, List<Node>> bookmarkedOrder = new HashMap<>();
    private Map<Integer, Pair[]> bookmarkedScores = new HashMap<>();
    private boolean cachingScores = true;
    private IKnowledge knowledge = new Knowledge2();
    private List<Set<Node>> prefixes;

    public TeyssierScorer(Score score) {
        this.score = score;
    }

    public TeyssierScorer(IndependenceTest test) {
        this.test = test;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double score(List<Node> order) {
        this.order = new ArrayList<>(order);

        this.variables = score != null ? score.getVariables() : test.getVariables();

        this.scores = new Pair[order.size()];
        this.prefixes = new ArrayList<>();
        for (int i = 0; i < order.size(); i++) prefixes.add(null);
        initializeScores();
        bookmark(1);
        return score();
    }

    public double score() {
        return sum();
    }

    public boolean moveLeft(Node v) {
        int index = order.indexOf(v);
        if (index <= 0) return false;
        if (index >= order.size()) return false;

        Node v1 = order.get(index - 1);
        Node v2 = order.get(index);

        order.set(index - 1, v2);
        order.set(index, v1);

        recalculate(index);
        recalculate(index - 1);

        return true;
    }

    public boolean moveRight(Node v) {
        int index = order.indexOf(v);
        if (index < 0) return false;
        if (index >= order.size() - 1) return false;

        Node v1 = order.get(index + 1);
        Node v2 = order.get(index);

        order.set(index + 1, v2);
        order.set(index, v1);

        recalculate(index);
        recalculate(index + 1);

        return true;
    }

    public void moveTo(Node v, int i) {
        int index = order.indexOf(v);
        if (index == i) return;

        if (i < index) {
            while (--index >= i) {
                if (!moveLeft(v)) break;
            }
        } else {
            while (++index <= i) {
                if (!moveRight(v)) break;
            }
        }
    }

    public void swap(Node m, Node n) {
        int i = order.indexOf(m);
        int j = order.indexOf(n);

        order.set(i, n);
        order.set(j, m);

        score(order);
    }

    public List<Node> getOrder() {
        return new ArrayList<>(order);
    }

    public int indexOf(Node v) {
        return order.indexOf(v);
    }

    private double sum() {
        double score = 0;

        for (int i = 0; i < order.size(); i++) {
            score += scores[i].getScore();
        }

        return score;
    }

    private void initializeScores() {
        movingScores = new Pair[order.size()];
        Arrays.fill(movingScores, null);

        for (int i = 0; i < order.size(); i++) {
            recalculate(i);
        }

        bookmarkedOrder = new HashMap<>();
        bookmarkedScores = new HashMap<>();
        nodeMoving = null;
    }

    private double score(Node n, Set<Node> pi) {
        if (cachingScores) {
            ScoreKey key = new ScoreKey(n, pi);
            if (cache.containsKey(key)) {
                return cache.get(key).getScore();
            }
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double v = this.score.localScore(variables.indexOf(n), parentIndices);

        if (cachingScores) {
            ScoreKey key = new ScoreKey(n, pi);
            cache.put(key, new Pair(pi, v));
        }

        return v;
    }

    private List<Node> getPrefix(int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    private void recalculate(int p) {

        Node n = order.get(p);

        if (score instanceof GraphScore) {
            ((GraphScore) score).setN(n);
            ((GraphScore) score).setPrefix(getPrefix(p));
        }


        if (n != nodeMoving) {
            nodeMoving = n;
            for (int k = p; k < movingScores.length; k++) {
                movingScores[k] = null;
            }
        } else if (movingScores[p] != null) {
            scores[p] = movingScores[p];
            return;
        }

        if (new HashSet<>(getPrefix(p)).equals(prefixes.get(p))) {
            prefixes.set(p, new HashSet<>(getPrefix(p)));
        } else {
            scores[p] = getGrowShrink(p);
        }

        if (n == nodeMoving) {
            movingScores[p] = scores[p];
        }
    }

    public Set<Node> getMb(int p) {
        return scores[p].getMb();
    }

    private Pair getGrowShrink(int p) {
        if (test != null) {
            return getGrowShrinkIndep(p);
        }

        Node n = order.get(p);

        Set<Node> mb = new HashSet<>();
        boolean changed = true;

        double sMax = score(n, new HashSet<>());
        Set<Node> prefix = new HashSet<>(getPrefix(p));

        // Grow-shrink
        while (changed) {
            changed = false;

            // Let z be the node that maximizes the score...
            Node z = null;

            for (Node z0 : prefix) {
                if (mb.contains(z0)) continue;

                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;
                mb.add(z0);
                double s2 = score(n, mb);

                if (s2 > sMax) {
                    sMax = s2;
                    z = z0;
                }

                mb.remove(z0);
            }

            if (z != null) {
                mb.add(z);
                changed = true;
                boolean changed2 = true;

                while (changed2) {
                    changed2 = false;

                    Node w = null;

                    for (Node z0 : new HashSet<>(mb)) {
                        mb.remove(z0);

                        double s2 = score(n, mb);

                        if (s2 > sMax) {
                            sMax = s2;
                            w = z0;
                        }

                        mb.add(z0);
                    }

                    if (w != null) {
                        mb.remove(w);
                        changed2 = true;
                    }
                }
            }
        }

        return new Pair(mb, -sMax);
//        return new Pair(mb, mb.size());
    }

    private Pair getGrowShrinkIndep(int p) {
        Node n = order.get(p);

        Set<Node> mb = new HashSet<>();
        boolean changed = true;

        Set<Node> prefix = new HashSet<>(getPrefix(p));

        // Grow-shrink
        while (changed) {
            changed = false;

            for (Node z0 : prefix) {
                if (mb.contains(z0)) continue;

                if (test.isDependent(n, z0, new ArrayList<>(mb))) {
                    mb.add(z0);
                    changed = true;
                }
            }

            boolean changed2 = true;

            while (changed2) {
                changed2 = false;

                for (Node z1 : new HashSet<>(mb)) {
                    Set<Node> _mb = new HashSet<>(mb);
                    _mb.remove(z1);

                    if (test.isIndependent(n, z1, new ArrayList<>(_mb))) {
                        mb.remove(z1);
                        changed2 = true;
                    }
                }
            }
        }

        return new Pair(mb, mb.size());
    }

    public void bookmark(int index) {
        bookmarkedOrder.computeIfAbsent(index, k -> new ArrayList<>());

        bookmarkedOrder.get(index).clear();
        bookmarkedOrder.get(index).addAll(order);

        bookmarkedScores.computeIfAbsent(index, k -> new Pair[scores.length]);

        System.arraycopy(scores, 0, bookmarkedScores.get(index), 0, scores.length);

        bookmarkedMovingScores.computeIfAbsent(index, k -> new Pair[scores.length]);

        bookmarkedMovingScores.put(index, Arrays.copyOf(movingScores, movingScores.length));
        bookmarkedNodeMoving.put(index, nodeMoving);

//        System.out.println("BOOKMARKING: " + order + " score = " + score());

    }

    public void restoreBookmark(int index) {
        if (bookmarkedOrder == null) {
            bookmarkedOrder = new HashMap<>();
            bookmarkedOrder.put(index, new ArrayList<>());
        }

        order.clear();
        order.addAll(bookmarkedOrder.get(index));

        System.arraycopy(bookmarkedScores.get(index), 0, scores, 0, bookmarkedScores.get(index).length);
        System.arraycopy(bookmarkedMovingScores.get(index), 0, movingScores, 0, bookmarkedMovingScores.get(index).length);

        nodeMoving = bookmarkedNodeMoving.get(index);
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
            numEdges += getMb(p).size();
        }

        return numEdges;
    }

    private static class Pair {
        private final Set<Node> mb;
        private final double score;

        private Pair(Set<Node> mb, double score) {
            this.mb = mb;
            this.score = score;
        }

        public Set<Node> getMb() {
            return mb;
        }

        public double getScore() {
            return score;
        }
    }

    public static class ScoreKey {
        private final Node y;
        private final Set<Node> prefix;

        public ScoreKey(Node y, Set<Node> prefix) {
            this.y = y;
            this.prefix = new HashSet<>(prefix);
        }

        public int hashCode() {
            return 17 * y.hashCode() + 3 * prefix.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof ScoreKey)) {
                return false;
            }

            ScoreKey spec = (ScoreKey) o;
            return spec.y.equals(this.y) && spec.getPrefix().equals(this.prefix);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPrefix() {
            return prefix;
        }
    }
}
