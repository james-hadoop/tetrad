package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;

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
    private List<Node> variables;
    private List<Node> order;
    private final Score score;
    private double[] scores;
    private Node nodeMoving = null;
    private double[] movingScores = null;
    private List<Node> bookmarkedOrder;
    private double[] bookmarkedScores;
    private double[] bookmarkedMovingScores;
    private Node bookmarkedNodeMoving;
    private boolean cachingScores = true;
    private final Map<ScoreKey, Double> cache = new HashMap<>();
    private IKnowledge knowledge = new Knowledge2();
    private final Map<ScoreKey, Double> precedentScores = new HashMap<>();

    public TeyssierScorer(Score score) {
        this.score = score;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double score(List<Node> order) {
        this.order = new ArrayList<>(order);
        this.variables = score.getVariables();
        this.scores = new double[order.size()];
        this.precedentScores.clear();
        initializeScores();
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
        if (index == order.size() - 1) return false;

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
            score += scores[i];
        }

        return score;
    }

    private void initializeScores() {
        movingScores = new double[order.size()];
        Arrays.fill(movingScores, Double.NaN);

        for (int i = 0; i < order.size(); i++) {
            double nodeScore = getGrowShrink(i);
            scores[i] = nodeScore;
        }

        bookmarkedOrder = null;
        bookmarkedScores = null;
        nodeMoving = null;
    }

    private double score(Node n, Set<Node> pi) {
        if (cachingScores && cache.containsKey(new ScoreKey(n, pi))) {
            return cache.get(new ScoreKey(n, pi));
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double v = -this.score.localScore(variables.indexOf(n), parentIndices);

        if (cachingScores) {
            cache.put(new ScoreKey(n, pi), v);
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

    private void recalculate(int i) {
        Set<Node> precedent = new HashSet<>(getPrefix(i));
        ScoreKey key = new ScoreKey(order.get(i), new HashSet<>(precedent));

        if (precedentScores.containsKey(key)) {
            scores[i] = precedentScores.get(key);
            return;
        }

        double score = getGrowShrink(i);

        precedentScores.put(key, score);
        scores[i] = score;
    }

    private double getGrowShrink(int p) {
        Node n = order.get(p);

        if (n != nodeMoving) {
            nodeMoving = n;
            for (int i = p; i < movingScores.length; i++) {
                movingScores[i] = Double.NaN;
            }
        } else if (!Double.isNaN(movingScores[p])) {
            return movingScores[p];
        }

        Set<Node> mb = new HashSet<>();
        boolean changed = true;
        double sMin = score(n, new HashSet<>());;
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

                if (s2 < sMin) {
                    sMin = s2;
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

                        if (s2 < sMin) {
                            sMin = s2;
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

        if (n == nodeMoving) {
            movingScores[p] = sMin;
        }

        return sMin;
    }

    public Set<Node> getMb(int p) {
        Node n = order.get(p);

        Set<Node> mb = new HashSet<>();
        boolean changed = true;
        double sMin = score(n, new HashSet<>());
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

                if (s2 < sMin) {
                    sMin = s2;
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

                        if (s2 < sMin) {
                            sMin = s2;
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

        return mb;
    }

    public void bookmark() {
        if (bookmarkedOrder == null) {
            bookmarkedOrder = new ArrayList<>();
        }

        bookmarkedOrder.clear();
        bookmarkedOrder.addAll(order);

        if (bookmarkedScores == null) {
            bookmarkedScores = new double[scores.length];
        }

        System.arraycopy(scores, 0, bookmarkedScores, 0, scores.length);

        if (bookmarkedMovingScores == null) {
            bookmarkedMovingScores = new double[scores.length];
        }

        bookmarkedMovingScores = Arrays.copyOf(movingScores, movingScores.length);
        bookmarkedNodeMoving = nodeMoving;
    }

    public void restoreBookmark() {
        if (bookmarkedOrder == null) {
            throw new IllegalArgumentException("The state has not yet been bookmarked.");
        }

        order.clear();
        order.addAll(bookmarkedOrder);

        System.arraycopy(bookmarkedScores, 0, scores, 0, bookmarkedScores.length);
        System.arraycopy(bookmarkedMovingScores, 0, movingScores, 0, bookmarkedMovingScores.length);

        nodeMoving = bookmarkedNodeMoving;
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public int size() {
        return order.size();
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
