package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;

public class TeyssierScorer {
    private final Score score;
    private final Map<ScoreKey, Double> cache = new HashMap<>();
    private final IndTestScore test;
    private List<Node> variables;
    private double[] scores;
    private LinkedList<Node> order;
    private IKnowledge knowledge = new Knowledge2();
    private Node nodeMoving = null;
    private double[] movingScores = null;
    private Map<Node, Set<Node>> associates;

    public TeyssierScorer(Score score) {
        this.score = score;
        this.test = new IndTestScore(score);
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double score(List<Node> order) {
        this.order = new LinkedList<>(order);
        this.variables = score.getVariables();
        this.associates = getAssociates(variables);
        this.scores = new double[variables.size()];
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

        order.remove(v);
        order.add(index - 1, v);

        recalculate(index);
        recalculate(index - 1);

        return true;
    }

    public boolean moveRight(Node v) {
        int index = order.indexOf(v);
        if (index < 0) return false;
        if (index == order.size() - 1) return false;

        order.remove(v);
        order.add(index + 1, v);

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

    public List<Node> getOrder() {
        return new ArrayList<>(order);
    }

    public int indexOf(Node v) {
        return order.indexOf(v);
    }

    public boolean isAssociated(Node x, Node y) {
        return test.isDependent(x, y);
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
    }

    private double score(Node n, Set<Node> pi) {
        if (cache.containsKey(new ScoreKey(n, pi))) {
            return cache.get(new ScoreKey(n, pi));
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double v = -this.score.localScore(variables.indexOf(n), parentIndices);

        cache.put(new ScoreKey(n, pi), v);

        return v;
    }

    private List<Node> getPrefix(int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    private void recalculate(int x) {
        double score = getGrowShrink(x);
        scores[x] = score;
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

        Set<Node> pi = new HashSet<>();
        boolean changed = true;
        double s_new = POSITIVE_INFINITY;
        double s_node = score(n, new HashSet<>());

        // Grow-shrink
        while (changed) {
            changed = false;

            // Let z be the node that maximizes the score...
            Node z = null;

            for (Node z0 : new HashSet<>(getPrefix(p))) {
                if (pi.contains(z0)) continue;
                if (!associates.get(n).contains(z0)) continue;
                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;
                pi.add(z0);
                double s2 = score(n, pi);

                if (s2 < s_new) {
                    s_new = s2;
                    z = z0;
                }

                pi.remove(z0);
            }

            if (s_new < s_node) {
                pi.add(z);
                s_node = s_new;
                changed = true;
            }

            if (z != null) {
                pi.add(z);
                s_node = s_new;
                changed = true;
                boolean changed2 = true;

                while (changed2) {
                    changed2 = false;

                    Node w = null;

                    for (Node z0 : new HashSet<>(pi)) {
                        pi.remove(z0);

                        double s2 = score(n, pi);

                        if (s2 < s_new) {
                            s_new = s2;
                            w = z0;
                        }

                        pi.add(z0);
                    }

                    if (w != null) {
                        pi.remove(w);
                        s_node = s_new;
                        changed2 = true;
                    }
                }
            }

        }

        if (n == nodeMoving) {
            movingScores[p] = s_node;
        }

        return s_node;
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

    public Map<Node, Set<Node>> getAssociates(List<Node> variables) {
        Map<Node, Set<Node>> associates = new HashMap<>();

        for (Node v : variables) {
            Set<Node> nodes = new HashSet<>();

            for (Node w : variables) {
                if (isAssociated(w, v)) {
                    nodes.add(w);
                }
            }

            associates.put(v, nodes);
        }

        return associates;
    }
}
