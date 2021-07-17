package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeEqualityMode;

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
    private final Map<Integer, List<Node>> bookmarkedOrder = new HashMap<>();
    private final Map<Integer, Pair[]> bookmarkedScores = new HashMap<>();
    private Score score;
    private IndependenceTest test;
    private List<Node> variables;
    private LinkedList<Node> order;
    private Pair[] scores;
    private boolean cachingScores = true;
    private IKnowledge knowledge = new Knowledge2();
    private List<Set<Node>> prefixes;

    public TeyssierScorer(Score score) {
        this.score = score;
        this.order = new LinkedList<>(score.getVariables());

//        NodeEqualityMode.setEqualityMode(NodeEqualityMode.Type.OBJECT);
    }

    public TeyssierScorer(IndependenceTest test) {
        this.test = test;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double score(List<Node> order) {
        this.order = new LinkedList<>(order);

        this.variables = score != null ? score.getVariables() : test.getVariables();

        this.scores = new Pair[order.size()];
        this.prefixes = new ArrayList<>();
        for (int i = 0; i < order.size(); i++) prefixes.add(null);
        initializeScores();
        return score();
    }

    public double score() {
        return sum();
    }

    public boolean promote(Node v) {
        int index = order.indexOf(v);
        if (!(index >= 1 && index <= order.size() - 1)) return false;

        Node v1 = order.get(index - 1);
        Node v2 = order.get(index);

        order.set(index - 1, v2);
        order.set(index, v1);

        recalculate(index - 1);
        recalculate(index);

//        System.out.println("Promoting " + v + " order = " + order + " # edges = " + getNumEdges());

        return true;
    }

    public boolean demote(Node v) {
        int index = order.indexOf(v);
        if (!(index >= 0 && index <= order.size() - 2)) return false;

        Node v1 = order.get(index);
        Node v2 = order.get(index + 1);

        order.set(index, v2);
        order.set(index + 1, v1);

        recalculate(index);
        recalculate(index + 1);

//        System.out.println("\tDeomoting " + v + " order = " + order + " # edges = " + getNumEdges());

        return true;
    }

    public void moveTo(Node v, int toIndex) {
        order.remove(v);
        order.add(toIndex, v);

        for (int i = 0; i < order.size(); i++) {
            recalculate(i);
        }

//        int fromIndex = order.indexOf(v);
//
//        while (toIndex < fromIndex) {
//            if (!promote(v)) break;
//            fromIndex--;
//        }
//
//        while (toIndex > fromIndex) {
//            if (!demote(v)) break;
//            fromIndex++;
//        }
    }

    public void swap(Node m, Node n) {
        int i = order.indexOf(m);
        int j = order.indexOf(n);

        moveTo(m, j);
        moveTo(n, i);
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
        scores = new Pair[order.size()];
        bookmarkedScores.clear();
        bookmarkedOrder.clear();

        for (int i = 0; i < order.size(); i++) {
            recalculate(i);
        }
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
        if (!new HashSet<>(getPrefix(p)).equals(prefixes.get(p))) {
            scores[p] = getGrowShrink(p);
        }
    }

    public Set<Node> getMb(int p) {
        return scores[p].getMb();
    }

    private Pair getGrowShrink(int p) {
        if (false) {
            if (test != null) {
                return getPearlIndep(p);
            } else {
                return getPearlScore(p);
            }
        }else {
            if (test != null) {
                return getGrowShrinkIndep(p);
            } else {

                Node n = order.get(p);

                Set<Node> mb = new HashSet<>();
                boolean changed = true;

                double sMax = score(n, new HashSet<>());
                List<Node> prefix = getPrefix(p);

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
            }
        }
    }

    private Pair getGrowShrinkIndep(int p) {
        Node n = order.get(p);

        List<Node> mb = new ArrayList<>();
        boolean changed = true;

        List<Node> prefix = getPrefix(p);

        // Grow-shrink
        while (changed) {
            changed = false;

            for (Node z0 : prefix) {
                if (mb.contains(z0)) continue;

                if (test.isDependent(n, z0, mb)) {
                    mb.add(z0);
                    changed = true;
                }

                boolean changed2 = true;

                while (changed2) {
                    changed2 = false;

                    for (Node z1 : new ArrayList<>(mb)) {
                        mb.remove(z1);

                        if (test.isIndependent(n, z1, mb)) {
                            changed2 = true;
                        } else {
                            mb.add(z1);
                        }
                    }
                }
            }
        }

        return new Pair(new HashSet<>(mb), mb.size());
    }

    private Pair getPearlScore(int p) {
        Node n = order.get(p);

        List<Node> prefix = getPrefix(p);
        prefix.remove(n);
        Set<Node> parents = new HashSet<>();

        for (Node z : prefix) {
            List<Node> w = new ArrayList<>(prefix);
            w.remove(z);

            double s2 = score(n, parents);

            w.remove(z);

            double s1 = score(n, parents);

            if (s2 <= s1) {
                parents.add(z);
            }
        }

        return new Pair(parents, parents.size());
    }

    private Pair getPearlIndep(int p) {
        Node n = order.get(p);

        List<Node> prefix = getPrefix(p);
        prefix.remove(n);
        Set<Node> parents = new HashSet<>();

        for (Node z : prefix) {
            List<Node> w = new ArrayList<>(prefix);
            w.remove(z);

            if (test.isDependent(n, z, w)) {
                parents.add(z);
            }
        }

        return new Pair(parents, parents.size());
    }

    public synchronized void bookmark(int index) {
        if (order == null) throw new IllegalArgumentException("Need to score an initial order first.");

        bookmarkedOrder.put(index, new ArrayList<>(order));
        bookmarkedScores.put(index, Arrays.copyOf(scores, scores.length));
    }

    public synchronized void flipToBookmark(int index) {
        order = new LinkedList<>(bookmarkedOrder.get(index));
        scores = Arrays.copyOf(bookmarkedScores.get(index), bookmarkedScores.get(index).length);
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

    public void shuffleVariables() {
        order = new LinkedList<>(order);
        shuffle(order);
        score(order);
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
        private final Set<Node> pi;

        public ScoreKey(Node y, Set<Node> pi) {
            this.y = y;
            this.pi = new HashSet<>(pi);
        }

        public int hashCode() {
            return y.hashCode() + 11 * pi.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof ScoreKey)) {
                return false;
            }

            ScoreKey spec = (ScoreKey) o;
            return spec.y.equals(this.y) && spec.pi.equals(this.pi);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPi() {
            return pi;
        }
    }
}
