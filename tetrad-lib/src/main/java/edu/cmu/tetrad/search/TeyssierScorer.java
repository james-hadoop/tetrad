package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Node;

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
    private Node[] bookmarkedOrder = null;
    private Pair[] bookmarkedScores = null;
    private Score score;
    private IndependenceTest test;
    private final List<Node> variables;
    private Node[] order;
    private Pair[] scores;
    private boolean cachingScores = true;
    private IKnowledge knowledge = new Knowledge2();
    private LinkedList<Set<Node>> prefixes;

    public TeyssierScorer(Score score) {
        this.variables = score.getVariables();
        this.score = score;
        score(variables);
    }

    public TeyssierScorer(IndependenceTest test) {
        this.variables = score.getVariables();
        this.test = test;
        score(variables);
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double score(List<Node> order) {
        this.order = new Node[order.size()];
        for (int i = 0; i < order.size(); i++) this.order[i] = order.get(i);
        initializeScores();
        return score();
    }

    private void initializeScores() {
        initializeArrays();

        for (int i = 0; i < order.length; i++) {
            recalculate(i);
        }
    }

    private void initializeArrays() {
        this.scores = new Pair[order.length];

        this.prefixes = new LinkedList<>();
        for (int i = 0; i < order.length; i++) this.prefixes.add(null);

        this.bookmarkedOrder = new Node[order.length];
        setNodes(this.bookmarkedOrder, this.order);

        this.bookmarkedScores = new Pair[order.length];
        setPairs(this.bookmarkedScores, this.scores);

        this.prefixes = new LinkedList<>();
        for (int i = 0; i < order.length; i++) this.prefixes.add(null);
    }

    public double score() {
        return sum();
    }

    public boolean promote(Node v) {
        int index = index(order, v);

        if (!(index >= 1 && index <= order.length - 1)) return false;

        Node v1 = order[index - 1];
        Node v2 = order[index];

        order[index - 1] = v2;
        order[index] = v1;

        recalculate(index - 1);
        recalculate(index);

        return true;
    }

    public boolean demote(Node v) {
        int index = index(order, v);

        if (!(index >= 1 && index <= order.length - 1)) return false;

        Node v1 = order[index];
        Node v2 = order[index + 1];

        order[index] = v2;
        order[index + 1] = v1;

        recalculate(index);
        recalculate(index + 1);

        return true;
    }

    public void moveTo(Node v, int toIndex) {
        int fromIndex = index(order, v);

        while (toIndex < fromIndex) {
            if (!promote(v)) break;
            fromIndex--;
        }

        while (toIndex > fromIndex) {
            if (!demote(v)) break;
            fromIndex++;
        }
    }

    private int index(Node[] order, Node v) {
        return Arrays.binarySearch(order, v);
    }

    public void moveToFirst(Node v) {
        moveTo(v, 0);
    }

    public void moveToLast(Node v) {
        moveTo(v, size() - 1);
    }

    public void swap(Node m, Node n) {
        int i = index(order, m);
        int j = index(order, n);

        moveTo(m, j);
        moveTo(n, i);
    }

    public List<Node> getOrder() {
        return Arrays.asList(this.order);
    }

    public int indexOf(Node v) {
        return index(order, v);
    }

    private double sum() {
        double score = 0;

        for (int i = 0; i < order.length; i++) {
            score += scores[i].getScore();
        }

        return score;
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
        return new ArrayList<>(Arrays.asList(order).subList(0, i));
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

                Node n = order[p];

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
//                return new Pair(mb, mb.size());
            }
        }
    }

    private Pair getGrowShrinkIndep(int p) {
        Node n = order[p];

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
        Node n = order[p];

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
        Node n = order[p];

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

    public synchronized void bookmark() {
        setNodes(bookmarkedOrder, order);
        setPairs(bookmarkedScores, scores);
    }

    public synchronized void goToBookmark() {
        setNodes(order, bookmarkedOrder);
        setPairs(scores, bookmarkedScores);
    }

    private void setNodes(Node[] order1, Node[] order2) {
        System.arraycopy(order2, 0, order1, 0, order1.length);
    }

    private void setPairs(Pair[] order1, Pair[] order2) {
        System.arraycopy(order2, 0, order1, 0, order1.length);
    }

    public void setCachingScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public int size() {
        return order.length;
    }

    public int getNumEdges() {
        int numEdges = 0;

        for (int p = 0; p < order.length; p++) {
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
