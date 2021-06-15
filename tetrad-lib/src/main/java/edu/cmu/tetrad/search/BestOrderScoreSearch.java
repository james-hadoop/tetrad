package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;
import static java.util.Collections.shuffle;


/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {
    private final Score score;

    public BestOrderScoreSearch(Score score) {
        this.score = score;
    }

    public Graph search(List<Node> initialOrder, int method) {
        List<Node> order;

        switch (method) {
            case 1:
                order = getBestOrdering1(initialOrder);
                break;
            case 2:
                order = getBestOrdering2(initialOrder);
                break;
            default:
                throw new IllegalArgumentException("Expecting 1 or 2.");
        }

        K3 k3 = new K3(score);

        List<Set<Node>> pis = k3.scoreResult(order, true).getPis();


        Graph G1 = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            for (Node z : pis.get(i)) {
                G1.addDirectedEdge(z, order.get(i));
            }
        }

        return SearchGraphUtils.patternForDag(G1);
    }

    // Using K3.
    public List<Node> getBestOrdering1(List<Node> variables) {
        K3 k3 = new K3(score);

        long start = System.currentTimeMillis();

        List<Node> b0 = new ArrayList<>(variables);
        double s0 = k3.score(b0);

        System.out.println("Original order = " + b0);

        boolean changed = true;

        int count = 0;

        List<Node> br = new ArrayList<>(b0);

        Map<Node, Set<Node>> associates = getAssociates(variables, k3);

        br = new ArrayList<>(br);
        double sr = k3.score(br);

        for (int r = 0; r < 1; r++) {

            // Until you can't do it anymore...
            while (changed) {
                changed = false;

                List<Node> b2 = new ArrayList<>(br);
                double s2 = sr;

                // ...pick a variable v...
                for (Node v : variables) {

                    // ...and move v to an optimal position.
                    for (Node w : associates.get(v)) {
                        int j = b2.indexOf(w);
                        List<Node> b1 = new ArrayList<>(b2);
                        b1.remove(v);
                        b1.add(j, v);

                        count++;

                        double s1 = k3.score(b1);

                        if (s1 < s2) {
                            s2 = s1;
                            b2 = b1;
                        }
                    }
                }

                // Output the best order you find.
                if (s2 < sr) {
                    br = b2;
                    sr = s2;
                    changed = true;
                }
            }

            if (sr < s0) {
                s0 = sr;
                b0 = new ArrayList<>(br);

                System.out.println("Updated order = " + b0 + " count = " + count);
            }

            shuffle(br);
            sr = k3.score(br);
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return b0;
    }

    // Using Teyssier and Kohler's neighor swaps.
    public List<Node> getBestOrdering2(List<Node> initialOrder) {
        long start = System.currentTimeMillis();

        List<Node> br = new ArrayList<>(initialOrder);
        Scorer scorer = new Scorer(score, br);

        double best = POSITIVE_INFINITY;

        // Pick a variable v, move it to the first position and then score all
        // subsequent positions and pick the best one. Pick another variable v and do the same.
        // Etc. until you can't move v's any longer.
        for (int r = 0; r < 10; r++) {
            for (Node node : initialOrder) {
                scorer.moveTo(node, 0);
                double bestScore = scorer.score();
                int bestIndex = scorer.indexOf(node);
                boolean moved;

                do {
                    moved = scorer.moveRight(node);

                    if (scorer.score() < bestScore) {
                        bestScore = scorer.score();
                        bestIndex = scorer.indexOf(node);
                    }
                } while (moved);

                scorer.moveTo(node, bestIndex);
            }

            if (scorer.score() < best) {
                best = scorer.score();
            } else {
                break;
            }
        }

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return scorer.getOrder();
    }

    private Map<Node, Set<Node>> getAssociates(List<Node> variables, K3 k3) {
        Map<Node, Set<Node>> associates = new HashMap<>();

        for (Node v : variables) {
            Set<Node> nodes = new HashSet<>();

            for (Node w : variables) {
                if (k3.isAssociated(w, v)) {
                    nodes.add(w);
                }
            }

            associates.put(v, nodes);
        }

        return associates;
    }


    public static class Subproblem {
        private final Node y;
        private final Set<Node> prefix;

        public Subproblem(Node y, Set<Node> prefix) {
            this.y = y;
            this.prefix = new HashSet<>(prefix);
        }

        public int hashCode() {
            return 17 * y.hashCode() + 3 * prefix.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof Subproblem)) {
                return false;
            }

            Subproblem spec = (Subproblem) o;
            return spec.y.equals(this.y) && spec.getPrefix().equals(this.prefix);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPrefix() {
            return prefix;
        }
    }

    public static class Scorer {
        private final Score score;
        private final K3 k3;
        private final Map<Subproblem, Double> allScores = new HashMap<>();
        List<Node> variables;
        double[] scores;
        Map<Node, Set<Node>> pis = new HashMap<>();
        private List<Node> order;
        private IKnowledge knowledge = new Knowledge2();
        private Node nodeMoving = null;
        private double[] movingScores = null;

        public Scorer(Score score, List<Node> order) {
            this.score = score;
            this.k3 = new K3(score);
            this.variables = score.getVariables();
            this.scores = new double[variables.size()];
            this.order = order;
        }

        public void setKnowledge(IKnowledge knowledge) {
            this.knowledge = knowledge;
        }

        public double score(List<Node> order) {
            setOrder(order);
            return score2();
        }

        public double score() {
            return score2();
        }

        public double score1() {
            return k3.score(order);
        }

        public double score2() {
            return sum();
        }

        private double sum() {
            double score = 0;

            for (int i = 0; i < order.size(); i++) {
                score += scores[i];
            }

            return score;
        }

        private void initializeScores() {
            for (int i = 0; i < order.size(); i++) {
                double nodeScore = getScore(order, i);
                scores[i] = nodeScore;
            }
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

        public void moveToFirst(Node v) {
            boolean moved;

            do {
                moved = moveLeft(v);
            } while (moved);
        }

        public boolean moveTo(Node v, int i) {
            int index = order.indexOf(v);
            if (index == i) return false;

            if (i < index) {
                while (--index >= i) {
                    if (!moveLeft(v)) break;
                }
            } else {
                while (++index <= i) {
                    if (!moveRight(v)) break;
                }
            }

            return true;
        }

        private void recalculate(int x) {
            double score = getScore(order, x);
            scores[x] = score;
        }

        private double getScore(List<Node> br, int p) {
            Node n = br.get(p);

            if (n != nodeMoving) {
                nodeMoving = n;
                movingScores = new double[order.size()];
                Arrays.fill(movingScores, Double.NaN);
            } else if (!Double.isNaN(movingScores[p])) {
                return movingScores[p];
            }

            Set<Node> pi = new HashSet<>();
            boolean changed = true;
            double s_new = POSITIVE_INFINITY;
            double s_node = score(n, new HashSet<>());

            while (changed) {
                changed = false;

                // Let z be the node that maximizes the score...
                Node z = null;

                for (Node z0 : new HashSet<>(getPrefix(br, p))) {
                    if (pi.contains(z0)) continue;
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
//                    boolean changed2 = true;
//
//                    while (changed2) {
//                        changed2 = false;
//
//                        for (Node z0 : new HashSet<>(pi)) {
//                            if (z0 == z) continue;
//                            pi.remove(z0);
//
//                            double s2 = score(n, pi);
//
//                            if (s2 < s_new) {
//                                s_new = s2;
//                                z = z0;
//                            }
//
//                            pi.add(z0);
//                        }
//
//                        if (s_new < s_node) {
//                            pi.remove(z);
//                            s_node = s_new;
//                            changed2 = true;
//                        }
//                    }
                }
            }

            pis.put(n, pi);

            if (n == nodeMoving) {
                movingScores[order.indexOf(n)] = s_node;
            }

            return s_node;
        }

        private double score(Node n, Set<Node> pi) {
            if (allScores.containsKey(new Subproblem(n, pi))) {
                return allScores.get(new Subproblem(n, pi));
            }

            int[] parentIndices = new int[pi.size()];

            int k = 0;

            for (Node p : pi) {
                parentIndices[k++] = variables.indexOf(p);
            }

            double v = -this.score.localScore(variables.indexOf(n), parentIndices);

            allScores.put(new Subproblem(n, pi), v);

            return v;
        }

        public List<Node> getOrder() {
            return new ArrayList<>(order);
        }

        public void setOrder(List<Node> order) {
            this.order = new ArrayList<>(order);
            initializeScores();
        }

        public List<Node> getPrefix(List<Node> order, int i) {
            List<Node> prefix = new ArrayList<>();
            for (int j = 0; j < i; j++) {
                prefix.add(order.get(j));
            }
            return prefix;
        }

        public int indexOf(Node v) {
            return order.indexOf(v);
        }
    }
}
