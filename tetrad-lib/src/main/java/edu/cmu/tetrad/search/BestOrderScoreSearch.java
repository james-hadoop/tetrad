package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;


/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {
    private final Score score;
    private final List<Node> variables;
    private IKnowledge knowledge = new Knowledge2();
    private final Map<Node, Set<Node>> pis = new HashMap<>();
    private final Map<Subproblem, Double> allScores = new HashMap<>();

    public BestOrderScoreSearch(Score score) {
        this.score = score;
        this.variables = score.getVariables();
    }

    public Graph search(List<Node> initialOrder) {
        List<Node> order = getBestOrder(initialOrder);

        Graph G1 = new EdgeListGraph(order);

        for (Node n : order) {
            for (Node z : pis.get(n)) {
                G1.addDirectedEdge(z, n);
            }
        }

        return SearchGraphUtils.patternForDag(G1);
    }

    public List<Node> getBestOrder(List<Node> initialOrder) {
        long start = System.currentTimeMillis();
        LinkedList<Node> order = new LinkedList<>(initialOrder);

        // Modeled on an insertion sort. Each insertions a BIC comparison of models over the same
        // variables.
        for (int i = 0; i < order.size(); i++) {

            // Pick a various nodes beyond i and move it up into the i.
            for (int j = i; j < order.size(); j++) {
                Node v = order.get(j);
                move(order, v, 0);
                double[] scores = new double[i];

                for (int k = 0; k < i; k++) {
                    scores[k] = getScore(order, k);
                }

                double min = sumScores(scores);
                int b = -0;

                for (int l = 1; l < i; l++) {
                    swap(order, l, scores);
                    double ss = sumScores(scores);

                    if (ss < min) {
                        b = l;
                        min = ss;
                    }
                }

                move(order, v, b);
            }
        }

        long stop = System.currentTimeMillis();
        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
        System.out.println("Final order = " + order);

        return order;
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }


    private void move(LinkedList<Node> br, Node v, int b) {
        br.remove(v);
        br.add(b, v);
    }

    private double sumScores(double[] scores) {
        double s2 = 0;
        for (double _s : scores) s2 += _s;
        return s2;
    }

    private double getScore(LinkedList<Node> br, int p) {
        Node n = br.get(p);
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
                boolean changed2 = true;

                while (changed2) {
                    changed2 = false;

                    for (Node z0 : new HashSet<>(pi)) {
                        if (z0 == z) continue;
                        pi.remove(z0);

                        double s2 = score(n, pi);

                        if (s2 < s_new) {
                            s_new = s2;
                            z = z0;
                        }

                        pi.add(z0);
                    }

                    if (s_new < s_node) {
                        pi.remove(z);
                        s_node = s_new;
                        changed2 = true;
                    }
                }
            }
        }

        pis.put(n, pi);

        return s_node;
    }

    private void swap(LinkedList<Node> br, int ww, double[] scores) {
        move(br, br.get(ww - 1), ww);
        recalculate(ww - 1, scores, br);
        recalculate(ww, scores, br);
    }

    private void recalculate(int x, double[] scores, LinkedList<Node> br) {
        double score = getScore(br, x);
        scores[x] = score;
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

    private double score(Node n, Set<Node> pi) {
        if (allScores.containsKey(new Subproblem(n, pi))) {
            return allScores.get(new Subproblem(n, pi));
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double score = -this.score.localScore(variables.indexOf(n), parentIndices);

        allScores.put(new Subproblem(n, pi), score);

        return score;
    }
}
