package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements an optimization of the K2 algorithm.
 *
 * @author josephramsey
 */
public class K3 implements ForwardScore {

    // The score used.
    private final Score _score;

    // Used to find dependent variables.
    private final IndTestScore test;

    // Map from subproblems to parent sets.
    private final Map<Subproblem, Set<Node>> subproblemPis = new HashMap<>();

    // Map from subproblem to scores.
    private final Map<Subproblem, Double> subproblemScores = new HashMap<>();

    /**
     * Constructs a K3 search
     *
     * @param score the score used. A score that works well with FGES (GES) will do fine.
     */
    public K3(Score score) {
        this._score = score;
        this.test = new IndTestScore(score);
    }

    /**
     * Searches for a DAG with the given topological order.
     */
    public Graph search(List<Node> order) {
        ScoreResult result = scoreResult(order);
        List<Set<Node>> pis = result.getPis();

        Graph G1 = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            for (Node z : pis.get(i)) {
                G1.addDirectedEdge(z, n);
            }
        }

        return G1;
    }

    /**
     * Returns the score of a graph found with the given topological order.
     */
    public double score(List<Node> order) {
        ScoreResult result = scoreResult(order);
        return result.getScore();
    }

    /**
     * Returns true iff x an y are associated.
     */
    @Override
    public boolean isAssociated(Node x, Node y) {
        return test.isDependent(x, y);
    }

    private ScoreResult scoreResult(List<Node> order) {
        List<Node> variables = _score.getVariables();
        double score = 0;
        List<Set<Node>> pis = new ArrayList<>();

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Subproblem subproblem = subproblem(order, n);

            double s_node = score(variables, n, new HashSet<>());
            Set<Node> pi;

            if (subproblemPis.containsKey(subproblem)) {
                s_node = subproblemScores.get(subproblem);
                pi = subproblemPis.get(subproblem);
            } else {
                pi = new HashSet<>();
                boolean changed = true;
                double s_new = POSITIVE_INFINITY;

                while (changed) {
                    changed = false;

                    // Let z be the node that maximizes the score...
                    Node z = null;

                    for (Node z0 : subproblem.getPredeceessors()) {
                        if (pi.contains(z0)) continue;
                        pi.add(z0);

                        double s2 = score(variables, n, pi);

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

                    for (Node z0 : new HashSet<>(pi)) {
                        pi.remove(z0);

                        double s2 = score(variables, n, pi);

                        if (s2 < s_new) {
                            s_new = s2;
                            z = z0;
                        }

                        pi.add(z0);
                    }

                    if (s_new < s_node) {
                        pi.remove(z);
                        s_node = s_new;
                        changed = true;
                    }
                }

                this.subproblemPis.put(subproblem, pi);
                this.subproblemScores.put(subproblem, s_node);
            }

            pis.add(pi);
            score += s_node;
        }

        return new ScoreResult(pis, score);
    }

    private double score(List<Node> variables, Node n, Set<Node> pi) {
        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        return -_score.localScore(variables.indexOf(n), parentIndices);
    }

    private Subproblem subproblem(List<Node> order, Node n) {
        Set<Node> _predecessors = new HashSet<>();
        for (int j = 0; j < order.indexOf(n); j++) _predecessors.add(order.get(j));
        return new Subproblem(n, _predecessors);
    }

    private static class Subproblem {
        private final Node y;
        private final Set<Node> predeceessors;

        public Subproblem(Node y, Set<Node> predecessors) {
            this.y = y;
            this.predeceessors = new HashSet<>(predecessors);
        }

        public int hashCode() {
            return 17 * y.hashCode() + 3 * predeceessors.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof Subproblem)) {
                return false;
            }

            Subproblem spec = (Subproblem) o;
            return spec.y.equals(this.y) && spec.predeceessors.equals(this.predeceessors);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPredeceessors() {
            return predeceessors;
        }
    }

    private static class ScoreResult {
        private final List<Set<Node>> pis;
        private final double score;

        public ScoreResult(List<Set<Node>> pis, double score) {
            this.pis = pis;
            this.score = score;
        }

        public List<Set<Node>> getPis() {
            return pis;
        }

        public double getScore() {
            return score;
        }
    }
}
