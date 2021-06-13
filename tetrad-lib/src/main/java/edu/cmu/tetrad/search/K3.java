package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

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
    private final Map<Subproblem, List<Node>> subproblemPis = new HashMap<>();

    // Map from subproblem to scores.
    private final Map<Subproblem, Double> subproblemScores = new HashMap<>();

    // All scores of variables given their putative parents.
    private final Map<Subproblem, Double> allScores = new HashMap<>();

    // The variables from the score.
    List<Node> variables;

    // The knowledge; we use forbidden edges from this.
    private IKnowledge knowledge = new Knowledge2();

    /**
     * Constructs a K3 search
     *
     * @param score the score used. A score that works well with FGES (GES) will do fine.
     */
    public K3(Score score) {
        this._score = score;
        this.test = new IndTestScore(score);
        this.variables = _score.getVariables();
    }

    /**
     * Searches for a DAG with the given topological order.
     */
    public Graph search(List<Node> order) {
        ScoreResult result = scoreResult(order);
        List<List<Node>> pis = result.getPis();

        Graph G1 = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            for (Node z : pis.get(i)) {
                G1.addDirectedEdge(z, n);
            }
        }

        return G1;
    }

    @Override
    public double score(List<Node> order) {
        return scoreResult(order).getScore();
    }

    /**
     * Returns the score of a graph found with the given topological order.
     */
    public ScoreResult scoreResult(List<Node> order) {
        double score = 0;

        List<List<Node>> pis = new ArrayList<>();

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Subproblem subproblem = subproblem(order, n);

            double s_node = score(n, new ArrayList<>());
            List<Node> pi;

            if (subproblemPis.containsKey(subproblem)) {
                s_node = subproblemScores.get(subproblem);
                pi = subproblemPis.get(subproblem);
            } else {
                pi = new ArrayList<>();
                boolean changed = true;
                double s_new = POSITIVE_INFINITY;

                while (changed) {
                    changed = false;

                    // Let z be the node that maximizes the score...
                    Node z = null;

                    for (Node z0 : subproblem.getPrefix()) {
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

                this.subproblemPis.put(subproblem, pi);
                this.subproblemScores.put(subproblem, s_node);
            }

            pis.add(pi);

            score += s_node;
        }

        return new ScoreResult(pis, score);
    }

    /**
     * Returns true iff x an y are associated.
     */
    @Override
    public boolean isAssociated(Node x, Node y) {
        return test.isDependent(x, y);
    }

    @Override
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    private double score(Node n, List<Node> pi) {
        if (allScores.containsKey(new Subproblem(n, pi))) {
            return allScores.get(new Subproblem(n, pi));
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double score = -_score.localScore(variables.indexOf(n), parentIndices);

        allScores.put(new Subproblem(n, pi), score);

        return score;
    }

    private Subproblem subproblem(List<Node> order, Node n) {
        List<Node> prefix = getPrefix(order, n);
        return new Subproblem(n, prefix);
    }

    @NotNull
    public List<Node> getPrefix(List<Node> order, Node n) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < order.indexOf(n); j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    private static class Subproblem {
        private final Node y;
        private final Set<Node> prefix;

        public Subproblem(Node y, List<Node> prefix) {
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
            return spec.y.equals(this.y) && spec.prefix.equals(this.prefix);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPrefix() {
            return prefix;
        }
    }

    private static class ScoreResult {
        private final List<List<Node>> pis;
        private final double score;

        public ScoreResult(List<List<Node>> pis, double score) {
            this.pis = pis;
            this.score = score;
        }

        public List<List<Node>> getPis() {
            return pis;
        }

        public double getScore() {
            return score;
        }
    }
}
