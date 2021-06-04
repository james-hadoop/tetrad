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
public class K2 implements FastForward {

    // The score used.
    private final Score _score;
    private final Map<Node, Set<Node>> families = new HashMap<>();
    private final Map<Node, Set<Node>> familyPis = new HashMap<>();
    private final Map<ScoreSpec, Double> scores = new WeakHashMap<>();
    private final IndTestScore test;

    /**
     * Constructs a FFS search
     *
     * @param score the score used. A score that works well with FGES (GES) will do fine.
     */
    public K2(Score score) {
        this._score = score;
        this.test = new IndTestScore(score);
    }

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

    public double score(List<Node> order) {
        ScoreResult result = scoreResult(order);
        return result.getScore();
    }

    @Override
    public boolean isAssociated(Node w, Node v) {
        return test.isDependent(w, v);
    }

    private ScoreResult scoreResult(List<Node> order) {
        List<Node> variables = _score.getVariables();
        double score = 0;

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Set<Node> families = family(order, n);
            Set<Node> previousFamilies = this.families.get(n);
            Set<Node> previousFamilyPis = familyPis.get(n);

            double s_node = score(variables, n, new HashSet<>());

            if (previousFamilies != null && families.equals(previousFamilyPis)
                    && previousFamilies.equals(families)) {
                double score1 = score(variables, n, previousFamilyPis);

                if (score1 < s_node) {
                    s_node = score1;
                }
            } else {

                Set<Node> pi = new HashSet<>();
                boolean changed = true;
                boolean add;

                if (previousFamilies != null) {
                    if (families.containsAll(previousFamilies)) {
                        add = true;
                    } else if (previousFamilies.containsAll(families)) {
                        add = false;
                    } else {
                        add = true;
                        familyPis.put(n, null);
                    }
                } else {
                    add = true;
                    familyPis.put(n, null);
                }

                double s_new = POSITIVE_INFINITY;

                if (familyPis.get(n) != null) {
                    familyPis.get(n).retainAll(families);

                    while (changed) {
                        changed = false;

                        Node z = null;

                        {
                            for (Node z0 : familyPis.get(n)) {
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
                        }
                    }
                }

                changed = true;

                while (changed) {
                    changed = false;

                    // Let z be the node that maximizes the score...
                    Node z = null;

                    if (add) {
                        for (Node z0 : families) {
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
                    }

                    if (!add) {
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
                }

                this.families.put(n, families);
                this.familyPis.put(n, pi);
            }

            score += s_node;
        }

        List<Set<Node>> pis = new ArrayList<>();

        for (Node n : order) {
            Set<Node> pi = familyPis.get(n);
            pis.add(pi);
        }

        return new ScoreResult(pis, score);
    }

    private double score(List<Node> variables, Node n, Set<Node> pi) {
        ScoreSpec key = new ScoreSpec(n, pi);

        Double score = scores.get(key);

        if (score != null) {
            return score;
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        score = -_score.localScore(variables.indexOf(n), parentIndices);

        key = new ScoreSpec(n, new HashSet<>(pi));
        scores.put(key, score);

        return score;
    }

    private Set<Node> family(List<Node> order, Node n) {
        Set<Node> _predecessors = new HashSet<>();
        for (int j = 0; j < order.indexOf(n); j++) _predecessors.add(order.get(j));
        return _predecessors;
    }

    private static class ScoreSpec {
        private final Node y;
        private final Set<Node> p;

        public ScoreSpec(Node y, Set<Node> p) {
            this.y = y;
            this.p = p;
        }

        public int hashCode() {
            return 17 * y.hashCode() + 3 * p.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof ScoreSpec)) {
                return false;
            }

            ScoreSpec spec = (ScoreSpec) o;
            return spec.y.equals(this.y) && spec.p.equals(this.p);
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
