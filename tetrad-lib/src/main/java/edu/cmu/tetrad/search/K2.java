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
    private final Map<Node, Set<Node>> predecessors = new HashMap<>();
    private final Map<Node, Set<Node>> previousPis = new HashMap<>();
    private final Map<List<Node>, Double> permutationScores = new WeakHashMap<>();
    private final Map<ScoreSpec, Double> predecessorScores = new WeakHashMap<>();
    private final Map<ScoreSpec, Set<Node>> predecessorPis = new WeakHashMap<>();
    private final Map<ScoreSpec, Double> scores = new WeakHashMap<>();

    /**
     * Constructs a FFS search
     *
     * @param score the score used. A score that works well with FGES (GES) will do fine.
     */
    public K2(Score score) {
        this._score = score;
    }

    public Graph search(List<Node> order) {
        List<Node> variables = _score.getVariables();
        double score = 0;

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Set<Node> predecessors = predecessors(order, n);
            Set<Node> previousPredecessors = this.predecessors.get(n);
            Set<Node> previousPi = previousPis.get(n);

            double s_node = score(variables, n, new HashSet<>());

            List<Node> _predecessors = new ArrayList<>(predecessors);
            _predecessors.add(0, n);

            if (predecessorScores.get(new ScoreSpec(n, new HashSet<>(_predecessors))) != null) {
                double _score1 = predecessorScores.get(new ScoreSpec(n, new HashSet<>(_predecessors)));

                if (_score1 < score) {
                    s_node = _score1;
                }
            } else {

//            if (previousPredecessors != null && predecessors.equals(previousPi)
//                    && previousPredecessors.equals(predecessors)) {
//                double score1 = score(variables, n, previousPi);
//
//                if (score1 < s_node) {
//                    s_node = score1;
//                }
//            } else {

                Set<Node> pi = new HashSet<>();
                boolean changed = true;
                boolean add;

                if (previousPredecessors != null) {
                    // If adding, add = true
                    // If removing, add = false
                    // If mixture, add = true but forget previous pi's (i.e. build from empty)

                    if (predecessors.containsAll(previousPredecessors)) {
                        add = true;
                    } else if (previousPredecessors.containsAll(predecessors)) {
                        add = false;
                    } else {
                        add = true;
                        previousPis.put(n, null);
                    }
                } else {
                    add = true;
                    previousPis.put(n, null);
                }

                double s_new = POSITIVE_INFINITY;

                if (previousPis.get(n) != null) {
                    previousPis.get(n).retainAll(predecessors);

                    while (changed) {
                        changed = false;

                        Node z = null;

                        {
                            for (Node z0 : previousPis.get(n)) {
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
                        for (Node z0 : predecessors) {
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

                this.predecessors.put(n, predecessors);
                this.previousPis.put(n, pi);

                this.predecessorScores.put(new ScoreSpec(n, new HashSet<>(_predecessors)), s_node);
                this.predecessorPis.put(new ScoreSpec(n, new HashSet<>(_predecessors)), pi);
            }

            score += s_node;
        }

        Graph G1 = new EdgeListGraph(variables);

        for (Node n : order) {
            for (Node z : previousPis.get(n)) {
                G1.addDirectedEdge(z, n);
            }
        }

        return G1;
    }

    private Set<Node> getNodes(Set<Node> predecessors) {
        return predecessors;
    }

    public double score(List<Node> order) {
        Double score0 = permutationScores.get(order);

        if (score0 != null) {
            return score0;
        }

        List<Node> variables = _score.getVariables();
        double score = 0;

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Set<Node> predecessors = predecessors(order, n);
            Set<Node> previousPredecessors = this.predecessors.get(n);
            Set<Node> previousPi = previousPis.get(n);

            double s_node = score(variables, n, new HashSet<>());

            if (previousPredecessors != null && predecessors.equals(previousPi)
                    && previousPredecessors.equals(predecessors)) {
                double score1 = score(variables, n, previousPi);

                if (score1 < s_node) {
                    s_node = score1;
                }
            } else {

                Set<Node> pi = new HashSet<>();
                boolean changed = true;
                boolean add;

                if (previousPredecessors != null) {
                    // If adding, add = true
                    // If removing, add = false
                    // If mixture, add = true but forget previous pi's (i.e. build from empty)

                    if (predecessors.containsAll(previousPredecessors)) {
                        add = true;
                    } else if (previousPredecessors.containsAll(predecessors)) {
                        add = false;
                    } else {
                        add = true;
                        previousPis.put(n, null);
                    }
                } else {
                    add = true;
                    previousPis.put(n, null);
                }

                double s_new = POSITIVE_INFINITY;

                if (previousPis.get(n) != null) {
                    previousPis.get(n).retainAll(predecessors);

                    while (changed) {
                        changed = false;

                        Node z = null;

                        {
                            for (Node z0 : previousPis.get(n)) {
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
                        for (Node z0 : predecessors) {
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

                this.predecessors.put(n, predecessors);
                this.previousPis.put(n, pi);
            }

            score += s_node;
        }

        this.permutationScores.put(new ArrayList<>(order), score);
        return score;
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

    private Set<Node> predecessors(List<Node> order, Node n) {
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
}
