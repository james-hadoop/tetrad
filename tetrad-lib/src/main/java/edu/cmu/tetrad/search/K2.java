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
    private double score = Double.NaN;
    private final Map<Record, Double> hashedScores = new WeakHashMap<>();



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
        this.score = 0;

        for (int i = 0; i < order.size(); i++) {
            Node n = order.get(i);
            Set<Node> predecessors = predecessors(order, n);

            Set<Node> pi = new HashSet<>();
            boolean changed = true;
            double s_node = score(variables, n, new HashSet<>());

            Set<Node> possibleParents = predecessors;
            boolean add = true;

            if (this.predecessors.get(n) != null) {
                if (predecessors.containsAll(this.predecessors.get(n))) {
//                    possibleParents = new HashSet<>(predecessors);
//                    possibleParents.retainAll(this.predecessors.get(n));
                    add = true;
                } else if (this.predecessors.get(n).containsAll(predecessors)) {
                    add = false;
                } else {
                    previousPis.put(n, null);
                }
            }

            double s_new = POSITIVE_INFINITY;

            if (previousPis.get(n) != null) {
                previousPis.get(n).retainAll(predecessors);

//                pi = previousPis.get(n);


                while (changed) {
                    changed = false;

                    Node z = null;

                    {
                        for (Node z0 : previousPis.get(n)) {
//                            if (!predecessors.contains(z0)) continue;
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
//            double s_new = POSITIVE_INFINITY;

            while (changed) {
                changed = false;

                // Let z be the node that maximizes the score...
                Node z = null;

                if (add) {
                    for (Node z0 : possibleParents) {
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

    private double score(List<Node> variables, Node n, Set<Node> pi) {
//        Record key = new Record(n, pi);
//        Double value = hashedScores.get(key);
//
//        if (value != null) {
//            return value;
//        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double score = -_score.localScore(variables.indexOf(n), parentIndices);

//        hashedScores.put(key, score);

        return score;
    }

    private Set<Node> predecessors(List<Node> order, Node n) {
        Set<Node> _predecessors = new HashSet<>();
        for (int j = 0; j < order.indexOf(n); j++) _predecessors.add(order.get(j));
        return _predecessors;
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }

    private static class Record {
        private final Node n;
        private final Set<Node> parents;

        public Record(Node n, Set<Node> parents) {
            this.n = n;
            this.parents = parents;
        }

        public int hashCode() {
            return n.hashCode() + 5 * parents.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof Record)) {
                return false;
            }

            Record r = (Record) o;

            return r.n.equals(this.n) && r.parents.equals(this.parents);
        }

        public Node getN() {
            return n;
        }

        public Set<Node> getParents() {
            return parents;
        }
    }
}
