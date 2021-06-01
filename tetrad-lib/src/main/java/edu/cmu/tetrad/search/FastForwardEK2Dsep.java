package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.Scorer;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.POSITIVE_INFINITY;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FFS) algorithm.
 *
 * @author josephramsey
 */
public class FastForwardEK2Dsep implements FastForward {

    // The score used (default FML BIC). Lower is better.
    private final Score _score;
    private final Scorer scorer;
    private double score = Double.NaN;

    /**
     * Constructs a FFS search
     *
     * @param score the scorer used, by default FML BIC (for linear models). The score
     *              in general should be lower for better models.
     */
    public FastForwardEK2Dsep(Score score, Scorer scorer) {
        this._score = score;
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        List<Node> variables = _score.getVariables();
        Graph G1 = new EdgeListGraph(variables);

        double s = 0;

        for (int i = 0; i < order.size(); i++) {
            Node y = order.get(i);

            List<Node> pi = new ArrayList<>();
            double s_node = 0;
            boolean changed = true;

            while (changed) {
                changed = false;

                // Let x be the node that maximizes the score...
                Node x = null;
                double s_new = POSITIVE_INFINITY;

                for (int j = 0; j < i; j++) {
                    Node x0 = order.get(j);
                    if (pi.contains(x0)) continue;
                    double s2 = s_node + diff(variables, x0, y, pi);

                    if (s2 < s_new) {
                        s_new = s2;
                        x = x0;
                    }
                }

                if (s_new < s_node) {
                    s_node = s_node + diff(variables, x, y, pi);
                    pi.add(x);
                    G1.addDirectedEdge(x, y);
                    changed = true;
                }

                for (Node x0 : new ArrayList<>(pi)) {
                    pi.remove(x0);
                    double s2 = s_node - diff(variables, x0, y, pi);

                    if (s2 < s_new) {
                        s_new = s2;
                        x = x0;
                    }

                    pi.add(x0);
                }

                if (s_new < s_node) {
                    pi.remove(x);
                    s_node = s_node - diff(variables, x, y, pi);
                    G1.addDirectedEdge(x, y);
                    changed = true;
                }
            }

            s += s_node;
        }

        score = scorer.score(G1);
        return G1;
    }

    public Graph search2(List<Node> order) {
        List<Node> variables = _score.getVariables();
        Graph G1 = new EdgeListGraph(variables);

        double s = 0;

        for (int i = 0; i < order.size(); i++) {
            Node y = order.get(i);

            List<Node> pi = new ArrayList<>();
            double s_node = 0;
            boolean changed = true;

            while (changed) {
                changed = false;

                // Let x be the node that maximizes the score...
                Node x = null;
                double s_new = POSITIVE_INFINITY;

                for (int j = 0; j < i; j++) {
                    Node x0 = order.get(j);
                    if (pi.contains(x0)) continue;
                    double s2 = s_node + diff(variables, x0, y, pi);

                    if (s2 < s_new) {
                        s_new = s2;
                        x = x0;
                    }
                }

                if (s_new < s_node) {
                    s_node = s_node + diff(variables, x, y, pi);
                    pi.add(x);
                    G1.addDirectedEdge(x, y);
                    changed = true;
                }

                for (Node x0 : new ArrayList<>(pi)) {
                    pi.remove(x0);
                    double s2 = s_node - diff(variables, x0, y, pi);

                    if (s2 < s_new) {
                        s_new = s2;
                        x = x0;
                    }

                    pi.add(x0);
                }

                if (s_new < s_node) {
                    pi.remove(x);
                    s_node = s_node - diff(variables, x, y, pi);
                    G1.addDirectedEdge(x, y);
                    changed = true;
                }
            }

            s += s_node;
        }

        score = s;
        return G1;
    }

    private double diff(List<Node> variables, Node x, Node y, List<Node> pi) {
        int[] parentIndices = new int[pi.size()];

        for (int k = 0; k < pi.size(); k++) {
            parentIndices[k] = variables.indexOf(pi.get(k));
        }

        return -_score.localScoreDiff(variables.indexOf(x), variables.indexOf(y), parentIndices);
    }

    private double score(List<Node> variables, Node n, List<Node> pi) {
        int[] parentIndices = new int[pi.size()];

        for (int k = 0; k < pi.size(); k++) {
            parentIndices[k] = variables.indexOf(pi.get(k));
        }

        return -_score.localScore(variables.indexOf(n), parentIndices);
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }
}
