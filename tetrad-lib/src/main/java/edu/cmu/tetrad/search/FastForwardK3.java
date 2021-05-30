package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.Scorer;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.NEGATIVE_INFINITY;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FFS) algorithm.
 *
 * @author josephramsey
 */
public class FastForwardK3 implements FastForward {

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
    public FastForwardK3(Score score, Scorer scorer) {
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

        for (Node n : _score.getVariables()) {
            List<Node> pi = new ArrayList<>();
            double s_old = _score.localScore(variables.indexOf(n), new int[0]);
            boolean changed = true;

            while (changed) {
                changed = false;

                // Let z be the node that maximizes the score...
                Node z = null;
                double s_new = NEGATIVE_INFINITY;

                for (Node z0 : variables) {
                    if (z0 == n) continue;

                    if (!(order.indexOf(z0) < order.indexOf(n))) {
                        continue;
                    }

                    if (pi.contains(z0)) continue;
                    pi.add(z0);

                    int[] parentIndices = new int[pi.size()];

                    for (int k = 0; k < pi.size(); k++) {
                        parentIndices[k] = variables.indexOf(pi.get(k));
                    }

                    double s2 = _score.localScore(variables.indexOf(n), parentIndices);

                    if (s2 > s_new) {
                        s_new = s2;
                        z = z0;
                    }

                    pi.remove(z0);
                }

                if (s_new > s_old) {
                    pi.add(z);
                    G1.addDirectedEdge(z, n);
                    s_old = s_new;
                    changed = true;
                }

                for (Node z0 : pi) {
                    if (z0 == n) continue;

                    if (pi.contains(z0)) continue;
                    pi.remove(z0);

                    int[] parentIndices = new int[pi.size()];

                    for (int k = 0; k < pi.size(); k++) {
                        parentIndices[k] = variables.indexOf(pi.get(k));
                    }

                    double s2 = _score.localScore(variables.indexOf(n), parentIndices);

                    if (s2 > s_new) {
                        s_new = s2;
                        z = z0;
                    }

                    pi.add(z0);
                }

                if (s_new > s_old) {
                    pi.remove(z);
                    G1.removeEdge(z, n);
                    s_old = s_new;
                    changed = true;
                }
            }
        }

        score = scorer.score(G1);
        return G1;
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }
}
