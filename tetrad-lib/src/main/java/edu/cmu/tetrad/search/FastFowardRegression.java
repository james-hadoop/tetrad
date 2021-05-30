package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionCovariance;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.sem.Scorer;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.abs;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FFS) algorithm.
 *
 * @author josephramsey
 */
public class FastFowardRegression implements FastForward {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;
    private double score = Double.NaN;

    /**
     * Constructs a FFS search
     *
     * @param scorer the scorer used, by default FML BIC (for linear models). The score
     *               in general should be lower for better models.
     */
    public FastFowardRegression(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        RegressionCovariance regression = new RegressionCovariance(new CovarianceMatrix(scorer.getDataSet()));
        Graph G2 = new EdgeListGraph(order);

        for (int i = 0; i < order.size(); i++) {
            Node y = order.get(i);

            List<Node> parents = new ArrayList<>();
            for (int j = 0; j < i; j++) {
                parents.add(order.get(j));
            }

            RegressionResult result = regression.regress(y, parents);

            for (int j = 0; j < i; j++) {
                if (abs(result.getCoef()[j + 1]) > 0.1) {
                    G2.addDirectedEdge(parents.get(j), order.get(i));
                }
            }
        }

        score = scorer.score(G2);
        return G2;
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }
}
