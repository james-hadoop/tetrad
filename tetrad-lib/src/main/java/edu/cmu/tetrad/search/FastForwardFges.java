package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.Scorer;

import java.util.List;

/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC), given variables
 * im causal order. Implements the Global Score Search (FFS) algorithm.
 *
 * <p>
 * public class FastFowardFges  {
 * <p>
 * // The score used (default FML BIC). Lower is better.
 * private final Scorer scorer;
 * private double score = Double.NaN;
 * <p>
 * /**
 * Constructs a FFS search
 */
public class FastForwardFges implements FastForward {
    private final Scorer scorer;
    private double score;

    public FastForwardFges(Scorer scorer) {
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        DataSet data = scorer.getDataSet();

        Score _score = new SemBicScore(scorer.getDataSet());

        IKnowledge knowledge = new Knowledge2();
        final List<Node> variables = data.getVariables();

        for (int i = 0; i < variables.size(); i++) {
            knowledge.addToTier(i, order.get(i).getName());
        }

        Fges fges = new Fges(_score);
        fges.setKnowledge(knowledge);

        final Graph G2 = fges.search();
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
