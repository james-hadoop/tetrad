package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
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
public class FastFowardDsep implements FastForward {

    // The score used (default FML BIC). Lower is better.
    private final Scorer scorer;
    private final Graph trueDag;
    private double score = Double.NaN;

    /**
     * Constructs a FFS search
     */
    public FastFowardDsep(Graph trueDag, Scorer scorer) {
        this.trueDag = trueDag;
        this.scorer = scorer;
    }

    /**
     * Does the search.
     *
     * @param order The variables in causal order.
     * @return The estimated DAG.
     */
    public Graph search(List<Node> order) {
        Score _score = new GraphScore(trueDag);

        IKnowledge knowledge = new Knowledge2();
        final List<Node> variables = _score.getVariables();

        for (int i = 0; i < variables.size(); i++) {
            knowledge.addToTier(i, order.get(i).getName());
        }

        Pc pc = new Pc(new IndTestDSep(trueDag));
        pc.setKnowledge(knowledge);
        Graph G2 = pc.search();

//        Fges fges = new Fges(_score);
//        fges.setKnowledge(knowledge);
//        final Graph G2 = fges.search();

//        score = scorer.score(G2);
        return G2;
    }

    /**
     * Returns the score of the most recent search.
     */
    public double score() {
        return score;
    }
}
