package edu.cmu.tetrad.algcomparison.graph;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Stores a single graph for use in simulations, etc.
 *
 * @author jdramsey
 */
public class SingleGraph implements RandomGraph {
    static final long serialVersionUID = 23L;

    private final Graph graph;

    public SingleGraph(Graph graph) {
        this.graph = graph;
    }

    @Override
    public Graph createGraph(Parameters parameters) {
        return this.graph;
    }

    @Override
    public String getDescription() {
        return "Graph supplied by user";
    }

    @Override
    public List<String> getParameters() {
        return new ArrayList<>();
    }
}
