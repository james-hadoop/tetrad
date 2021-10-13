package edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * PC.
 *
 * @author jdramsey
 */
public class SingleGraphAlg implements Algorithm {

    static final long serialVersionUID = 23L;
    private final Graph graph;

    public SingleGraphAlg(Graph graph) {
        this.graph = graph;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        return graph;
    }

    @Override
    public String getDescription() {
        return "Given graph";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> getParameters() {
        return new ArrayList<>();
    }
}
