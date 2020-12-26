package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.IndTestDSep;
import edu.cmu.tetrad.search.IndTestTSep;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.TConnection;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for D-separation test. Requires a true DAG as input.
 *
 * @author jdramsey
 */
@TestOfIndependence(
        name = "T-Separation Test",
        command = "t-sep-test",
        dataType = DataType.Graph
)
public class TSeparationTest implements IndependenceWrapper {

    static final long serialVersionUID = 23L;
    private Graph graph;

    /**
     * Use this empty constructor to satisfy the java reflection
     */
    public TSeparationTest() {

    }

    public TSeparationTest(Graph graph) {
        this.graph = graph;
    }

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        if (dataSet == null) {
            IndTestTSep indTestTSep = new IndTestTSep(graph);
            indTestTSep.setTimeLimit(parameters.getDouble(Params.TIME_LIMIT));
            return indTestTSep;
        } else {
            throw new IllegalArgumentException("Expecting no data for a d-separation test.");
        }
    }

    @Override
    public String getDescription() {
        return "T-Separation Test";
    }

    @Override
    public DataType getDataType() {
        return DataType.Graph;
    }

    @Override
    public List<String> getParameters() {
        ArrayList<String> params = new ArrayList<>();
        params.add(Params.TIME_LIMIT);
        return params;
    }

    public void setGraph(Graph graph) {
        this.graph = graph;
    }
    
}
