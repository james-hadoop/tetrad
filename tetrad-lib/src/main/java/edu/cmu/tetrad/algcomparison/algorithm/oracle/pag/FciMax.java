package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * FCI.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FCI-Max",
        command = "fcimax",
        algoType = AlgType.allow_latent_common_causes
)
@Bootstrapping
public class FciMax implements Algorithm, TakesInitialGraph, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;

    public FciMax() {
    }

    public FciMax(IndependenceWrapper test) {
        this.test = test;
    }

    public FciMax(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            if (algorithm != null) {
                initialGraph = algorithm.search(dataSet, parameters, trueGraph);
            }

            edu.cmu.tetrad.search.FciMax search = new edu.cmu.tetrad.search.FciMax(test.getTest(dataSet, parameters, trueGraph));
            search.setDepth(parameters.getInt(Params.DEPTH));
            search.setKnowledge(dataSet.getKnowledge());
            search.setMaxPathLength(parameters.getInt(Params.MAX_PATH_LENGTH));
            search.setCompleteRuleSetUsed(parameters.getBoolean(Params.COMPLETE_RULE_SET_USED));
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));

            return search.search();
        } else {
            FciMax algorithm = new FciMax(test);

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, algorithm, parameters.getInt(Params.NUMBER_RESAMPLING));
            search.setKnowledge(data.getKnowledge());

            search.setPercentResampleSize(parameters.getDouble(Params.PERCENT_RESAMPLE_SIZE));
            search.setResamplingWithReplacement(parameters.getBoolean(Params.RESAMPLING_WITH_REPLACEMENT));

            ResamplingEdgeEnsemble edgeEnsemble = ResamplingEdgeEnsemble.Highest;
            switch (parameters.getInt(Params.RESAMPLING_ENSEMBLE, 1)) {
                case 0:
                    edgeEnsemble = ResamplingEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = ResamplingEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = ResamplingEdgeEnsemble.Majority;
            }
            search.setEdgeEnsemble(edgeEnsemble);
            search.setAddOriginalDataset(parameters.getBoolean(Params.ADD_ORIGINAL_DATASET));

            search.setParameters(parameters);
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));
            return search.search();
        }
    }

    public String getDescription() {
        return "FCI-Max (Fast Causal Inference) using " + test.getDescription()
                + (algorithm != null ? " with initial graph from "
                + algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        parameters.add(Params.DEPTH);
        parameters.add(Params.MAX_PATH_LENGTH);
        parameters.add(Params.COMPLETE_RULE_SET_USED);

        parameters.add(Params.VERBOSE);
        return parameters;
    }

    @Override
    public Graph getInitialGraph() {
        return initialGraph;
    }

    @Override
    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    @Override
    public void setInitialGraph(Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    @Override
    public IndependenceWrapper getIndependenceWrapper() {
        return test;
    }

    @Override
    public void setIndependenceWrapper(IndependenceWrapper test) {
        this.test = test;
    }

}
