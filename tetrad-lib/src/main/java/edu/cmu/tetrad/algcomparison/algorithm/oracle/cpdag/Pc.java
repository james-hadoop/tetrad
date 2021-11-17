package edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag;

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
import edu.cmu.tetrad.search.PcAll;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

@edu.cmu.tetrad.annotation.Algorithm(
        name = "PC",
        command = "pc",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class Pc implements Algorithm, TakesInitialGraph, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;

    public Pc(IndependenceWrapper test) {
        this.test = test;
    }

    public Pc(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    public Pc() {
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            if (algorithm != null) {
//            	initialGraph = algorithm.search(dataSet, parameters);
            }
            edu.cmu.tetrad.search.PcAll search = new edu.cmu.tetrad.search.PcAll(test.getTest(dataSet, parameters, trueGraph), initialGraph);
            search.setDepth(parameters.getInt(Params.DEPTH));
            search.setKnowledge(dataSet.getKnowledge());
            search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.REGULAR);
            search.setConcurrent(edu.cmu.tetrad.search.PcAll.Concurrent.NO);
            search.setColliderDiscovery(PcAll.ColliderDiscovery.FAS_SEPSETS);
            search.setConflictRule(PcAll.ConflictRule.PRIORITY);
            search.setHeuristic(parameters.getInt(Params.FAS_HEURISTIC));
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));
            return search.search();
        } else {
            Pc algorithm = new Pc(test);

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

    @Override
    public String getDescription() {
        return "PC (\"Peter and Clark\"), Priority Rule, using " + test.getDescription()
                + (algorithm != null ? " with initial graph from " +
                algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
//        parameters.add(Params.DEPTH);

        parameters.add(Params.STABLE_FAS);
        parameters.add(Params.CONCURRENT_FAS);
//        parameters.add(Params.COLLIDER_DISCOVERY_RULE);
        parameters.add(Params.CONFLICT_RULE);
        parameters.add(Params.DEPTH);
        parameters.add(Params.FAS_HEURISTIC);
//        parameters.add(Params.USE_MAX_P_ORIENTATION_HEURISTIC);
//        parameters.add(Params.MAX_P_ORIENTATION_MAX_PATH_LENGTH);


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
    public void setIndependenceWrapper(IndependenceWrapper independenceWrapper) {
        this.test = independenceWrapper;
    }
}