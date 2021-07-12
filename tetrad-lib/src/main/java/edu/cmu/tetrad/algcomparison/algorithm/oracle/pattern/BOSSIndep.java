package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.Boss;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * BOSS (Greedy Sparest Permutation, Simplified).
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "BOSSIndep",
        command = "boss-indep",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class BOSSIndep implements Algorithm, HasKnowledge, TakesIndependenceWrapper, TakesInitialGraph {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test = null;
    //    private ScoreWrapper score = null;
    private IKnowledge knowledge = new Knowledge2();
    private Graph initialGraph;
    private Algorithm algorithm;

    public BOSSIndep() {

    }

    public BOSSIndep(IndependenceWrapper test) {
        this.test = test;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (algorithm != null) {
            this.initialGraph = algorithm.search(dataSet, parameters);
        }

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            IndependenceTest test = this.test.getTest(dataSet, parameters);
            Boss boss = new Boss(test);
            boss.setCachingScores(parameters.getBoolean(Params.CACHE_SCORES));
            boss.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            boss.setReturnCpdag(true);

            if (parameters.getInt(Params.BOSS_METHOD) == 1) {
                boss.setMethod(Boss.Method.BOSS_ALL_INDICES);
            } else if (parameters.getInt(Params.BOSS_METHOD) == 2) {
                boss.setMethod(Boss.Method.BOSS_PROMOTION);
            } else if (parameters.getInt(Params.BOSS_METHOD) == 3) {
                boss.setMethod(Boss.Method.SP);
            } else {
                throw new IllegalArgumentException("Unexpected method: " + parameters.getInt(Params.BOSS_METHOD));
            }

            return boss.search(test.getVariables());
        } else {
            BOSSIndep fges = new BOSSIndep();

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, fges, parameters.getInt(Params.NUMBER_RESAMPLING));
            search.setKnowledge(knowledge);

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
    public Graph getComparisonGraph(Graph graph) {
        return new EdgeListGraph(graph);
    }

    @Override
    public String getDescription() {
        return "BOSS from test";
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        ArrayList<String> params = new ArrayList<>();
        params.add(Params.CACHE_SCORES);
        params.add(Params.NUM_STARTS);
        params.add(Params.BOSS_METHOD);
        return params;
    }

    @Override
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    @Override
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
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
