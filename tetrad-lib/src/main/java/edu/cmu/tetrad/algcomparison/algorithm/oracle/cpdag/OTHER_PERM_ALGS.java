package edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.OtherPermAlgs;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * BOSS (Best Order Score Search).
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "OTHER_PERM_ALGS",
        command = "other_perm_algs",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class OTHER_PERM_ALGS implements Algorithm, UsesScoreWrapper, TakesIndependenceWrapper {
    static final long serialVersionUID = 23L;
    private ScoreWrapper score = null;
    private IndependenceWrapper test;

    public OTHER_PERM_ALGS() {
        // Used in reflection; do not delete.
    }

    public OTHER_PERM_ALGS(ScoreWrapper score, IndependenceWrapper test) {
        this.score = score;
        this.test = test;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (trueGraph != null && dataSet.getVariables().size() == trueGraph.getNumNodes()) {
            trueGraph = GraphUtils.replaceNodes(trueGraph, dataSet.getVariables());
        }

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            Score score = this.score.getScore(dataSet, parameters);

            IndependenceTest test = this.test.getTest(dataSet, parameters, trueGraph);

            test.setVerbose(parameters.getBoolean(Params.VERBOSE));

            OtherPermAlgs otherPermAlgs;

            otherPermAlgs = new OtherPermAlgs(test, score);

            OtherPermAlgs.Method method;

            switch (parameters.getInt(Params.OTHER_PERM_METHOD)) {
                case 1:
                    method = OtherPermAlgs.Method.RCG;
                    break;
                case 2:
                    method = OtherPermAlgs.Method.GSP;
                    break;
                case 3:
                    method = OtherPermAlgs.Method.ESP;
                    break;
                case 4:
                    method = OtherPermAlgs.Method.SP;
                    break;
                default:
                    throw new IllegalStateException("Pick a number from 1 to 4: " +
                            "1 = RCG, 2 = GSP, 3 = ESP, 4 = SP");
            }

            System.out.println("Picked " + method);

            otherPermAlgs.setMethod(method);
            otherPermAlgs.setUseScore(parameters.getBoolean(Params.GRASP_USE_SCORE));

            otherPermAlgs.setCacheScores(parameters.getBoolean(Params.CACHE_SCORES));
            otherPermAlgs.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            otherPermAlgs.setVerbose(parameters.getBoolean(Params.VERBOSE));
            otherPermAlgs.setKnowledge(dataSet.getKnowledge());
            otherPermAlgs.setUsePearl(parameters.getBoolean(Params.GRASP_USE_VERMA_PEARL));
            otherPermAlgs.setDepth(parameters.getInt(Params.GRASP_DEPTH));
            otherPermAlgs.setNumRounds(parameters.getInt(Params.NUM_ROUNDS));

            otherPermAlgs.bestOrder(score.getVariables());
            return otherPermAlgs.getGraph(parameters.getBoolean(Params.OUTPUT_CPDAG));
        } else {
            OTHER_PERM_ALGS algorithm = new OTHER_PERM_ALGS(score, test);

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
        return "BOSS (Best Order Score Search) using " + test.getDescription()
                + " or " + score.getDescription();
    }

    @Override
    public DataType getDataType() {
        return score.getDataType();
    }

    @Override
    public List<String> getParameters() {
        ArrayList<String> params = new ArrayList<>();
        params.add(Params.CACHE_SCORES);
        params.add(Params.NUM_STARTS);
        params.add(Params.BOSS_SCORE_TYPE);
        params.add(Params.GRASP_USE_SCORE);
//        params.add(Params.OUTPUT_CPDAG);
        params.add(Params.GRASP_DEPTH);
        params.add(Params.OTHER_PERM_METHOD);
        params.add(Params.NUM_ROUNDS);
        params.add(Params.VERBOSE);
        return params;
    }

    @Override
    public ScoreWrapper getScoreWrapper() {
        return score;
    }

    @Override
    public void setScoreWrapper(ScoreWrapper score) {
        this.score = score;
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
