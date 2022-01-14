package edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * GRaSP.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "GRaSP",
        command = "grasp",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class Grasp implements Algorithm, UsesScoreWrapper {
    static final long serialVersionUID = 23L;
    private ScoreWrapper score = null;

    public Grasp() {
        // Used in reflection; do not delete.
    }

    public Grasp(ScoreWrapper score) {
        this.score = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            Score score = this.score.getScore(dataSet, parameters);

            GRaSP grasp = new GRaSP(score);

            grasp.setCacheScores(parameters.getBoolean(Params.CACHE_SCORES));
            grasp.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            grasp.setVerbose(parameters.getBoolean(Params.VERBOSE));
            grasp.setKnowledge(dataSet.getKnowledge());
            grasp.setParentCalculation(TeyssierScorer.ParentCalculation.GrowShrinkMb);
            grasp.setDepth(parameters.getInt(Params.GSP_DEPTH));
            grasp.setUseTuck(parameters.getBoolean(Params.USE_TUCK));

            if (parameters.getBoolean(Params.BOSS_SCORE_TYPE)) {
                grasp.setScoreType(TeyssierScorer.ScoreType.Edge);
            } else {
                grasp.setScoreType(TeyssierScorer.ScoreType.SCORE);
            }

            grasp.bestOrder(score.getVariables());
            return grasp.getGraph(parameters.getBoolean(Params.OUTPUT_CPDAG));
        } else {
            Grasp algorithm = new Grasp(score);

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
        return "GRaSP using " + score.getDescription();
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
        params.add(Params.OUTPUT_CPDAG);
        params.add(Params.GSP_DEPTH);
        params.add(Params.USE_TUCK);
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

}
