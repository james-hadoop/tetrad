package edu.cmu.tetrad.algcomparison.algorithm.multi;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
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

import static edu.cmu.tetrad.util.Params.*;

/**
 * Wraps the IMaGES algorithm for continuous variables.
 * </p>
 * Requires that the parameter 'randomSelectionSize' be set to indicate how many
 * datasets should be taken at a time (randomly). This cannot given multiple values.
 *
 * @author jdramsey
 */
@Bootstrapping
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FASK2",
        command = "fask2",
        algoType = AlgType.forbid_latent_common_causes,
        dataType = DataType.Continuous
)
public class Fask2 implements Algorithm, UsesScoreWrapper {
    static final long serialVersionUID = 23L;
    private ScoreWrapper score = null;

    // Don't delete.
    public Fask2() {

    }

    public Fask2(ScoreWrapper score) {
        this.score = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        DataSet _data = (DataSet) dataSet;

        for (int j = 0; j < _data.getNumColumns(); j++) {
            for (int i = 0; i < _data.getNumRows(); i++) {
                if (Double.isNaN(_data.getDouble(i, j))) {
                    throw new IllegalArgumentException("Please remove or impute missing values.");
                }
            }
        }

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            edu.cmu.tetrad.search.Fask2 search;

            search = new edu.cmu.tetrad.search.Fask2(score.getScore(dataSet, parameters), (DataSet) dataSet);

            search.setDepth(parameters.getInt(GSP_DEPTH));
            search.setNumRounds(parameters.getInt(NUM_ROUNDS));
            search.setSkewEdgeThreshold(parameters.getDouble(SKEW_EDGE_THRESHOLD));
            search.setOrientationAlpha(parameters.getDouble(ORIENTATION_ALPHA));
            search.setTwoCycleScreeningCutoff(parameters.getDouble(TWO_CYCLE_SCREENING_THRESHOLD));
            search.setDelta(parameters.getDouble(FASK_DELTA));
            search.setEmpirical(!parameters.getBoolean(FASK_NONEMPIRICAL));
            search.setKnowledge(dataSet.getKnowledge());
            search.setVerbose(parameters.getBoolean(VERBOSE));

            int lrRule = parameters.getInt(FASK_LEFT_RIGHT_RULE);

            if (lrRule == 1) {
                search.setLeftRight(edu.cmu.tetrad.search.Fask2.LeftRight.FASK1);
            } else if (lrRule == 2) {
                search.setLeftRight(edu.cmu.tetrad.search.Fask2.LeftRight.FASK2);
            } else if (lrRule == 3) {
                search.setLeftRight(edu.cmu.tetrad.search.Fask2.LeftRight.RSKEW);
            } else if (lrRule == 4) {
                search.setLeftRight(edu.cmu.tetrad.search.Fask2.LeftRight.SKEW);
            } else if (lrRule == 5) {
                search.setLeftRight(edu.cmu.tetrad.search.Fask2.LeftRight.TANH);
            } else {
                throw new IllegalStateException("Unconfigured left right rule index: " + lrRule);
            }

            return search.search();
        } else {
            Fask2 fask = new Fask2(score);

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, fask, parameters.getInt(Params.NUMBER_RESAMPLING));
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
            search.setVerbose(parameters.getBoolean(VERBOSE));
            return search.search();
        }
    }


    @Override
    public String getDescription() {
        if (score != null) {
            return "FASK2 using " + score.getDescription();
        }

        throw new IllegalStateException("Score not available here.");
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        if (score != null) {
            parameters.addAll(score.getParameters());
        }

        parameters.add(GSP_DEPTH);
        parameters.add(NUM_ROUNDS);
        parameters.add(SKEW_EDGE_THRESHOLD);
        parameters.add(TWO_CYCLE_SCREENING_THRESHOLD);
        parameters.add(ORIENTATION_ALPHA);
        parameters.add(FASK_DELTA);
        parameters.add(FASK_LEFT_RIGHT_RULE);
        parameters.add(FASK_NONEMPIRICAL);
        parameters.add(VERBOSE);
        return parameters;
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