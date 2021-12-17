package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

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
import edu.cmu.tetrad.search.Bfci;
import edu.cmu.tetrad.search.TeyssierScorer;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


/**
 * BFCI.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "BFCI",
        command = "bfci",
        algoType = AlgType.allow_latent_common_causes
)
@Bootstrapping
public class BFCI implements Algorithm, UsesScoreWrapper, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private ScoreWrapper score;

    public BFCI() {
    }

    public BFCI(ScoreWrapper score, IndependenceWrapper test) {
        this.test = test;
        this.score = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            Bfci search = new Bfci(test.getTest(dataSet, parameters, trueGraph), score.getScore(dataSet, parameters));
            search.setKnowledge(dataSet.getKnowledge());
            search.setMaxPathLength(parameters.getInt(Params.MAX_PATH_LENGTH));
            search.setCompleteRuleSetUsed(parameters.getBoolean(Params.COMPLETE_RULE_SET_USED));
            search.setTriangleDepth(parameters.getInt(Params.MAX_PERM_SIZE));

            search.setCacheScores(parameters.getBoolean(Params.CACHE_SCORES));
            search.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));
            search.setUseScore(parameters.getBoolean(Params.USE_SCORE));
            search.setKnowledge(dataSet.getKnowledge());

            if (parameters.getBoolean(Params.BOSS_SCORE_TYPE)) {
                search.setScoreType(TeyssierScorer.ScoreType.Edge);
            } else {
                search.setScoreType(TeyssierScorer.ScoreType.SCORE);
            }

            Object obj = parameters.get(Params.PRINT_STREAM);

            if (obj instanceof PrintStream) {
                search.setOut((PrintStream) obj);
            }

            return search.search();
        } else {
            BFCI algorithm = new BFCI(score, test);
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
        return "BFCI (BOSS-FCI) using " + test.getDescription()
                + " or " + score.getDescription();
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();

        params.add(Params.MAX_PATH_LENGTH);
        params.add(Params.COMPLETE_RULE_SET_USED);

        params.add(Params.CACHE_SCORES);
        params.add(Params.NUM_STARTS);
        params.add(Params.BOSS_SCORE_TYPE);
        params.add(Params.USE_SCORE);
        params.add(Params.MAX_PERM_SIZE);
        params.add(Params.VERBOSE);

        params.add(Params.VERBOSE);
        return params;
    }

    @Override
    public IndependenceWrapper getIndependenceWrapper() {
        return test;
    }

    @Override
    public void setIndependenceWrapper(IndependenceWrapper test) {
        this.test = test;
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
