package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;
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
public class BFCI implements Algorithm, HasKnowledge, UsesScoreWrapper, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private ScoreWrapper score;
    private IKnowledge knowledge = new Knowledge2();

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
//            search.setMaxDegree(parameters.getInt(Params.MAX_DEGREE));
            search.setKnowledge(knowledge);
            search.setMaxPathLength(parameters.getInt(Params.MAX_PATH_LENGTH));
            search.setCompleteRuleSetUsed(parameters.getBoolean(Params.COMPLETE_RULE_SET_USED));

            search.setBreakTies(parameters.getBoolean(Params.BREAK_TIES));
            search.setCacheScores(parameters.getBoolean(Params.CACHE_SCORES));
            search.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));

            search.setUseScore(parameters.getBoolean(Params.USE_SCORE));
//            boss.setGspDepth(parameters.getInt(Params.DEPTH));

//            search.setMethod(Boss.Method.BOSS);

//            if (parameters.getInt(Params.BOSS_METHOD) == 1) {
//                search.setMethod(Boss.Method.BOSS);
//            } else if (parameters.getInt(Params.BOSS_METHOD) == 2) {
//                search.setMethod(Boss.Method.SP);
//            } else {
//                throw new IllegalArgumentException("Unexpected method: " + parameters.getInt(Params.BOSS_METHOD));
//            }

            if (parameters.getBoolean(Params.BOSS_SCORE_TYPE) ) {
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

//    @Override
//    public Graph getComparisonGraph(Graph graph) {
//        return new DagToPag2(graph).convert();
//    }

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

//        params.add(Params.FAITHFULNESS_ASSUMED);
        params.add(Params.MAX_DEGREE);
        params.add(Params.MAX_PATH_LENGTH);
        params.add(Params.COMPLETE_RULE_SET_USED);

        params.add(Params.CACHE_SCORES);
        params.add(Params.NUM_STARTS);
//        params.add(Params.BOSS_METHOD);
        params.add(Params.BOSS_SCORE_TYPE);
        params.add(Params.BREAK_TIES);
        params.add(Params.USE_SCORE);
        params.add(Params.VERBOSE);

        params.add(Params.VERBOSE);
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
