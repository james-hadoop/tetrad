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
import edu.cmu.tetrad.search.Boss;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.TeyssierScorer;
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
        name = "BOSSALL",
        command = "bossall",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class BOSSALL implements Algorithm, UsesScoreWrapper, TakesIndependenceWrapper {
    static final long serialVersionUID = 23L;
    private ScoreWrapper score = null;
    private IndependenceWrapper test;

    public BOSSALL() {

    }

    public BOSSALL(ScoreWrapper score, IndependenceWrapper test) {
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

            Boss boss;

            boss = new Boss(test, score);

            Boss.Method method;

            switch (parameters.getInt(Params.BOSS_METHOD)) {
                case 1:
                    method = Boss.Method.BOSS;
                    break;
                case 2:
                    method = Boss.Method.GASP;
                    break;
                case 3:
                    method = Boss.Method.QUICK_GASP;
                    break;
                case 4:
                    method = Boss.Method.ESP;
                    break;
                case 5:
                    method = Boss.Method.TSP;
                    break;
                case 6:
                    method = Boss.Method.SP;
                    break;
                default:
                    throw new IllegalStateException("Pick a number from 1 to 6: 1 = BOSS, 2 = GASP, " +
                            "3 = QUICK_GASP, 4 = ESP, 5 = TSP, 6 = SP");
            }

            System.out.println("Picked " + method);

            boss.setMethod(method);
            boss.setUseScore(parameters.getBoolean(Params.USE_SCORE));

            boss.setCacheScores(parameters.getBoolean(Params.CACHE_SCORES));
            boss.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            boss.setVerbose(parameters.getBoolean(Params.VERBOSE));
            boss.setKnowledge(dataSet.getKnowledge());
            boss.setMaxPermSize(parameters.getInt(Params.MAX_PERM_SIZE));
            boss.setParentCalculation(TeyssierScorer.ParentCalculation.GrowShrinkMb);
            boss.setDepth(parameters.getInt(Params.GSP_DEPTH));

            if (parameters.getBoolean(Params.BOSS_SCORE_TYPE)) {
                boss.setScoreType(TeyssierScorer.ScoreType.Edge);
            } else {
                boss.setScoreType(TeyssierScorer.ScoreType.SCORE);
            }

            boss.bestOrder(score.getVariables());
            return boss.getGraph(parameters.getBoolean(Params.OUTPUT_CPDAG));
        } else {
            BOSSALL algorithm = new BOSSALL(score, test);

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
        return "BOSSALL (Best Order Score Search) using " + test.getDescription()
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
        params.add(Params.USE_SCORE);
        params.add(Params.OUTPUT_CPDAG);
        params.add(Params.MAX_PERM_SIZE);
        params.add(Params.GSP_DEPTH);
        params.add(Params.BOSS_METHOD);
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
