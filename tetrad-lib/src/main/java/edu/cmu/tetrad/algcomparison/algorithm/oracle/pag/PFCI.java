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
import edu.cmu.tetrad.search.Pfci;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


/**
 * Adjusts GFCI to use a permutation algorithm (such as GRaSP) to do the initial
 * steps of finding adjacencies and unshielded colliders.
 *
 * GFCI reference is this:
 *
 * J.M. Ogarrio and P. Spirtes and J. Ramsey, "A Hybrid Causal Search Algorithm
 * for Latent Variable Models," JMLR 2016.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "PFCI",
        command = "pfci",
        algoType = AlgType.allow_latent_common_causes
)
@Bootstrapping
public class PFCI implements Algorithm, UsesScoreWrapper, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private ScoreWrapper score;

    public PFCI() {
        // Used for reflection; do not delete.
    }

    public PFCI(ScoreWrapper score, IndependenceWrapper test) {
        this.test = test;
        this.score = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            Pfci search = new Pfci(test.getTest(dataSet, parameters, trueGraph), score.getScore(dataSet, parameters));
            search.setKnowledge(dataSet.getKnowledge());
            search.setMaxPathLength(parameters.getInt(Params.MAX_PATH_LENGTH));
            search.setCompleteRuleSetUsed(parameters.getBoolean(Params.COMPLETE_RULE_SET_USED));

            search.setDepth(parameters.getInt(Params.GRASP_DEPTH));
            search.setUncoveredDepth(parameters.getInt(Params.GRASP_UNCOVERED_DEPTH));
            search.setUseForwardTuckOnly(parameters.getBoolean(Params.GRASP_FORWARD_TUCK_ONLY));
            search.setUsePearl(parameters.getBoolean(Params.GRASP_USE_PEARL));
            search.setTimeout(parameters.getInt(Params.TIMEOUT));
            search.setVerbose(parameters.getBoolean(Params.VERBOSE));
            search.setNumStarts(parameters.getInt(Params.NUM_STARTS));
            search.setGraspAlg(parameters.getBoolean(Params.GRASP_ALG));

            Object obj = parameters.get(Params.PRINT_STREAM);

            if (obj instanceof PrintStream) {
                search.setOut((PrintStream) obj);
            }

            return search.search();
        } else {
            PFCI algorithm = new PFCI(score, test);
            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, algorithm, parameters.getInt(Params.NUMBER_RESAMPLING));
            search.setKnowledge(data.getKnowledge());

            search.setPercentResampleSize(parameters.getDouble(Params.PERCENT_RESAMPLE_SIZE));
            search.setResamplingWithReplacement(parameters.getBoolean(Params.RESAMPLING_WITH_REPLACEMENT));

            ResamplingEdgeEnsemble edgeEnsemble;

            switch (parameters.getInt(Params.RESAMPLING_ENSEMBLE, 1)) {
                case 0:
                    edgeEnsemble = ResamplingEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = ResamplingEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = ResamplingEdgeEnsemble.Majority;
                    break;
                default:
                    throw new IllegalStateException("Unexpected value: " + parameters.getInt(Params.RESAMPLING_ENSEMBLE, 1));
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
        return "PFCI (Permutation-step FCI) using " + test.getDescription()
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
        params.add(Params.MAX_PATH_LENGTH);

        // Flags
        params.add(Params.GRASP_DEPTH);
        params.add(Params.GRASP_UNCOVERED_DEPTH);
        params.add(Params.GRASP_FORWARD_TUCK_ONLY);
        params.add(Params.GRASP_USE_PEARL);
        params.add(Params.TIMEOUT);
        params.add(Params.VERBOSE);
        params.add(Params.NUM_STARTS);
        params.add(Params.GRASP_ALG);

        // Parameters
        params.add(Params.NUM_STARTS);

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
