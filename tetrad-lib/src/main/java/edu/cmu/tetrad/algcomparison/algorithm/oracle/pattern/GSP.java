package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.GreedySparsestPermutation;
import edu.cmu.tetrad.search.K2;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * BOSS (Best Order Scoring Search).
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "GSP",
        command = "gsp",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class GSP implements Algorithm, HasKnowledge, UsesScoreWrapper {

    static final long serialVersionUID = 23L;

    private IKnowledge knowledge = new Knowledge2();
    private ScoreWrapper scoreWrapper;

    public GSP() {

    }

    public GSP(ScoreWrapper score) {
        this.scoreWrapper = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {

            Score score = getScoreWrapper().getScore(dataSet, parameters);

            GreedySparsestPermutation search = new GreedySparsestPermutation(new K2(score));
            search.setMethod(GreedySparsestPermutation.Method.NORECURSIVE);
            search.setNumRestarts(5);

            List<Node> variables = new ArrayList<>(score.getVariables());
            Graph graph = search.search(variables);

            System.out.println("Score for original order = " + search.getScoreOriginalOrder());
            System.out.println("Score for learned order = " + search.getScoreLearnedOrder());

            return SearchGraphUtils.patternForDag(graph);
        } else {
            GSP fges = new GSP();

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
        return "GSP";
    }

    @Override
    public DataType getDataType() {
        return scoreWrapper.getDataType();
    }

    @Override
    public List<String> getParameters() {
        return new ArrayList<>();
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
    public ScoreWrapper getScoreWrapper() {
        return scoreWrapper;
    }

    @Override
    public void setScoreWrapper(ScoreWrapper score) {
        this.scoreWrapper = score;
    }
}
