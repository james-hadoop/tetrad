package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.FgesBicScorer;
import edu.cmu.tetrad.sem.Scorer;
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
        name = "BOSS",
        command = "boss",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class BOSS implements Algorithm, HasKnowledge, UsesScoreWrapper, TakesInitialGraph {

    static final long serialVersionUID = 23L;

    private IKnowledge knowledge = new Knowledge2();
    private ScoreWrapper scoreWrapper;
    private Graph initialGraph;
    private Algorithm algorithm;

    public BOSS() {

    }

    public BOSS(ScoreWrapper score) {
        this.scoreWrapper = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (algorithm != null) {
            this.initialGraph = algorithm.search(dataSet, parameters);
        }

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            Score score = getScoreWrapper().getScore(dataSet, parameters);

            BestOrderScoreSearch search = new BestOrderScoreSearch(new K2(score));
            search.setAlgorithm(parameters.getInt(Params.DEPTH));
//            BestOrderScoreSearch search = new BestOrderScoreSearch(new FastFowardDsep(initialGraph, scorer));
            List<Node> variables = new ArrayList<>(score.getVariables());
            Graph graph = search.search(variables);

            System.out.println("Score for original order = " + search.getScoreOriginalOrder());
            System.out.println("Score for learned order = " + search.getScoreLearnedOrder());

//            return graph;
            return SearchGraphUtils.patternForDag(graph);
        } else {
            BOSS fges = new BOSS();

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
        return "BOSS K2";
    }

    @Override
    public DataType getDataType() {
        return scoreWrapper.getDataType();
    }

    @Override
    public List<String> getParameters() {
        ArrayList<String> strings = new ArrayList<>();
        strings.add(Params.DEPTH);
        return strings;
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

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public void setInitialGraph(Algorithm algorithm) {
        this.algorithm = algorithm;
    }
}
