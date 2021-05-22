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
import edu.cmu.tetrad.search.GlobalScoreSearch;
import edu.cmu.tetrad.sem.FgesBicScorer;
import edu.cmu.tetrad.sem.FmlBicScorer;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

/**
 * GSS (Global Scoring Search).
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FSS",
        command = "fss",
        algoType = AlgType.forbid_latent_common_causes
)
@Bootstrapping
public class Fss implements Algorithm, HasKnowledge, UsesScoreWrapper {

    static final long serialVersionUID = 23L;

    private IKnowledge knowledge = new Knowledge2();
    private ScoreWrapper scoreWrapper;

    public Fss() {

    }

    public Fss(ScoreWrapper score) {
        this.scoreWrapper = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {

            Scorer scorer;
            if (dataSet.isContinuous()) {
                scorer = new FmlBicScorer((DataSet) dataSet,
                        parameters.getDouble(Params.PENALTY_DISCOUNT));
            } else {
                scorer = new FgesBicScorer(getScoreWrapper().getScore(dataSet, parameters), ((DataSet) dataSet));
            }

            GlobalScoreSearch search
                    = new GlobalScoreSearch(scorer);
//            search.setKnowledge(knowledge);

            List<Node> variables = new ArrayList<>(scorer.getVariables());
//            Collections.sort(variables);

            return search.search(variables);
        } else {
            Fss fges = new Fss();

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
        return "Global Scoring Search";
    }

    @Override
    public DataType getDataType() {
        return scoreWrapper.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.SYMMETRIC_FIRST_STEP);
        parameters.add(Params.MAX_DEGREE);
        parameters.add(Params.PARALLELISM);
        parameters.add(Params.FAITHFULNESS_ASSUMED);

        parameters.add(Params.VERBOSE);
        parameters.add(Params.MEEK_VERBOSE);

        return parameters;
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
    public void setScoreWrapper(ScoreWrapper score) {
        this.scoreWrapper = score;
    }

    @Override
    public ScoreWrapper getScoreWrapper() {
        return scoreWrapper;
    }
}
