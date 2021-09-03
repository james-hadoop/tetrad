package edu.cmu.tetrad.algcomparison.algorithm.multi;

import edu.cmu.tetrad.algcomparison.algorithm.MultiDataSetAlgorithm;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Fges;
import edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.annotation.Bootstrapping;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SemBicScoreMultiFas;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Wraps the MultiFask algorithm for continuous variables.
 * </p>
 * Requires that the parameter 'randomSelectionSize' be set to indicate how many
 * datasets should be taken at a time (randomly). This cannot given multiple
 * values.
 *
 * @author jdramsey
 */

// Using FaskVote now in place of MultiFask. Keeping the name "MultiFask" in the interface.
//@edu.cmu.tetrad.annotation.Algorithm(
//        name = "MultiFask",
//        command = "multi-fask",
//        algoType = AlgType.forbid_latent_common_causes,
//        dataType = DataType.Continuous
//)
@Bootstrapping
public class MultiFaskV1 implements MultiDataSetAlgorithm, HasKnowledge {

    static final long serialVersionUID = 23L;
    private IKnowledge knowledge = new Knowledge2();

    public MultiFaskV1() {

    }

    @Override
    public Graph search(List<DataModel> dataSets, Parameters parameters) {
        for (DataModel d : dataSets) {
            DataSet _data = (DataSet) d;

            for (int j = 0; j < _data.getNumColumns(); j++) {
                for (int i = 0; i < _data.getNumRows(); i++) {
                    if (Double.isNaN(_data.getDouble(i, j))) {
                        throw new IllegalArgumentException("Please remove or impute missing values.");
                    }
                }
            }
        }

        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            List<DataSet> _dataSets = new ArrayList<>();
            for (DataModel d : dataSets) {
                _dataSets.add((DataSet) d);
            }
            final SemBicScoreMultiFas score = new SemBicScoreMultiFas(dataSets);
            score.setPenaltyDiscount(parameters.getDouble(Params.PENALTY_DISCOUNT));
            edu.cmu.tetrad.search.MultiFaskV1 search = new edu.cmu.tetrad.search.MultiFaskV1(_dataSets, score);
            search.setKnowledge(knowledge);
            return search.search();
        } else {
            MultiFaskV1 imagesSemBic = new MultiFaskV1();

            List<DataSet> datasets = new ArrayList<>();

            for (DataModel dataModel : dataSets) {
                datasets.add((DataSet) dataModel);
            }
            GeneralResamplingTest search = new GeneralResamplingTest(datasets, imagesSemBic, parameters.getInt(Params.NUMBER_RESAMPLING));
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
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        if (parameters.getInt(Params.NUMBER_RESAMPLING) < 1) {
            return search(Collections.singletonList((DataModel) DataUtils.getContinuousDataSet(dataSet)), parameters);
        } else {
            MultiFaskV1 imagesSemBic = new MultiFaskV1();

            List<DataSet> dataSets = Collections.singletonList(DataUtils.getContinuousDataSet(dataSet));
            GeneralResamplingTest search = new GeneralResamplingTest(dataSets, imagesSemBic, parameters.getInt(Params.NUMBER_RESAMPLING));
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
//        return new EdgeListGraph(graph);
////        return SearchGraphUtils.cpdagForDag(graph);
////        return new TsDagToPag(new EdgeListGraph(graph)).convert();
//    }

    @Override
    public String getDescription() {
        return "IMaGES for continuous variables (using the SEM BIC score)";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        // MultiFask uses SemBicScore internally, so we'll need to add the score parameters too - Zhou
        List<String> parameters = new LinkedList<>();
        parameters.addAll((new Fges()).getParameters());
        parameters.addAll((new LinearGaussianBicScore()).getParameters());
        parameters.add(Params.NUM_RUNS);
        parameters.add(Params.RANDOM_SELECTION_SIZE);

        parameters.add(Params.VERBOSE);

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
}
