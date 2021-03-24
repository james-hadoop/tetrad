package edu.cmu.tetrad.algcomparison.score;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for linear, Gaussian SEM BIC score.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Score(
        name = "EBIC Score",
        command = "ebic-score",
        dataType = {DataType.Continuous, DataType.Covariance}
)
public class EBicScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;

        edu.cmu.tetrad.search.EBic eBicScore;

        if (dataSet instanceof DataSet) {
            eBicScore = new edu.cmu.tetrad.search.EBic((DataSet) this.dataSet);
        } else if (dataSet instanceof ICovarianceMatrix) {
            eBicScore = new edu.cmu.tetrad.search.EBic((ICovarianceMatrix) this.dataSet);
        } else {
            throw new IllegalArgumentException("Expecting either a dataset or a covariance matrix.");
        }

        eBicScore.setGamma(parameters.getDouble(Params.PENALTY_DISCOUNT));

        return eBicScore;
    }

    @Override
    public String getDescription() {
        return "EBIC Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.PENALTY_DISCOUNT);
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }

}
