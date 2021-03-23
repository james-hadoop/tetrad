package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestScore;
import edu.cmu.tetrad.search.IndependenceTest;
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
@TestOfIndependence(
        name = "Zhang-Shen Bound Test",
        command = "zs-test",
        dataType = {DataType.Continuous, DataType.Covariance}
)
public class ZhangShenBoundTest implements IndependenceWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;

        edu.cmu.tetrad.search.ZhangShenBoundScore score;

        if (dataSet instanceof DataSet) {
            score = new edu.cmu.tetrad.search.ZhangShenBoundScore((DataSet) this.dataSet);
        } else if (dataSet instanceof ICovarianceMatrix) {
            score = new edu.cmu.tetrad.search.ZhangShenBoundScore((ICovarianceMatrix) this.dataSet);
        } else {
            throw new IllegalArgumentException("Expecting either a dataset or a covariance matrix.");
        }

        score.setCalculateSquaredEuclideanNorms(parameters.getBoolean(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED));
        score.setRiskBound(parameters.getDouble(Params.ZS_RISK_BOUND));
        score.setCorrelationThreshold(parameters.getDouble(Params.CORRELATION_THRESHOLD));
        score.setTakeLog(parameters.getBoolean(Params.TAKE_LOGS));
//        score.setPenaltyDiscount(parameters.getDouble(Params.PENALTY_DISCOUNT));
//        score.setTrueErrorVariance(parameters.getDouble(Params.TRUE_ERROR_VARIANCE));
        return new IndTestScore(score, dataSet);
    }

    @Override
    public String getDescription() {
        return "Zhang-Shen Bound Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.ZS_RISK_BOUND);
        parameters.add(Params.CORRELATION_THRESHOLD);
//        parameters.add(Params.PENALTY_DISCOUNT);
        parameters.add(Params.TAKE_LOGS);
        parameters.add(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED);
//        parameters.add(Params.TRUE_ERROR_VARIANCE);
        return parameters;
    }

    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }
}
