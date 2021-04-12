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
        name = "Kim et al. Scores",
        command = "kim-scores",
        dataType = {DataType.Continuous, DataType.Covariance}
)
public class KimEtAlScores implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;

        edu.cmu.tetrad.search.KimEtAlScores kimEtAlScores;

        if (dataSet instanceof DataSet) {
            kimEtAlScores = new edu.cmu.tetrad.search.KimEtAlScores((DataSet) this.dataSet,
                    parameters.getDouble(Params.CORRELATION_THRESHOLD));
        } else if (dataSet instanceof ICovarianceMatrix) {
            kimEtAlScores = new edu.cmu.tetrad.search.KimEtAlScores((ICovarianceMatrix) this.dataSet,
                    parameters.getDouble(Params.CORRELATION_THRESHOLD));
        } else {
            throw new IllegalArgumentException("Expecting either a dataset or a covariance matrix.");
        }

        kimEtAlScores.setTrueErrorVariance(parameters.getDouble(Params.TRUE_ERROR_VARIANCE));
        kimEtAlScores.setLambda(parameters.getDouble(Params.MANUAL_LAMBDA));
        kimEtAlScores.setPenaltyDiscount(parameters.getDouble(Params.PENALTY_DISCOUNT));
//        kimEtAlScores.setTakeLog(parameters.getBoolean(Params.TAKE_LOGS));
//        kimEtAlScores.setCalculateSquareEuclideanNorms(parameters.getBoolean(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED));

        switch (parameters.getInt(Params.SEM_GIC_RULE)) {
            case 1:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.BIC);
                break;
            case 2:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.GIC2);
                break;
            case 3:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.RIC);
                break;
            case 4:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.RICc);
                break;
            case 5:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.GIC5);
                break;
            case 6:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.GIC6);
                break;
            case 7:
                kimEtAlScores.setRuleType(edu.cmu.tetrad.search.KimEtAlScores.RuleType.MANUAL);
                break;
            default:
                throw new IllegalStateException("Expecting one of the available options: " + parameters.getInt(Params.SEM_BIC_RULE));
        }

        return kimEtAlScores;
    }

    @Override
    public String getDescription() {
        return "Kim et al. Scores";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.TRUE_ERROR_VARIANCE);
        parameters.add(Params.MANUAL_LAMBDA);
        parameters.add(Params.SEM_GIC_RULE);
        parameters.add(Params.PENALTY_DISCOUNT);
        parameters.add(Params.CORRELATION_THRESHOLD);
//        parameters.add(Params.TAKE_LOGS);
//        parameters.add(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED);
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }

}
