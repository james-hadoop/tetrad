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
        name = "Sem GIC Score",
        command = "sem-gic-score",
        dataType = {DataType.Continuous, DataType.Covariance}
)
public class SemGicScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;

        edu.cmu.tetrad.search.SemGicScore semGicScore;

        if (dataSet instanceof DataSet) {
            semGicScore = new edu.cmu.tetrad.search.SemGicScore((DataSet) this.dataSet);
        } else if (dataSet instanceof ICovarianceMatrix) {
            semGicScore = new edu.cmu.tetrad.search.SemGicScore((ICovarianceMatrix) this.dataSet);
        } else {
            throw new IllegalArgumentException("Expecting either a dataset or a covariance matrix.");
        }

        semGicScore.setTrueErrorVariance(parameters.getDouble(Params.PENALTY_DISCOUNT));

        switch (parameters.getInt(Params.SEM_GIC_RULE)) {
            case 1:
                semGicScore.setRuleType(edu.cmu.tetrad.search.SemGicScore.RuleType.BIC);
                break;
            case 2:
                semGicScore.setRuleType(edu.cmu.tetrad.search.SemGicScore.RuleType.RICc);
                break;
            case 3:
                semGicScore.setRuleType(edu.cmu.tetrad.search.SemGicScore.RuleType.GIC5);
                break;
            case 4:
                semGicScore.setRuleType(edu.cmu.tetrad.search.SemGicScore.RuleType.GIC6);
                break;
            default:
                throw new IllegalStateException("Expecting 1, 2, 3 or 4: " + parameters.getInt(Params.SEM_BIC_RULE));
        }

        return semGicScore;
    }

    @Override
    public String getDescription() {
        return "Sem GIC Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.TRUE_ERROR_VARIANCE);
        parameters.add(Params.SEM_GIC_RULE);
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }

}
