package edu.cmu.tetrad.algcomparison.score;

import edu.cmu.tetrad.data.*;
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
        name = "Sem BIC Score",
        command = "sem-bic-score",
        dataType = {DataType.Continuous, DataType.Covariance}
)
public class LinearGaussianBicScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;

        edu.cmu.tetrad.search.SemBicScore semBicScore;

        if (dataSet instanceof DataSet) {
            semBicScore = new edu.cmu.tetrad.search.SemBicScore((DataSet) this.dataSet);
        } else if (dataSet instanceof ICovarianceMatrix) {
            semBicScore = new edu.cmu.tetrad.search.SemBicScore((ICovarianceMatrix) this.dataSet);
        } else {
            throw new IllegalArgumentException("Expecting either a dataset or a covariance matrix.");
        }

        semBicScore.setPenaltyDiscount(parameters.getDouble(Params.PENALTY_DISCOUNT));
        semBicScore.setStructurePrior(parameters.getDouble(Params.SEM_BIC_STRUCTURE_PRIOR));

        switch (parameters.getInt(Params.SEM_BIC_RULE)) {
            case 1:
                semBicScore.setRuleType(edu.cmu.tetrad.search.SemBicScore.RuleType.CHICKERING);
                break;
            case 2:
                semBicScore.setRuleType(edu.cmu.tetrad.search.SemBicScore.RuleType.NANDY);
                break;
            default:
                throw new IllegalStateException("Expecting 1 or 2: " + parameters.getInt(Params.SEM_BIC_RULE));
        }

        return semBicScore;
    }

    @Override
    public String getDescription() {
        return "Linear, Gaussian BIC Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add(Params.PENALTY_DISCOUNT);
        parameters.add(Params.SEM_BIC_STRUCTURE_PRIOR);
        parameters.add(Params.SEM_BIC_RULE);
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }

}
