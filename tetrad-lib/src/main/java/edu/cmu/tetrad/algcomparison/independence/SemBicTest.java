package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.Experimental;
import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;

import java.util.*;

/**
 * Wrapper for Fisher Z test.
 *
 * @author jdramsey
 */
@TestOfIndependence(
        name = "SEM BIC Test",
        command = "sem-bic-test",
        dataType = {DataType.Continuous, DataType.Covariance}
)
//@Experimental
public class SemBicTest implements IndependenceWrapper {

    static final long serialVersionUID = 23L;

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        SemBicScore score;

        if (dataSet instanceof ICovarianceMatrix) {
            score = new SemBicScore((ICovarianceMatrix) dataSet);
        } else {
            score = new SemBicScore((DataSet) dataSet);
        }

        score.setPenaltyDiscount(parameters.getDouble(Params.PENALTY_DISCOUNT));
//        score.setStructurePrior(parameters.getDouble(Params.SEM_BIC_STRUCTURE_PRIOR));

        switch (parameters.getInt(Params.SEM_BIC_RULE)) {
            case 1:
                score.setRuleType(SemBicScore.RuleType.CHICKERING);
                break;
            case 2:
                score.setRuleType(edu.cmu.tetrad.search.SemBicScore.RuleType.NANDY);
                break;
            case 3:
                score.setRuleType(SemBicScore.RuleType.HIGH_DIMENSIONAL);
                break;
            default:
                throw new IllegalStateException("Expecting 1, 2, 3: " + parameters.getInt(Params.SEM_BIC_RULE));
        }

        return new IndTestScore(score, dataSet);
    }

    @Override
    public String getDescription() {
        return "SEM BIC Test";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add(Params.PENALTY_DISCOUNT);
        params.add(Params.SEM_BIC_STRUCTURE_PRIOR);
        params.add(Params.SEM_BIC_RULE);
        return params;
    }
}
