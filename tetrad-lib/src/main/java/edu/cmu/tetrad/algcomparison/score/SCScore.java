package edu.cmu.tetrad.algcomparison.score;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.ConditionalGaussianScore;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.StochasticComplexityScore;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for the Stochastic Complexity score.
 * Encompasses the fNML and qNML scores.
 *
 * @author Bryan Andrews
 */
@edu.cmu.tetrad.annotation.Score(
        name = "Stochastic Complexity (SC) Score",
        command = "stoch-comp-bic",
        dataType = DataType.Discrete
)
public class SCScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;
        final StochasticComplexityScore stochasticComplexityScore =
                new StochasticComplexityScore(DataUtils.getDiscreteDataSet(dataSet),
                        parameters.getBoolean("scoreEq"));
        return stochasticComplexityScore;
    }

    @Override
    public String getDescription() {
        return "Stochastic Complexity Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Discrete;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("scoreEq");
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return  dataSet.getVariable(name);
    }
}
