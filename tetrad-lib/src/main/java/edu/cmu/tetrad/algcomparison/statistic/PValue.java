package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.sem.FmlBicScorer;

import static java.lang.Math.tanh;

/**
 * Returns the p-value of a linear, Gaussian model..
 *
 * @author jdramsey
 */
public class PValue implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "PValue";
    }

    @Override
    public String getDescription() {
        return "PValue";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        FmlBicScorer scorer = new FmlBicScorer((DataSet) dataModel, 1);
        scorer.score(SearchGraphUtils.dagFromPattern(estGraph));
        return scorer.getPValue();
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value / 10);
    }
}

