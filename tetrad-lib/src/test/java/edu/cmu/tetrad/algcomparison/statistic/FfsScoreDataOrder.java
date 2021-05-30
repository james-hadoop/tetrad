package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FastForwardGlobalScore;
import edu.cmu.tetrad.sem.FmlBicScorer;

import static java.lang.Math.tanh;

/**
 * Returns the p-value of a linear, Gaussian model..
 *
 * @author jdramsey
 */
public class FfsScoreDataOrder implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "FfsScoreDO";
    }

    @Override
    public String getDescription() {
        return "FFS score for variables in data order";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        FastForwardGlobalScore fss = new FastForwardGlobalScore(new FmlBicScorer(((DataSet) dataModel), 1));
        fss.search(dataModel.getVariables());
        System.out.println("FFS score for variables in data order = " + fss.score());
        System.out.println("Data order = " + dataModel.getVariables());
        return fss.score();
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value / 10);
    }
}

