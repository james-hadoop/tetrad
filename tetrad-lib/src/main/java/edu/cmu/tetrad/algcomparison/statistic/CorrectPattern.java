package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.AdjacencyConfusion;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class CorrectPattern implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "CorrPatt";
    }

    @Override
    public String getDescription() {
        return "Correct Pattern = 1, incorrect pattern = 0";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph pattern1 = SearchGraphUtils.patternForDag(trueGraph);
        Graph pattern2 = SearchGraphUtils.patternForDag(estGraph);
        if (pattern1.equals(pattern2)) return 1; else return 0;
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
