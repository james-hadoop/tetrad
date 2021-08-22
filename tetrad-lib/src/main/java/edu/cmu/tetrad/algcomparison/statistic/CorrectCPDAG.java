package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.SearchGraphUtils;

/**
 * Outputs 1 if the skeleton is correct, 0 if not..
 *
 * @author jdramsey
 */
public class CorrectCPDAG implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "CorrCpdag";
    }

    @Override
    public String getDescription() {
        return "Correct CPDAG";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        return SearchGraphUtils.patternForDag(trueGraph).equals(SearchGraphUtils.patternForDag(estGraph)) ?
                1 : 0;
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
