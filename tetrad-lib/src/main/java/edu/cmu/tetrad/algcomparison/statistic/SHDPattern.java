package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.SearchGraphUtils;

/**
 * Calculates the structural Hamming distance (SHD) between the estimated graph and
 * the true graph.
 *
 * @author jdramsey
 */
public class SHDPattern implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "SHDPat";
    }

    @Override
    public String getDescription() {
        return "Structural Hamming Distance comparing CPDAGs";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        GraphUtils.GraphComparison comparison
                = SearchGraphUtils.getGraphComparison3(SearchGraphUtils.patternFromDag(estGraph),
                SearchGraphUtils.patternFromDag(trueGraph), System.out);
        return comparison.getShd();
    }

    @Override
    /**
     * This will be given the index of the SHD stat.
     */
    public double getNormValue(double value) {
        return 1.0 - Math.tanh(0.001 * value);
    }
}
