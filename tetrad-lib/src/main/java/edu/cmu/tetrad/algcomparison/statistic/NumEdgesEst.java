package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;

import static java.lang.Math.tanh;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class NumEdgesEst implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "EdgE";
    }

    @Override
    public String getDescription() {
        return "Number of edges in the estimated graph";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        return estGraph.getNumEdges();
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value);
    }
}
