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
public class AvgDegreeTrue implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AvDT";
    }

    @Override
    public String getDescription() {
        return "Average degree of the true graph";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        int V = trueGraph.getNumNodes();
        int E = trueGraph.getNumEdges();
        return (2.0 * E) / (double) (V);
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value);
    }
}
