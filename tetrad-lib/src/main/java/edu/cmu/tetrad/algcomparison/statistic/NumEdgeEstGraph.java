package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.graph.Graph;

import static java.lang.Math.tanh;

/**
 * The arrow precision. This counts arrowheads maniacally, wherever they occur in the graphs.
 * The true positives are the number of arrowheads in both the true and estimated graphs.
 * Thus, if the true contains X*->Y and estimated graph either does not contain an edge from
 * X to Y or else does not contain an arrowhead at X for an edge from X to Y, one false
 * positive is counted. Similarly for false negatives.
 *
 * @author jdramsey
 */
public class NumEdgeEstGraph implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "NumEdgesEst";
    }

    @Override
    public String getDescription() {
        return "Number of Edges in the Estimated Graph";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        return estGraph.getNumEdges();
    }

    @Override
    public double getNormValue(double value) {
        return 1.0 / value;
    }
}
