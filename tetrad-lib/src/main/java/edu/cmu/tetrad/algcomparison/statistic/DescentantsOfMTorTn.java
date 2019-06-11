package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.List;

/**
 * The arrow precision. This counts arrowheads maniacally, wherever they occur in the graphs.
 * The true positives are the number of arrowheads in both the true and estimated graphs.
 * Thus, if the true contains X*->Y and estimated graph either does not contain an edge from
 * X to Y or else does not contain an arrowhead at X for an edge from X to Y, one false
 * positive is counted. Similarly for false negatives.
 *
 * @author jdramsey
 */
public class DescentantsOfMTorTn implements Statistic {
    static final long serialVersionUID = 23L;



    @Override
    public String getAbbreviation() {
        return "mtorDescTn";
    }

    @Override
    public String getDescription() {
        return "Descendants of mTor TN";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        HillStats hillStats = new HillStats(trueGraph, estGraph);
        List<Node> tn = hillStats.getTn();

        return tn.size();
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
