package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

/**
 * The arrow precision. This counts arrowheads maniacally, wherever they occur in the graphs.
 * The true positives are the number of arrowheads in both the true and estimated graphs.
 * Thus, if the true contains X*->Y and estimated graph either does not contain an edge from
 * X to Y or else does not contain an arrowhead at X for an edge from X to Y, one false
 * positive is counted. Similarly for false negatives.
 *
 * @author jdramsey
 */
public class TrekPrecision implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "TrekP";
    }

    @Override
    public String getDescription() {
        return "Trek Precision";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph _trueGraph = GraphUtils.replaceNodes(trueGraph, estGraph.getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(estGraph, _trueGraph.getNodes());



        int tp = 0;
        int fp = 0;

        for (Edge edge : _estGraph.getEdges()) {
            if (trueGraph.getNode(edge.getNode1().getName()) == null || trueGraph.getNode(edge.getNode2().getName()) == null) {
                continue;
            }
            if (estGraph.getNode(edge.getNode1().getName()) == null || estGraph.getNode(edge.getNode2().getName()) == null) {
                continue;
            }

            if (_trueGraph.existsTrek(edge.getNode1(), edge.getNode2())) {
                tp++;
            } else {
                fp++;
            }
        }

        return tp / (double) (tp + fp);
//        return estGraph.getNumEdges();

}

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
