package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.ArrowConfusion;
import edu.cmu.tetrad.graph.*;

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
public class NonancestorPrecision implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "NonAncP";
    }

    @Override
    public String getDescription() {
        return "Nonancestor Precision";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        Graph _trueGraph = GraphUtils.replaceNodes(trueGraph, estGraph.getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(estGraph, _trueGraph.getNodes());

        int tp = 0;
        int fp = 0;

        for (Edge edge : estGraph.getEdges()) {
            if (trueGraph.getNode(edge.getNode1().getName()) == null || trueGraph.getNode(edge.getNode2().getName()) == null) {
                continue;
            }
            if (estGraph.getNode(edge.getNode1().getName()) == null || estGraph.getNode(edge.getNode2().getName()) == null) {
                continue;
            }

            if (!_trueGraph.existsTrek(edge.getNode1(), edge.getNode2())) {
                System.out.println("Not a trek: " + edge);
                continue;
            } else {
                List<List<Node>> treks = GraphUtils.treks(trueGraph, edge.getNode1(), edge.getNode2(), 20);
//                List<Node> trek = treks.get(0);
//
//                int length = Integer.MAX_VALUE;
//                for (List<Node> _trek : treks) {
//                    if (_trek.size() < length) {
//                        trek = _trek;
//                        length = _trek.size();
//                    }
//                }

//                for (List<Node> _trek : treks) {
//                    Edge _edge = trueGraph.getEdge(_trek.get(0), _trek.get(1));
//                    if (_edge.pointsTowards(_trek.get(1))) {
//                        trek = _trek;
//                        break;
//                    }
//                }

                System.out.println("Edge in PAG: " + edge);

                for (List<Node> trek : treks) {
                    System.out.println("\tTrek: " + GraphUtils.pathString(trek, trueGraph));
                }
            }

            if (edge.getEndpoint1() == Endpoint.ARROW) {
                if (!_trueGraph.isAncestorOf(edge.getNode1(), edge.getNode2())) {
                    tp++;
                } else {
                    fp++;
                }
            }

            if (edge.getEndpoint2() == Endpoint.ARROW) {
                if (!_trueGraph.isAncestorOf(edge.getNode2(), edge.getNode1())) {
                    tp++;
                } else {
                    fp++;
                }
            }
        }

        return tp / (double) (tp + fp);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
