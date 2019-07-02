package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
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
public class AncestorRecall implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AncR";
    }

    @Override
    public String getDescription() {
        return "Ancestor Recall";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph _trueGraph = GraphUtils.replaceNodes(trueGraph, estGraph.getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(estGraph, _trueGraph.getNodes());

        int tp = 0;
        int fp = 0;

        List<Node> nodes = _trueGraph.getNodes();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) continue;

                Node x = nodes.get(i);
                Node y = nodes.get(j);

//                if (_trueGraph.getNode(x.getName()) == null || _trueGraph.getNode(y.getName()) == null) {
//                    continue;
//                }
//
//                if (_estGraph.getNode(x.getName()) == null || _estGraph.getNode(y.getName()) == null) {
//                    continue;
//                }

//                if (!estGraph.existsTrek(x, y) && !estGraph.existsTrek(y, x)) continue;
//                if (!trueGraph.existsTrek(x, y) && !trueGraph.existsTrek(y, x)) continue;

                if (_trueGraph.isAncestorOf(x, y)) {
                    if (_estGraph.isAncestorOf(x, y)) {
                        tp++;
                    } else {
                        fp++;
                    }
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
