package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.RandomUtil;

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
public class PercentNonancestor implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "%NonAnc";
    }

    @Override
    public String getDescription() {
        return "Percent Monancestor";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        List<Node> nodes = trueGraph.getNodes();
        int tp = 0;
        int all = 0;

//        for (int i = 0; i < nodes.size(); i++) {
//            for (int j = 0; j < nodes.size(); j++) {
//                if (i == j) {
//                    continue;
//                }
//
//                final Node x = nodes.get(i);
//                final Node y = nodes.get(j);
//
//                if (!trueGraph.isAncestorOf(x, y)) {
//                    tp++;
//                }
//
////                if (trueGraph.existsTrek(nodes.get(i), nodes.get(j))) {
////                    tp++;
////                }
//
//                all++;
//
//            }
//        }

        for (int k = 0; k < 1000; k++) {
            int i = RandomUtil.getInstance().nextInt(nodes.size());
            int j = RandomUtil.getInstance().nextInt(nodes.size());
            if (i == j) {
                k--;
                continue;
            }

            final Node x = nodes.get(i);
            final Node y = nodes.get(j);

            if (!trueGraph.existsTrek(x, y)) {
                k--;
                continue;
            }

            if (!trueGraph.isAncestorOf(x, y)) {
                tp++;
            }

            all++;
        }

        return tp / (double) (all);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
