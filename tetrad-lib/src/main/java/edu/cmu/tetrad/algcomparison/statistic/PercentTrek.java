package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
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
public class PercentTrek implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "%Trek";
    }

    @Override
    public String getDescription() {
        return "Percent Trek";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        List<Node> nodes = trueGraph.getNodes();
        int tp = 0;
        int all = 0;

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (trueGraph.existsTrek(nodes.get(i), nodes.get(j))) {
                    tp++;
                }

                all++;

            }
        }
//        for (int k = 0; k < 1000; k++) {
//            int i = RandomUtil.getInstance().nextInt(nodes.size());
//            int j = RandomUtil.getInstance().nextInt(nodes.size());
//            if (i == j) {
//                k--;
//                continue;
//            }
//
//            if (trueGraph.existsTrek(nodes.get(i), nodes.get(j))) {
//                tp++;
//            }
//
//            all++;
//        }

        return tp / (double) (all);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
