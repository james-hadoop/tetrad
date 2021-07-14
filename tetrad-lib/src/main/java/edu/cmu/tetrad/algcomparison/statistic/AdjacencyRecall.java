package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.AdjacencyConfusion;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

/**
 * The adjacency recall. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class AdjacencyRecall implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AR";
    }

    @Override
    public String getDescription() {
        return "Adjacency Recall";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());
        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(trueGraph, estGraph);
        int adjTp = adjConfusion.getAdjTp();
        int adjFn = adjConfusion.getAdjFn();
        return adjTp / (double) (adjTp + adjFn);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
