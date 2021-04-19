package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.AdjacencyConfusion;
import edu.cmu.tetrad.algcomparison.statistic.utils.ArrowConfusion;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class AhpcCriterion implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AHCPC-Cr";
    }

    @Override
    public String getDescription() {
        return "AHPC Criterion";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        double ap = new AdjacencyPrecision().getValue(trueGraph, estGraph, dataModel);
        double e = estGraph.getNumEdges();

        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());

        if (estGraph == null) throw new NullPointerException();

        int misorientationm = 0;

        for (Edge edge : trueGraph.getEdges()) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            Edge edge2 = estGraph.getEdge(x, y);

            if (edge2 == null) continue;

            if (edge.pointsTowards(x) && edge2.pointsTowards(y)) misorientationm++;
            else if (edge.pointsTowards(y) && edge2.pointsTowards(x)) misorientationm++;
        }

        int m = misorientationm;

        return ap * (1.0 - ap) * e;// ? 1.0 : 0.0;
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
