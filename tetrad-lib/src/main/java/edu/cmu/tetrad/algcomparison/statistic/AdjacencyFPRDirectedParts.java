package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.AdjacencyConfusion;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;

/**
 * The adjacency true positive rate. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class AdjacencyFPRDirectedParts implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AFPRDP";
    }

    @Override
    public String getDescription() {
        return "Adjacency False Positive Rate Directed Parts";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        EdgeListGraph true2 = new EdgeListGraph(trueGraph.getNodes());
        EdgeListGraph est2 = new EdgeListGraph(estGraph.getNodes());

        Graph truePattern = SearchGraphUtils.patternForDag(trueGraph);

        for (Edge edge : truePattern.getEdges()) {
            if (edge.isDirected()) {
                true2.addEdge(edge);
                if (estGraph.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                    est2.addEdge(trueGraph.getEdge(edge.getNode1(), edge.getNode2()));
                }
            }
        }

        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(true2, est2);
        int adjFp = adjConfusion.getAdjFp();
        int adjTn = adjConfusion.getAdjTn();
        return adjFp / (double) (adjFp + adjTn);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
