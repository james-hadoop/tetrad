package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class AhpBoundFewHubs implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AHP_Bound_Few_Hubs";
    }

    @Override
    public String getDescription() {
        return "Arrowhead Precision Bound Few Hubs";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {

        Set<Edge> U = new HashSet<>();
        Set<Edge> L = new HashSet<>();

        for (Node node : estGraph.getNodes()) {
            List<Node> adj = estGraph.getAdjacentNodes(node);
            if (adj.size() < 2) continue;

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> _adj = GraphUtils.asList(choice, adj);
                Node x = _adj.get(0);
                Node y = _adj.get(1);

                if (!estGraph.isAdjacentTo(x, y)) {
                    U.add(estGraph.getEdge(node, x));
                    U.add(estGraph.getEdge(node, y));

                    if (trueGraph.isAdjacentTo(x, y)) {
                        L.add(estGraph.getEdge(node, x));
                        L.add(estGraph.getEdge(node, y));
                    }
                }
            }
        }

        int f = 0;

        for (Edge e : L) {
            Edge e2 = trueGraph.getEdge(e.getNode1(), e.getNode2());
            if (!e.equals(e2)) {
                f++;
            }
        }

        double r = f / (double) L.size();

        // Count the number of arrows in estGraph
        int A = 0;

        for (Edge edge : estGraph.getEdges()) {
            if (edge.getEndpoint1() == Endpoint.ARROW) A++;
            if (edge.getEndpoint2() == Endpoint.ARROW) A++;
        }

        // Calculate the density of trueGraph;
        int e = trueGraph.getNumEdges();
        int v = trueGraph.getNumNodes();
        double d = 2 * e / (double) (v * (v - 1));

        // Calculate bound.
        return 1.0 - 2 * r * (U.size() / (double) A) * d;
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
