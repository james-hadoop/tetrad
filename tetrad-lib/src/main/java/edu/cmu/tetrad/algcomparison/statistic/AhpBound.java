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
public class AhpBound implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AHP_Bound";
    }

    @Override
    public String getDescription() {
        return "Arrowhead Precision Bound";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {

        AhpR ahpR = new AhpR();
        double r = ahpR.getValue(trueGraph, estGraph, dataModel);

        Set<Edge> U = new HashSet<>();
        Set<Edge> L = new HashSet<>();

        for (Node node : estGraph.getNodes()) {
            List<Node> adj = estGraph.getAdjacentNodes(node);

            for (int i = 0; i < adj.size(); i++) {
                for (int j = i + 1; j < adj.size(); j++) {
                    Node x = adj.get(i);
                    Node y = adj.get(j);

                    // Unshielded triple
                    if (!estGraph.isAdjacentTo(x, y)) {
                        U.add(estGraph.getEdge(node, x));
                        U.add(estGraph.getEdge(node, y));

                        // In L
                        if (trueGraph.isAdjacentTo(x, y)) {// && trueGraph.isAdjacentTo(node, x) && trueGraph.isAdjacentTo(node, y)) {
                            L.add(estGraph.getEdge(node, x));
                            L.add(estGraph.getEdge(node, y));
                        }
                    }
                }
            }
        }

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

//        d = L.size() / (double) U.size();

        int f = 0;

        for (Edge e1 : L) {
            Node x = e1.getNode1();
            Node y = e1.getNode2();

            Edge e2 = trueGraph.getEdge(x, y);

            if (e2 == null) f++;
            else if (e1.pointsTowards(x) && e2.pointsTowards(y)) f++;
            else if (e1.pointsTowards(y) && e2.pointsTowards(x)) f++;
        }

        double bound;

        if (L.isEmpty() || U.isEmpty() || A == 0) bound = Double.NaN;
        else bound = 1.0 - (r) * (L.size() / (double) U.size()) * (U.size() / (double) A);

        System.out.println("f = " + f + " L = " + L.size() + " U = " + U.size() + " A = " + A + " bound = " + bound
            + " U / L = " + L.size() / (double) U.size() + " d = " + d);

        // Calculate bound.
        return bound;
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
