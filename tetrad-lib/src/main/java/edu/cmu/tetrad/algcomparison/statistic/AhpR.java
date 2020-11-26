package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Returns the "r" statistic for calculating the bound for AHPC.
 *
 * @author jdramsey
 */
public class AhpR implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AHP_Rho";
    }

    @Override
    public String getDescription() {
        return "Arrowhead Precision R";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {

        // Count the number of unshielded triples in estGraph.
        Set<Edge> L = new HashSet<>();

        for (Node node : estGraph.getNodes()) {
            List<Node> adj = estGraph.getAdjacentNodes(node);

            for (int i = 0; i < adj.size(); i++) {
                for (int j = i + 1; j < adj.size(); j++) {
                    Node x = adj.get(i);
                    Node y = adj.get(j);

                    // In U
                    if (!estGraph.isAdjacentTo(x, y)) {

                        // In L
                        if (trueGraph.isAdjacentTo(x, y)) {// && trueGraph.isAdjacentTo(node, x) && trueGraph.isAdjacentTo(node, y)) {
                            L.add(estGraph.getEdge(node, x));
                            L.add(estGraph.getEdge(node, y));
                        }
                    }
                }
            }
        }

        // f = edges in L that are misoriented or missing in trueGraph
        int f = 0;

        for (Edge e1 : L) {
            Node x = e1.getNode1();
            Node y = e1.getNode2();

            Edge e2 = trueGraph.getEdge(x, y);

            if (e2 == null) f++;
            else if (e1.pointsTowards(x) && e2.pointsTowards(y)) f++;
            else if (e1.pointsTowards(y) && e2.pointsTowards(x)) f++;
        }

        return f / (double) L.size();
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
