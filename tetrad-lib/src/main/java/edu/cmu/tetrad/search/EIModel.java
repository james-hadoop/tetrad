package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 *
 */
public class EIModel {
    private final Graph graph;
    private final Map<Edge, Double> times;
    private final Map<Edge, Integer> excitements;
    private final TConnection tconn;

    public EIModel(Graph graph, Map<Edge, Double> times, Map<Edge, Integer> excitatory, double timeLimit) {
        this.graph = graph;
        this.times = times;
        this.excitements = excitatory;
        this.tconn = new TConnection();
        this.tconn.setPathType(TConnection.PathType.DIRECT);
        this.tconn.setTimeLimit(timeLimit);
    }

    public int getWellenPrediction(Node x, Node y, List<Node> z) {
        List<LinkedList<Node>> paths = tconn.findTemporalPaths(graph, x, y, z, times);

        boolean allE = true;
        boolean allI = true;

        for (LinkedList<Node> path : paths) {
            int e = 1;

            for (int i = 0; i < path.size() - 1; i++) {
                Edge edge = graph.getEdge(path.get(i), path.get(i + 1));
                if (excitements.get(edge) == 1) {
                    e *= 1;
                } else if (excitements.get(edge) == 2) {
                    e *= -1;
                } else if (excitements.get(edge) == 3) {
                    allI = false;
                    allE = false;
                }
            }

            if (e == 1) {
                allI = false;
            }

            if (e == -1) {
                allE = false;
            }
        }

        if (allE) {
            return 1;
        } else if (allI) {
            return 2;
        } else {
            return 3;
        }
    }
}
