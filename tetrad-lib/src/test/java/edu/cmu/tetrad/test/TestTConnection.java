package edu.cmu.tetrad.test;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.TConnection;
import edu.cmu.tetrad.util.RandomUtil;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.io.File;
import java.util.*;

public class TestTConnection {

    private static String pathString(Graph graph, List<Node> path, List<Node> conditioningVars) {
        StringBuilder buf = new StringBuilder();

        if (path.size() < 2) {
            return "";
        }

        buf.append(path.get(0).toString());

        if (conditioningVars.contains(path.get(0))) {
            buf.append("(C)");
        }

        for (int m = 1; m < path.size(); m++) {
            Node n0 = path.get(m - 1);
            Node n1 = path.get(m);
            Node n2 = null;

            if (m + 1 < path.size()) {
                n2 = path.get(m + 1);
            }

            Edge edge = graph.getEdge(n0, n1);

            if (edge == null) {
                buf.append("(-)");
            } else {
                if (n0 == n1) {
                    buf.append("-[").append(n0).append("-->").append(n1).append("]");
                    continue;
                }

                Endpoint endpoint0 = edge.getProximalEndpoint(n0);
                Endpoint endpoint1 = edge.getProximalEndpoint(n1);

                if (endpoint0 == Endpoint.ARROW) {
                    buf.append("<");
                } else if (endpoint0 == Endpoint.TAIL) {
                    buf.append("-");
                } else if (endpoint0 == Endpoint.CIRCLE) {
                    buf.append("o");
                }

                buf.append("-");

                if (endpoint1 == Endpoint.ARROW) {
                    buf.append(">");
                } else if (endpoint1 == Endpoint.TAIL) {
                    buf.append("-");
                } else if (endpoint1 == Endpoint.CIRCLE) {
                    buf.append("o");
                }
            }

            buf.append(n1.toString());

            if (conditioningVars.contains(n1)) {
                buf.append("(C)");
            }

            if (n2 != null && fork(n0, n1, n2, graph)) {
                buf.append("(F)");
            }
        }

        return buf.toString();
    }

    private void printPaths(Graph graph, Node x, Node y, List<Node> cond, List<LinkedList<Node>> temporalPaths,
                            double timeLimit) {
        System.out.println("\n==========================================");
        System.out.println();

        System.out.println("Paths from " + x + " to " + y + " conditioning on " + (cond.isEmpty() ? "NOTHING" : cond));
        System.out.println("Time limit = " + timeLimit + " ms\n");

        if (temporalPaths.isEmpty()) {
            System.out.println("NO PATHS\n");
        } else {
            for (List<Node> path : temporalPaths) {
                System.out.println(pathString(graph, path, cond));
            }

            System.out.println();
        }
    }

    @NotNull
    private Map<Edge, Double> getRandomTimes(Graph graph) {
        Map<Edge, Double> times = new HashMap<>();

        for (Edge edge : graph.getEdges()) {
            times.put(edge, RandomUtil.getInstance().nextInt(60) + 20.);
        }
        return times;
    }

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());
        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        Map<Edge, Double> times = getRandomTimes(graph);

        Node x5 = node(graph, "X5");
        Node x15 = node(graph, "X15");
        Node x11 = node(graph, "X11");
        Node x13 = node(graph, "X13");
        Node x7 = node(graph, "X7");
        Node x8 = node(graph, "X8");
        Node x9 = node(graph, "X9");

        List<Node> cond = new ArrayList<>();

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.DIRECT);
        tconn.setTimeLimit(300);

        List<LinkedList<Node>> temporalPaths = tconn.findTemporalPaths(graph, x5, x15, cond, times);
        printPaths(graph, x5, x15, cond, temporalPaths, tconn.getTimeLimit());

        assert(temporalPaths.contains(path(x5, x11, x13, x15)));
        assert(temporalPaths.contains(path(x5, x7, x8, x15)));
        assert(temporalPaths.contains(path(x5, x7, x9, x11, x13, x15)));
    }

    @Test
    public void test2() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());

        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        Map<Edge, Double> times = getRandomTimes(graph);

        Node x6 = node(graph, "X6");
        Node x3 = node(graph, "X3");
        Node x12 = node(graph, "X12");

        List<Node> cond = new ArrayList<>();
        cond.add(node(graph, "X12"));

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.PATH);
        tconn.setTimeLimit(1000);

        List<LinkedList<Node>> temporalPaths = tconn.findTemporalPaths(graph, x6, x3, cond, times);
        printPaths(graph, x6, x3, cond, temporalPaths, tconn.getTimeLimit());

        assert(temporalPaths.contains(path(x6, x12, x3)));

    }

    @Test
    public void test3() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());

        Graph g = GraphUtils.loadGraphTxt(new File("src/test/resources/graph8.txt"));

        Map<Edge, Double> times = getRandomTimes(g);

        Node x1 = node(g, "X1");
        Node x2 = node(g, "X2");
        Node x3 = node(g, "X3");
        Node x4 = node(g, "X4");

        List<Node> cond = new ArrayList<>();

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.DIRECT);
        tconn.setTimeLimit(300);

        List<LinkedList<Node>> temporalPaths = tconn.findTemporalPaths(g, x1, x2, cond, times);
        printPaths(g, x1, x2, cond, temporalPaths, tconn.getTimeLimit());

        assert(temporalPaths.contains(path(x1, x2)));
        assert(temporalPaths.contains(path(x1, x2, x3, x4, x1, x2)));

    }

    public Node node(Graph graph, String x1) {
        return graph.getNode(x1);
    }

    @Test
    public void test4() {
        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");
        Node x5 = new GraphNode("X5");
        Node x6 = new GraphNode("X6");
        Node x7 = new GraphNode("X7");
        Node x8 = new GraphNode("X8");

        List<Node> nodes = new ArrayList<>();
        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);
        nodes.add(x5);
        nodes.add(x6);
        nodes.add(x7);
        nodes.add(x8);

        Graph graph = new EdgeListGraph(nodes);

        graph.addDirectedEdge(x1, x3);
        graph.addDirectedEdge(x2, x3);

        graph.addDirectedEdge(x3, x4);
        graph.addDirectedEdge(x4, x5);
        graph.addDirectedEdge(x5, x6);
        graph.addDirectedEdge(x6, x7);
        graph.addDirectedEdge(x7, x8);


        Map<Edge, Double> times = getRandomTimes(graph);

        List<Node> cond = new ArrayList<>();
        cond.add(x8);

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.PATH);
        tconn.setTimeLimit(600);

        List<LinkedList<Node>> temporalPaths = tconn.findTemporalPaths(graph, x1, x2, cond, times);
        printPaths(graph, x1, x2, cond, temporalPaths, tconn.getTimeLimit());

        assert(temporalPaths.contains(path(x1, x3, x2)));
    }

    private LinkedList<Node> path(Node...nodes) {
        LinkedList<Node> path = new LinkedList<>();
        for (Node node : nodes) path.addLast(node);
        return path;
    }

    private static boolean fork(Node prev, Node x, Node c, Graph graph) {
        if (prev == null) return false;

        Edge e1 = graph.getEdge(prev, x);
        Edge e2 = graph.getEdge(x, c);

        if (e1 == null || e2 == null) {
            throw new IllegalArgumentException("Not an edge in the graph: <" + prev + ", " + x + ", " + c + ">");
        }

        return e1.pointsTowards(prev) && e2.pointsTowards(c);
    }
}
