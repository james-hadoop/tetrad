package edu.cmu.tetrad.test;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.TConnection;
import edu.cmu.tetrad.util.RandomUtil;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Tests the TConnection model.
 *
 * @author jdramsey@andrew.cmu.edu 2021.1.4
 */
public class TestTConnection {

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
//        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));
        Graph graph = GraphUtils.randomGraph(15, 0, 15, 100, 100, 100, true);

        // Helper class to build up the necessary graph and time map from data. Here we pick random times.
        TConnection.Records records = new TConnection.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    RandomUtil.getInstance().nextUniform(20, 80));
        }

        try {
            records.toFile(new File("/Users/josephramsey/Downloads/records-t-connection.txt").toPath());
            records = TConnection.Records.fromFile(new File("/Users/josephramsey/Downloads/records-t-connection.txt").toPath());
            System.out.println(records);
        } catch (IOException e) {
            e.printStackTrace();
        }

        graph = records.getGraph();
        Map<Edge, Double> times = records.getTimesAvg();

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

        List<LinkedList<Node>> temporalPaths = tconn.findPaths(graph, x5, x15, cond, times);
        TConnection.printPaths(graph, x5, x15, cond, temporalPaths, tconn.getTimeLimit(), times);

//        assert (temporalPaths.contains(path(x5, x11, x13, x15)));
//        assert (temporalPaths.contains(path(x5, x7, x8, x15)));
//        assert (temporalPaths.contains(path(x5, x7, x9, x11, x13, x15)));
    }

    @Test
    public void test2() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());

        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        TConnection.Records records = new TConnection.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    RandomUtil.getInstance().nextUniform(20, 80));
        }

        graph = records.getGraph();
        Map<Edge, Double> times = records.getTimesAvg();

        Node x6 = node(graph, "X6");
        Node x3 = node(graph, "X3");
        Node x12 = node(graph, "X12");

        List<Node> cond = new ArrayList<>();
        cond.add(node(graph, "X12"));

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.ALL_PATHS);
        tconn.setTimeLimit(1000);

        List<LinkedList<Node>> temporalPaths = tconn.findPaths(graph, x6, x3, cond, times);
        TConnection.printPaths(graph, x6, x3, cond, temporalPaths, tconn.getTimeLimit(), times);

        assert (temporalPaths.contains(path(x6, x12, x3)));

    }

    @Test
    public void test3() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());

        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph8.txt"));

        TConnection.Records records = new TConnection.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    RandomUtil.getInstance().nextUniform(20, 80));
        }

        graph = records.getGraph();
        Map<Edge, Double> times = records.getTimesAvg();

        Node x1 = node(graph, "X1");
        Node x2 = node(graph, "X2");
        Node x3 = node(graph, "X3");
        Node x4 = node(graph, "X4");

        List<Node> cond = new ArrayList<>();

        TConnection tconn = new TConnection();

        tconn.setPathType(TConnection.PathType.DIRECT);
        tconn.setTimeLimit(300);

        List<LinkedList<Node>> temporalPaths = tconn.findPaths(graph, x1, x2, cond, times);
        TConnection.printPaths(graph, x1, x2, cond, temporalPaths, tconn.getTimeLimit(), times);

        assert (temporalPaths.contains(path(x1, x2)));
        assert (temporalPaths.contains(path(x1, x2, x3, x4, x1, x2)));

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

        tconn.setPathType(TConnection.PathType.ALL_PATHS);
        tconn.setTimeLimit(600);

        List<LinkedList<Node>> temporalPaths = tconn.findPaths(graph, x1, x2, cond, times);
        TConnection.printPaths(graph, x1, x2, cond, temporalPaths, tconn.getTimeLimit(), times);

        assert (temporalPaths.contains(path(x1, x3, x2)));
    }

    @Test
    public void test5() {
        try {
            TConnection.Records records = TConnection.Records.fromFile(new File("/Users/josephramsey/Downloads/data.latency.txt").toPath());
            System.out.println(records);

            TConnection tconn = new TConnection();

            tconn.setPathType(TConnection.PathType.DIRECT);
            tconn.setTimeLimit(300);

            Graph graph = records.getGraph();
            List<Node> nodes = graph.getNodes();

            for (int i = 0; i < nodes.size(); i++) {
                for (int j = 0; j < nodes.size(); j++) {
                    Node x1 = nodes.get(i);
                    Node x2 = nodes.get(j);
                    List<Node> cond = new ArrayList<>();
                    List<LinkedList<Node>> temporalPaths = tconn.findPaths(graph, x1, x2, cond, records.getTimesAvg());

                    System.out.println(x1 + " " + x2 + " " + temporalPaths);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Node node(Graph graph, String x1) {
        return graph.getNode(x1);
    }

    private LinkedList<Node> path(Node... nodes) {
        LinkedList<Node> path = new LinkedList<>();
        for (Node node : nodes) path.addLast(node);
        return path;
    }

}
