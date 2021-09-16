package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.RandomUtil;
import org.jetbrains.annotations.NotNull;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * Calculates the d-connection of X and Y conditional on variables Z, but allowing only directed paths the sum
 * of the time delays of edges clock in at under some given time limit T. The data are a directed (possibly cyclic)
 * graph, time delays > 0 for each edge in the graph, and a time limit > 0. Time delays for edges are given as a map
 * from Edge to Double.
 * <p>
 * The method for finding paths that clock in under T takes as input a graph g, nodes x, y, and set of nodes z, and
 * a time delay map, as above, and output all paths that clock in in under T. Three types of edges may be found
 * The first is directed paths from x to y conditional on z--that is causal influences from x to y conditional on z
 * that clock in under T. The second is treks--that is, paths from some w to both x and y, where each of the legs w to
 * x and w to y must clock in under the time limit. Unless w = x, these are not causal paths from x to y, but they are
 * d-connecting paths. The third allows for both such forks and also for colliders along paths from r to s and t to s
 * where s or some descendant of s is conditioned on. Note that for colliders where a descendant is conditioned on,
 * both directed paths from their sources to the collider must clock in under the time limit.
 * <p>
 * Methods for TConnection are provided, one in which times are taken to be random numbers between 20 and 80, and
 * another in which times may be provided by the user. The semantics is that x is t-connected to y just in case
 * there is a t-connecting path from x to y for which all directed paths clock in under the given time limit.
 * <p>
 * There are two properties for this class, one a PathType, set to DIRECT, TREK, or PATH, and the other a time limit.
 * The PathType is by default set to PATH, and the time limit is by default set to 1000. Note that as the time limit
 * increases without bound, if the path type is set to PATH, and if the times for edges are fixed, the isTConnectioned
 * method always reverts to the d-connection relation.
 *
 * @author jdramsey@andrew.cmu.edu 2021.1.4
 */
public class TConnection {
    private PathType pathType = PathType.ALL_PATHS;
    private double timeLimit = 1000;

    public List<LinkedList<Node>> findPaths(Graph g, Node x, Node y, List<Node> z, Map<Edge, Double> times) {
        if (!g.containsNode(x)) {
            throw new IllegalArgumentException("That graph doesn't contain the 'x' node.");
        }

        if (!g.containsNode(y)) {
            throw new IllegalArgumentException("That graph doesn't contain the 'y' node.");
        }

        for (Node _z : z) {
            if (!g.containsNode(_z)) {
                throw new IllegalArgumentException("That graph doesn't contain one of the conditioning nodes.");
            }
        }

        for (Edge e : g.getEdges()) {
            if (times.get(e) == null) {
                throw new IllegalArgumentException("Edge " + e + " is not assigned a time delay.");
            }

            if (times.get(e) <= 0) {
                throw new IllegalArgumentException("Edge " + e + " is assigned a time delay that is <= 0 ms: "
                        + times.get(e));
            }
        }

        LinkedList<Node> path = new LinkedList<>();
        path.add(x);
        List<LinkedList<Node>> paths = new ArrayList<>();
        Set<Triple> triples = new HashSet<>();
        findTemporalPathsVisit(g, null, x, y, z, path, paths, times, 0, triples, new double[]{0});
        return paths;
    }

    public boolean isTConnectedTo(Graph graph, Node x, Node y, List<Node> z) {
        return isTConnectedTo(graph, x, y, z, getRandomTimes(graph));
    }

    public boolean isTConnectedTo(Graph graph, Node x, Node y, List<Node> z, Map<Edge, Double> times) {
        List<LinkedList<Node>> temporalPaths = findPaths(graph, x, y, z, times);
        return !temporalPaths.isEmpty();
    }

    public void setPathType(PathType pathType) {
        this.pathType = pathType;
    }

    public double getTimeLimit() {
        return timeLimit;
    }

    public void setTimeLimit(double timeLimit) {
        if (timeLimit < 0) {
            throw new IllegalArgumentException("Time limit must be >= 0 ms: " + timeLimit);
        }

        this.timeLimit = timeLimit;
    }

    /**
     * Build the graph and times for TConnection record by record.
     */
    public static class Records {
        private final Graph graph;
        private final Map<Edge, Double> times;

        public Records() {
            graph = new EdgeListGraph();
            times = new HashMap<>();
        }

        public void addRecord(String x, String y, double time) {
            if (x == null) {
                throw new IllegalArgumentException("x is null.");
            }

            if (y == null) {
                throw new IllegalArgumentException("y is null.");
            }

            if (time <= 0) {
                throw new IllegalArgumentException("Time is <= 0 ms: " + time);
            }

            if (graph.getNode(x) == null) {
                graph.addNode(new GraphNode(x));
            }

            if (graph.getNode(y) == null) {
                graph.addNode(new GraphNode(y));
            }

            graph.addDirectedEdge(graph.getNode(x), graph.getNode(y));

            Edge e = graph.getEdge(graph.getNode(x), graph.getNode(y));

            times.put(e, time);
        }

        public Graph getGraph() {
            return graph;
        }

        public Map<Edge, Double> getTimes() {
            return times;
        }

        public String toString() {
            StringBuilder s = new StringBuilder();
            NumberFormat nf = new DecimalFormat("0.0000");

            List<Edge> edges = new ArrayList<>(graph.getEdges());

            s.append("Var1\tVar2\tDelay\n");

            for (int i = 0; i < edges.size(); i++) {
                Edge e = edges.get(i);
                s.append(e.getNode1().toString()).append("\t");
                s.append(e.getNode2().toString()).append("\t");
                s.append(nf.format(times.get(e))).append("\t");

                if (i < edges.size() - 1) {
                    s.append("\n");
                }
            }

            return s.toString();
        }

        public void toFile(Path path) throws IOException {
            FileWriter writer = new FileWriter(path.toFile());
            writer.write(toString());
            writer.close();
        }

        public static TConnection.Records fromFile(Path path) throws IOException {
            TConnection.Records records = new TConnection.Records();

            BufferedReader r = new BufferedReader(new FileReader(path.toFile()));
            String line;

            System.out.println("Skipping first line of file as header.");
            r.readLine();

            while ((line = r.readLine()) != null) {
                String[] tokens = line.split("\t");

                if (tokens.length != 3) {
                    throw new IllegalArgumentException("Expecting 3 tokens: " + line);
                }

                String v1 = tokens[0];
                String v2 = tokens[1];
                double time = Double.parseDouble(tokens[2]);

                records.addRecord(v1, v2, time);
            }

            return records;
        }
    }

    public static String pathString(Graph graph, List<Node> path, List<Node> conditioningVars) {
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

            if (n2 != null && TConnection.fork(n0, n1, n2, graph)) {
                buf.append("(F)");
            }
        }

        return buf.toString();
    }

    public static boolean directed(List<Node> path, Graph graph) {
        for (int i = 0; i < path.size() - 1; i++) {
            Edge e = graph.getEdge(path.get(i), path.get(i + 1));
            if (!e.pointsTowards(path.get(i + 1))) {
                return false;
            }
        }

        return true;
    }

    public static double sumTimes(List<Node> path, Graph graph, Map<Edge, Double> times) {
        double sum = 0.0;

        for (int i = 0; i < path.size() - 1; i++) {
            Edge e = graph.getEdge(path.get(i), path.get(i + 1));

            if (e.pointsTowards(path.get(i))) {
                sum -= times.get(e);
            } else {
                sum += times.get(e);
            }
        }

        return sum;
    }

    public static void printPaths(Graph graph, Node x, Node y, List<Node> cond, List<LinkedList<Node>> temporalPaths,
                                  double timeLimit, Map<Edge, Double> times) {
        System.out.println("Paths from " + x + " to " + y + (cond.isEmpty() ? (" conditioning on " + cond) : ""));
        System.out.println("Time limit = " + timeLimit + " ms\n");

        if (temporalPaths.isEmpty()) {
            System.out.println("NO PATHS");
        } else {
            NumberFormat nf = new DecimalFormat("0.00");

            for (List<Node> path : temporalPaths) {
                if (directed(path, graph)) {
                    System.out.print(nf.format(sumTimes(path, graph, times)) + ":\t");
                }
                System.out.println(pathString(graph, path, cond));
            }

            System.out.println();
        }
    }

    private boolean findTemporalPathsVisit(Graph g, Node prev, Node x, Node to, List<Node> z,
                                           LinkedList<Node> path, List<LinkedList<Node>> paths,
                                           Map<Edge, Double> times,
                                           double sum, Set<Triple> triples, double[] descTime) {
        List<Node> adjacentNodes = g.getAdjacentNodes(x);

        boolean foundCond = false;

        for (Node c : adjacentNodes) {
            foundCond = foundCond || z.contains(x);

            if (prev != null && prev == c) continue;

            boolean fork = fork(prev, x, c, g);
            boolean collider = collider(prev, x, c, g);

            if (fork && triples.contains(new Triple(prev, x, c))) continue;
            if (collider && triples.contains(new Triple(prev, x, c))) continue;

            if (fork) triples.add(new Triple(prev, x, c));
            if (collider) triples.add(new Triple(prev, x, c));

            if (c == prev) continue;
            Edge edge = g.getEdge(x, c);
            double time = times.get(g.getEdge(x, c));

            if (prev != null) {
                if (!((collider && (foundCond)) || (!collider && !z.contains(x)))) {
//                    sum -= time;
                    continue;
                }
            }

            if (pathType == PathType.DIRECT && (x != c && !edge.pointsTowards(c))) {
                continue;
            }

            sum += time;
            path.addLast(c);

            if (sum > timeLimit || sum <= 0) {
                sum -= time;
                path.removeLast();
                continue;
            }

            if (c == to) {
                paths.add(new LinkedList<>(path));
            }

            if (pathType == PathType.TREK && fork) {
                foundCond = findTemporalPathsVisit(g, x, c, to, z, path, paths, times, 0, triples, descTime);
            } else if (pathType == PathType.ALL_PATHS && fork) {
                foundCond = findTemporalPathsVisit(g, x, c, to, z, path, paths, times, 0, triples, descTime);
            } else if (pathType == PathType.ALL_PATHS && collider) {
                foundCond = findTemporalPathsVisit(g, x, c, to, z, path, paths, times, descTime[0], triples, descTime);
            } else {
                foundCond = findTemporalPathsVisit(g, x, c, to, z, path, paths, times, sum, triples, descTime);
            }

            if (pathType == PathType.ALL_PATHS && prev != null && !(fork || collider)) {
                if (foundCond) {
                    descTime[0] += time;
                    sum -= time;
                    path.removeLast();
                    continue;
                }
            }

            sum -= time;
            path.removeLast();
        }

        return foundCond;
    }

    private static boolean collider(Node prev, Node x, Node c, Graph graph) {
        if (prev == null) return false;

        Edge e1 = graph.getEdge(prev, x);
        Edge e2 = graph.getEdge(x, c);

        if (e1 == null || e2 == null) {
            throw new IllegalArgumentException("Not an edge in the graph.");
        }

        return e1.pointsTowards(x) && e2.pointsTowards(x);
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

    @NotNull
    private Map<Edge, Double> getRandomTimes(Graph graph) {
        Map<Edge, Double> times = new HashMap<>();

        for (Edge edge : graph.getEdges()) {
            times.put(edge, RandomUtil.getInstance().nextInt(60) + 20.);
        }

        return times;
    }

    public enum PathType {DIRECT, TREK, ALL_PATHS}
}
