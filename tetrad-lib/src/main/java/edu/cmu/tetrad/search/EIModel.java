package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;

import java.io.*;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * An excitatory/inhibitory model, for TMS data.
 *
 * @author jdramsey@andrew.cmu.edu 2021.1.4
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
        this.tconn.setTimeLimit(timeLimit);
    }

    public Return getQualitativePrediction(Node x, Node y, List<Node> z) {
        tconn.setPathType(TConnection.PathType.DIRECT);
        List<LinkedList<Node>> paths = new LinkedList<>();

        if (x.getName().equals("?") || y.getName().equals("?")) {
            if (x.getName().equals("?") && !y.getName().equals("?")) {
                for (Node _x : graph.getNodes()) {
                    paths.addAll(tconn.findPaths(graph, _x, y, z, times));
                }
            } else if (!x.getName().equals("?") && y.getName().equals("?")) {
                for (Node _y : graph.getNodes()) {
                    paths.addAll(tconn.findPaths(graph, x, _y, z, times));
                }
            } else if (x.getName().equals("?") && y.getName().equals("?")) {
                for (Node _x : graph.getNodes()) {
                    for (Node _y : graph.getNodes()) {
                        paths.addAll(tconn.findPaths(graph, _x, _y, z, times));
                    }
                }
            }
        } else {
            paths.addAll(tconn.findPaths(graph, x, y, z, times));
        }

        LinkedList<Integer> pathExcitements = new LinkedList<>();

        if (paths.isEmpty()) return new Return(paths, pathExcitements, 0);

        boolean allI = true;
        boolean allE = true;

        for (LinkedList<Node> path : paths) {
            int prod = 1;

            for (int i = 0; i < path.size() - 1; i++) {
                Edge edge = graph.getEdge(path.get(i), path.get(i + 1));

                if (excitements.get(edge) == 1) {
                    prod *= 1;
                } else if (excitements.get(edge) == 2) {
                    prod *= -1;
                } else if (excitements.get(edge) == 3) {
                    prod = 0;
                } else {
                    throw new IllegalStateException();
                }
            }

            if (prod == 1) {
                pathExcitements.add(1);
                allE = false;
            }

            if (prod == -1) {
                pathExcitements.add(2);
                allI = false;
            }

            if (prod == 0) {
                pathExcitements.add(3);
                allE = false;
                allI = false;
            }
        }

        if (allI && !allE) {
            return new Return(paths, pathExcitements, 1);
        } else if (!allI && allE) {
            return new Return(paths, pathExcitements, 2);
        } else if (!allI /*&& !allE*/) {
            return new Return(paths, pathExcitements, 3);
        } else {
            throw new IllegalStateException();
        }
    }

    /**
     * Build the graph, times, and excitements for the EI model record by record.
     */
    public static class Records {
        private final Graph graph;
        private final Map<Edge, Double> times;
        private final Map<Edge, Integer> excitements;

        public Records() {
            graph = new EdgeListGraph();
            times = new HashMap<>();
            excitements = new HashMap<>();
        }

        public void addRecord(String var1, String var2, double time, Integer excitement) {
            if (var1 == null) {
                throw new IllegalArgumentException("var1 is null.");
            }

            if (var2 == null) {
                throw new IllegalArgumentException("var2 is null.");
            }

            if (time <= 0) {
                throw new IllegalArgumentException("time is <= 0: " + time);
            }

            if (excitement < 1 || excitement > 3) {
                throw new IllegalArgumentException("excitement must be an integer in [1, 3].");
            }

            if (graph.getNode(var1) == null) {
                graph.addNode(new GraphNode(var1));
            }

            if (graph.getNode(var2) == null) {
                graph.addNode(new GraphNode(var2));
            }

            graph.addDirectedEdge(graph.getNode(var1), graph.getNode(var2));

            Edge e = graph.getEdge(graph.getNode(var1), graph.getNode(var2));

            times.put(e, time);
            excitements.put(e, excitement);
        }

        public Graph getGraph() {
            return graph;
        }

        public Map<Edge, Double> getTimes() {
            return times;
        }

        public Map<Edge, Integer> getExcitements() {
            return excitements;
        }

        public String toString() {
            StringBuilder s = new StringBuilder();
            NumberFormat nf = new DecimalFormat("0");

            List<Edge> edges = new ArrayList<>(graph.getEdges());

            s.append("Source\tTarget\tDelay\tExcitatory\n");

            for (int i = 0; i < edges.size(); i++) {
                Edge e = edges.get(i);
                s.append(e.getNode1().toString()).append("\t");
                s.append(e.getNode2().toString()).append("\t");
                s.append(nf.format(times.get(e))).append("\t");

                switch (excitements.get(e)) {
                    case 0:
                        s.append("0");
                        break;
                    case 1:
                        s.append("+");
                        break;
                    case 2:
                        s.append("-");
                        break;
                    case 3:
                        s.append("+/-");
                        break;
                    default:
                        throw new IllegalArgumentException("Expecting 0, -, +, or +/- for excitement: " + excitements.get(e));
                }

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

        public static Records fromFile(Path path) throws IOException {
            Records records = new Records();

            BufferedReader r = new BufferedReader(new FileReader(path.toFile()));

//            System.out.println("Skipping first line of file as header.");

            r.readLine();

            String line;

            while ((line = r.readLine()) != null) {
                String[] tokens = line.split("\t");

                if (tokens.length != 4) {
                    throw new IllegalArgumentException("Expecting 4 tokens: " + line);
                }

                String v1 = tokens[0];
                String v2 = tokens[1];
                double time = Double.parseDouble(tokens[2]);

                String excitement = tokens[3];
                int _excitement;

                switch (excitement) {
                    case "0":
                        _excitement = 0;
                        break;
                    case "+":
                        _excitement = 1;
                        break;
                    case "-":
                        _excitement = 2;
                        break;
                    case "+/-":
                        _excitement = 3;
                        break;
                    default:
                        throw new IllegalArgumentException("Expecting 0, -, +, or +/- for excitement: " + excitement);
                }

                records.addRecord(v1, v2, time, _excitement);
            }

            return records;
        }
    }

    public static void printPaths(Graph graph, Node x, Node y, List<Node> z, Return wellmanPrediction,
                                  double timeLimit, Map<Edge, Double> times) {
//        System.out.println("\n==========================================");
        System.out.println();
        NumberFormat nf = new DecimalFormat("0");

        TConnection tconn = new TConnection();
        tconn.setPathType(TConnection.PathType.DIRECT);
        tconn.setTimeLimit(timeLimit);

        List<LinkedList<Node>> paths = new LinkedList<>();

        if (x.getName().equals("?") || y.getName().equals("?")) {
            if (x.getName().equals("?") && !y.getName().equals("?")) {
                for (Node _x : graph.getNodes()) {
                    paths.addAll(tconn.findPaths(graph, _x, y, z, times));
                }
            } else if (!x.getName().equals("?") && y.getName().equals("?")) {
                for (Node _y : graph.getNodes()) {
                    paths.addAll(tconn.findPaths(graph, x, _y, z, times));
                }
            } else if (x.getName().equals("?") && y.getName().equals("?")) {
                for (Node _x : graph.getNodes()) {
                    for (Node _y : graph.getNodes()) {
                        paths.addAll(tconn.findPaths(graph, _x, _y, z, times));
                    }
                }
            }
        } else {
            paths.addAll(tconn.findPaths(graph, x, y, z, times));
        }

        System.out.println("Paths from " + x + " to " + y + (z.isEmpty() ? (" conditioning on " + z) : ""));
        System.out.println("Time limit = " + nf.format(timeLimit) + " ms\n");

        List<LinkedList<Node>> temporalPaths = wellmanPrediction.getPaths();

        if (temporalPaths.isEmpty()) {
            System.out.println("NO PATHS");
        } else {
            for (LinkedList<Node> path : paths) {
                System.out.print("(" + nf.format(TConnection.sumTimes(path, graph, times)) + " ms) ");
                System.out.print(TConnection.pathString(graph, path, z));

                if (temporalPaths.contains(path)) {
                    System.out.print("  ");

                    switch (wellmanPrediction.getPathExcitements().get(temporalPaths.indexOf(path))) {
                        case 0:
                            System.out.print("0");
                            break;
                        case 1:
                            System.out.print("+");
                            break;
                        case 2:
                            System.out.print("-");
                            break;
                        case 3:
                            System.out.print("+/-");
                            break;
                        default:
                            throw new IllegalStateException("Path prediction should be 0, 1, 2, or 3: " + wellmanPrediction);
                    }
                }

                System.out.println();
            }
        }
    }

    private static Node node(Graph graph, String x1) {
        return graph.getNode(x1);
    }

    public static void printEiResult(String file, String timeLimit, String x, String y,
                                     String... cond) {
        Records records = null;

        try {
            records = Records.fromFile(new File(file).toPath());
//            System.out.println(records);
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (records == null) {
            throw new NullPointerException("Parsing did not recover records.");
        }

        Graph graph = records.getGraph();
        Map<Edge, Double> times = records.getTimes();
        Map<Edge, Integer> excitements = records.getExcitements();

        Node _x = x.equals("?") ? new GraphNode("?") : node(graph, x);
        Node _y = y.equals("?") ? new GraphNode("?") : node(graph, y);

        if (_x == null || _y == null) {
            throw new IllegalArgumentException("Source or target not a node in this graph or a wildcard.");
        }

        double _timeLimit = Double.parseDouble(timeLimit);

        EIModel wellman = new EIModel(graph, times, excitements, _timeLimit);

        List<Node> _cond = new ArrayList<>();
        for (String s : cond) _cond.add(graph.getNode(s));

        Return wellmanPrediction = wellman.getQualitativePrediction(_x, _y, _cond);

        printPaths(graph, _x, _y, _cond, wellmanPrediction, _timeLimit, times);

        System.out.print("\nQualitative excitatory/inhibitory prediction for " + x + " >>> " + y + ": ");

        switch (wellmanPrediction.getPrediction()) {
            case 0:
                System.out.println("0");
                break;
            case 1:
                System.out.println("+");
                break;
            case 2:
                System.out.println("-");
                break;
            case 3:
                System.out.println("+/-");
                break;
            default:
                throw new IllegalStateException("Prediction should be 0, 1, 2, or 3: " + wellmanPrediction);
        }
    }

    public static class Return {
        private final List<LinkedList<Node>> paths;
        private final List<Integer> pathExcitements;
        private final int excitement;

        public Return(List<LinkedList<Node>> paths, List<Integer> pathExcitements, int predictedExcitement) {
            this.paths = paths;
            this.pathExcitements = pathExcitements;
            this.excitement = predictedExcitement;
        }

        public List<LinkedList<Node>> getPaths() {
            return paths;
        }

        public List<Integer> getPathExcitements() {
            return pathExcitements;
        }

        public int getPrediction() {
            return excitement;
        }
    }

    public static void main(String... args) {
        if (args.length < 4) throw new IllegalArgumentException("Need at least 4 arguments");
        String file = args[0];
        String timeLimit = args[1];
        String x = args[2];
        String y = args[3];
        String[] cond = new String[args.length - 4];
        System.arraycopy(args, 4, cond, 0, cond.length);
        try {
            printEiResult(file, timeLimit, x, y, cond);
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
