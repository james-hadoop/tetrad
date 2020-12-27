package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;

import java.io.*;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 *
 */
public class EIModel {
    private final Graph graph;
    private final Map<Edge, Double> times;
    private final Map<Edge, Integer> excitements;
    private final TConnection tconn;

    private int[][] prod = {{1, 2, 0, 3}, {2, 1, 0, 3}, {0, 0, 0, 0}, {3, 3, 0, 3}};
    private int[][] sum = {{1, 2, 0, 3}, {2, 1, 0, 3}, {0, 0, 0, 0}, {3, 3, 0, 3}};

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

        if (paths.isEmpty()) return 0;

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

        if (allE && !allI) {
            return 1;
        } else if (allI && !allE) {
            return 2;
        } else {
            return 3;
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
            NumberFormat nf = new DecimalFormat("0.0000");

            List<Edge> edges = new ArrayList<>(graph.getEdges());

            s.append("Var1\tVar2\tDelay\tExcitatory\n");

            for (int i = 0; i < edges.size(); i++) {
                Edge e= edges.get(i);
                s.append(e.getNode1().toString()).append("\t");
                s.append(e.getNode2().toString()).append("\t");
                s.append(nf.format(times.get(e))).append("\t");
                s.append(excitements.get(e));

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

            System.out.println("Skipping first line of file as header.");

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
                int excitement = Integer.parseInt(tokens[3]);

                records.addRecord(v1, v2, time, excitement);
            }

            return records;
        }
    }
}
