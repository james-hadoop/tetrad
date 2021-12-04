package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.algorithm.mixed.Mgm;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.*;
import edu.cmu.tetrad.algcomparison.independence.ConditionalGaussianLRT;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ConditionalGaussianBicScore;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import edu.pitt.dbmi.data.reader.tabular.MixedTabularDatasetFileReader;
import org.junit.Test;

import java.io.*;
import java.util.*;

public class TestIrvine {


    @Test
    public void test1() {
        File dir = new File("/Users/josephramsey/Downloads/GitHub/tetrad-example-analyses");

        File[] files = dir.listFiles();

        for (File dataDir : files) {


            if (dataDir.getName().startsWith("README")) continue;
            if (dataDir.getName().startsWith("LICENSE")) continue;
            if (dataDir.getName().startsWith(".")) continue;

            if (!dataDir.getName().contains("auto-mpg")) continue;
            ;

            File dataDir2 = new File(dataDir, "tetrad.formatted.data");

            File[] obj = dataDir2.listFiles();

            if (obj == null) {
                System.out.println("Couldn't fine the directory for " + dataDir.getName());
            }

            if (obj != null) {
                for (File f : obj) {
                    if (f.getName().startsWith(".")) continue;

                    String name = f.getName();

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.PC, Algorithm.CPC, Algorithm.PCMAX};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FGES, Algorithm.BOSS};

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FCI, Algorithm.FCIMAX, Algorithm.RFCI};
                    Algorithm[] algorithms = new Algorithm[]{Algorithm.BOSS, Algorithm.BFCI};

                    IndependenceWrapper test;
                    ScoreWrapper score;

                    Parameters params = new Parameters();
                    params.set(Params.ALPHA, 0.001);
                    params.set(Params.PENALTY_DISCOUNT, 2);
                    params.set(Params.NUM_STARTS, 10);
                    params.set(Params.COMPLETE_RULE_SET_USED, true);
                    params.set(Params.BOSS_SCORE_TYPE, false);

                    for (Algorithm algorithm : algorithms) {
                        File outfile = new File(dataDir, "results1/"
                                + algorithm.toString().toLowerCase(Locale.ROOT) + "." + name);

                        if (outfile.exists()) continue;

                        DataSet _data = null;

                        if (name.contains("continuous")) {

                            try {
                                Data data = new ContinuousTabularDatasetFileReader(f.toPath(), Delimiter.TAB).readInData();
                                _data = (DataSet) DataConvertUtils.toDataModel(data);

                                if (_data == null) {
                                    System.out.println("Couldn't find that dataset: " + f);
                                } else {
                                    _data = DataUtils.shuffleColumns(_data);
                                }

                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            assert _data != null;
                            test = new FisherZ();
                            score = new edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore();
                        } else {

                            String[] tokens = name.split("\\.");

                            int numCategories = 0;

                            for (int i = 0; i < tokens.length; i++) {
                                if ("mixed".equals(tokens[i])) {
                                    numCategories = Integer.parseInt(tokens[i + 1]);
                                }
                            }

                            try {
                                Data data = new MixedTabularDatasetFileReader(f.toPath(), Delimiter.TAB,
                                        numCategories).readInData();
                                _data = (DataSet) DataConvertUtils.toDataModel(data);

                                if (_data == null) {
                                    System.out.println("Couldn't find that dataset: " + f);
                                } else {
                                    _data = DataUtils.shuffleColumns(_data);
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            assert _data != null;
                            test = new ConditionalGaussianLRT();
                            score = new ConditionalGaussianBicScore();
                        }

                        edu.cmu.tetrad.algcomparison.algorithm.Algorithm search = null;

                        if (algorithm == Algorithm.PC) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc(test);
                        } else if (algorithm == Algorithm.CPC) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Cpc(test);
                        } else if (algorithm == Algorithm.PCMAX) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.PcMax(test);
                        } else if (algorithm == Algorithm.FGES) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Fges(score);
                        } else if (algorithm == Algorithm.BOSS) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.BOSS(score, test);
                        } else if (algorithm == Algorithm.FCI) {
                            search = new Fci(test);
                        } else if (algorithm == Algorithm.FCIMAX) {
                            search = new FciMax(test);
                        } else if (algorithm == Algorithm.RFCI) {
                            search = new Rfci(test);
                        } else if (algorithm == Algorithm.GFCI) {
                            search = new Gfci(score, test);
                        } else if (algorithm == Algorithm.BFCI) {
                            search = new BFCI(score, test);
                        } else if (algorithm == Algorithm.MGM) {
                            search = new Mgm();
                        }

                        assert search != null;
                        Graph graph = search.search(_data, params, null);

                        try {
                            if (!outfile.getParentFile().mkdirs()) {
                                System.out.println("Couldn't make directories");
                            }
                            PrintWriter writer = new PrintWriter(outfile);
                            writer.println(graph);
                            writer.close();
                        } catch (FileNotFoundException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
        }
    }

    @Test
    public void test2() {
        File dir = new File("/Users/josephramsey/Downloads/GitHub/tetrad-example-analyses");

        File[] files = dir.listFiles();

        for (File dataDir : files) {
            if (dataDir.getName().startsWith("README")) continue;
            if (dataDir.getName().startsWith("LICENSE")) continue;
            if (dataDir.getName().startsWith(".")) continue;

            File dataDir2 = new File(dataDir, "tetrad.formatted.data");

            File[] obj = dataDir2.listFiles();

            if (obj == null) {
                System.out.println("Couldn't fine the directory for " + dataDir.getName());
            }

            if (obj != null) {
                for (File f : obj) {
                    if (f.getName().startsWith(".")) continue;

                    String name = f.getName();

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.PC, Algorithm.CPC, Algorithm.PCMAX};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FGES, Algorithm.BOSS};

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FCI, Algorithm.FCIMAX, Algorithm.RFCI};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.GFCI, Algorithm.BFCI};

                    Parameters params = new Parameters();
                    params.set(Params.ALPHA, 0.001);
                    params.set(Params.PENALTY_DISCOUNT, 2);
                    params.set(Params.NUM_STARTS, 5);
                    params.set(Params.COMPLETE_RULE_SET_USED, true);

                    for (Algorithm algorithm : Algorithm.values()) {
                        File outfile = new File(dataDir, "results1/"
                                + algorithm.toString().toLowerCase(Locale.ROOT) + "." + name);
                        File knowledgeFile = new File(dataDir, "ground.truth.knowledge/knowledge." + name);

                        if (!knowledgeFile.exists()) {
                            continue;
                        }

                        IKnowledge knowledge = null;

                        try {
                            knowledge = DataUtils.parseKnowledge(knowledgeFile, DelimiterType.TAB, "//");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        if (knowledge == null) {
                            continue;
                        }

//                        System.out.println(knowledge);

                        Graph graph = GraphUtils.loadGraphTxt(outfile);
                        boolean satisfies = true;

                        for (Edge edge : graph.getEdges()) {
                            Node x = edge.getNode1();
                            Node y = edge.getNode2();

                            String xn = x.getName();
                            String yn = y.getName();

                            if (Objects.equals(graph.getEdge(x, y), Edges.directedEdge(x, y)) && knowledge.isForbidden(xn, yn)) {
//                                System.out.println("Violates: " + edge + " " + outfile.getName());
                                satisfies = false;
                            } else if (Objects.equals(graph.getEdge(y, x), Edges.directedEdge(y, x)) && knowledge.isForbidden(yn, xn)) {
//                                System.out.println("Violates: " + edge + " " + outfile.getName());
                                satisfies = false;
                            }
//                            else if (graph.isDirectedFromTo(y, x) && knowledge.isRequired(xn, yn)) {
//                                satisfies = false;
//                            } else if (graph.isDirectedFromTo(x, y) && knowledge.isRequired(yn, xn)) {
//                                satisfies = false;
//                            }
                        }

                        if (satisfies) {
                            System.out.println("Satisfies ground truth: " + outfile.getName() + "(" + graph.getNumEdges() + " edge)");
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testMeasles() {
        File dir = new File("/Users/josephramsey/Downloads/GitHub/tetrad-example-analyses/hungarian-measles");
        File lagFile = new File(dir, "results.time.lag.runs/boss.lag2.txt");
        Graph lag = GraphUtils.loadGraphTxt(lagFile);

        // Find nodes for adj.
        Graph adj = new EdgeListGraph();

        for (Edge edge : lag.getEdges()) {
            Node fromNode = edge.getNode1();
            Node toNode = edge.getNode2();

            String from = fromNode.getName();
            String to = toNode.getName();

            String[] tokens1 = from.split(":");
            String[] tokens2 = to.split(":");

            if (adj.getNode(tokens1[0]) == null) {
                adj.addNode(new ContinuousVariable(tokens1[0]));
            }

            if (adj.getNode(tokens2[0]) == null) {
                adj.addNode(new ContinuousVariable(tokens2[0]));
            }
        }

        // Read in adjacencies and put them in adj as undirected edges.
        File adjacentFile = new File(dir, "data/hungary_county_edges.txt");

        System.out.println("FROM\tLAGFROM\tTO\tLAG2,\tDISTANCE");

        try {
            BufferedReader in = new BufferedReader(new FileReader(adjacentFile));

            in.readLine();

            String line;

            while ((line = in.readLine()) != null) {
                String[] tokens = line.split("\t");

                String name1 = tokens[0];
                String name2 = tokens[1];

                Node node1 = adj.getNode(name1);
                Node node2 = adj.getNode(name2);

                if (node1 == node2) continue;

                adj.addUndirectedEdge(node1, node2);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        int sum0 = 0;
        int sum1 = 0;
        int sum2 = 0;

        int count0 = 0;
        int count1 = 0;
        int count2 = 0;

        int max0 = 0;
        int max1 = 0;
        int max2 = 0;

        for (Edge edge : lag.getEdges()) {
            if (!edge.isDirected()) {
                continue;
            }

            Node fromNode = Edges.getDirectedEdgeTail(edge);
            Node toNode = Edges.getDirectedEdgeHead(edge);

            String from = fromNode.getName();
            String to = toNode.getName();

            String[] tokens1 = from.split(":");
            String[] tokens2 = to.split(":");

            if (tokens1[0].equals("BUDAPEST")) {
                tokens1[0] = "PEST";
            }

            if (tokens2[0].equals("BUDAPEST")) {
                tokens2[0] = "PEST";
            }

            String index1 = tokens1.length > 1 ? tokens1[tokens1.length - 1] : "0";
            String index2 = tokens2.length > 1 ? tokens2[tokens2.length - 1] : "0";

//            if (Integer.parseInt(index1) < Integer.parseInt(index2)) {
//                continue;
//            }

            if (Integer.parseInt(index2) != 0) {
                continue;
            }

            Node adj1 = adj.getNode(tokens1[0]);
            Node adj2 = adj.getNode(tokens2[0]);

            int distance = distance(adj1, adj2, adj);

            System.out.println(tokens1[0]
                    + "\t" + index1
                    + "\t" + tokens2[0]
                    + "\t" + index2
                    + "\t" + distance);

            int _lag = Integer.parseInt(index1) - Integer.parseInt(index2);

            switch (_lag) {
                case 0:
                    sum0 += distance;
                    count0++;
                    if (distance > max0) max0 = distance;
                    break;
                case 1:
                    sum1 += distance;
                    count1++;
                    if (distance > max1) max1 = distance;
                    break;
                case 2:
                    sum2 += distance;
                    count2++;
                    if (distance > max2) max2 = distance;
                    break;
                default:
                    throw new IllegalStateException();
            }
        }

        List<Node> adjNodes = adj.getNodes();
        double sum = 0.0;
        int count = 0;
        int max = 0;

        for (int i = 0; i < adjNodes.size(); i++) {
            for (int j = i + 1; j < adjNodes.size(); j++) {
                if (adjNodes.get(i).getName().equals("BUDAPEST")
                        || adjNodes.get(j).getName().equals("BUDAPEST")) {
                    continue;
                }

                int d = distance(adjNodes.get(i), adjNodes.get(j), adj);
                sum += d;
                count += 1;

                if (d > max) {
                    max = d;
                }
            }
        }

        double avg = sum / count;

        System.out.println();

        System.out.println("Average distance between counties = " + avg);
        System.out.println("Maximum distance between counties = " + max);

        System.out.println();

        System.out.println("Average distance lag 0 = " + (sum0 / (double) count0));
        System.out.println("Average distance lag 1 = " + (sum1 / (double) count1));
        System.out.println("Average distance lag 2 = " + (sum2 / (double) count2));

        System.out.println();

        System.out.println("Max distance lag 0 = " + (max0));
        System.out.println("Max distance lag 1 = " + (max1));
        System.out.println("Max distance lag 2 = " + (max2));

        System.out.println();

        System.out.println("# lag 0 = " + (count0));
        System.out.println("# lag 1 = " + (count1));
        System.out.println("# lag 2 = " + (count2));
    }

    private int distance(Node adj1, Node adj2, Graph adj) {
        int distance = 0;
        Set<Node> all = new HashSet<>();
        all.add(adj1);

        while (!all.contains(adj2)) {
            all.addAll(adjacent(adj, all));
            distance++;
        }

        return distance;
    }

    private Set<Node> adjacent(Graph adj, Node adj1) {
        return new HashSet<>(adj.getAdjacentNodes(adj1));
    }

    private Set<Node> adjacent(Graph adj, Set<Node> adj1) {
        Set<Node> all = new HashSet<>();

        for (Node _adj : adj1) {
            all.addAll(adj.getAdjacentNodes(_adj));
        }

        return all;
    }

    private enum Algorithm {PC, CPC, PCMAX, FGES, BOSS, FCI, FCIMAX, RFCI, GFCI, BFCI, MGM}
}

