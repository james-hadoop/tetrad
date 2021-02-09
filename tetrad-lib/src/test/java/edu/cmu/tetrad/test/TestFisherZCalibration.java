package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.graph.GraphUtils.existsSemiDirectedPath;
import static edu.cmu.tetrad.graph.GraphUtils.traverseSemiDirected;
import static java.lang.StrictMath.abs;
import static org.junit.Assert.assertEquals;

public class TestFisherZCalibration {

//    public static void main(String... args) {
//        test1();
//    }

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(105034020L);
        doTest();
    }

    @Test
    public void doTest() {
        Parameters parameters = new Parameters();
        parameters.set(Params.ALPHA, 0.01);
        parameters.set(Params.DEPTH, 2);
        parameters.set(Params.STRUCTURE_PRIOR, 0);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.SAMPLE_SIZE, 10000);
        parameters.set(Params.SEM_BIC_RULE, 1);
        parameters.set(Params.PENALTY_DISCOUNT, 1);

        parameters.set(Params.NUM_MEASURES, 20);
        parameters.set(Params.AVG_DEGREE, 4);
        int numDraws = 1000;

        int maxNumEdges = parameters.getInt(Params.NUM_MEASURES) * parameters.getInt(Params.AVG_DEGREE) / 2;

        Graph graph = GraphUtils.randomDag(parameters.getInt(Params.NUM_MEASURES),
                0, maxNumEdges, 100,
                100, 100, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet data = im.simulateData(parameters.getInt(Params.SAMPLE_SIZE), false);

        IndependenceTest test1 = new FisherZ().getTest(new CovarianceMatrix(data), parameters);
        IndependenceTest test2 = new SemBicTest().getTest(data, parameters);

        List<Node> variables = data.getVariables();
        graph = GraphUtils.replaceNodes(graph, variables);

        IndependenceTest dsep = new IndTestDSep(graph);

        for (int depth : new int[]{0, 1, 2, 3, 4, 5, 7, 8, 9, 10}) {
            testOneDepth(parameters, numDraws, test1, test2, variables, dsep, depth);
        }
    }

    private void testOneDepth(Parameters parameters, int numDraws, IndependenceTest test1, IndependenceTest test2, List<Node> variables, IndependenceTest dsep, int depth) {
        int countSame = 0;
        int fn1 = 0;
        int fn2 = 0;
        int fp1 = 0;
        int fp2 = 0;
        int ds = 0;
        int dc = 0;

        for (int i = 0; i < numDraws; i++) {
            Collections.shuffle(variables);
//            Collections.shuffle(variables);
//            Collections.shuffle(variables);

            Node x = variables.get(0);
            Node y = variables.get(1);

            List<Node> z = new ArrayList<>();
            for (int j = 0; j < depth; j++) {
                z.add(variables.get(j + 2));
            }

            boolean fzInd = test1.isIndependent(x, y, z);
            boolean sembInd = test2.isIndependent(x, y, z);
            boolean _dsep = dsep.isIndependent(x, y, z);

            if (fzInd == sembInd) countSame++;

            if (fzInd && !_dsep) fn1++;
            if (!fzInd && _dsep) fp1++;
            if (sembInd && !_dsep) fn2++;
            if (!sembInd && _dsep) fp2++;
            if (_dsep) ds++;
            if (!_dsep) dc++;
        }

        TextTable table = new TextTable(3, 3);
        table.setToken(0, 1, "FP");
        table.setToken(0, 2, "FN");
        table.setToken(1, 0, "Fisher Z");
        table.setToken(2, 0, "Local Consistency Criterion");

        table.setToken(1, 1, "" + fp1);
        table.setToken(1, 2, "" + fn1);
        table.setToken(2, 1, "" + fp2);
        table.setToken(2, 2, "" + fn2);

        System.out.println();
        System.out.println("Depth = " + depth);
        System.out.println();
        System.out.println("Same = " + countSame + " out of " + numDraws);
        System.out.println("# dsep = " + ds);
        System.out.println("# dconn = " + dc);
        System.out.println();
        System.out.println(table);

        System.out.println();

        double alpha = parameters.getDouble(Params.ALPHA);
        System.out.println("alpha = " + alpha);
        double alphaHat = fp1 / (double) ds;
        double betaHat = fn1 / (double) dc;
        double alphaHat2 = fp2 / (double) ds;
        double betaHat2 = fn2 / (double) dc;
        System.out.println("alpha^ for Fisher Z = " + alphaHat);
        System.out.println("beta^ for Fisher Z = " + betaHat);
        System.out.println("alpha^ for Local Scoring Consistency Criterion = " + alphaHat2);
        System.out.println("beta^ for Local Scoring Consistency Criterion = " + betaHat2);

//        Assert.assertTrue(abs(alpha - alphaHat) < 2 * alpha);
    }

    //    @Test
    public void test2() {
        Parameters parameters = new Parameters();
        parameters.set(Params.PENALTY_DISCOUNT, 1);
        parameters.set(Params.STRUCTURE_PRIOR, 0);
        parameters.set(Params.COEF_LOW, .2);
        parameters.set(Params.COEF_HIGH, 1);
        int sampleSize = 20000;

        Graph graph = GraphUtils.randomDag(100, 0, 300, 100,
                100, 100, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet data = im.simulateData(sampleSize, false);

        Fges fges = new Fges(new SemBicScore(data));
        fges.setVerbose(true);

        fges.search();

//        fges = new Fges(new GraphScore(graph));
//        fges.setVerbose(true);
//        fges.search();


    }

    @Test
    public void test3() {

//        RandomUtil.getInstance().setSeed(92883342449L);

        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 5);
        parameters.set(Params.NUM_MEASURES, 20);
        parameters.set(Params.AVG_DEGREE, 4);
        parameters.set(Params.SAMPLE_SIZE, 100000);
//        parameters.set(Params.SAMPLE_SIZE, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000);
        parameters.set(Params.SEM_BIC_RULE, 1);
        parameters.set(Params.SEM_BIC_STRUCTURE_PRIOR, 0);
        parameters.set(Params.PENALTY_DISCOUNT, 1);
//        parameters.set(Params.TDEPTH, 1);
        parameters.set(Params.FAITHFULNESS_ASSUMED, true);
        parameters.set(Params.SYMMETRIC_FIRST_STEP, false);
        parameters.set(Params.ADJUST_ORIENTATIONS, false);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1.0);
        parameters.set(Params.VAR_HIGH, 3.0);
        parameters.set(Params.VERBOSE, false);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);
        parameters.set(Params.USE_EQUIVALENT_SAMPLE_SIZE, false);


        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.NUM_RUNS));
        statistics.add(new ParameterColumn(Params.SEM_BIC_RULE));
//        statistics.add(new ParameterColumn(Params.NUM_MEASURES));
//        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.PENALTY_DISCOUNT));
        statistics.add(new ParameterColumn(Params.SEM_BIC_STRUCTURE_PRIOR));

        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
//        statistics.add(new ArrowheadPrecision());
//        statistics.add(new ArrowheadRecall());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
//        statistics.add(new AhpBound());

        statistics.add(new F1Adj());
        statistics.add(new F1Arrow());
        statistics.add(new SHD());

        statistics.add(new CyclicEst());

        statistics.add(new ElapsedTime());

        statistics.setWeight("SHD", .5);

        Algorithms algorithms = new Algorithms();

        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.SemBicScore()));
//
        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setSaveGraphs(false);
        comparison.setSaveData(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);

        comparison.compareFromSimulations(
                "/Users/josephramsey/tetrad/comparison2040", simulations, algorithms, statistics, parameters);
    }

    @Test
    public void test4() {
        Parameters parameters = new Parameters();

        parameters.set("numRuns", 1);
        parameters.set("differentGraphs", true);
        parameters.set("sampleSize", 1000);

        parameters.set("numMeasures", 50);
        parameters.set("numLatents", 0);
        parameters.set("avgDegree", 6);
        parameters.set("maxDegree", 100);
        parameters.set("maxIndegree", 100);
        parameters.set("maxOutdegree", 100);
        parameters.set("connected", false);
//        parameters.set("depth", 1);
        parameters.set(Params.FAITHFULNESS_ASSUMED, false);

        parameters.set("coefLow", 0.2);
        parameters.set("coefHigh", 1.0);
        parameters.set("varLow", 1);
        parameters.set("varHigh", 3);
        parameters.set("verbose", false);
        parameters.set("coefSymmetric", true);
        parameters.set("randomizeColumns", true);

        parameters.set("penaltyDiscount", 1);
        parameters.set(Params.STRUCTURE_PRIOR, 0);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn("thresholdAlpha"));
        statistics.add(new ParameterColumn("penaltyDiscount"));

        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.SemBicScore()
        ));

        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(false);
        comparison.setShowSimulationIndices(false);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveGraphs(true);

        comparison.compareFromSimulations("comparisonWayne", simulations, algorithms, statistics, parameters);
    }

    @Test
    public void test5() {
        Graph graph = GraphUtils.randomGraph(7, 0, 10, 100, 100, 100, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet data = im.simulateData(1000, true);
        Score score = new SemBicScore(data);

        Fges fges = new Fges(score);

        Graph pattern = fges.search();

//        MeekRules rules = new MeekRules();
//        rules.setUndirectUnforcedEdges(true);

//        System.out.println("First");
//        rules.orientImplied(pattern);
//        rules.orientImplied(pattern);
//
//        System.out.println("Second");
//        rules.orientImplied(pattern);

        if (pattern == null) throw new IllegalArgumentException();

        for (Node node : pattern.getNodes()) {
            if (undirectUnforcedEdges(node, pattern)) {
                System.out.println("Could undirect " + node);
            }
        }
    }

    private boolean undirectUnforcedEdges(Node node, Graph graph) {
        System.out.println("Undirecting unforced edges for " + node);

        Set<Node> parentsToUndirect = new HashSet<>();
        List<Node> parents = graph.getParents(node);

        NEXT_EDGE:
        for (Node x : parents) {
            for (Node parent : parents) {
                if (parent != x) {
                    if (!graph.isAdjacentTo(parent, x)) {
                        continue NEXT_EDGE;
                    }
                }
            }

            parentsToUndirect.add(x);
        }

        System.out.println("   Parents to undirect for " + node + " = " + parentsToUndirect);

        return (!parentsToUndirect.isEmpty());
    }

    @Test
    public void test6() {
        List<Node> nodes = new ArrayList<>();
        Node x = new GraphNode("x");
        Node y = new GraphNode("y");
        Node z = new GraphNode("z");

        nodes.add(x);
        nodes.add(y);
        nodes.add(z);

        Graph graph = new EdgeListGraph(nodes);

        graph.addUndirectedEdge(x, y);
        graph.addDirectedEdge(x, z);
        graph.addDirectedEdge(y, z);

        MeekRules rules = new MeekRules();

        rules.revertToUnshieldedColliders(nodes, graph);

        System.out.println(graph);
    }

    @Test
    public void test7() {
        NumberFormat nf = new DecimalFormat("00");

        for (int j = 3; j <= 8; j++) {
            System.out.println("\nj = " + j);


            double apsum = 0.0;
            double arsum = 0.0;
            double ahpsum = 0.0;
            double ahrsum = 0.0;
            double f1adjsum = 0.0;
            double f1arrowsum = 0.0;

            for (int i = 0; i < 50; i++) {
                File dataFile = new File("/Users/josephramsey/Downloads/Triplet_A_Star_Data"
                        + "/n60_prob0." + nf.format(j) + "_N500_pow1"
                        + "/raw_data_13" + nf.format(i) + ".csv");
                File graphFile = new File("/Users/josephramsey/Downloads/Triplet_A_Star_Data"
                        + "/n60_prob0." + nf.format(j) + "_N500_pow1"
                        + "/true_dag_13" + nf.format(i) + ".txt");

//                System.out.println(dataFile);

                DataSet rawdata = null;

                try {
                    ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(dataFile.toPath(), Delimiter.COMMA);
                    reader.setHasHeader(false);
                    Data data = reader.readInData();
                    rawdata = (DataSet) DataConvertUtils.toDataModel(data);
                } catch (IOException ioException) {
                    ioException.printStackTrace();
                }

//                System.out.println(rawdata);

                SemBicScore score = new SemBicScore(rawdata);
                score.setUseEquivalentSampleSize(false);
                score.setRuleType(SemBicScore.RuleType.CHICKERING);
                score.setPenaltyDiscount(1);

                Fges fges = new Fges(score);
                fges.setFaithfulnessAssumed(false);
                fges.setTurning(false);
//                fges.setVerbose(true);

                Graph out = fges.search();

//                System.out.println(out);


//                Graph graph = GraphUtils.loadRSpecial(graphFile);
                Graph graph = GraphUtils.loadGraphTxt(graphFile);


                List<Node> nodes = graph.getNodes();

                for (int k = 1; k <= nodes.size(); k++) {
                    String name = "X" + k;
                    Node node = graph.getNode(name);

                    if (node == null) {
                        System.out.println(name);
                    }

                    node.setName("C" + k);
                }

                graph = GraphUtils.replaceNodes(graph, out.getNodes());


//                graph = SearchGraphUtils.patternForDag(graph);

//                System.out.println(graph);

//                System.out.println("AP = " + new AdjacencyPrecision().getValue(graph, out, rawdata));
//                System.out.println("AR = " + new AdjacencyRecall().getValue(graph, out, rawdata));
//                System.out.println("AHP = " + new ArrowheadPrecision().getValue(graph, out, rawdata));
//                System.out.println("AHR = " + new AdjacencyRecall().getValue(graph, out, rawdata));
//                System.out.println("F1Adj = " + new F1Adj().getValue(graph, out, rawdata));
//                System.out.println("F1Arrow = " + new F1Arrow().getValue(graph, out, rawdata));
//                System.out.println();

                apsum += new AdjacencyPrecision().getValue(graph, out, rawdata);
                arsum += new AdjacencyRecall().getValue(graph, out, rawdata);

                ahpsum += new ArrowheadPrecision().getValue(graph, out, rawdata);
                ahrsum += new ArrowheadRecall().getValue(graph, out, rawdata);

                f1adjsum += new F1Adj().getValue(graph, out, rawdata);
                f1arrowsum += new F1Arrow().getValue(graph, out, rawdata);
            }

            System.out.println("AP = " + apsum / 50.);
            System.out.println("AR = " + arsum / 50.);
            System.out.println("AHP = " + ahpsum / 50.);
            System.out.println("AHR = " + ahrsum / 50.);
            System.out.println("F1Adj = " + f1adjsum / 50.);
            System.out.println("F1Arrow = " + f1arrowsum / 50.);
        }
    }

    @Test
    public void test8() {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < 20; i++) {
            nodes.add(new GraphNode("X" + i));
        }

        EdgeListGraph graph = new EdgeListGraph(nodes);

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                graph.addDirectedEdge(nodes.get(i), nodes.get(j));
                graph.addDirectedEdge(nodes.get(j), nodes.get(i));
            }
        }

        List<Edge> edges = new ArrayList<>(graph.getEdges());

        EdgeListGraph out = new EdgeListGraph(nodes);

        while (!edges.isEmpty()) {
            Edge edge = edges.get(RandomUtil.getInstance().nextInt(edges.size()));
            edges.remove(edge);

            Node x = edge.getNode1();
            Node y = edge.getNode2();

            if (out.isAdjacentTo(x, y)) continue;

            List<Node> tNeighbors = getTNeighbors(x, y, out);

            List<Node> T = new ArrayList<>();

            for (Node tNeighbor : tNeighbors) {
                if (RandomUtil.getInstance().nextInt(2) == 1) {
                    T.add(tNeighbor);
                }
            }

            Set<Node> union = new HashSet<>(T);
            union.addAll(getNaYX(x, y, out));

            if (!isClique(union, out)) {
                continue;
            }

            if (!semidirectedPathCondition(y, x, union, out)) {
                continue;
            }

            out.addEdge(edge);

            for (Node t : T) {
                if (out.getEdge(t, y).isDirected()) {
                    throw new IllegalArgumentException();
                }

                out.removeEdge(t, y);
                out.addDirectedEdge(t, y);
            }

            boolean cycle3 = out.existsDirectedCycle();

            MeekRules rules = new MeekRules();
            rules.orientImplied(out);

            boolean cycle4 = out.existsDirectedCycle();

            if (!cycle3 && cycle4) {
                System.out.println("cycle4 cycle");

                for (Node n : nodes) {
                    List<Node> path = out.findDirectedPath(n, n);
                    if (!path.isEmpty()) {
                        System.out.println(out.getNumEdges() + "." + GraphUtils.pathString(path, out));
                    }
                }

            }
        }
    }

    // Find all adj that are connected to Y by an undirected edge that are adjacent to X (that is, by undirected or
    // directed edge).
    private Set<Node> getNaYX(Node x, Node y, Graph graph) {
        List<Node> adj = graph.getAdjacentNodes(y);
        Set<Node> nayx = new HashSet<>();

        for (Node z : adj) {
            if (z == x) {
                continue;
            }
            Edge yz = graph.getEdge(y, z);
            if (!Edges.isUndirectedEdge(yz)) {
                continue;
            }
            if (!graph.isAdjacentTo(z, x)) {
                continue;
            }
            nayx.add(z);
        }

        return nayx;
    }

    // Returns true iif the given set forms a clique in the given graph.
    private boolean isClique(Set<Node> nodes, Graph graph) {
        List<Node> _nodes = new ArrayList<>(nodes);
        for (int i = 0; i < _nodes.size(); i++) {
            for (int j = i + 1; j < _nodes.size(); j++) {
                if (!graph.isAdjacentTo(_nodes.get(i), _nodes.get(j))) {
                    return false;
                }
            }
        }

        return true;
    }

    private boolean semidirectedPathCondition(Node from, Node to, Set<Node> cond, Graph graph) {
        if (from == to) throw new IllegalArgumentException();

        Queue<Node> Q = new LinkedList<>();
        Set<Node> V = new HashSet<>();

        for (Node u : graph.getAdjacentNodes(from)) {
            Edge edge = graph.getEdge(from, u);
            Node c = traverseSemiDirected(from, edge);

            if (c == null) {
                continue;
            }

            if (!V.contains(c)) {
                V.add(c);
                Q.offer(c);
            }
        }

        while (!Q.isEmpty()) {
            Node t = Q.remove();

            if (cond.contains(t)) {
                continue;
            }

            if (t == to) {
                return false;
            }

            for (Node u : graph.getAdjacentNodes(t)) {
                Edge edge = graph.getEdge(t, u);
                Node c = traverseSemiDirected(t, edge);

                if (c == null) {
                    continue;
                }

                if (!V.contains(c)) {
                    V.add(c);
                    Q.offer(c);
                }
            }
        }

        return true;
    }

    // Get all adj that are connected to Y by an undirected edge and not adjacent to X.
    private List<Node> getTNeighbors(Node x, Node y, Graph graph) {
        List<Edge> yEdges = graph.getEdges(y);
        List<Node> tNeighbors = new ArrayList<>();

        for (Edge edge : yEdges) {
            if (!Edges.isUndirectedEdge(edge)) {
                continue;
            }

            Node z = edge.getDistalNode(y);

            if (graph.isAdjacentTo(z, x)) {
                continue;
            }

            tNeighbors.add(z);
        }

        return tNeighbors;
    }

    @Test
    public void test9() {
        Graph dag = GraphUtils.randomDag(20, 0, 40, 100,
                100, 100, false);

        MeekRules rules = new MeekRules();

        Graph pattern = new EdgeListGraph(dag);

        rules.orientImplied(pattern);

        Graph pattern2 = new EdgeListGraph(pattern);

        rules.orientImplied(pattern2);

        if (!pattern.equals(pattern2)) {
            System.out.println("Not equal");
        }

        System.out.println(pattern);

        Graph dag2 = SearchGraphUtils.dagFromPattern(pattern2);

        rules.orientImplied(dag2);

        if (!pattern.equals(dag2)) {
            System.out.println("Not equal");
        }




    }

    @Test
    public void test10() {
        int numNodes = 8;
        int aveDegree = 4;
        int numIterations = 1;

        for (int i = 0; i < numIterations; i++) {
            Graph dag = GraphUtils.randomDag(numNodes, 0, aveDegree * numNodes / 2,
                    10, 10, 10, false);
            GraphScore score = new GraphScore(dag);
//            score.setVerbose(true);
            Fges fges = new Fges(score);
            fges.setFaithfulnessAssumed(false);
            fges.setSymmetricFirstStep(true);
            fges.setVerbose(true);
//            fges.setTDepth(1);
            Graph pattern1 = fges.search();
            Graph pattern2 = new Pc(new IndTestDSep(dag)).search();
            assertEquals(pattern2, pattern1);
        }
    }
}
