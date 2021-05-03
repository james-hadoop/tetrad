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
import edu.cmu.tetrad.sem.StandardizedSemIm;
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

import static edu.cmu.tetrad.graph.GraphUtils.loadRSpecial;
import static edu.cmu.tetrad.graph.GraphUtils.traverseSemiDirected;
import static java.lang.Math.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestFisherZCalibration {

    public static Graph loadMyGraph(File file, boolean dir) {
        DataSet eg = null;

        try {
            ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(file.toPath(), Delimiter.COMMA);
            reader.setHasHeader(false);
            Data data = reader.readInData();
            eg = (DataSet) DataConvertUtils.toDataModel(data);
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }

        if (eg == null) throw new NullPointerException();

        List<Node> vars = eg.getVariables();

        Graph graph = new EdgeListGraph(vars);

        for (int i = 0; i < vars.size(); i++) {
            for (int j = 0; j < vars.size(); j++) {
                if (i == j) continue;
                if (eg.getDouble(i, j) == 1 && eg.getDouble(j, i) == 1) {
                    if (!graph.isAdjacentTo(vars.get(i), vars.get(j))) {
                        graph.addUndirectedEdge(vars.get(i), vars.get(j));
                    }
                } else if (dir && eg.getDouble(i, j) == 0 && eg.getDouble(j, i) == 1) {
                    graph.addDirectedEdge(vars.get(i), vars.get(j));
                } else if (!dir && eg.getDouble(i, j) == 1 && eg.getDouble(j, i) == 0) {
                    graph.addDirectedEdge(vars.get(i), vars.get(j));
                }
            }
        }

        return graph;
    }

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
    }

    @Test
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
    }

    @Test
    public void test3() {

//        RandomUtil.getInstance().setSeed(92883342449L);

        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 10);

        parameters.set(Params.AVG_DEGREE, 4);

        int N = 1000;
        parameters.set(Params.SAMPLE_SIZE, N);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1);
        parameters.set(Params.VAR_HIGH, 2);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);
        parameters.set(Params.SEM_IM_SIMULATION_TYPE, false);
        parameters.set(Params.CORRELATION_THRESHOLD, 1);

        parameters.set(Params.SYMMETRIC_FIRST_STEP, false);
        parameters.set(Params.VERBOSE, false);
        parameters.set(Params.FAITHFULNESS_ASSUMED, true);
        parameters.set(Params.PARALLELISM, 20);

        parameters.set(Params.SEM_BIC_STRUCTURE_PRIOR, 0);

        parameters.set(Params.INTERVAL_BETWEEN_SHOCKS, 50);
        parameters.set(Params.INTERVAL_BETWEEN_RECORDINGS, 50);
        parameters.set(Params.SELF_LOOP_COEF, 0);
        parameters.set(Params.FISHER_EPSILON, 0.001);

        parameters.set(Params.PENALTY_DISCOUNT, 2);
        parameters.set(Params.TAKE_LOGS, true);
        parameters.set(Params.TRUE_ERROR_VARIANCE, 1);
        parameters.set(Params.DISCRETIZE, false);

        parameters.set(Params.STABLE_FAS, true);
        parameters.set(Params.CONCURRENT_FAS, true);
        parameters.set(Params.COLLIDER_DISCOVERY_RULE, 3);
        parameters.set(Params.CONFLICT_RULE, 3);
        parameters.set(Params.DEPTH, -1);
        parameters.set(Params.USE_MAX_P_ORIENTATION_HEURISTIC, true);
        parameters.set(Params.MAX_P_ORIENTATION_MAX_PATH_LENGTH, -1);

        parameters.set(Params.SEM_GIC_RULE, 1, 2, 3, 4, 5, 6);
//        parameters.set(Params.ZS_RISK_BOUND, 5. / vars);

//        double gamma = 1- 1.0 / (2 * (log(vars) / log(N)));

        parameters.set(Params.EBIC_GAMMA, 1);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.NUM_RUNS));
        statistics.add(new ParameterColumn(Params.NUM_MEASURES));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));

        statistics.add(new ParameterColumn(Params.SEM_GIC_RULE));
        statistics.add(new ParameterColumn(Params.ZS_RISK_BOUND));
        statistics.add(new ParameterColumn(Params.EBIC_GAMMA));
//
        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new AdjacencyFDR());
        statistics.add(new ArrowheadFDR());
        statistics.add(new AhpcCriterion());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                new edu.cmu.tetrad.algcomparison.score.SemBicScore()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.KimEtAlScores()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.ZhangShenBoundScore()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.EbicScore()));

        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setSaveGraphs(false);
        comparison.setSaveData(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        int[] vars = {20};//, 100, 200};
//        int[] vars = {1000, 2000};

        for (int _vars : vars) {
            parameters.set(Params.NUM_MEASURES, _vars);
//            parameters.set(Params.SEM_GIC_RULE, /*1, 2,*/ /*3,*/ 4, /*5,*/ 6);
            comparison.compareFromSimulations(
                    "/Users/josephramsey/tetrad/comparison-" + _vars, simulations, algorithms, statistics, parameters);
        }
    }

    @Test
    public void test3a() {

        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 1);

        parameters.set(Params.NUM_MEASURES, 600);
        parameters.set(Params.AVG_DEGREE, 2.8);

        parameters.set(Params.SAMPLE_SIZE, 200);
        parameters.set(Params.COEF_LOW, 0.1);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1);
        parameters.set(Params.VAR_HIGH, 2);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);
        parameters.set(Params.SEM_IM_SIMULATION_TYPE, false);

        parameters.set(Params.SYMMETRIC_FIRST_STEP, false);
        parameters.set(Params.VERBOSE, true);
        parameters.set(Params.FAITHFULNESS_ASSUMED, true);
        parameters.set(Params.PARALLELISM, 20);

        parameters.set(Params.SEM_GIC_RULE, 1, 2, 3, 4, 5, 6);
//        parameters.set(Params.PENALTY_DISCOUNT, .05);
        parameters.set(Params.TAKE_LOGS, true);
        parameters.set(Params.TRUE_ERROR_VARIANCE, 1);
        parameters.set(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED, false);
        parameters.set(Params.ZS_RISK_BOUND, 0, .1, .2, .3, .4, .5, .6, .7);


        parameters.set(Params.USE_MAX_P_ORIENTATION_HEURISTIC, true);
        parameters.set(Params.MAX_P_ORIENTATION_MAX_PATH_LENGTH, -1);

//        parameters.set(Params.ZS_RISK_BOUND, 0, 0.001, 0.01, 0.199, 0.3);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.SEM_GIC_RULE));
        statistics.add(new ParameterColumn(Params.ZS_RISK_BOUND));
        statistics.add(new NumberOfNodesTrue());
        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new AdjacencyTPR());
        statistics.add(new AdjacencyFPR());
        statistics.add(new AdjacencyTPRDirectedParts());
        statistics.add(new AdjacencyFPRDirectedParts());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                new edu.cmu.tetrad.algcomparison.score.KimEtAlScores()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.ZhangShenBoundScore()));
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                new edu.cmu.tetrad.algcomparison.score.EBicScore()));
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                new edu.cmu.tetrad.algcomparison.score.SemBicScore()));

//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll(new KimEtAlTests(), null));
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll(new SemBicTest(), null));

        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));
//        simulations.add(new LinearFisherModel(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setSaveGraphs(false);
        comparison.setSaveData(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        comparison.compareFromSimulations(
                "/Users/josephramsey/tetrad/comparison2040", simulations, algorithms, statistics, parameters);
    }

    @Test
    public void test3b() {
        int[] numMeasures = {600};//{300, 600, 1200, 2400, 500};
        double[] avgDegree = {6};//{2.0, 2.8, 3.5, 4.0, 10.0};
        int[] sampleSize = {400};//{100, 200, 300, 400, 300};

        boolean faithfulness = true;
        boolean zhangshen = false;

        int[] indices = {0};

        zhangshen = true;

        for (int i : indices) {
            visit3b(numMeasures[i], avgDegree[i], sampleSize[i], faithfulness, zhangshen);
        }

    }

    private void visit3b(int numMeasures, double value, int value1, boolean faithfulness, boolean zhangshen) {
        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 1);

        parameters.set(Params.NUM_MEASURES, numMeasures);
        parameters.set(Params.AVG_DEGREE, value);

//        parameters.set(Params.SAMPLE_SIZE, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000);
        parameters.set(Params.SAMPLE_SIZE, value1);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1);
        parameters.set(Params.VAR_HIGH, 3);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);
        parameters.set(Params.SEM_IM_SIMULATION_TYPE, false);

        parameters.set(Params.SYMMETRIC_FIRST_STEP, false);
        parameters.set(Params.VERBOSE, true);
        parameters.set(Params.FAITHFULNESS_ASSUMED, faithfulness);
        parameters.set(Params.PARALLELISM, 20);

        parameters.set(Params.SEM_GIC_RULE, 4, 6);
        parameters.set(Params.EBIC_GAMMA, .8, .9, .95, 1);
        parameters.set(Params.TAKE_LOGS, true);
        parameters.set(Params.TRUE_ERROR_VARIANCE, 1);
        parameters.set(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED, false);
        parameters.set(Params.ZS_RISK_BOUND, 0, 0.01, 0.1, .2, .3);

        parameters.set(Params.USE_MAX_P_ORIENTATION_HEURISTIC, true);
        parameters.set(Params.MAX_P_ORIENTATION_MAX_PATH_LENGTH, -1);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.SEM_GIC_RULE));
        statistics.add(new ParameterColumn(Params.ZS_RISK_BOUND));
        statistics.add(new ParameterColumn(Params.EBIC_GAMMA));
        statistics.add(new NumberOfNodesTrue());
        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new AdjacencyTPR());
        statistics.add(new AdjacencyFPR());
        statistics.add(new AdjacencyTPRDirectedParts());
        statistics.add(new AdjacencyFPRDirectedParts());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

//        if (zhangshen) {
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.ZhangShenBoundScore()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.EbicScore()));
//        } else {
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.KimEtAlScores()));
//            algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                    new edu.cmu.tetrad.algcomparison.score.SemBicScore()));
//        }

        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));
//        simulations.add(new LinearFisherModel(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setSaveGraphs(false);
        comparison.setSaveData(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        String file = "Comparison." + numMeasures;
        if (faithfulness) file += ".fa";
        if (zhangshen) file += ".zs";
        file += ".txt";

        comparison.compareFromSimulations(
                "/Users/josephramsey/tetrad/comparison3b", simulations, file, algorithms, statistics, parameters);
    }

    @Test
    public void test3c() {
        int numMeasures = 50;
        double avgDegree = 4;
        Integer[] sampleSizes = new Integer[]{200, 500, 1000, 2000, 5000, 10000};//, 20000,
//            50000};
//        Integer[] sampleSizes = new Integer[]{1000, 10000, 100000};

        boolean faithfulness = true;
        boolean zhangshen = true;

        visit3c(numMeasures, avgDegree, sampleSizes, faithfulness, zhangshen);

    }

    private void visit3c(int numMeasures, double avgDegree, Integer[] sampleSizes, boolean faithfulness, boolean zhangshen) {
        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 2);

        parameters.set(Params.NUM_MEASURES, numMeasures);
        parameters.set(Params.AVG_DEGREE, avgDegree);

        parameters.set(Params.SAMPLE_SIZE, sampleSizes);
        parameters.set(Params.COEF_LOW, 0.1);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1);
        parameters.set(Params.VAR_HIGH, 2);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);
        parameters.set(Params.SEM_IM_SIMULATION_TYPE, false);

        parameters.set(Params.SYMMETRIC_FIRST_STEP, false);
        parameters.set(Params.VERBOSE, true);
        parameters.set(Params.FAITHFULNESS_ASSUMED, faithfulness);
        parameters.set(Params.PARALLELISM, 20);

        parameters.set(Params.SEM_GIC_RULE, 1);
        parameters.set(Params.PENALTY_DISCOUNT, 2);
        parameters.set(Params.TAKE_LOGS, true);
        parameters.set(Params.TRUE_ERROR_VARIANCE, 1);
        parameters.set(Params.CALCULATE_EUCLIDEAN_NORM_SQUARED, false);
        parameters.set(Params.ZS_RISK_BOUND, 0.01);//0, .1, .2, .3, .4, .5, .6, .7, .8, .9);

        parameters.set(Params.USE_MAX_P_ORIENTATION_HEURISTIC, true);
        parameters.set(Params.MAX_P_ORIENTATION_MAX_PATH_LENGTH, -1);

        parameters.set(Params.COLLIDER_DISCOVERY_RULE, 3);
        parameters.set(Params.ALPHA, 0.01);
//        parameters.set(Params.CONFLICT_RULE);

        parameters.set(Params.SEM_GIC_RULE, 1, 2, 3, 4, 5, 6);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.SEM_GIC_RULE));
        statistics.add(new ParameterColumn(Params.ZS_RISK_BOUND));
        statistics.add(new NumberOfNodesTrue());
        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
//        statistics.add(new AdjacencyRecall());
//        statistics.add(new AdjacencyTPR());
//        statistics.add(new AdjacencyFPR());
//        statistics.add(new AdjacencyTPRDirectedParts());
//        statistics.add(new AdjacencyFPRDirectedParts());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new AdjacencyTP());
        statistics.add(new AdjacencyFP());
//        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new Misorientations());
        statistics.add(new AhpcCriterion());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.ZhangShenBoundScore()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.EbicScore()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.KimEtAlScores()));
//            algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
//                    new edu.cmu.tetrad.algcomparison.score.SemBicScore()));


//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll(new FisherZ()));

        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));
//        simulations.add(new LinearFisherModel(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setSaveGraphs(false);
        comparison.setSaveData(false);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        String file = "Comparison." + numMeasures;
        file += "." + avgDegree;
        if (faithfulness) file += ".fa";
        if (zhangshen) file += ".zs";
        file += ".txt";

        comparison.compareFromSimulations(
                "/Users/josephramsey/tetrad/comparison3c", simulations, file, algorithms, statistics, parameters);
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

        parameters.set("coefLow", 0.199);
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
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);
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

        rules.setRevertToUnshieldedColliders(false);
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
                score.setRuleType(SemBicScore.RuleType.CHICKERING);
                score.setPenaltyDiscount(2);

                Fges fges = new Fges(score, 4);
                fges.setMaxDegree(100);
                fges.setSymmetricFirstStep(false);

                Graph out = fges.search();
                Graph graph = loadRSpecial(graphFile);

                List<Node> nodes = graph.getNodes();

                for (int k = 1; k <= nodes.size(); k++) {
                    String name = "X" + k;
                    Node node = graph.getNode(name);

                    if (node == null) {
                        System.out.println(name);
                    } else {
                        node.setName("C" + k);
                    }
                }

                graph = GraphUtils.replaceNodes(graph, out.getNodes());

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
    public void test7a() {
        NumberFormat nf = new DecimalFormat("0000");

        String prefix = "tges1";

        System.out.println("\nComparison of graph files with prefix '" + prefix + "' to true graph");

        for (int i = 9200; i <= 9209; i++) {
            File dataFile = new File("/Users/josephramsey/Downloads/ges_turning"
                    + "/raw_data_" + nf.format(i) + ".csv");
            File graphFile = new File("/Users/josephramsey/Downloads/ges_turning"
                    + "/true_model_" + nf.format(i));

//            File gesFile = new File("/Users/josephramsey/Downloads/ges_turning/" + prefix + "_N10000_dag_" + nf.format(i) + ".csv");

            DataSet rawdata = null;

            try {
                ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(dataFile.toPath(), Delimiter.COMMA);
                reader.setHasHeader(false);
                Data data = reader.readInData();
                rawdata = (DataSet) DataConvertUtils.toDataModel(data);
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }

            SemBicScore score = new SemBicScore(rawdata);
            score.setRuleType(SemBicScore.RuleType.CHICKERING);
            score.setPenaltyDiscount(2);

            Fges fges = new Fges(score, 4);
            fges.setMaxDegree(100);
            fges.setSymmetricFirstStep(false);

            Graph trueDag = loadMyGraph(graphFile, true);
//            Graph gesGraph = loadMyGraph(gesFile, false);

            Graph est = fges.search();
            trueDag = GraphUtils.replaceNodes(trueDag, est.getNodes());

            Graph truedag = SearchGraphUtils.dagFromPattern(trueDag);
            Graph dagout = SearchGraphUtils.dagFromPattern(est);

            System.out.println("\nModel " + i + "\n");
            NumberFormat nf2 = new DecimalFormat("0.####");

            System.out.println("True # edges = " + nf2.format(new NumberOfEdgesTrue().getValue(truedag, dagout, rawdata)));
            System.out.println("Est # edges = " + nf2.format(new NumberOfEdgesEst().getValue(truedag, dagout, rawdata)));
            System.out.println("AP = " + nf2.format(new AdjacencyPrecision().getValue(truedag, dagout, rawdata)));
            System.out.println("AR = " + nf2.format(new AdjacencyRecall().getValue(truedag, dagout, rawdata)));
            System.out.println("AHP = " + nf2.format(new ArrowheadPrecision().getValue(truedag, dagout, rawdata)));
            System.out.println("AHR = " + nf2.format(new ArrowheadRecall().getValue(truedag, dagout, rawdata)));
            System.out.println("F1Adj = " + nf2.format(new F1Adj().getValue(truedag, dagout, rawdata)));
            System.out.println("F1Arrow = " + nf2.format(new F1Arrow().getValue(truedag, dagout, rawdata)));
            System.out.println("BIC true = " + nf2.format(new BicTrue().getValue(truedag, dagout, rawdata)));
            System.out.println("BIC est = " + nf2.format(new BicEst().getValue(truedag, dagout, rawdata)));
            System.out.println("BIC diff = " + nf2.format(new BicDiff().getValue(truedag, dagout, rawdata)));
            System.out.println("Cycle in est DAG? = " + nf2.format(new CyclicEst().getValue(truedag, dagout, rawdata)));

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
        int numNodes = 20;
        int aveDegree = 3;
        int numIterations = 1;

        for (int i = 0; i < numIterations; i++) {
            Graph dag = GraphUtils.randomDag(numNodes, 0, aveDegree * numNodes / 2,
                    100, 100, 100, false);
            GraphScore score = new GraphScore(dag);
//            score.setVerbose(true);
            Fges fges = new Fges(score, 1);
//            fges.setSymmetricFirstStep(true);
            fges.setVerbose(true);
//            fges.setTDepth(1);
            Graph pattern1 = fges.search();
            Graph pattern2 = new Pc(new IndTestDSep(dag)).search();
            assertEquals(pattern2, pattern1);
        }
    }

    @Test
    public void test11() {
        int numNodes = 5;
        int aveDegree = 2;
        int numIterations = 5;

        int N = 100000;

        long seed = 1614264484974L;
//
        RandomUtil.getInstance().setSeed(seed);

        for (int i = 0; i < numIterations; i++) {
            System.out.println("seed = " + RandomUtil.getInstance().getSeed());

            Graph trueDag = GraphUtils.randomDag(numNodes, 0, aveDegree * numNodes / 2,
                    100, 100, 100, false);

            Parameters p = new Parameters();
            p.set(Params.VAR_LOW, 1);
            p.set(Params.VAR_HIGH, 1);
            p.set(Params.COEF_LOW, 0);
            p.set(Params.COEF_HIGH, 1);

            SemPm pm = new SemPm(trueDag);
            SemIm im = new SemIm(pm, p);
            DataSet D = im.simulateData(N, false);

            SemBicScore score = new SemBicScore(D);
            score.setRuleType(SemBicScore.RuleType.CHICKERING);
            score.setPenaltyDiscount(1);
            score.setVerbose(false);
            Fges fges = new Fges(score, 1);
            fges.setTrueDag(trueDag);
            fges.setVerbose(true);
            fges.setMeekVerbose(true);
            Graph estPattern = fges.search();

            estPattern = GraphUtils.replaceNodes(estPattern, trueDag.getNodes());
            Graph estDag = SearchGraphUtils.dagFromPattern(estPattern);

            System.out.println("True DAG = " + trueDag);
            System.out.println("FGES pattern = " + estPattern);
            System.out.println("FGES sample DAG = " + estDag);

            Graph truePattern = SearchGraphUtils.patternFromDag(trueDag);

            double ap = new AdjacencyPrecision().getValue(truePattern, estPattern, D);
            double ar = new AdjacencyRecall().getValue(truePattern, estPattern, D);
            double ahpc = new ArrowheadPrecisionCommonEdges().getValue(truePattern, estPattern, D);
            double ahrc = new ArrowheadRecallCommonEdges().getValue(truePattern, estPattern, D);

            NumberFormat nf = new DecimalFormat("0.00");

            System.out.println("AP = " + nf.format(ap) + " AR = " + nf.format(ar) + " AHPC = " + nf.format(ahpc) + " AHRC = " + nf.format(ahrc));

            assert estPattern != null;
            assertTrue("Markov failed", markovTo(estDag, trueDag));
        }
    }

    private boolean markovTo(Graph estDag, Graph trueDag) {
        List<Node> variables = estDag.getNodes();
        IndTestDSep dsepTrue = new IndTestDSep(trueDag);

        boolean ret = true;

        for (Node x : variables) {
            List<Node> parents = estDag.getParents(x);

            for (Node y : variables) {
                if (parents.contains(y)) continue;
                if (x == y) continue;
                if (estDag.isAncestorOf(x, y)) continue;

                boolean ind1 = dsepTrue.isIndependent(x, y, parents);

                System.out.println("ind1 = " + ind1);

                if (!ind1) {
                    System.out.println("Not Markov to the true graph "
                            + " because " + SearchLogUtils.independenceFact(x, y, parents)
                            + " not d-separated in true graph.");
                    ret = false;
                }
            }
        }

        return ret;
    }

    @Test
    public void test12() {
        int pn = 100;

        for (int m0 : new int[]{0, 1, 2, 3, 4}) {
            for (int n : new int[]{1000}) {
                for (double lambda = 1; lambda <= 20; lambda += 0.5) {
//                        for (double gamma = 1; gamma <= 1; gamma++) {
//                            double lambda = 2 * log(pn) + (1 + gamma) * log(log(pn));
                    double pmin = 2 - pow(1 + exp(-(lambda - 1) / 2.) * sqrt(lambda), pn - m0);

                    double pd = lambda / log(n);

                    System.out.println("pn = " + pn + " m0 = " + m0 + " lambda = "
                            + lambda + " pmin = " + pmin + " pd = " + pd);
                }
            }
        }
    }

    @Test
    public void test13() {
        int pn = 500;
        int n = 300;
        double riskBound = 0.01;

        NumberFormat nf = new DecimalFormat("0.00000");

        for (int m0 = 0; m0 < pn; m0++) {
            double lambda = ZhangShenBoundScore.zhangShenLambda(m0, pn, riskBound);

            double p = ZhangShenBoundScore.getP(pn, m0, lambda);


            double pd = lambda / log(n);

            System.out.println(" m0 = " + m0 + " pn = " + pn + " riskBound = " + nf.format(riskBound)
                    + " lambda = " + nf.format(lambda) + " pd = " + nf.format(pd) + " recovered risk bound = " + nf.format(1 - p));
        }
    }


    @Test
    public void test14() {
        int pn = 8;
        int m0 = 4;
        int n = 1000;

        NumberFormat nf = new DecimalFormat("0.00");

        double maxLambda = Double.NEGATIVE_INFINITY;

        for (double pmin = 1.0; pmin >= 0.0; pmin -= 0.01) {
            double lambda = ZhangShenBoundScore.zhangShenLambda(m0, pn, 1 - pmin);
            double pd = lambda / log(n);

            System.out.println(" m0 = " + m0 + " pn = " + pn + " pmin = " + nf.format(pmin)
                    + " lambda = " + nf.format(lambda) + " pd = " + nf.format(pd));
        }
    }

    @Test
    public void test15() {
        int pn = 15000;
        int n = 300;

        double[] lambda = new double[7];

        lambda[1] = log(n);
        lambda[2] = pow(pn, 1. / 3.);
        lambda[3] = 2 * (log(pn));
        lambda[4] = 2 * (log(pn) + log(log(pn)));
        lambda[5] = log(log(n)) * (log(pn));
        lambda[6] = log(n) * log(pn);

        NumberFormat nf = new DecimalFormat("0.00");

        for (int i = 1; i <= 6; i++) {
//            System.out.println("pn = " + pn + " n = " + n + " lambda" + i + " = " + nf.format(lambda[i])
//                    + " effective penalty discount = " + nf.format((2 * lambda[i]) / log(n)));
            System.out.println("GIC" + i + " lambda = " + nf.format(lambda[i])
                    + " penalty discount = " + nf.format((lambda[i]) / log(n)));
        }
    }

    @Test
    public void test16() {

        Graph gStar = new EdgeListGraph();
        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");

        gStar.addNode(x1);
        gStar.addNode(x2);
        gStar.addNode(x3);
        gStar.addNode(x4);

        gStar.addDirectedEdge(x1, x2);
        gStar.addDirectedEdge(x2, x3);
        gStar.addDirectedEdge(x3, x4);
        gStar.addDirectedEdge(x1, x4);

        SemPm semPm = new SemPm(gStar);
        SemIm semIm = new SemIm(semPm);

        Parameters parameters = new Parameters();
        parameters.set(Params.SAMPLE_SIZE, 1000000);

        StandardizedSemIm sem = new StandardizedSemIm(semIm, parameters);

        for (double d2 = -1; d2 <= 1; d2 += 0.1) {
            for (double d3 = -1; d3 <= 1; d3 += 0.1) {
                for (double d4 = -1; d4 <= 1; d4 += 0.1) {
                    for (double d1 = -1; d1 <= 1; d1 += 0.1) {
                        if ((d1 > -0.801 && d1 < -0.199) || (d1 > 0.199 && d1 < 0.801)) {
                            if ((d2 > -0.801 && d2 < -0.199) || (d2 > 0.199 && d2 < 0.801)) {
                                if ((d3 > -0.801 && d3 < -0.199) || (d3 > 0.199 && d3 < 0.801)) {
                                    if ((d4 > -0.801 && d4 < -0.199) || (d4 > 0.199 && d4 < 0.801)) {
                                        d1 = Math.round(d1 * 10000.00) / 10000.00;
                                        d2 = Math.round(d2 * 10000.00) / 10000.00;
                                        d3 = Math.round(d3 * 10000.00) / 10000.00;
                                        d4 = Math.round(d4 * 10000.00) / 10000.00;

//                                        System.out.println("d1 = " + d1 + " d2 = " + d2 + " d3 = " + d3 + " d4 = " + d4);

                                        try {
                                            assertTrue(sem.setEdgeCoefficient(x1, x2, d1));
                                            assertTrue(sem.setEdgeCoefficient(x2, x3, d2));
                                            assertTrue(sem.setEdgeCoefficient(x3, x4, d3));
                                            assertTrue(sem.setEdgeCoefficient(x1, x4, d4));

                                            DataSet dataSet = sem.simulateDataReducedForm(50000, false);

                                            CovarianceMatrix covarianceMatrix = new CovarianceMatrix(dataSet);
                                            IndependenceTest test = new IndTestFisherZ(covarianceMatrix, 0.0001);

                                            Node _x1 = dataSet.getVariable("X1");
                                            Node _x2 = dataSet.getVariable("X2");
                                            Node _x3 = dataSet.getVariable("X3");
                                            Node _x4 = dataSet.getVariable("X4");

                                            if (!test.isIndependent(_x1, _x2, Collections.singletonList(_x4))) continue;

                                            Fges fges = new Fges(new SemBicScore(dataSet));

                                            Graph out = fges.search();

                                            if (out.getNumEdges() == 5) {
                                                System.out.println("\nd1 = " + d1 + " d2 = " + d2 + " d3 = " + d3 + " d4 = " + d4);
                                                System.out.println(out);
                                                System.out.println(sem.getImplCovar());
                                            }
                                        } catch (AssertionError e) {
                                            // continue;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


    }
}
