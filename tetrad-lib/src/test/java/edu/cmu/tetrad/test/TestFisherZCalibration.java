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
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TextTable;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.log;
import static java.lang.StrictMath.abs;

public class TestFisherZCalibration {

//    public static void main(String... args) {
//        test1();
//    }

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(105034020L);
        toTest();
    }

    private void toTest() {
        Parameters parameters = new Parameters();
        parameters.set(Params.ALPHA, 0.01);
        parameters.set(Params.DEPTH, 2);
        parameters.set(Params.STRUCTURE_PRIOR, 0);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.SAMPLE_SIZE, 1000000);
        parameters.set(Params.SEM_BIC_RULE, 3);
        parameters.set(Params.PENALTY_DISCOUNT, 2);

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
        fges.setTrueGraph(graph);

        fges.search();

//        fges = new Fges(new GraphScore(graph));
//        fges.setVerbose(true);
//        fges.search();


    }

    @Test
    public void test3() {

//        RandomUtil.getInstance().setSeed(92883342449L);

        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 30);
        parameters.set(Params.NUM_MEASURES, 20);
        parameters.set(Params.AVG_DEGREE, 4.);
        parameters.set(Params.SAMPLE_SIZE, 5000);//, 50000, 500000);//, 200000, 500000);// 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000);
        parameters.set(Params.SEM_BIC_RULE, 3);
        parameters.set(Params.PENALTY_DISCOUNT, 1);
        parameters.set(Params.STRUCTURE_PRIOR, 0);
        parameters.set(Params.TDEPTH, -1);
        parameters.set(Params.FAITHFULNESS_ASSUMED, false);
        parameters.set(Params.TURNING, true);
        parameters.set(Params.COEF_LOW, 0);
        parameters.set(Params.COEF_HIGH, 1);
        parameters.set(Params.VAR_LOW, 1.0);
        parameters.set(Params.VAR_HIGH, 3.0);
        parameters.set(Params.VERBOSE, false);
        parameters.set(Params.RANDOMIZE_COLUMNS, true);

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn(Params.NUM_RUNS));
        statistics.add(new ParameterColumn(Params.SEM_BIC_RULE));
        statistics.add(new ParameterColumn(Params.NUM_MEASURES));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.PENALTY_DISCOUNT));

        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
//        statistics.add(new ArrowheadPrecision());
//        statistics.add(new ArrowheadRecall());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new AhpBound());

        statistics.add(new F1Adj());
//        statistics.add(new F1Arrow());
        statistics.add(new SHD());
        statistics.add(new ElapsedTime());

        statistics.setWeight("SHD", .5);

        Algorithms algorithms = new Algorithms();

        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges(
                new edu.cmu.tetrad.algcomparison.score.SemBicScore()));
//
        Simulations simulations = new Simulations();

        simulations.add(new SemSimulation(new RandomForward()));

        Comparison comparison = new Comparison();

//        comparison.setShowAlgorithmIndices(true);
//        comparison.setShowSimulationIndices(true);
//        comparison.setShowUtilities(true);
//        comparison.setSortByUtility(true);
        comparison.setSaveGraphs(false);
//        comparison.setSavePags(true);
        comparison.setSaveData(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);

        comparison.compareFromSimulations(
                "/Users/user/tetrad4/comparison2040", simulations, algorithms, statistics, parameters);
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
        parameters.set("depth", 1);
        parameters.set(Params.FAITHFULNESS_ASSUMED, true);

        parameters.set("coefLow", 0.2);
        parameters.set("coefHigh", 0.7);
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
}
