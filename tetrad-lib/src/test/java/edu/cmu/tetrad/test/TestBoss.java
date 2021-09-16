///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.BOSS;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.PcMax;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Fci;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.FciMax;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.*;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Rfci;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.graph.SingleGraph;
import edu.cmu.tetrad.algcomparison.independence.DSeparationTest;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore;
import edu.cmu.tetrad.algcomparison.score.ZhangShenBoundScore;
import edu.cmu.tetrad.algcomparison.simulation.LinearSemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.LinearSemIm;
import edu.cmu.tetrad.sem.LinearSemPm;
import edu.cmu.tetrad.sem.StandardizedLinearSemIm;
import edu.cmu.tetrad.util.*;
import org.apache.commons.collections4.OrderedMap;
import org.apache.commons.collections4.map.ListOrderedMap;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.search.Boss.Method.GSP;
import static java.util.Collections.shuffle;
import static java.util.Collections.sort;
import static org.junit.Assert.assertEquals;

/**
 * Tests to make sure the DelimiterType enumeration hasn't been tampered with.
 *
 * @author Joseph Ramsey
 */
@SuppressWarnings("ALL")
public final class TestBoss {

    @NotNull
    private static Boss getBoss(Score score, IndependenceTest test) {
        Boss boss;

        if (true) {
            boss = new Boss(test);
        } else {
            boss = new Boss(score);
        }

        return boss;
    }

    private static boolean printFailed(Graph g, Graph dag, String alg) {
        double ap = new AdjacencyPrecision().getValue(g, dag, null);
        double ar = new AdjacencyRecall().getValue(g, dag, null);
        double ahp = new ArrowheadPrecision().getValue(g, dag, null);
        double ahr = new ArrowheadRecall().getValue(g, dag, null);

        NumberFormat nf = new DecimalFormat("0.00");

        if (dag.getNumEdges() != g.getNumEdges()) {
            System.out.println("Failed " + alg +
                    " ap = " + nf.format(ap) + " ar = " + nf.format(ar)
                    + " ahp = " + nf.format(ahp) + " ahr = " + nf.format(ahr));
            return true;
        }

        return false;
    }

    private static void runTestLoop(Graph g, List<Node> order, Score score, IndependenceTest test, boolean useTest) {
        g = new EdgeListGraph(g);
        order = new ArrayList<>(order);

        Boss.Method[] methods = {Boss.Method.BOSS, Boss.Method.SP};

        for (Boss.Method method : methods) {
            Boss boss = getBoss(score, test);

            boss.setCacheScores(true);
            boss.setMethod(method);
            boss.setNumStarts(1);
            boss.setBreakTies(true);
            boss.setGspDepth(-1);
            boss.setVerbose(true);
            List<Node> perm = boss.bestOrder(order);
            Graph dag = boss.getGraph(perm, false);

            printFailed(g, dag, method + " " + order + " \n" + dag);
        }
    }

    @Test
    public void testBoss() {
        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 100);
        params.set(Params.AVG_DEGREE, 2, 4, 6);
        params.set(Params.SAMPLE_SIZE, 10000);
        params.set(Params.NUM_RUNS, 30);
        params.set(Params.RANDOMIZE_COLUMNS, false);
        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.COEF_LOW, 0.1);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.BOSS_METHOD, 1);
        params.set(Params.BREAK_TIES, true);


        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new ZhangShenBoundScore(), new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.DEPTH));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/testBoss",
                simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss2() {
        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 20);
        params.set(Params.AVG_DEGREE, 4);
        params.set(Params.SAMPLE_SIZE, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000,
                50000, 100000, 200000);
        params.set(Params.NUM_RUNS, 5);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.COEF_LOW, 0.2);
        params.set(Params.COEF_HIGH, 0.8);
        params.set(Params.VERBOSE, false);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.USE_SCORE, true);


        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new LinearGaussianBicScore(), new FisherZ()));
//        algorithms.add(new GSP(new SemBicScore()));
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new PcAll(new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new F1Adj());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/testBoss2",
                simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss3() {
        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 100);
        params.set(Params.AVG_DEGREE, 10);
        params.set(Params.SAMPLE_SIZE, 1000);
        params.set(Params.NUM_RUNS, 1);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VERBOSE, true);
        params.set(Params.COLLIDER_DISCOVERY_RULE, 2, 3);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.CACHE_SCORES, false);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new ZhangShenBoundScore(), new FisherZ()));
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new PcAll(new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/testBoss3",
                simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss4() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 200);
        params.set(Params.NUM_MEASURES, 100);
        params.set(Params.AVG_DEGREE, 2);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, 0);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 1);

        params.set(Params.BOSS_METHOD, 1);
        params.set(Params.BOSS_SCORE_TYPE, false);
        params.set(Params.BREAK_TIES, false);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);

        params.set(Params.PENALTY_DISCOUNT, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0);
//
        params.set(Params.EBIC_GAMMA, 0.5);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new edu.cmu.tetrad.algcomparison.score.FmlBicScore(), new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.PENALTY_DISCOUNT));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new AdjacencyTPR());
        statistics.add(new AdjacencyFPR());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/testBoss4", simulations, algorithms, statistics, params);
    }

    @Test
    public void testLuFigure3() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000);
        params.set(Params.NUM_MEASURES, 20);
        params.set(Params.AVG_DEGREE, 4);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, 0.2);
        params.set(Params.COEF_HIGH, 0.8);
        params.set(Params.COV_LOW, 0);
        params.set(Params.COV_HIGH, 0);
        params.set(Params.VAR_LOW, 1);
        params.set(Params.VAR_HIGH, 3);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 30);

        params.set(Params.BOSS_SCORE_TYPE, false);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);

        params.set(Params.PENALTY_DISCOUNT, 2.0);
        params.set(Params.ALPHA, 0.001);
//
        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new LinearGaussianBicScore(), new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new F1Adj());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
//        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/Lu.figure.3", simulations,
                algorithms, statistics, params);
    }

    @Test
    public void testLuFigure6() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 500);
        params.set(Params.NUM_MEASURES, 60);
        params.set(Params.AVG_DEGREE, 2, 4, 6, 8, 10, 12);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, 0.2);
        params.set(Params.COEF_HIGH, 0.8);
        params.set(Params.COV_LOW, 0);
        params.set(Params.COV_HIGH, 0);
        params.set(Params.VAR_LOW, 1);
        params.set(Params.VAR_HIGH, 3);
//        params.set(Params.COEF_SYMMETRIC, false);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 5);

        params.set(Params.BOSS_SCORE_TYPE, false);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.CACHE_SCORES, false);
        params.set(Params.NUM_STARTS, 1);

        params.set(Params.PENALTY_DISCOUNT, 2.0);
        params.set(Params.ALPHA, 0.001);
//
        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new LinearGaussianBicScore(), new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Cpc(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.PcMax(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison
                .algorithm.oracle.cpdag.Fges(new LinearGaussianBicScore()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
//        statistics.add(new ParameterColumn(Params.BOSS_METHOD));
//        statistics.add(new ParameterColumn(Params.BOSS_SCORE_TYPE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/Lu.figure.6", simulations,
                algorithms, statistics, params);
    }

    @Test
    public void testClark() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 500);
//        params.set(Params.NUM_MEASURES, 20);
//        params.set(Params.AVG_DEGREE, 2, 4, 6, 8);//, 10, 12);//, 14);//, 1 6, 18, 20);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, .2);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 100);

        params.set(Params.BOSS_METHOD, 1);
        params.set(Params.BOSS_SCORE_TYPE, false);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.CACHE_SCORES, false);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.USE_SCORE, true);

        params.set(Params.COLLIDER_DISCOVERY_RULE, 2);
        params.set(Params.PENALTY_DISCOUNT, 1.0);
        params.set(Params.ALPHA, 0.001);
//
        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new edu.cmu.tetrad.algcomparison.score.FmlBicScore(), new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison
                .algorithm.oracle.cpdag.PcAll(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison
                .algorithm.oracle.cpdag.Fges(new LinearGaussianBicScore()));

        Node x1 = new ContinuousVariable("X1");
        Node x2 = new ContinuousVariable("X2");
        Node x3 = new ContinuousVariable("X3");
        Node x4 = new ContinuousVariable("X4");
        Node x5 = new ContinuousVariable("X5");
        Node x6 = new ContinuousVariable("X6");
        Node x7 = new ContinuousVariable("X7");

        List<Node> nodes = new ArrayList<>();

        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);
        nodes.add(x5);
        nodes.add(x6);
        nodes.add(x7);

        Graph graph = new EdgeListGraph(nodes);

        graph.addDirectedEdge(x1, x2);
        graph.addDirectedEdge(x1, x3);
        graph.addDirectedEdge(x1, x4);
        graph.addDirectedEdge(x1, x5);
        graph.addDirectedEdge(x1, x6);

        graph.addDirectedEdge(x2, x7);
        graph.addDirectedEdge(x3, x7);
        graph.addDirectedEdge(x4, x7);
        graph.addDirectedEdge(x5, x7);
        graph.addDirectedEdge(x6, x7);

        System.out.println(graph);

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new SingleGraph(graph)));

        Statistics statistics = new Statistics();
//        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
//        statistics.add(new ParameterColumn(Params.BOSS_METHOD));
//        statistics.add(new ParameterColumn(Params.BOSS_SCORE_TYPE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/clark", simulations,
                algorithms, statistics, params);
    }


    @Test
    public void testBoss8() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 10000);
        params.set(Params.NUM_MEASURES, 5);
        params.set(Params.AVG_DEGREE, 1, 2, 3, 4);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, 0.1);
        params.set(Params.COEF_HIGH, 1.5);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 100);

        params.set(Params.BOSS_METHOD, 1);
        params.set(Params.BOSS_SCORE_TYPE, true);
        params.set(Params.BREAK_TIES, false);
        params.set(Params.CACHE_SCORES, false);
        params.set(Params.NUM_STARTS, 1);

        params.set(Params.PENALTY_DISCOUNT, 1);//.5, 1, 2, 3, 4, 5);
        params.set(Params.ALPHA, 0.001);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new LinearGaussianBicScore(), new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Cpc(new FisherZ()));
        algorithms.add(new PcMax(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Fges(new ZhangShenBoundScore()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.NUM_MEASURES));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
//        statistics.add(new ParameterColumn(Params.BOSS_METHOD));
//        statistics.add(new ParameterColumn(Params.BOSS_SCORE_TYPE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.PENALTY_DISCOUNT));
//        statistics.add(new CorrectCPDAG());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_CPDAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/testBoss8", simulations,
                algorithms, statistics, params);
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    //@Test
    public void testAllFacts() {
        List<Ret> allFacts = new ArrayList<>();
        allFacts.add(getFactsSimpleCanceling());
        allFacts.add(getFactsRaskutti());
        allFacts.add(getFigure6());
        allFacts.add(getFigure7());
        allFacts.add(getFigure8());
        allFacts.add(getFigure12());

        int count = 0;

        boolean printCpdag = false;


        Boss.Method[] methods = {GSP, Boss.Method.BOSS, Boss.Method.SP};

        for (Ret facts : allFacts) {
            count++;

            OrderedMap<String, Set<Graph>> graphs = new ListOrderedMap<>();
            OrderedMap<String, Set<String>> labels = new ListOrderedMap<>();

            System.out.println();
            System.out.println(facts.getLabel());
            System.out.println(facts.getFacts());

            List<Node> order = facts.facts.getVariables();

            int numRuns = 100;

            for (int t = 0; t < numRuns; t++) {
                shuffle(order);
                IndTestDSep test = new IndTestDSep(facts.getFacts());

                for (Boss.Method method : methods) {
                    Boss boss = new Boss(test);
                    boss.setCacheScores(true);
                    boss.setMethod(method);
                    boss.setNumStarts(1);

                    List<Node> perm = boss.bestOrder(test.getVariables());
                    Graph cpdag = boss.getGraph(perm, true);

                    if (graphs.get(method.toString()) == null) {
                        graphs.put(method.toString(), new HashSet<>());
                        labels.put(method.toString(), new HashSet<>());
                    }

                    graphs.get(method.toString()).add(cpdag);
                    labels.get(method.toString()).add(facts.getLabel());

                    if (printCpdag) {
                        System.out.println(cpdag);
                    }
                }
            }

            for (String label : labels.keySet()) {
                printGraphs(label, graphs);
            }
        }
    }

    //@Test
    public void testRaskutti() {
        Ret facts = getFactsRaskutti();

        IndTestDSep test = new IndTestDSep(facts.getFacts());

        Boss boss = new Boss(test);
        boss.setCacheScores(true);
        boss.setMethod(Boss.Method.BOSS);
        boss.setNumStarts(1);

        List<Node> variables = test.getVariables();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);

            List<Node> p2 = boss.bestOrder(test.getVariables());
            Graph cpdag = boss.getGraph(p2, true);

            System.out.println(p + " " + cpdag.getNumEdges());
        }
    }

    //@Test
    public void testFromData() {
        for (int i = 0; i < 10; i++) {
            System.out.println("\nRun " + (i + 1));

            Graph g = GraphUtils.randomGraph(15, 0, 30, 100,
                    100, 100, false);
            LinearSemPm pm = new LinearSemPm(g);
            Parameters params = new Parameters();
            params.set(Params.COEF_LOW, 0.1);
            params.set(Params.COEF_HIGH, 1);
            LinearSemIm im = new LinearSemIm(pm, params);
            DataSet data = im.simulateData(100000, false);

            data = DataUtils.shuffleColumns(data);
            List<Node> order = data.getVariables();

            edu.cmu.tetrad.search.EbicScore score = new edu.cmu.tetrad.search.EbicScore(data);
//            score.setPenaltyDiscount(2);
            edu.cmu.tetrad.search.IndependenceTest test = new edu.cmu.tetrad.search.IndTestFisherZ(data, 0.01);

            runTestLoop(g, order, score, test, false);
        }
    }

    //@Test
    public void testFromDsep() {
        for (int i = 0; i < 10; i++) {
            System.out.println("\nRun " + (i + 1));

            List<Node> nodes = new ArrayList<>();
            for (int t = 0; t < 4; t++) {
                nodes.add(new GraphNode("X" + (t + 1)));
            }

            Graph g = GraphUtils.randomGraph(nodes, 0, 4, 100,
                    100, 100, false);
            g = new EdgeListGraph(g);

            List<Node> order = g.getNodes();
            shuffle(order);

            IndependenceTest test = new IndTestDSep(g);

            runTestLoop(g, order, null, test, true);
        }
    }

    @Test
    public void testFromDsepSolus() {
        Parameters parameters = new Parameters();

        parameters.set(Params.NUM_MEASURES, 10);
        parameters.set(Params.AVG_DEGREE, 1, 2, 3, 4, 5, 6, 7, 8, 9);
        parameters.set(Params.NUM_RUNS, 1);
        parameters.set(Params.BOSS_METHOD, 1);
        parameters.set(Params.BOSS_SCORE_TYPE, true);
        parameters.set(Params.BREAK_TIES, true);
        parameters.set(Params.SAMPLE_SIZE, 10000);
        parameters.set(Params.ALPHA, 0.001);
        parameters.set(Params.DEPTH, -1);
        parameters.set(Params.USE_SCORE, false);

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));
//        simulations.add(new SemSimulationTrueModel(new RandomForward()));

        Algorithms algorithms = new Algorithms();
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.PcAll(new DSeparationTest()));
        algorithms.add(new BOSS(new LinearGaussianBicScore(), new DSeparationTest()));
//        algorithms.add(new BOSS(new edu.cmu.tetrad.algcomparison.score.FmlBicScore()));
//        algorithms.add(new BOSSIndep(new FisherZ()));
//        algorithms.add(new GSPIndep(new DSeparationTest()));
//        algorithms.add(new GSPIndep(new FisherZ()));

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/soluscomparison",
                simulations, algorithms, statistics, parameters);
    }

    //@Test
    public void testFromDsep1() {

        List<Node> nodes = new ArrayList<>();

        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");

        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);

        List<Edge> edges = new ArrayList<>();
        edges.add(Edges.directedEdge(x1, x2));
        edges.add(Edges.directedEdge(x1, x3));
        edges.add(Edges.directedEdge(x2, x4));
        edges.add(Edges.directedEdge(x3, x4));

        Graph g = new EdgeListGraph(nodes);
        for (Edge e : edges) g.addEdge(e);

        List<Node> order = new ArrayList<>();

        order.add(x4);
        order.add(x2);
        order.add(x3);
        order.add(x1);

        IndependenceTest test = new IndTestDSep(g);

        runTestLoop(g, order, null, test, true);
    }

    //@Test
    public void testFromDsep2() {
        List<Node> nodes = new ArrayList<>();

        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");

        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);

        List<Edge> edges = new ArrayList<>();
        edges.add(Edges.directedEdge(x1, x2));
        edges.add(Edges.directedEdge(x1, x3));
        edges.add(Edges.directedEdge(x2, x4));
        edges.add(Edges.directedEdge(x3, x4));

        for (int i = 0; i < 100; i++) {
            Graph g = new EdgeListGraph(nodes);

            shuffle(edges);

            for (Edge e : edges) {
                g.addEdge(e);
            }

            IndependenceTest test = new IndTestDSep(g);

            List<Node> order = new ArrayList<>();

            order.add(x1);
            order.add(x2);
            order.add(x3);
            order.add(x4);

            Collections.shuffle(order);
            runTestLoop(g, order, null, test, true);
//            runTestLoop2(g, order, null, test, true);
        }
    }

    private boolean isCpdagForDag(Graph cpdag, Graph dag) {
        if (!GraphUtils.undirectedGraph(cpdag).equals(GraphUtils.undirectedGraph(dag))) {
            return false;
        }

        return true;
    }

    private void printGraphs(String label, Map<String, Set<Graph>> graphs) {
        if (!graphs.containsKey(label)) return;

        List<Graph> _graphs = new ArrayList<>(graphs.get(label));

        System.out.println("===== " + label + "\n");

        for (int i = 0; i < _graphs.size(); i++) {
            System.out.println("Found this CPDAG (" + (i + 1) + "):\n\n" + _graphs.get(i) + "\n");
        }
    }

    //@Test
    public void test12345() {
        IndependenceFacts facts = getFigure6().facts;
        List<Node> nodes = facts.getVariables();

        sort(nodes);
        TeyssierScorer scorer = new TeyssierScorer(new IndTestDSep(facts));
        assert (scorer.score(nodes) == 7);
    }

    public Ret getFactsSimpleCanceling() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x3, list(x2)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));
        facts.add(new IndependenceFact(x3, x1, list(x2)));
        facts.add(new IndependenceFact(x1, x4, list())); // unfaithful.

        return new Ret("Simple 4-node path canceling model", facts);
    }

    public Ret getFigure8() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");
        Node x5 = new GraphNode("5");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x3, list(x2)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));
        facts.add(new IndependenceFact(x4, x5, list()));

        return new Ret("Solus Theorem 11, SMR !==> ESP (Figure 8)", facts);
    }

    public Ret getFigure12() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");
        Node x5 = new GraphNode("5");
        Node x6 = new GraphNode("6");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x3));
        facts.add(new IndependenceFact(x1, x5, list(x2, x3, x4)));
        facts.add(new IndependenceFact(x4, x6, list(x1, x2, x3, x5)));
        facts.add(new IndependenceFact(x1, x3, list(x2, x4, x5, x6)));

        return new Ret("Solus Theorem 12, TSP !==> Orientation Faithfulness (Figure 11)", facts);
    }

    public Ret getFigure7() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x2, list(x4)));
        facts.add(new IndependenceFact(x1, x3, list(x2)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));

        return new Ret("Solus Theorem 12, ESP !==> TSP (Figure 7)", facts);
    }

    public Ret getFigure6() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");
        Node x5 = new GraphNode("5");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x5, list(x2, x3)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));
        facts.add(new IndependenceFact(x3, x5, list(x1, x2, x4)));
        facts.add(new IndependenceFact(x1, x4, list(x2, x3, x5)));
        facts.add(new IndependenceFact(x1, x4, list(x2, x3)));

        return new Ret("Solus Theorem 11, TSP !==> Faithfulness counterexample (Figure 6)", facts);
    }

    public Ret getFacts6() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x2, list(x4)));
        facts.add(new IndependenceFact(x1, x3, list(x2)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));

        return new Ret("", facts);
    }

    public Ret getFactsRaskutti() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x1, x3, list(x2)));
        facts.add(new IndependenceFact(x2, x4, list(x1, x3)));
        facts.add(new IndependenceFact(x1, x2, list(x4)));

        return new Ret("Raskutti Theorem 2.4 SMR !==> Restricted Faithfulness", facts);
    }

    //@Test
    public void testFromIndependgetenceFactsa() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        List<Node> nodes = new ArrayList<>();
        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);

        Graph graph = new EdgeListGraph(nodes);
        graph.addDirectedEdge(x1, x2);
        graph.addDirectedEdge(x2, x3);
        graph.addDirectedEdge(x3, x4);
        graph.addDirectedEdge(x1, x4);

        GraphScore score = new GraphScore(graph);

        Boss boss = new Boss(score);
        boss.setCacheScores(false);
        System.out.println(boss.bestOrder(nodes));
    }

    //@Test
    public void testFromIndependenceFacts3() {
        Node x1 = new GraphNode("1");
        Node x2 = new GraphNode("2");
        Node x3 = new GraphNode("3");
        Node x4 = new GraphNode("4");

        List<Node> nodes = new ArrayList<>();
        nodes.add(x1);
        nodes.add(x2);
        nodes.add(x3);
        nodes.add(x4);

        Graph graph = new EdgeListGraph(nodes);
        graph.addDirectedEdge(x1, x2);
        graph.addDirectedEdge(x2, x3);
        graph.addDirectedEdge(x3, x4);
        graph.addDirectedEdge(x1, x4);

        LinearSemPm pm = new LinearSemPm(graph);
        LinearSemIm im = new LinearSemIm(pm);
        StandardizedLinearSemIm imsd = new StandardizedLinearSemIm(im, new Parameters());

        List<List<Node>> existingPaths = new ArrayList<>();

        boolean canceled = setPathsCanceling(x1, x4, imsd, existingPaths);

        if (!canceled) {
            System.out.println("Can't cancel those paths; either there's just one path or the paths overlap.");
        } else {
            System.out.println("Canceled paths from " + x1 + " to " + x4);
        }

        if (MatrixUtils.isPositiveDefinite(imsd.getImplCovar())) {
            System.out.println("Positive definite");
        }

        System.out.println(imsd.getImplCovar());

        DataSet data = imsd.simulateData(1000, false);

        Pc pc = new Pc(new IndTestFisherZ(data, 0.01));
        pc.setVerbose(true);
        System.out.println("PC " + pc.search());

        edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
        fges.setVerbose(true);
        System.out.println("FGES " + fges.search());

        Boss boss = new Boss(new edu.cmu.tetrad.search.SemBicScore(data));
        boss.setCacheScores(false);
        System.out.println("BOSS " + boss.bestOrder(data.getVariables()));
    }

    //@Test
    public void testFromIndependenceFacts4() {
        Graph graph = GraphUtils.randomGraph(10, 0, 20, 100, 100, 100, false);
        List<Node> nodes = graph.getNodes();

        LinearSemPm pm = new LinearSemPm(graph);
        LinearSemIm im = new LinearSemIm(pm);
        StandardizedLinearSemIm imsd = new StandardizedLinearSemIm(im, new Parameters());

        List<List<Node>> existingPaths = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) continue;
                boolean canceled = setPathsCanceling(nodes.get(i), nodes.get(j), imsd, existingPaths);

                if (canceled) {
                    System.out.println("Canceled " + nodes.get(i) + " -----> " + nodes.get(j));
                }
            }
        }

        System.out.println(imsd.getImplCovar());

        DataSet data = imsd.simulateData(1000, false);

        Pc pc = new Pc(new IndTestFisherZ(data, 0.01));
        pc.setVerbose(true);
        System.out.println("PC " + pc.search());

        edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
        fges.setVerbose(true);
        System.out.println("FGES " + fges.search());

        Boss boss = new Boss(new edu.cmu.tetrad.search.SemBicScore(data));
        boss.setCacheScores(false);
        System.out.println("BOSS " + boss.bestOrder(data.getVariables()));
    }

    private boolean setPathsCanceling(Node x1, Node x4, StandardizedLinearSemIm imsd, List<List<Node>> existingPaths) {
        SemGraph graph = imsd.getSemPm().getGraph();
        graph.setShowErrorTerms(false);

        List<List<Node>> paths = GraphUtils.allDirectedPathsFromTo(graph, x1, x4, -1);

        if (paths.size() < 2) return false;

        List<List<Node>> paths2 = new ArrayList<>();
        paths2.addAll(paths);
        paths2.addAll(existingPaths);

        for (int i = 0; i < paths2.size(); i++) {
            for (int j = i + 1; j < paths2.size(); j++) {
                List<Node> path1 = new ArrayList<>(paths2.get(i));
                List<Node> path2 = new ArrayList<>(paths2.get(j));
                path1.retainAll(path2);
                path1.remove(x1);
                path1.remove(x4);
                if (!path1.isEmpty()) return false;
            }
        }

        existingPaths.addAll(paths);

        List<Double> products = new ArrayList<>();

        for (List<Node> path : paths) {
            double p = 1.0;

            for (int i = 0; i < path.size() - 1; i++) {
                Node x = path.get(i);
                Node y = path.get(i + 1);
                p *= imsd.getEdgeCoef(x, y);
            }

            products.add(p);
        }

        double sum = 0;

        for (int i = 1; i < products.size(); i++) {
            sum += products.get(i);
        }

        double factor = 1;
        Set<NodePair> pairs = new HashSet<>();
        boolean changed = true;

        while (changed) {
            changed = false;

            for (int j = 1; j < paths.size(); j++) {
                List<Node> path = paths.get(j);
                double p = 1.0;

                for (int i = 0; i < path.size() - 1; i++) {
                    Node x = path.get(i);
                    Node y = path.get(i + 1);
                    if (pairs.contains(new NodePair(x, y))) continue;
                    boolean set = imsd.setEdgeCoefficient(x, y, -factor * imsd.getEdgeCoef(x, y) * (products.get(0)) / sum);
                    if (set) {
                        pairs.add(new NodePair(x, y));
                        changed = true;
                    }
                }

                products.add(p);
            }
        }

        return true;
    }

    private List<Node> list(Node... nodes) {
        List<Node> list = new ArrayList<>();
        Collections.addAll(list, nodes);
        return list;
    }

    //@Test
    public void testTeyssierIndices() {
        Graph graph = GraphUtils.randomGraph(10, 0, 10, 100, 100, 100, false);
        LinearSemPm pm = new LinearSemPm(graph);
        LinearSemIm im = new LinearSemIm(pm);
        DataSet dataSet = im.simulateData(1000, false);
        edu.cmu.tetrad.search.SemBicScore score = new edu.cmu.tetrad.search.SemBicScore(dataSet);
        TeyssierScorer scorer = new TeyssierScorer(score);

        List<Node> order = scorer.getOrder();
        scorer.score(order);

        Node v = order.get(5);

        for (int i = 0; i < 10; i++) {
            int r = RandomUtil.getInstance().nextInt(order.size());
            scorer.moveTo(v, r);
//            System.out.println(r + " " + scorer.indexOf(v));

            assertEquals(r, scorer.indexOf(v));
        }
    }

    @Test
    public void testBfci() {
        Parameters params = new Parameters();
        params.set(Params.SAMPLE_SIZE, 1000);
        params.set(Params.NUM_MEASURES, 25);
        params.set(Params.NUM_LATENTS, 8);
        params.set(Params.AVG_DEGREE, 6);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.COEF_LOW, 0);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VAR_LOW, 1);
        params.set(Params.VAR_HIGH, 3);
        params.set(Params.VERBOSE, true);

        params.set(Params.NUM_RUNS, 30);

        params.set(Params.BOSS_SCORE_TYPE, false);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 3);
        params.set(Params.USE_SCORE, true);

        params.set(Params.MAX_PATH_LENGTH, -1);
        params.set(Params.COMPLETE_RULE_SET_USED, true);

        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.ALPHA, 0.001);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new Fci(new FisherZ()));
        algorithms.add(new FciMax(new FisherZ()));
        algorithms.add(new Rfci(new FisherZ()));
        algorithms.add(new Gfci(new edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore(), new FisherZ()));
        algorithms.add(new BFCI(new edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore(), new FisherZ()));
//        algorithms.add(new BOSS(new SemBicScore(), new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
//        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
//        statistics.add(new ParameterColumn(Params.BOSS_METHOD));
//        statistics.add(new ParameterColumn(Params.BOSS_SCORE_TYPE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
//        statistics.add(new SHD_CPDAG());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
//        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_PAG);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss/bfci", simulations,
                algorithms, statistics, params);
    }

    private static class Ret {
        private final String label;
        private final IndependenceFacts facts;

        public Ret(String label, IndependenceFacts facts) {
            this.label = label;
            this.facts = facts;
        }

        public String getLabel() {
            return label;
        }

        public IndependenceFacts getFacts() {
            return facts;
        }
    }

    @Test
    public void testWayne() {
        List<Ret> allFacts = new ArrayList<>();
        allFacts.add(getFactsSimpleCanceling());
//        allFacts.add(getFactsRaskutti());
//        allFacts.add(getFigure6());
//        allFacts.add(getFigure7());
//        allFacts.add(getFigure8());
//        allFacts.add(getFigure12());

        int count = 0;

        boolean printCpdag = false;


        Boss.Method[] methods = {GSP, Boss.Method.BOSS, Boss.Method.SP};

        for (Ret facts : allFacts) {
            count++;

            TeyssierScorer scorer = new TeyssierScorer(new IndTestDSep(facts.getFacts()));

            OrderedMap<String, Set<Graph>> graphs = new ListOrderedMap<>();
            OrderedMap<String, Set<String>> labels = new ListOrderedMap<>();

            System.out.println();
            System.out.println(facts.getLabel());
            System.out.println(facts.getFacts());

            List<Node> variables = facts.facts.getVariables();

            PermutationGenerator gen = new PermutationGenerator(variables.size());
            int[] perm;

            while ((perm = gen.next()) != null) {
                List<Node> p = GraphUtils.asList(perm, variables);

                System.out.println("Initial permutation = " + p);

//                scorer.score(p);
//
//                System.out.println(scorer.getGraph(false));
//
                Boss2 boss = new Boss2(new IndTestDSep(facts.getFacts()));

                List<Node> order = boss.bestOrder(p);
                System.out.println(boss.getGraph(order, true));

            }

         }
    }
}





