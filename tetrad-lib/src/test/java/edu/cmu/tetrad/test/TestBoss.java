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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.BOSS;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.EbicScore;
import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.sem.StandardizedSemIm;
import edu.cmu.tetrad.util.*;
import org.apache.commons.collections4.OrderedMap;
import org.apache.commons.collections4.map.ListOrderedMap;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static java.util.Collections.shuffle;
import static java.util.Collections.sort;

/**
 * Tests to make sure the DelimiterType enumeration hasn't been tampered with.
 *
 * @author Joseph Ramsey
 */
@SuppressWarnings("ALL")
public final class TestBoss {


    private static void runTestLoop(Graph g, List<Node> order, Score score, IndependenceTest test, boolean useTest) {
        System.out.println("order = " + order);

        //        g = new EdgeListGraph(g);
//        order = new ArrayList<>(order);

//        {
//            Gsp gsp;
//
//            if (true) {
//                gsp = new Gsp(test);
//            } else {
//                gsp = new Gsp(score);
//            }
//
//            gsp.setCachingScores(true);
//            gsp.setNumStarts(1);
//            gsp.setGspDepth(5);
//            gsp.setReturnCpdag(true);
//            Graph pattern = gsp.search(order);
//            pattern = GraphUtils.replaceNodes(pattern, variables);
//
//            printFailed(SearchGraphUtils.patternForDag(g), pattern, "GSP");
//        }

//        {
//            Boss boss = getBoss(score, test);
//
//            boss.setCacheScores(true);
//            boss.setMethod(Boss.Method.BOSS_PROMOTION);
//            boss.setNumStarts(1);
//            List<Node> perm = boss.bestOrder(order);
//            Graph dag = boss.getGraph(perm, false);
//
//            printFailed(g, dag, "BOSS Promotion " + order + " " + perm
//                    + " " + g
//                    + " " + dag);
//        }

        {
            Boss boss = getBoss(score, test);

            boss.setCacheScores(true);
            boss.setMethod(Boss.Method.BOSS_ALL_INDICES);
            boss.setNumStarts(3);
            boss.setVerbose(true);
            List<Node> perm = boss.bestOrder(order);
            Graph dag = boss.getGraph(perm, false);

            printFailed(g, dag, "BOSS All Indices " + order// + " " + perm
//                    + " \n" + g
                    + " \n" + dag);
        }

//        {
//            Boss boss = getBoss(score, test);
//
//            boss.setCacheScores(true);
//            boss.setMethod(Boss.Method.SP);
//            boss.setNumStarts(1);
//            List<Node> perm = boss.bestOrder(order);
//            Graph dag = boss.getGraph(perm, false);
//
//            printFailed(g, dag, "SP " + order + " " + perm
//                    + " " + g
//                    + " " + dag);
//        }

//        {
//            Boss boss;
//
//            if (true) {
//                boss = new Boss(test);
//            } else {
//                boss = new Boss(score);
//            }
//
//            boss.setCacheScores(true);
//            boss.setMethod(Boss.Method.SP);
//            boss.setNumStarts(1);
//            List<Node> perm = boss.bestOrder(order);
//            Graph pattern = boss.getGraph(perm, true);
//            printFailed(SearchGraphUtils.patternForDag(g), pattern, "SP");
//        }

    }

    private static void runTestLoop2(Graph g, List<Node> order, Score score, IndependenceTest test, boolean useTest) {

        {
            Boss boss = getBoss(score, test);

            boss.setCacheScores(true);
            boss.setMethod(Boss.Method.SP);
            boss.setNumStarts(1);
            List<Node> perm = boss.bestOrder(order);
            Graph dag = boss.getGraph(perm, false);

//            if (dag.getNumEdges() != g.getNumEdges()) {
//                System.out.println("order = " + order + " # edges = " + dag.getNumEdges());
//            }
        }
    }

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

    @Test
    public void testBoss() {
        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 10);
        params.set(Params.AVG_DEGREE, 0, 1, 2, 4, 5, 6, 7, 8);
        params.set(Params.SAMPLE_SIZE, 100000);
        params.set(Params.NUM_RUNS, 20);
        params.set(Params.RANDOMIZE_COLUMNS, false);
        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.COEF_LOW, 0.25);
        params.set(Params.COEF_HIGH, 1.0);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.BOSS_METHOD, 1);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.DEPTH));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveData(false);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss", simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss2() {
        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 20);
        params.set(Params.AVG_DEGREE, 4);
        params.set(Params.SAMPLE_SIZE, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000,
                50000, 100000, 200000);
        params.set(Params.NUM_RUNS, 30);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0.2);
        params.set(Params.COEF_HIGH, 0.8);
        params.set(Params.VERBOSE, false);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
        algorithms.add(new Fges(new SemBicScore()));
        algorithms.add(new PcAll(new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveData(false);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss", simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss3() {
//        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 100);
        params.set(Params.AVG_DEGREE, 10);
        params.set(Params.SAMPLE_SIZE, 1000);
        params.set(Params.NUM_RUNS, 1);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VERBOSE, false);
        params.set(Params.COLLIDER_DISCOVERY_RULE, 2, 3);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.CACHE_SCORES, false);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new EbicScore()));
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new PcAll(new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveData(false);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss", simulations, algorithms, statistics, params);
    }

    private void swap(int vv, int ww, double[] scores, LinkedList<Node> br, int prefixSize, K3 k3) {
        if (vv < 0) throw new IllegalArgumentException("vv < 0");
        if (vv >= prefixSize) throw new IllegalArgumentException("vv > prefixSize");
        if (!(ww >= vv)) throw new IllegalArgumentException("w < prefixSize");
        if (ww != vv + 1) throw new IllegalArgumentException("ww is not directly subsequent to vv");

        Node v = br.get(vv);

        br.remove(v);
        br.add(ww, v);

        for (int r = 0; r < prefixSize; r++) {
            if (r == vv || r == ww) {
                double score1 = k3.getScore(new K3.Subproblem(br.get(r), new HashSet<>(getPrefix(br, r))));
                scores[r] = score1;
            }
        }
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    @Test
    public void testAllFacts() {
        List<Ret> allFacts = new ArrayList<>();
        allFacts.add(getFactsSimpleCanceling());
        allFacts.add(getFactsRaskutti());
        allFacts.add(getFigure6());
        allFacts.add(getFigure7());
        allFacts.add(getFigure8());
        allFacts.add(getFigure12());

        int count = 0;

        boolean printPattern = false;

        List<Boss.Method> bossMethods = new ArrayList<>();
        bossMethods.add(Boss.Method.SP);
        bossMethods.add(Boss.Method.BOSS_PROMOTION);
        bossMethods.add(Boss.Method.BOSS_ALL_INDICES);
//        bossMethods.add(BestOrderScoreSearch.Method.ESP);
//        bossMethods.add(BestOrderScoreSearch.Method.GSP);

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

                {
                    Gsp boss = new Gsp(test);
                    boss.setCachingScores(true);
                    boss.setNumStarts(1);
                    boss.setGspDepth(5);
                    Graph pattern = SearchGraphUtils.patternForDag(boss.search(order));

                    if (graphs.get("GSP") == null) {
                        graphs.put("GSP", new HashSet<>());
                        labels.put("GSP", new HashSet<>());
                    }

                    graphs.get("GSP").add(pattern);
                    labels.get("GSP").add(facts.getLabel());

                    if (printPattern) {
                        System.out.println(pattern);
                    }
                }

                {
                    Boss.Method method = Boss.Method.BOSS_PROMOTION;

                    Boss boss = new Boss(test);
                    boss.setCacheScores(true);
                    boss.setMethod(method);
                    boss.setNumStarts(1);

                    List<Node> perm = boss.bestOrder(test.getVariables());
                    Graph pattern = boss.getGraph(perm, false);

                    if (graphs.get(method.toString()) == null) {
                        graphs.put(method.toString(), new HashSet<>());
                        labels.put(method.toString(), new HashSet<>());
                    }

                    graphs.get(method.toString()).add(pattern);
                    labels.get(method.toString()).add(facts.getLabel());

                    if (printPattern) {
                        System.out.println(pattern);
                    }
                }

                {
                    Boss.Method method = Boss.Method.BOSS_ALL_INDICES;

                    Boss boss = new Boss(test);
                    boss.setCacheScores(true);
                    boss.setMethod(method);
                    boss.setNumStarts(1);
                    List<Node> perm = boss.bestOrder(test.getVariables());
                    Graph pattern = boss.getGraph(perm, true);

                    if (graphs.get(method.toString()) == null) {
                        graphs.put(method.toString(), new HashSet<>());
                        labels.put(method.toString(), new HashSet<>());
                    }

                    graphs.get(method.toString()).add(pattern);
                    labels.get(method.toString()).add(facts.getLabel());

                    if (printPattern) {
                        System.out.println(pattern);
                    }
                }

                {
                    Boss.Method method = Boss.Method.SP;

                    Boss boss = new Boss(test);
                    boss.setCacheScores(true);
                    boss.setMethod(method);
                    boss.setNumStarts(1);
                    List<Node> perm = boss.bestOrder(test.getVariables());
                    Graph pattern = boss.getGraph(perm, true);

                    if (graphs.get(method.toString()) == null) {
                        graphs.put(method.toString(), new HashSet<>());
                        labels.put(method.toString(), new HashSet<>());
                    }

                    graphs.get(method.toString()).add(pattern);
                    labels.get(method.toString()).add(facts.getLabel());

                    if (printPattern) {
                        System.out.println(pattern);
                    }
                }
            }

            printGraphs("GSP", graphs);
            printGraphs("BOSS_PROMOTION", graphs);
            printGraphs("BOSS_ALL_INDICES", graphs);
            printGraphs("SP", graphs);
        }
    }

    @Test
    public void testRaskutti() {
        Ret facts = getFactsRaskutti();

        IndTestDSep test = new IndTestDSep(facts.getFacts());

        Boss boss = new Boss(test);
        boss.setCacheScores(true);
        boss.setMethod(Boss.Method.BOSS_ALL_INDICES);
        boss.setNumStarts(1);

        List<Node> variables = test.getVariables();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);

            List<Node> p2 = boss.bestOrder(test.getVariables());
            Graph pattern = boss.getGraph(p2, true);

            System.out.println(p + " " + pattern.getNumEdges());
        }
    }

    @Test
    public void testFromData() {
        for (int i = 0; i < 10; i++) {
            System.out.println("\nRun " + (i + 1));

            Graph g = GraphUtils.randomGraph(15, 0, 30, 100,
                    100, 100, false);
            SemPm pm = new SemPm(g);
            Parameters params = new Parameters();
            params.set(Params.COEF_LOW, 0.1);
            params.set(Params.COEF_HIGH, 1);
            SemIm im = new SemIm(pm, params);
            DataSet data = im.simulateData(100000, false);

            data = DataUtils.shuffleColumns(data);
            List<Node> order = data.getVariables();

            edu.cmu.tetrad.search.EbicScore score = new edu.cmu.tetrad.search.EbicScore(data);
//            score.setPenaltyDiscount(2);
            edu.cmu.tetrad.search.IndependenceTest test = new edu.cmu.tetrad.search.IndTestFisherZ(data, 0.01);

            runTestLoop(g, order, score, test, false);
        }
    }

    @Test
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

    @Test
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

    private boolean isPatternForDag(Graph pattern, Graph dag) {
        if (!GraphUtils.undirectedGraph(pattern).equals(GraphUtils.undirectedGraph(dag))) {
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

    @Test
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

        return new Ret("Simple 4-node path canceling model that GES should get right", facts);
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

    @Test
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

    @Test
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

        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        StandardizedSemIm imsd = new StandardizedSemIm(im, new Parameters());

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

    @Test
    public void testFromIndependenceFacts4() {
        Graph graph = GraphUtils.randomGraph(10, 0, 20, 100, 100, 100, false);
        List<Node> nodes = graph.getNodes();

        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        StandardizedSemIm imsd = new StandardizedSemIm(im, new Parameters());

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

    private boolean setPathsCanceling(Node x1, Node x4, StandardizedSemIm imsd, List<List<Node>> existingPaths) {
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
}





