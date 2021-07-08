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
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.sem.StandardizedSemIm;
import edu.cmu.tetrad.util.MatrixUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;

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
    public void testAllFacts() {
        List<Ret> allFacts = new ArrayList<>();
        allFacts.add(getFactsSimpleCanceling());
        allFacts.add(getFactsRaskutti());
        allFacts.add(getFigure6());
        allFacts.add(getFigure7());
        allFacts.add(getFigure8());
        allFacts.add(getFigure12());

        int count = 0;

        List<BestOrderScoreSearch.Method> bossMethods = new ArrayList<>();
        bossMethods.add(BestOrderScoreSearch.Method.SP);
        bossMethods.add(BestOrderScoreSearch.Method.PROMOTION);
        bossMethods.add(BestOrderScoreSearch.Method.ALL_INDICES);
//        bossMethods.add(BestOrderScoreSearch.Method.ESP);
        bossMethods.add(BestOrderScoreSearch.Method.GSP);

        for (Ret facts : allFacts) {
            count++;

            List<Node> order = facts.facts.getVariables();

            System.out.println("Model " + count + ": " + facts.label + "\n");
            System.out.println("\nFacts: \n\n" + facts.facts);

            int numRuns = 5;

            for (BestOrderScoreSearch.Method method : bossMethods) {
                System.out.println("Method = " + method + " " + facts.label + "\n");

                Set<Graph> graphs = new HashSet<>();

                for (int i = 0; i < numRuns; i++) {
                    shuffle(order);

                    BestOrderScoreSearch boss = new BestOrderScoreSearch(new IndTestDSep(facts.getFacts()));
                    boss.setCachingScores(true);
                    boss.setMethod(method);
                    boss.setNumStarts(10);
                    Graph pattern = SearchGraphUtils.patternForDag(boss.search(order));
                    graphs.add(pattern);
                }

                printGraphs(graphs);

            }

//            {
//                {
//                    System.out.println("Method = " + "FGES" + "\n");
//                    Set<Graph> graphs = new HashSet<>();
//
//                    for (int i = 0; i < numRuns; i++) {
//                        shuffle(order);
//                        facts.getFacts().setOrder(order);
//                        edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new GraphScore(facts.getFacts()));
//                        Graph pattern = fges.search();
//                        graphs.add(pattern);
//                    }
//
//                    printGraphs(graphs);
//                }
//
//                {
//                    System.out.println("Method = " + "PC" + "\n");
//                    Set<Graph> graphs = new HashSet<>();
//
//                    for (int i = 0; i < numRuns; i++) {
//                        shuffle(order);
//                        facts.getFacts().setOrder(order);
//                        edu.cmu.tetrad.search.PcAll cpc = new edu.cmu.tetrad.search.PcAll(new IndTestDSep(facts.getFacts()), null);
//                        cpc.setColliderDiscovery(edu.cmu.tetrad.search.PcAll.ColliderDiscovery.FAS_SEPSETS);
//                        Graph pattern = cpc.search();
//                        graphs.add(pattern);
//                    }
//
//                    printGraphs(graphs);
//                }
//
//                {
//                    System.out.println("Method = " + "CPC" + "\n");
//                    Set<Graph> graphs = new HashSet<>();
//
//                    for (int i = 0; i < numRuns; i++) {
//                        shuffle(order);
//                        facts.facts.setOrder(order);
//                        edu.cmu.tetrad.search.PcAll cpc = new edu.cmu.tetrad.search.PcAll(new IndTestDSep(facts.facts), null);
//                        cpc.setColliderDiscovery(edu.cmu.tetrad.search.PcAll.ColliderDiscovery.CONSERVATIVE);
//                        Graph pattern = cpc.search();
//                        graphs.add(pattern);
//                    }
//
//                    printGraphs(graphs);
//                }
//            }

            System.out.println(("\n\n--------------\n"));
        }
    }

    private void printGraphs(Set<Graph> graphs) {
        List<Graph> _graphs = new ArrayList<>(graphs);

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
//        facts.add(new IndependenceFact(x1, x4, list(x2))); // unfaithful.

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

        BestOrderScoreSearch boss = new BestOrderScoreSearch(score);
        boss.setCachingScores(false);
        System.out.println(boss.search(score.getVariables()));
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

        BestOrderScoreSearch boss = new BestOrderScoreSearch(new edu.cmu.tetrad.search.SemBicScore(data));
        boss.setCachingScores(false);
        System.out.println("BOSS " + boss.search(data.getVariables()));
    }

    @Test
    public void testFromIndependenceFacts4() {
        Graph graph = GraphUtils.randomGraph(10, 0, 20, 100, 100, 100, false);
        List<Node> nodes = graph.getNodes();

        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        StandardizedSemIm imsd = new StandardizedSemIm(im, new Parameters());

        List<List<Node>> existingPaths = new ArrayList<>();

//        LOOK:
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) continue;
                boolean canceled = setPathsCanceling(nodes.get(i), nodes.get(j), imsd, existingPaths);

                if (canceled) {
                    System.out.println("Canceled " + nodes.get(i) + " -----> " + nodes.get(j));
//                    break LOOK;
                }
            }
        }

//        if (MatrixUtils.isPositiveDefinite(imsd.getImplCovar())) {
//            System.out.println("Positive definite");
//        }

        System.out.println(imsd.getImplCovar());

        DataSet data = imsd.simulateData(1000, false);

        Pc pc = new Pc(new IndTestFisherZ(data, 0.01));
        pc.setVerbose(true);
        System.out.println("PC " + pc.search());

        edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
        fges.setVerbose(true);
        System.out.println("FGES " + fges.search());

        BestOrderScoreSearch boss = new BestOrderScoreSearch(new edu.cmu.tetrad.search.SemBicScore(data));
        boss.setCachingScores(false);
        System.out.println("BOSS " + boss.search(data.getVariables()));
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

//        for (int i = 0; i < paths.get(0).size() - 1; i++) {
//            Node x = paths.get(0).get(i);
//            Node y = paths.get(0).get(i + 1);
//            imsd.setEdgeCoefficientUnchecked(x, y, factor * imsd.getEdgeCoef(x, y));
//        }


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

//        for (int i = 0; i < paths.get(0).size() - 1; i++) {
//            Node x = paths.get(0).get(i);
//            Node y = paths.get(0).get(i + 1);
//            imsd.setEdgeCoefficientUnchecked(x, y, factor * imsd.getEdgeCoef(x, y));
//        }


//        for (List<Node> path : paths) {
//            for (int i = 0; i < path.size() - 1; i++) {
//                Node x = path.get(i);
//                Node y = path.get(i + 1);
//                if (pairs.contains(new NodePair(x, y))) continue;
//                imsd.setEdgeCoefficientUnchecked(x, y, factor * imsd.getEdgeCoef(x, y));
//                pairs.add(new NodePair(x, y));
//            }
//        }

        return true;
    }


    private List<Node> list(Node... nodes) {
        List<Node> list = new ArrayList<>();
        Collections.addAll(list, nodes);
        return list;
    }
}





