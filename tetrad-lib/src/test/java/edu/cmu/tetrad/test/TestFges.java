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

import edu.cmu.tetrad.bayes.BayesIm;
import edu.cmu.tetrad.bayes.BayesPm;
import edu.cmu.tetrad.bayes.MlBayesIm;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.LargeScaleSimulation;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TextTable;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static junit.framework.TestCase.assertFalse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Joseph Ramsey
 */
public class TestFges {


    private PrintStream out = System.out;
//    private OutputStream out =

    //    @Test
    public void explore1() {
        RandomUtil.getInstance().setSeed(1450184147770L);

        int numVars = 10;
        double edgesPerNode = 1.0;
        int numCases = 1000;
        double penaltyDiscount = 2.0;

        final int numEdges = (int) (numVars * edgesPerNode);

        List<Node> vars = new ArrayList<>();

        for (int i = 0; i < numVars; i++) {
            vars.add(new ContinuousVariable("X" + i));
        }

        Graph dag = GraphUtils.randomGraphRandomForwardEdges(vars, 0, numEdges, 30, 15, 15, false, true);
//        printDegreeDistribution(dag, System.out);

        int[] causalOrdering = new int[vars.size()];

        for (int i = 0; i < vars.size(); i++) {
            causalOrdering[i] = i;
        }

        LargeScaleSimulation simulator = new LargeScaleSimulation(dag, vars, causalOrdering);
        simulator.setOut(out);
        DataSet data = simulator.simulateDataFisher(numCases);

//        ICovarianceMatrix cov = new CovarianceMatrix(data);
        ICovarianceMatrix cov = new CovarianceMatrix(data);
        SemBicScore score = new SemBicScore(cov);
        score.setPenaltyDiscount(penaltyDiscount);

        Fges fges = new Fges(score);
        fges.setVerbose(false);
        fges.setOut(out);
        fges.setFaithfulnessAssumed(true);
//        fges.setMaxIndegree(1);
        fges.setCycleBound(5);

        Graph estPattern = fges.search();

//        printDegreeDistribution(estPattern, out);

        final Graph truePattern = SearchGraphUtils.patternForDag(dag);

        int[][] counts = SearchGraphUtils.graphComparison(estPattern, truePattern, null);

        int[][] expectedCounts = {
                {2, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 8, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
        };

        for (int i = 0; i < counts.length; i++) {
            assertTrue(Arrays.equals(counts[i], expectedCounts[i]));
        }
//

//        System.out.println(MatrixUtils.toString(expectedCounts));
//        System.out.println(MatrixUtils.toString(counts));

    }

    @Test
    public void explore2() {
        RandomUtil.getInstance().setSeed(1457220623122L);

        int numVars = 20;
        double edgeFactor = 1.0;
        int numCases = 1000;
        double structurePrior = 1;
        double samplePrior = 1;

        List<Node> vars = new ArrayList<>();

        for (int i = 0; i < numVars; i++) {
            vars.add(new ContinuousVariable("X" + i));
        }

        Graph dag = GraphUtils.randomGraphRandomForwardEdges(vars, 0, (int) (numVars * edgeFactor),
                30, 15, 15, false, true);
//        printDegreeDistribution(dag, out);

        BayesPm pm = new BayesPm(dag, 2, 3);
        BayesIm im = new MlBayesIm(pm, MlBayesIm.RANDOM);
        DataSet data = im.simulateData(numCases, false);

//        out.println("Finishing simulation");

        BDeScore score = new BDeScore(data);
        score.setSamplePrior(samplePrior);
        score.setStructurePrior(structurePrior);

        Fges ges = new Fges(score);
        ges.setVerbose(false);
        ges.setFaithfulnessAssumed(false);

        Graph estPattern = ges.search();

        final Graph truePattern = SearchGraphUtils.patternForDag(dag);

        int[][] counts = SearchGraphUtils.graphComparison(estPattern, truePattern, null);

        int[][] expectedCounts = {
                {2, 0, 0, 0, 0, 1},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {2, 0, 0, 13, 0, 3},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
        };

//        for (int i = 0; i < counts.length; i++) {
//            assertTrue(Arrays.equals(counts[i], expectedCounts[i]));
//        }

//        System.out.println(MatrixUtils.toString(expectedCounts));
//        System.out.println(MatrixUtils.toString(counts));
//        System.out.println(RandomUtil.getInstance().getSeed());
    }


    @Test
    public void testExplore3() {
        Graph graph = GraphConverter.convert("A-->B,A-->C,B-->D,C-->D");
        Fges fges = new Fges(new GraphScore(graph));
        Graph pattern = fges.search();
        assertEquals(SearchGraphUtils.patternForDag(graph), pattern);
    }

    @Test
    public void testExplore4() {
        Graph graph = GraphConverter.convert("A-->B,A-->C,A-->D,B-->E,C-->E,D-->E");
        Fges fges = new Fges(new GraphScore(graph));
        Graph pattern = fges.search();
        assertEquals(SearchGraphUtils.patternForDag(graph), pattern);
    }

    @Test
    public void testExplore5() {
        Graph graph = GraphConverter.convert("A-->B,A-->C,A-->D,A->E,B-->F,C-->F,D-->F,E-->F");
        Fges fges = new Fges(new GraphScore(graph));
        fges.setFaithfulnessAssumed(false);
        Graph pattern = fges.search();
        assertEquals(SearchGraphUtils.patternForDag(graph), pattern);
    }


    @Test
    public void testFromGraphSimpleFges() {

        // This may fail if faithfulness is assumed but should pass if not.

        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");

        Graph g = new EdgeListGraph();
        g.addNode(x1);
        g.addNode(x2);
        g.addNode(x3);
        g.addNode(x4);

        g.addDirectedEdge(x1, x2);
        g.addDirectedEdge(x1, x3);
        g.addDirectedEdge(x4, x2);
        g.addDirectedEdge(x4, x3);

        Graph pattern1 = new Pc(new IndTestDSep(g)).search();
        Fges fges = new Fges(new GraphScore(g));
        fges.setFaithfulnessAssumed(true);
        Graph pattern2 = fges.search();

//        System.out.println(pattern1);
//        System.out.println(pattern2);

        assertEquals(pattern1, pattern2);
    }

    @Test
    public void testFromGraphSimpleFgesMb() {

        // This may fail if faithfulness is assumed but should pass if not.

        Node x1 = new GraphNode("X1");
        Node x2 = new GraphNode("X2");
        Node x3 = new GraphNode("X3");
        Node x4 = new GraphNode("X4");

        Graph dag = new EdgeListGraph();
        dag.addNode(x1);
        dag.addNode(x2);
        dag.addNode(x3);
        dag.addNode(x4);

        dag.addDirectedEdge(x1, x2);
        dag.addDirectedEdge(x1, x3);
        dag.addDirectedEdge(x4, x2);
        dag.addDirectedEdge(x4, x3);

        GraphScore fgesScore = new GraphScore(dag);

        Fges fges = new Fges(fgesScore);
        Graph pattern1 = fges.search();

        Set<Node> mb = new HashSet<>();
        mb.add(x1);

        mb.addAll(pattern1.getAdjacentNodes(x1));

        for (Node child : pattern1.getChildren(x1)) {
            mb.addAll(pattern1.getParents(child));
        }

        Graph mb1 = pattern1.subgraph(new ArrayList<>(mb));

        FgesMb fgesMb = new FgesMb(fgesScore);
        Graph mb2 = fgesMb.search(x1);

        assertEquals(mb1, mb2);
    }

    @Test
    public void testFgesMbFromGraph() {
        RandomUtil.getInstance().setSeed(1450184147770L);

        int numNodes = 20;
        int numIterations = 1;

        for (int i = 0; i < numIterations; i++) {
//            System.out.println("Iteration " + (i + 1));
            Graph dag = GraphUtils.randomDag(numNodes, 0, numNodes, 10, 10, 10, false);
            GraphScore fgesScore = new GraphScore(dag);

            Fges fges = new Fges(fgesScore);
            Graph pattern1 = fges.search();

            Node x1 = fgesScore.getVariable("X1");

            Set<Node> mb = new HashSet<>();
            mb.add(x1);

            mb.addAll(pattern1.getAdjacentNodes(x1));

            for (Node child : pattern1.getChildren(x1)) {
                mb.addAll(pattern1.getParents(child));
            }

            Graph mb1 = pattern1.subgraph(new ArrayList<>(mb));

            FgesMb fgesMb = new FgesMb(fgesScore);
            Graph mb2 = fgesMb.search(x1);

            assertEquals(mb1, mb2);
        }
    }


    private void printDegreeDistribution(Graph dag, PrintStream out) {
        int max = 0;

        for (Node node : dag.getNodes()) {
            int degree = dag.getAdjacentNodes(node).size();
            if (degree > max) max = degree;
        }

        int[] counts = new int[max + 1];
        Map<Integer, List<Node>> names = new HashMap<>();

        for (int i = 0; i <= max; i++) {
            names.put(i, new ArrayList<Node>());
        }

        for (Node node : dag.getNodes()) {
            int degree = dag.getAdjacentNodes(node).size();
            counts[degree]++;
            names.get(degree).add(node);
        }

        for (int k = 0; k < counts.length; k++) {
            if (counts[k] == 0) continue;

            out.print(k + " " + counts[k]);

            for (Node node : names.get(k)) {
                out.print(" " + node.getName());
            }

            out.println();
        }
    }

    private boolean ancestral(Node x, Node y, Graph graph) {
        return graph.isAncestorOf(x, y) || graph.isAncestorOf(y, x);
    }

    /**
     * Runs the PC algorithm on the graph X1 --> X2, X1 --> X3, X2 --> X4, X3 --> X4. Should produce X1 -- X2, X1 -- X3,
     * X2 --> X4, X3 --> X4.
     */
    @Test
    public void testSearch1() {
        checkSearch("X1-->X2,X1-->X3,X2-->X4,X3-->X4",
                "X1---X2,X1---X3,X2-->X4,X3-->X4");
    }

    /**
     * Runs the PC algorithm on the graph X1 --> X2, X1 --> X3, X2 --> X4, X3 --> X4. Should produce X1 -- X2, X1 -- X3,
     * X2 --> X4, X3 --> X4.
     */
    @Test
    public void testSearch2() {
        checkSearch("X1-->X2,X1-->X3,X2-->X4,X3-->X4",
                "X1---X2,X1---X3,X2-->X4,X3-->X4");
    }

    /**
     * This will fail if the orientation loop doesn't continue after the first orientation.
     */
    @Test
    public void testSearch3() {
        checkSearch("A-->D,A-->B,B-->D,C-->D,D-->E",
                "A-->D,A---B,B-->D,C-->D,D-->E");
    }

    /**
     * This will fail if the orientation loop doesn't continue after the first orientation.
     */
    @Test
    public void testSearch4() {
        IKnowledge knowledge = new Knowledge2();
        knowledge.setForbidden("B", "D");
        knowledge.setForbidden("D", "B");
        knowledge.setForbidden("C", "B");

        checkWithKnowledge("A-->B,C-->B,B-->D", /*"A---B,B-->C,D",*/"A---B,B-->C,A---C,D---A",
                knowledge);
    }

    @Test
    public void testSearch5() {
        IKnowledge knowledge = new Knowledge2();
        knowledge.setTier(1, Collections.singletonList("A"));
        knowledge.setTier(2, Collections.singletonList("B"));

        checkWithKnowledge("A-->B", "A-->B", knowledge);
    }

    @Test
    public void testCites() {
        String citesString = "164\n" +
                "ABILITY\tGPQ\tPREPROD\tQFJ\tSEX\tCITES\tPUBS\n" +
                "1.0\n" +
                ".62\t1.0\n" +
                ".25\t.09\t1.0\n" +
                ".16\t.28\t.07\t1.0\n" +
                "-.10\t.00\t.03\t.10\t1.0\n" +
                ".29\t.25\t.34\t.37\t.13\t1.0\n" +
                ".18\t.15\t.19\t.41\t.43\t.55\t1.0";

        char[] citesChars = citesString.toCharArray();
        DataReader reader = new DataReader();
        ICovarianceMatrix cov = reader.parseCovariance(citesChars);

        IKnowledge knowledge = new Knowledge2();

        knowledge.addToTier(1, "ABILITY");
        knowledge.addToTier(2, "GPQ");
        knowledge.addToTier(3, "QFJ");
        knowledge.addToTier(3, "PREPROD");
        knowledge.addToTier(4, "SEX");
        knowledge.addToTier(5, "PUBS");
        knowledge.addToTier(6, "CITES");

        SemBicScore score = new SemBicScore(cov);
        score.setPenaltyDiscount(1);
        Fges fges = new Fges(score);
        fges.setKnowledge(knowledge);

        fges.setVerbose(true);

        Graph pattern = fges.search();

        System.out.println(pattern);

        String trueString = "Graph Nodes:\n" +
                "ABILITY;GPQ;PREPROD;QFJ;SEX;CITES;PUBS\n" +
                "\n" +
                "Graph Edges:\n" +
                "1. ABILITY --> GPQ\n" +
                "2. ABILITY --> PREPROD\n" +
                "3. ABILITY --> PUBS\n" +
                "4. GPQ --> QFJ\n" +
                "5. PREPROD --> CITES\n" +
                "6. PUBS --> CITES\n" +
                "7. QFJ --> CITES\n" +
                "8. QFJ --> PUBS\n" +
                "9. SEX --> PUBS";



        Graph trueGraph = null;


        try {
            trueGraph = GraphUtils.readerToGraphTxt(trueString);
        } catch (IOException e) {
            e.printStackTrace();
        }

        pattern = GraphUtils.replaceNodes(pattern, trueGraph.getNodes());

        assertEquals(trueGraph, pattern);
    }

    /**
     * Presents the input graph to FCI and checks to make sure the output of FCI is equivalent to the given output
     * graph.
     */
    private void checkSearch(String inputGraph, String outputGraph) {

        // Set up graph and node objects.
        Graph graph = GraphConverter.convert(inputGraph);

        // Set up search.
        Fges fges = new Fges(new GraphScore(graph));

        // Run search
        Graph resultGraph = fges.search();

        // Build comparison graph.
        Graph trueGraph = GraphConverter.convert(outputGraph);

        // PrintUtil out problem and graphs.
//        System.out.println("\nInput graph:");
//        System.out.println(graph);
//        System.out.println("\nResult graph:");
//        System.out.println(resultGraph);
//        System.out.println("\nTrue graph:");
//        System.out.println(trueGraph);

        resultGraph = GraphUtils.replaceNodes(resultGraph, trueGraph.getNodes());

        // Do test.
        assertTrue(resultGraph.equals(trueGraph));
    }

    /**
     * Presents the input graph to FCI and checks to make sure the output of FCI is equivalent to the given output
     * graph.
     */
    private void checkWithKnowledge(String inputGraph, String answerGraph,
                                    IKnowledge knowledge) {
        // Set up graph and node objects.
        Graph input = GraphConverter.convert(inputGraph);

        // Set up search.
        Fges fges = new Fges(new GraphScore(input));

        // Set up search.
        fges.setKnowledge(knowledge);

        // Run search
        Graph result = fges.search();

        // Build comparison graph.
        Graph answer = GraphConverter.convert(answerGraph);
//        Graph answer = new PC(new IndTestDSep(input)).search();

//        System.out.println("Input = " + input);
//        System.out.println("Knowledge = " + knowledge);
//        System.out.println("Answer = " + answer);
//        System.out.println("Result graph = " + result);

        // Do test.
        assertEquals(answer, result);
    }

    @Test
    public void testPcStable2() {
        RandomUtil.getInstance().setSeed(1450030184196L);
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < 10; i++) {
            nodes.add(new ContinuousVariable("X" + (i + 1)));
        }

        Graph graph = GraphUtils.randomGraph(nodes, 0, 10, 30, 15, 15, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet data = im.simulateData(200, false);

        TetradLogger.getInstance().setForceLog(false);
        IndependenceTest test = new IndTestFisherZ(data, 0.05);

        PcStableMax pc = new PcStableMax(test);
        pc.setVerbose(false);
        Graph pattern = pc.search();

        for (int i = 0; i < 1; i++) {
            DataSet data2 = DataUtils.reorderColumns(data);
            IndependenceTest test2 = new IndTestFisherZ(data2, 0.05);
            PcStableMax pc2 = new PcStableMax(test2);
            pc2.setVerbose(false);
            Graph pattern2 = pc2.search();
            assertTrue(pattern.equals(pattern2));
        }
    }


    @Test
    public void testFromGraph() {
        int numNodes = 10;
        int numIterations = 1;

        for (int i = 0; i < numIterations; i++) {
//            System.out.println("Iteration " + (i + 1));
            Graph dag = GraphUtils.randomDag(numNodes, 0, 2 * numNodes, 10, 10, 10, false);
            Fges fges = new Fges(new GraphScore(dag));
            fges.setFaithfulnessAssumed(false);
            Graph pattern1 = fges.search();
            Graph pattern2 = new Pc(new IndTestDSep(dag)).search();
//            System.out.println(pattern2);
            assertEquals(pattern2, pattern1);
        }
    }


    @Test
    public void testFromGraphWithForbiddenKnowledge() {
        int numNodes = 20;
        int numIterations = 20;

        for (int i = 0; i < numIterations; i++) {
            System.out.println("Iteration " + (i + 1));
            Graph dag = GraphUtils.randomDag(numNodes, 0, numNodes, 10, 10, 10, false);
            Graph knowledgeGraph = GraphUtils.randomDag(numNodes, 0, numNodes, 10, 10, 10, false);
            knowledgeGraph = GraphUtils.replaceNodes(knowledgeGraph, dag.getNodes());

            IKnowledge knowledge = forbiddenKnowledge(knowledgeGraph);

            Fges fges = new Fges(new GraphScore(dag));
            fges.setFaithfulnessAssumed(true);
            fges.setKnowledge(knowledge);
            Graph pattern1 = fges.search();

            for (Edge edge : knowledgeGraph.getEdges()) {
                Node x = Edges.getDirectedEdgeTail(edge);
                Node y = Edges.getDirectedEdgeHead(edge);

                if (pattern1.isParentOf(x, y)) {
                    System.out.println("Knowledge violated: " + edge + " x = " + x + " y = " + y);
                }

                assertFalse(pattern1.isParentOf(x, y));
            }
        }
    }

    @Test
    public void testFromGraphWithRequiredKnowledge() {
        int numNodes = 20;
        int numIterations = 20;

        for (int i = 0; i < numIterations; i++) {
            System.out.println("Iteration " + (i + 1));
            Graph dag = GraphUtils.randomDag(numNodes, 0, numNodes, 10, 10, 10, false);
            Graph knowledgeGraph = GraphUtils.randomDag(numNodes, 0, numNodes, 10, 10, 10, false);
            knowledgeGraph = GraphUtils.replaceNodes(knowledgeGraph, dag.getNodes());

            IKnowledge knowledge = requiredKnowledge(knowledgeGraph);

            Fges fges = new Fges(new GraphScore(dag));
            fges.setFaithfulnessAssumed(true);
            fges.setKnowledge(knowledge);
            Graph pattern1 = fges.search();

            for (Edge edge : knowledgeGraph.getEdges()) {
                Node x = Edges.getDirectedEdgeTail(edge);
                Node y = Edges.getDirectedEdgeHead(edge);

                if (!pattern1.isParentOf(x, y)) {
                    System.out.println("Knowledge violated: " + edge + " x = " + x + " y = " + y);
                }

                assertTrue (pattern1.isParentOf(x, y));
            }
        }
    }


    private IKnowledge forbiddenKnowledge(Graph graph) {
        IKnowledge knowledge = new Knowledge2(graph.getNodeNames());

        List<Node> nodes = graph.getNodes();

        for (Edge edge : graph.getEdges()) {
            Node n1 = Edges.getDirectedEdgeTail(edge);
            Node n2 = Edges.getDirectedEdgeHead(edge);

            if (n1.getName().startsWith("E_") || n2.getName().startsWith("E_")) {
                continue;
            }

            knowledge.setForbidden(n1.getName(), n2.getName());
        }

        return knowledge;
    }

    private IKnowledge requiredKnowledge(Graph graph) {

        IKnowledge knowledge = new Knowledge2(graph.getNodeNames());

        for (Edge edge : graph.getEdges()) {
            Node n1 = Edges.getDirectedEdgeTail(edge);
            Node n2 = Edges.getDirectedEdgeHead(edge);

            if (n1.getName().startsWith("E_") || n2.getName().startsWith("E_")) {
                continue;
            }

            knowledge.setRequired(n1.getName(), n2.getName());
        }

        return knowledge;
    }

    private Graph getSubgraph(Graph graph, boolean discrete1, boolean discrete2, DataSet dataSet) {
        Graph newGraph = new EdgeListGraph(graph.getNodes());

        for (Edge edge : graph.getEdges()) {
            Node node1 = dataSet.getVariable(edge.getNode1().getName());
            Node node2 = dataSet.getVariable(edge.getNode2().getName());

            if (discrete1 && node1 instanceof DiscreteVariable) {
                if (discrete2 && node2 instanceof DiscreteVariable) {
                    newGraph.addEdge(edge);
                }
            } else if (!discrete1 && node1 instanceof ContinuousVariable) {
                if (!discrete2 && node2 instanceof ContinuousVariable) {
                    newGraph.addEdge(edge);
                }
            } else if ((discrete1 && !discrete2) || (!discrete1 && discrete2)) {
                if (node1 instanceof ContinuousVariable && node2 instanceof DiscreteVariable) {
                    newGraph.addEdge(edge);
                } else if (node1 instanceof DiscreteVariable && node2 instanceof ContinuousVariable) {
                    newGraph.addEdge(edge);
                }
            }
        }

        return newGraph;
    }

    //    @Test
    public void testAjData() {
        double penalty = 4;

        try {

            for (int i = 0; i < 50; i++) {
                File dataPath = new File("/Users/jdramsey/Documents/LAB_NOTEBOOK.2012.04.20/2016.05.25/" +
                        "Simulated_data_for_Madelyn/simulation/data/DAG_" + i + "_data.txt");
                DataReader reader = new DataReader();
                DataSet Dk = reader.parseTabular(dataPath);

                File graphPath = new File("/Users/jdramsey/Documents/LAB_NOTEBOOK.2012.04.20/2016.05.25/" +
                        "Simulated_data_for_Madelyn/simulation/networks/DAG_" + i + "_graph.txt");

                Graph dag = GraphUtils.loadGraphTxt(graphPath);

                long start = System.currentTimeMillis();

//            Graph pattern = searchSemFges(Dk);
//            Graph pattern = searchBdeuFges(Dk, k);
                Graph pattern = searchMixedFges(Dk, penalty);

                long stop = System.currentTimeMillis();

                long elapsed = stop - start;
                long elapsedSeconds = elapsed / 1000;

                Graph truePattern = SearchGraphUtils.patternForDag(dag);

                GraphUtils.GraphComparison comparison = SearchGraphUtils.getGraphComparison3(pattern, truePattern, System.out);
                NumberFormat nf = new DecimalFormat("0.00");

                System.out.println(i +
                        "\t" + nf.format(comparison.getAdjPrec()) +
                        "\t" + nf.format(comparison.getAdjRec()) +
                        "\t" + nf.format(comparison.getAhdPrec()) +
                        "\t" + nf.format(comparison.getAhdRec()) +
                        "\t" + elapsedSeconds);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Graph searchSemFges(DataSet Dk, double penalty) {
        Dk = DataUtils.convertNumericalDiscreteToContinuous(Dk);
        SemBicScore score = new SemBicScore(new CovarianceMatrix(Dk));
        score.setPenaltyDiscount(penalty);
        Fges fges = new Fges(score);
        return fges.search();
    }

    private Graph searchBdeuFges(DataSet Dk, int k) {
        Discretizer discretizer = new Discretizer(Dk);
        List<Node> nodes = Dk.getVariables();

        for (Node node : nodes) {
            if (node instanceof ContinuousVariable) {
                discretizer.equalIntervals(node, k);
            }
        }

        Dk = discretizer.discretize();

        BDeuScore score = new BDeuScore(Dk);
        score.setSamplePrior(1.0);
        score.setStructurePrior(1.0);
        Fges fges = new Fges(score);
        return fges.search();
    }

    private Graph searchMixedFges(DataSet dk, double penalty) {
        MixedBicScore score = new MixedBicScore(dk);
        score.setPenaltyDiscount(penalty);
        Fges fges = new Fges(score);
        return fges.search();
    }

    private String getHeader(int u) {
        String header;

        switch (u) {
            case 0:
                header = "All edges";
                break;
            case 1:
                header = "Discrete-discrete";
                break;
            case 2:
                header = "Discrete-continuous";
                break;
            case 3:
                header = "Continuous-continuous";
                break;
            default:
                throw new IllegalStateException();
        }
        return header;
    }

    private void printBestStats(double[][][][] allAllRet, String[] algorithms, String[] statLabels,
                                int maxCategories, double ofInterestCutoff) {
        TextTable table = new TextTable(allAllRet.length + 1, allAllRet[0][0].length + 1);


        class Pair {
            private String algorithm;
            private double stat;

            public Pair(String algorithm, double stat) {
                this.algorithm = algorithm;
                this.stat = stat;
            }

            public String getAlgorithm() {
                return algorithm;
            }

            public double getStat() {
                return stat;
            }
        }


        System.out.println();
        System.out.println("And the winners are... !");

        for (int u = 0; u < 4; u++) {
            for (int numCategories = 2; numCategories <= maxCategories; numCategories++) {

                table.setToken(numCategories - 1, 0, numCategories + "");

                for (int statIndex = 0; statIndex < allAllRet[numCategories - 2][0].length; statIndex++) {
//                double maxStat = Double.NaN;
                    String maxAlg = "-";

                    List<Pair> algStats = new ArrayList<>();

                    for (int t = 0; t < algorithms.length; t++) {
                        double stat = allAllRet[numCategories - 2][t][statIndex][u];
                        if (!Double.isNaN(stat)) {
                            algStats.add(new Pair(algorithms[t], stat));
                        }
                    }

                    if (algStats.isEmpty()) {
                        maxAlg = "-";
                    } else {
                        Collections.sort(algStats, new Comparator<Pair>() {

                            @Override
                            public int compare(Pair o1, Pair o2) {
                                return -Double.compare(o1.getStat(), o2.getStat());
                            }
                        });

                        double maxStat = algStats.get(0).getStat();
                        maxAlg = algStats.get(0).getAlgorithm();

                        double minStat = algStats.get(algStats.size() - 1).getStat();

                        double diff = maxStat - minStat;
                        double ofInterest = maxStat - ofInterestCutoff * (diff);

                        for (int i = 1; i < algStats.size(); i++) {
                            if (algStats.get(i).getStat() >= ofInterest) {
                                maxAlg += "," + algStats.get(i).getAlgorithm();
                            }
                        }
                    }

                    table.setToken(numCategories - 1, statIndex + 1, maxAlg);
                }
            }

            for (int j = 0; j < statLabels.length; j++) {
                table.setToken(0, j + 1, statLabels[j]);
            }

            System.out.println();
            System.out.println(getHeader(u));
            System.out.println();

            System.out.println(table.toString());
        }


        NumberFormat nf = new DecimalFormat("0.00");

        System.out.println();
        System.out.println("Details:");
        System.out.println();
        System.out.println("Average statistics");

        for (int u = 0; u < 4; u++) {
            System.out.println();
            System.out.println(getHeader(u));
            System.out.println();

            for (int numCategories = 2; numCategories <= maxCategories; numCategories++) {
                System.out.println("\n# categories = " + numCategories);

                for (int t = 0; t < algorithms.length; t++) {
                    String algorithm = algorithms[t];

                    System.out.println("\nAlgorithm = " + algorithm);
                    System.out.println();

                    for (int statIndex = 0; statIndex < allAllRet[numCategories - 2][0].length; statIndex++) {
                        String statLabel = statLabels[statIndex];
                        double stat = allAllRet[numCategories - 2][t][statIndex][u];
                        System.out.println("\tAverage" + statLabel + " = " + nf.format(stat));
                    }
                }
            }
        }

    }

    @Test
    public void test7() {
        for (int i = 0; i < 10; i++) {

            Graph graph = GraphUtils.randomGraph(10, 0,
                    10, 10, 10, 10, false);
            SemPm semPm = new SemPm(graph);
            SemIm semIm = new SemIm(semPm);
            DataSet dataSet = semIm.simulateData(1000, false);

            Fges fges = new Fges(new SemBicScore(new CovarianceMatrix(dataSet)));
            Graph pattern = fges.search();

            Graph dag = dagFromPattern(pattern);

            assertFalse(dag.existsDirectedCycle());
        }
    }

    private Graph dagFromPattern(Graph pattern) {
        Graph dag = new EdgeListGraph(pattern);

        MeekRules rules = new MeekRules();

        WHILE:
        while (true) {
            List<Edge> edges = new ArrayList<>(dag.getEdges());

            for (Edge edge : edges) {
                if (Edges.isUndirectedEdge(edge)) {
                    Node x = edge.getNode1();
                    Node y = edge.getNode2();

                    List<Node> okx = dag.getAdjacentNodes(x);
                    okx.removeAll(dag.getChildren(x));
                    okx.remove(y);

                    List<Node> oky = dag.getAdjacentNodes(y);
                    oky.removeAll(dag.getChildren(y));
                    oky.remove(x);

                    if (!okx.isEmpty()) {
                        Node other = okx.get(0);
                        dag.removeEdge(other, x);
                        dag.removeEdge(y, x);
                        dag.addDirectedEdge(other, x);
                        dag.addDirectedEdge(y, x);
                    } else if (!oky.isEmpty()) {
                        Node other = oky.get(0);
                        dag.removeEdge(other, y);
                        dag.removeEdge(x, y);
                        dag.addDirectedEdge(other, y);
                        dag.addDirectedEdge(x, y);
                    } else {
                        dag.removeEdge(x, y);
                        dag.addDirectedEdge(x, y);
                    }

                    rules.orientImplied(dag);
                    continue WHILE;
                }
            }

            break;
        }

        return dag;
    }

    @Test
    public void testSemBicDiffs() {
        final int N = 1000;
        int numCond = 3;

        Graph graph = GraphUtils.randomGraph(10,0, 20, 100,
                100, 100, false);
        final List<Node> nodes = graph.getNodes();
        buildIndexing(nodes);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet dataSet = im.simulateData(N, false);
        SemBicScore score = new SemBicScore(dataSet);

        IndTestDSep dsep = new IndTestDSep(graph);
        int count = 1;

        for (int i = 0; i < 10000; i++) {
            Collections.shuffle(nodes);

            Node x = nodes.get(0);
            Node y = nodes.get(1);
            Set<Node> z = new HashSet<>();

            for (int c = 3; c <= 2 + numCond; c++) {
                z.add(nodes.get(c));
            }

            final boolean _dsep = dsep.isIndependent(x, y, new ArrayList<>(z));
            final double diff = scoreGraphChange(x, y, z, hashIndices, score) ;
            final boolean diffNegative = diff < 0;

            if (!_dsep && _dsep != diffNegative) {
                System.out.println(count++ + "\t" + (_dsep ? "dsep" : "dconn") + "\t" + (diffNegative ? "indep" : "dep") + "\tdiff = " + diff);
            }
        }

    }

    private double scoreGraphChange(Node x, Node y, Set<Node> parents,
                                    Map<Node, Integer> hashIndices, SemBicScore score) {
        int yIndex = hashIndices.get(y);

        if (x == y) {
            throw new IllegalArgumentException();
        }
        if (parents.contains(y)) {
            throw new IllegalArgumentException();
        }

        if (parents.contains(x)) {
            throw new IllegalArgumentException();
        }

        int[] parentIndices = new int[parents.size()];

        int count = 0;
        for (Node parent : parents) {
            parentIndices[count++] = hashIndices.get(parent);
        }

        return score.localScoreDiff(hashIndices.get(x), yIndex, parentIndices);
    }

    private HashMap<Node, Integer> hashIndices;

    // Maps adj to their indices for quick lookup.
    private void buildIndexing(List<Node> nodes) {
        this.hashIndices = new HashMap<>();

        int i = -1;

        for (Node n : nodes) {
            this.hashIndices.put(n, ++i);
        }
    }

    public static void main(String... args) {
    }

}




