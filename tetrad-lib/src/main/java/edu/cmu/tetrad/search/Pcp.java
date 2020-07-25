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

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;

/**
 * Implements the PCP algorithm. The original idea fior this was due to Eric Strobl and significantly revised by '
 * Wayne Lam and Peter Spirtes.
 *
 * @author Joseph Ramsey.
 */
public class Pcp implements GraphSearch {

    /**
     * The independence test used for the PC search.g
     */
    private final IndependenceTest independenceTest;

    /**
     * Sepset information accumulated in the search.
     */
    private SepsetMap sepsets;

    /**
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    private int depth = 1000;

    /**
     * The graph that's constructed during the search.
     */
    private Graph graph;

    /**
     * Elapsed time of the most recent search.
     */
    private long elapsedTime;

    /**
     * True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
     * graphs).
     */
    private boolean aggressivelyPreventCycles = false;

    // P-values from adjacency search
    private Map<NodePair, List<Double>> p1 = null;//new HashMap<>();

    // P-values from collider orientation
    private final Map<NodePair, Double> p2 = new HashMap<>();

    // Number of p-values used.
    private Map<NodePair, Set<Object>> I = new HashMap<>();

    // Set of ambiguous edges. Not used yet.
    private final Set<NodePair> Amb = new HashSet<>();

    // Set of nodes.
    private List<Node> nodes;

    // A p-value associated with each edge.
    private final Map<NodePair, Double> pp = new HashMap<>();

    //=============================CONSTRUCTORS==========================//

    /**
     * Constructs a new PC search using the given independence test as oracle.
     *
     * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
     *                         independence test, for fear of duplicating the data set!
     */
    public Pcp(IndependenceTest independenceTest) {
        if (independenceTest == null) {
            throw new NullPointerException();
        }

        this.independenceTest = independenceTest;
    }

    //==============================PUBLIC METHODS========================//

    /**
     * Runs PC starting with a complete graph over all nodes of the given conditional independence test, using the given
     * independence test and knowledge and returns the resultant graph. The returned graph will be a pattern if the
     * independence information is consistent with the hypothesis that there are no latent common causes. It may,
     * however, contain cycles or bidirected edges if this assumption is not born out, either due to the actual presence
     * of latent common causes, or due to statistical errors in conditional independence judgments.
     */
    @Override
    public Graph search() {
        return search(independenceTest.getVariables());
    }

    /**
     * Runs PC starting with a commplete graph over the given list of nodes, using the given independence test and
     * knowledge and returns the resultant graph. The returned graph will be a pattern if the independence information
     * is consistent with the hypothesis that there are no latent common causes. It may, however, contain cycles or
     * bidirected edges if this assumption is not born out, either due to the actual presence of latent common causes,
     * or due to statistical errors in conditional independence judgments.
     * <p>
     * All of the given nodes must be in the domatein of the given conditional independence test.
     */
    public Graph search(List<Node> nodes) {
        FasPcp fas = new FasPcp(getIndependenceTest());
        return search(fas, nodes);
    }

    public Graph search(FasPcp fas, List<Node> nodes) {
        this.nodes = nodes;

        long startTime = System.currentTimeMillis();

        if (getIndependenceTest() == null) {
            throw new NullPointerException();
        }

        List<Node> allNodes = getIndependenceTest().getVariables();

        if (!allNodes.containsAll(nodes)) {
            throw new IllegalArgumentException("All of the given nodes must " +
                    "be in the domain of the independence test provided.");
        }

        // Fast adjacency search, recording some stats. Algorithm 1.

        fas.setDepth(getDepth());

        graph = fas.algorithm1();
        sepsets = fas.getSepsets();

        this.p1 = fas.getP1();
        I = fas.getI();

        // Algorithm 2, orient colliders using sepsets, recording statistics.
        // The arguments to this are class fields.
        algorithm2(this.sepsets);

        // Algorithm 3
        // The arguments to this are class fields.
        algorithm3();

        // Make all circles into tails finally.
        for (Edge edge : graph.getEdges()) {
            Edge edge2 = new Edge(edge);

            if (edge2.getEndpoint1() == Endpoint.CIRCLE) {
                edge2.setEndpoint1(Endpoint.TAIL);
            }

            if (edge2.getEndpoint2() == Endpoint.CIRCLE) {
                edge2.setEndpoint2(Endpoint.TAIL);
            }

            graph.removeEdge(edge);
            graph.addEdge(edge2);
        }

        this.elapsedTime = System.currentTimeMillis() - startTime;

        return graph;
    }


    /**
     * @return true iff edges will not be added if they would create cycles.
     */
    public boolean isAggressivelyPreventCycles() {
        return this.aggressivelyPreventCycles;
    }

    /**
     * @param aggressivelyPreventCycles Set to true just in case edges will not be addeds if they would create cycles.
     */
    public void setAggressivelyPreventCycles(boolean aggressivelyPreventCycles) {
        this.aggressivelyPreventCycles = aggressivelyPreventCycles;
    }

    /**
     * @return the independence test being used in the search.
     */
    public IndependenceTest getIndependenceTest() {
        return independenceTest;
    }

    /**
     * @return the sepset map from the most recent search. Non-null after the first call to <code>search()</code>.
     */
    public SepsetMap getSepsets() {
        return this.sepsets;
    }

    /**
     * @return the current depth of search--that is, the maximum number of conditioning nodes for any conditional
     * independence checked.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
     * checked.
     *
     * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
     *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
     *              machines.
     */
    public void setDepth(int depth) {
        if (depth < -1) {
            throw new IllegalArgumentException("Depth must be -1 or >= 0: " + depth);
        }

        if (depth > 1000) {
            throw new IllegalArgumentException("Depth must be <= 1000.");
        }

        this.depth = depth;
    }

    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    public long getElapsedTime() {
        return elapsedTime;
    }

    public Set<Edge> getAdjacencies() {
        return new HashSet<>(graph.getEdges());
    }

    public Set<Edge> getNonadjacencies() {
        Graph complete = GraphUtils.completeGraph(graph);
        Set<Edge> nonAdjacencies = complete.getEdges();
        Graph undirected = GraphUtils.undirectedGraph(graph);
        nonAdjacencies.removeAll(undirected.getEdges());
        return new HashSet<>(nonAdjacencies);
    }

    //===============================PRIVATE METHODS=======================//

    public List<Node> getNodes() {
        return graph.getNodes();
    }

    /**
     * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
     * case y is in Sepset({x, z}).
     */
    public void algorithm2(SepsetMap set) {

        ArrayList<Object> r0 = new ArrayList<>();

        TetradLogger.getInstance().log("details", "Starting Collider Orientation:");

        List<Node> nodes = graph.getNodes();

        List<Double> Ppp = new ArrayList<>();

        for (Node y : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(y);

            if (adjacentNodes.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node x = adjacentNodes.get(combination[0]);
                Node z = adjacentNodes.get(combination[1]);

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(x, z)) {
                    continue;
                }

                if (r0.contains(new Triple(x, y, z))) continue;

                List<Node> sepset = set.get(x, z);

                if (!sepset.contains(y)) {
//                    graph.removeEdge(x, y);
//                    graph.removeEdge(z, y);

                    graph.setEndpoint(x, y, Endpoint.ARROW);
                    graph.setEndpoint(z, y, Endpoint.ARROW);

//                    graph.addDirectedEdge(x, y);
//                    graph.addDirectedEdge(z, y);

                    r0.add(new Triple(x, y, z));

                    TetradLogger.getInstance().log("colliderOrientations", SearchLogUtils.colliderOrientedMsg(x, y, z, sepset));

                    List<Node> adj = graph.getAdjacentNodes(x);

                    for (int i = 0; i <= depth; i++) {
                        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), depth);
                        int[] choice;

                        while ((choice = gen.next()) != null) {
                            List<Node> c = GraphUtils.asList(choice, adj);

                            // We synchronize so that the p-value of the result of the test is necessarily returned.
                            synchronized (independenceTest) {
                                independenceTest.isIndependent(x, z, c);
                                Ppp.add(independenceTest.getPValue());
                            }
                        }
                    }

                    adj = graph.getAdjacentNodes(z);

                    for (int i = 0; i <= depth; i++) {
                        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), depth);
                        int[] choice;

                        while ((choice = gen.next()) != null) {
                            List<Node> c = GraphUtils.asList(choice, adj);

                            // We synchronize so that the p-value of the result of the test is necessarily returned.
                            synchronized (independenceTest) {
                                independenceTest.isIndependent(x, z, c);
                                Ppp.add(independenceTest.getPValue());
                            }
                        }
                    }

                    // Now add all p-values to Ppp that would have been considered by CPC
                    // Node "NodePair" is an unordered pairs of nodes, {X, Y}.
                    pp.put(new NodePair(z, y), max(p1.get(new NodePair(x, y)), Ppp));
                    pp.put(new NodePair(x, y), max(p1.get(new NodePair(z, y)), Ppp));

                    Ppp.clear();
                }
            }
        }

        for (Edge edge : graph.getEdges()) {
            final Node x = edge.getNode1();
            final Node y = edge.getNode2();

            if (Edges.isBidirectedEdge(edge)) {
                graph.removeEdge(edge);
                graph.addUndirectedEdge(x, y);
                Amb.add(new NodePair(x, y));

                for (Node parent : graph.getParents(x)) {
                    graph.removeEdge(x, parent);
                    graph.addUndirectedEdge(x, parent);
                    Amb.add(new NodePair(x, parent));
                }

                for (Node parent : graph.getParents(y)) {
                    graph.removeEdge(parent, y);
                    graph.addUndirectedEdge(parent, y);
                    Amb.add(new NodePair(parent, y));
                }
            }
        }

        // Note "Amb" is not used yet.

        for (Edge edge : graph.getEdges()) {
            if (edge.isDirected()) {
                Node x = Edges.getDirectedEdgeHead(edge);
                Node y = Edges.getDirectedEdgeTail(edge);

                p2.put(new NodePair(x, y), sum(Collections.singletonList(pp.get(new NodePair(x, y)))));
                insertI(new NodePair(x, y), new Object());
            }
        }
    }


    private void algorithm3() {
        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node node : nodes) {
                changed = changed || meekR1(node, graph);
                changed = changed || meekR2(node, graph);
                changed = changed || meekR3(node, graph);
            }
        }
    }

//    private double fdr(double alpha, int m, int k) {
//        double m1 = 0;
//        for (int i = 1; i <= m; i++) m1 += 1.0 / i;
//        return m * alpha * m1 / (double) k;
//    }

    private void insertI(NodePair key, Object o) {
        I.get(key).add(o);
    }

    public Map<NodePair, Set<Object>> getI() {
        return I;
    }

    private double max(List<Double> p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double d : p) {
            if (d > max) max = d;
        }

        return max;
    }

    private double max(List<Double> p1, List<Double> p2) {
        List<Double> p3 = new ArrayList<>(p1);
        p3.addAll(p2);
        return max(p3);
    }

    private double sum(List<Double> p) {
        double sum = 0;

        for (double d : p) {
            sum += d;
        }

        return sum;
    }

    /**
     * @return the string in nodelist which matches string in BK.
     */
    public static Node translate(String a, List<Node> nodes) {
        for (Node node : nodes) {
            if ((node.getName()).equals(a)) {
                return node;
            }
        }

        return null;
    }

    /**
     * Meek's rule R1: if a-->b, b---c, and a not adj to c, then a-->c
     */
    private boolean meekR1(Node b, Graph graph) {
        boolean changed = false;

        List<Node> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            return false;
        }

        ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
        int[] choice;

        while ((choice = cg.next()) != null) {
            List<Node> nodes = GraphUtils.asList(choice, adjacentNodes);
            Node a = nodes.get(0);
            Node c = nodes.get(1);

            changed = changed || r1Helper(a, b, c, graph);
            changed = changed || r1Helper(c, b, a, graph);
        }

        return changed;
    }

    private boolean r1Helper(Node a, Node b, Node c, Graph graph) {
        boolean changed = false;

        if (!graph.isAdjacentTo(a, c) && graph.isDirectedFromTo(a, b)
                && graph.isAdjacentTo(b, c) && doesntPointToward(b, c, graph)) {
            if (isNotshieldedNoncollider(a, b, c, graph)) {
                return false;
            }

            for (Node d : graph.getParents(b)) {
                increment(b, d, p2.get(new NodePair(d, b)));
            }

            changed = direct(b, c, graph);
        }

        return changed;
    }

    private void increment(Node x, Node y, double pp) {
        this.pp.putIfAbsent(new NodePair(x, y), 0.0);
        this.pp.put(new NodePair(x, y), pp);
    }

    /**
     * If a-->b-->c, a--c, then b-->c.
     */
    private boolean meekR2(Node c, Graph graph) {
        boolean changed = false;

        List<Node> adjacentNodes = graph.getAdjacentNodes(c);

        if (adjacentNodes.size() < 2) {
            return false;
        }

        ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
        int[] choice;

        while ((choice = cg.next()) != null) {
            List<Node> nodes = GraphUtils.asList(choice, adjacentNodes);
            Node a = nodes.get(0);
            Node b = nodes.get(1);

            changed = changed || r2Helper(a, b, c, graph);
            changed = changed || r2Helper(b, a, c, graph);
            changed = changed || r2Helper(a, c, b, graph);
            changed = changed || r2Helper(c, a, b, graph);
        }

        return changed;
    }

    private boolean r2Helper(Node a, Node b, Node c, Graph graph) {
        boolean changed = false;

        if (graph.isDirectedFromTo(a, b) &&
                graph.isDirectedFromTo(b, c) &&
                graph.isAdjacentTo(a, c) && doesntPointToward(a, c, graph)) {

            for (Node d : graph.getParents(c)) {
                if (graph.isDirectedFromTo(b, d)) {
                    double p1 = p2.get(new NodePair(a, d));
                    double p2 = this.p2.get(new NodePair(d, c));

                    increment(a, c, Math.max(p1, p2));
                }
            }

            changed = direct(a, c, graph);
        }

        return changed;
    }

    private boolean isUndirectedFromTo(Node a, Node c, Graph graph) {
        Edge edge = graph.getEdge(a, c);
        return edge != null && Edges.isNondirectedEdge(edge);
    }

    /**
     * Meek's rule R3. If a--b, a--c, a--d, c-->b, d-->b, then orient a-->b.
     */
    private boolean meekR3(Node a, Graph graph) {
        boolean changed = false;

        List<Node> adjacentNodes = graph.getAdjacentNodes(a);

        if (adjacentNodes.size() < 3) {
            return false;
        }

        for (Node d : adjacentNodes) {
            if (Edges.isUndirectedEdge(graph.getEdge(a, d))) {
                List<Node> otherAdjacents = new ArrayList<>(adjacentNodes);
                otherAdjacents.remove(d);

                ChoiceGenerator cg = new ChoiceGenerator(otherAdjacents.size(), 2);
                int[] choice;

                while ((choice = cg.next()) != null) {
                    List<Node> nodes = GraphUtils.asList(choice, otherAdjacents);
                    Node b = nodes.get(0);
                    Node c = nodes.get(1);

                    boolean isKite = isKite(a, d, b, c, graph);

                    if (isKite) {
                        if (isNotshieldedNoncollider(c, d, b, graph)) {
                            continue;
                        }

                        for (Node e : graph.getParents(a)) {
                            if (graph.isAdjacentTo(b, e)) {
                                double p1 = max(this.p1.get(new NodePair(d, e)));
                                double p2 = this.p2.get(new NodePair(e, a));

                                increment(a, c, Math.max(p1, p2));
                            }
                        }

                        changed = changed || direct(d, a, graph);
                    }
                }
            }
        }

        return changed;
    }

    private boolean isKite(Node a, Node d, Node b, Node c, Graph graph) {
        boolean b4 = isUndirectedFromTo(d, c, graph);
        boolean b5 = isUndirectedFromTo(d, b, graph);
        boolean b6 = graph.isDirectedFromTo(b, a);
        boolean b7 = graph.isDirectedFromTo(c, a);
        boolean b8 = graph.isAdjacentTo(d, a) && doesntPointToward(d, a, graph);

        return b4 && b5 && b6 && b7 && b8;
    }

    private boolean doesntPointToward(Node a, Node d, Graph graph) {
        final Edge edge = graph.getEdge(a, d);
        return edge.getProximalEndpoint(d) != Endpoint.ARROW;
    }

    private boolean direct(Node a, Node b, Graph graph) {
        Edge before = graph.getEdge(a, b);

        Edge after = new Edge(before);

        after = new Edge(a, b, after.getProximalEndpoint(a), Endpoint.ARROW);

        graph.removeEdge(before);
        graph.addEdge(after);

        return true;
    }

    private static boolean isNotshieldedNoncollider(Node a, Node b, Node c,
                                                    Graph graph) {
        if (!graph.isAdjacentTo(a, b)) {
            return true;
        }

        if (!graph.isAdjacentTo(c, b)) {
            return true;
        }

        if (graph.isAdjacentTo(a, c)) {
            return true;
        }

        if (graph.isAmbiguousTriple(a, b, c)) {
            return true;
        }

        return graph.getEndpoint(a, b) == Endpoint.ARROW &&
                graph.getEndpoint(c, b) == Endpoint.ARROW;

    }
}





