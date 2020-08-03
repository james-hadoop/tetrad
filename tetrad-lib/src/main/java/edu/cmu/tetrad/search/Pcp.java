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

import java.util.*;

import static java.util.Collections.addAll;

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
    public Graph search() {
        if (getIndependenceTest() == null) {
            throw new NullPointerException();
        }

        List<Node> nodes = getIndependenceTest().getVariables();
        double alpha = independenceTest.getAlpha();

        // Algorithm 1

        Graph G1 = completeGraph(nodes);

        Map<OrderedPair<Node>, List<Node>> S = new HashMap<>();
        Map<OrderedPair<Node>, Set<Double>> Vp = new HashMap<>();
        Map<OrderedPair<Node>, Double> P1 = new HashMap<>();

        int l = -1;

        do {
            l = l + 1;

            Map<Node, List<Node>> a = new HashMap<>();

            for (Node X : nodes) {
                List<Node> adj = G1.getAdjacentNodes(X);
                a.put(X, adj);
            }

            for (Edge e : new HashSet<>(G1.getEdges())) {
                Node x = e.getNode1();
                Node y = e.getNode2();

                List<Node> aa = new ArrayList<>(a.get(x));
                aa.remove(y);

                if (aa.size() < l) continue;

                ChoiceGenerator gen = new ChoiceGenerator(aa.size(), l);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    List<Node> SS = GraphUtils.asList(choice, aa);

                    double p = pvalue(x, y, SS);

                    if (p <= alpha) {
                        addP(Vp, x, y, p);
                        addP(Vp, y, x, p);
                    } else {
                        G1.removeEdge(e);
                        setList(S, x, y, SS);
                        setList(S, y, x, SS);
                        clear(Vp, x, y);
                        clear(Vp, y, x);
                        break;
                    }
                }
            }
        } while (degree(G1) >= l + 2);

        for (OrderedPair<Node> pair : new HashSet<>(Vp.keySet())) {
            if (!Vp.get(pair).isEmpty()) {
                setP(P1, pair, max(Vp.get(pair)));
            } else {
                Vp.remove(pair);
            }
        }

        this.graph = G1;

        // algorithm 2

        Set<Triple> ut = new HashSet<>();

        for (Node y : nodes) {
            List<Node> adj = G1.getAdjacentNodes(y);

            if (adj.size() < 2) continue;

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                Node x = adj.get(choice[0]);
                Node z = adj.get(choice[1]);

                if (G1.isAdjacentTo(x, z)) continue;
                ut.add(new Triple(x, y, z));
            }
        }

        Map<OrderedPair<Node>, Double> P2 = new HashMap<>();
        List<List<Node>> R0 = new ArrayList<>();

        for (Triple triple : ut) {
            Node x = triple.getX();
            Node y = triple.getY();
            Node z = triple.getZ();

            List<Double> _V = new ArrayList<>();

            if (!S.get(new OrderedPair<>(x, z)).contains(y)) {
                graph.setEndpoint(x, y, Endpoint.ARROW);
                graph.setEndpoint(z, y, Endpoint.ARROW);

                addRecord(R0, x, y, z);
                addRecord(R0, z, y, x);

                List<Node> adjx = graph.getAdjacentNodes(x);

                DepthChoiceGenerator genx = new DepthChoiceGenerator(adjx.size(), adjx.size());
                int[] choicex;

                while ((choicex = genx.next()) != null) {
                    List<Node> cond = GraphUtils.asList(choicex, adjx);
                    double px = pvalue(x, z, cond);
                    _V.add(px);
                }

                List<Node> adjy = graph.getAdjacentNodes(x);

                DepthChoiceGenerator geny = new DepthChoiceGenerator(adjy.size(), adjy.size());
                int[] choicey;

                while ((choicey = geny.next()) != null) {
                    List<Node> cond = GraphUtils.asList(choicey, adjy);
                    double px = pvalue(x, z, cond);
                    _V.add(px);
                }

                addP(Vp, z, y, max(P1.get(new OrderedPair<>(x, y)), max(_V)));
                addP(Vp, x, y, max(P1.get(new OrderedPair<>(z, y)), max(_V)));
            }
        }

        // unorientation procedure for A2
        Graph G2 = new EdgeListGraph(G1);

        for (Edge edge : G1.getEdges()) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            List<Node> intox = G1.getNodesInTo(x, Endpoint.ARROW);
            for (Node w : intox) {
                if (G1.getEndpoint(x, w) == Endpoint.CIRCLE) {
                    G2.removeEdge(x, w);
                    G2.addUndirectedEdge(x, w);
                }
            }

            List<Node> intoy = G1.getNodesInTo(y, Endpoint.ARROW);
            for (Node w : intoy) {
                if (G1.getEndpoint(y, w) == Endpoint.CIRCLE) {
                    G2.removeEdge(y, w);
                    G2.addUndirectedEdge(y, w);
                }
            }

            for (Edge edgexy : graph.getEdges()) {
                if (!Edges.isDirectedEdge(edgexy)) continue;

                if (edgexy.pointsTowards(y)) {
                    setP(P2, new OrderedPair<>(x, y), sum(Vp.get(new OrderedPair<>(x, y))));
                }
            }
        }

        return G2;
    }

    private void addRecord(List<List<Node>> R, Node...x) {
        List<Node> l = new ArrayList<>();
        addAll(l, x);
        R.add(l);
    }

    private boolean existsRecord(List<List<Node>> R, Node...x) {
        List<Node> l = new ArrayList<>();
        addAll(l, x);
        return R.contains(l);
    }

    private Graph completeGraph(List<Node> nodes) {
        return GraphUtils.completeGraph(new EdgeListGraph(nodes));
    }

    private void clear(Map<OrderedPair<Node>, Set<Double>> v, Node x, Node y) {
        v.computeIfAbsent(new OrderedPair<>(x, y), k -> new HashSet<>());
        v.get(new OrderedPair<>(x, y)).clear();
    }

    private void setList(Map<OrderedPair<Node>, List<Node>> s, Node x, Node y, List<Node> SS) {
        s.put(new OrderedPair<>(x, y), SS);
    }

    private void addP(Map<OrderedPair<Node>, Set<Double>> v, Node x, Node y, double p) {
        v.computeIfAbsent(new OrderedPair<>(x, y), k -> new HashSet<>());
        v.get(new OrderedPair<>(x, y)).add(p);
    }

    private void setP(Map<OrderedPair<Node>, Double> V, OrderedPair<Node> pair, double p) {
        V.put(pair, p);
    }

    private double max(List<Double> p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double d : p) {
            if (d > max) max = d;
        }

        return max;
    }

    private double max(Set<Double> p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double d : p) {
            if (d > max) max = d;
        }

        return max;
    }

    private double max(double p1, double p2) {
        return Math.max(p1, p2);
    }

    private double sum(Set<Double> p) {
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

    private boolean isUndirectedFromTo(Node a, Node c, Graph graph) {
        Edge edge = graph.getEdge(a, c);
        return edge != null && Edges.isNondirectedEdge(edge);
    }

    private int degree(Graph graph) {
        int max = 0;

        for (Node x : graph.getNodes()) {
            int degree = graph.getDegree(x);

            if (degree > max) {
                max = degree;
            }
        }

        return max;
    }

    private double pvalue(Node x, Node y, List<Node> aa) {
        synchronized (independenceTest) {
            independenceTest.isIndependent(x, y, aa);
            return independenceTest.getPValue();
        }
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
}





