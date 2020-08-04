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
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    private int depth = 1000;

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
        Map<OrderedPair<Node>, Set<Double>> V = new HashMap<>();
        Map<OrderedPair<Node>, Double> P1 = new HashMap<>();

        int l = -1;

        while (degree(G1) - 1 >= l) {
            l = l + 1;

            Map<Node, List<Node>> a = new HashMap<>();

            for (Node X : nodes) {
                List<Node> adj = G1.getAdjacentNodes(X);
                a.put(X, adj);
            }

            for (Node x : nodes) {
                for (Node y : a.get(x)) {
                    List<Node> aa = new ArrayList<>(a.get(x));
                    aa.remove(y);

                    if (aa.size() < l) continue;

                    ChoiceGenerator gen = new ChoiceGenerator(aa.size(), l);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        List<Node> SS = GraphUtils.asList(choice, aa);

                        double p = pvalue(x, y, SS);

                        if (p <= alpha) {
                            addP(V, x, y, p);
                            addP(V, y, x, p);
                        } else {
                            G1.removeEdge(x, y);
                            setList(S, x, y, SS);
                            setList(S, y, x, SS);
                            clear(V, x, y);
                            clear(V, y, x);
                            break;
                        }
                    }
                }
            }
        }

        for (OrderedPair<Node> pair : new HashSet<>(V.keySet())) {
            if (!V.get(pair).isEmpty()) {
                setP(P1, pair, max(V.get(pair)));
            } else {
                P1.remove(pair);
            }
        }

        // algorithm 2

        Set<List<Node>> R0 = new HashSet<>();
        Set<Double> T = new HashSet<>();
        Map<OrderedPair<Node>, Set<Double>> Tp = new HashMap<>();
        Map<OrderedPair<Node>, Double> P2 = new HashMap<>();
        Set<List<Node>> amb = new HashSet<>();
        List<Triple> ut = getUT(G1);

        for (Triple triple : ut) {
            Node x = triple.getX();
            Node y = triple.getY();
            Node z = triple.getZ();

            if (G1.isAdjacentTo(x, z)) continue;

            if (!S.get(new OrderedPair<>(x, z)).contains(y)) {
                G1.setEndpoint(x, y, Endpoint.ARROW);
                G1.setEndpoint(z, y, Endpoint.ARROW);

                addRecord(R0, x, y, z);
                addRecord(R0, z, y, x);

                List<List<Node>> c = getC(x, z, y, G1);

                for (List<Node> cond : c) {
                    double p = pvalue(x, z, cond);
                    T.add(p);
                }

                addP(Tp, z, y, max(P1.get(new OrderedPair<>(x, y)), max(T)));
                addP(Tp, x, y, max(P1.get(new OrderedPair<>(z, y)), max(T)));
            }
        }

        // unorientation procedure for A2
        Graph G2 = new EdgeListGraph(G1);

        for (Edge edge : G1.getEdges()) {
            if (!Edges.isBidirectedEdge(edge)) continue;

            Node x = edge.getNode1();
            Node y = edge.getNode2();

            List<Node> intox = G1.getNodesInTo(x, Endpoint.ARROW);
            for (Node w : intox) {
                G2.removeEdge(x, w);
                G2.addUndirectedEdge(x, w);
                addRecord(amb, w, x);
                addRecord(amb, x, w);
            }

            List<Node> intoy = G1.getNodesInTo(y, Endpoint.ARROW);
            for (Node w : intoy) {
                G2.removeEdge(y, w);
                G2.addUndirectedEdge(y, w);
                addRecord(amb, w, y);
                addRecord(amb, y, w);
            }

            for (Edge edgexy : G1.getEdges()) {
                if (!Edges.isDirectedEdge(edgexy)) continue;

                if (edgexy.pointsTowards(y)) {
                    setP(P2, new OrderedPair<>(x, y), sum(Tp.get(new OrderedPair<>(x, y))));
                }
            }
        }

        // algorithm 3
        Set<List<Node>> R1 = new HashSet<>();
        Set<List<Node>> R2 = new HashSet<>();
        Set<List<Node>> R3 = new HashSet<>();

        Set<List<Node>> tri = getTri(G1);
        Set<List<Node>> kite = getKite(G1);

        boolean loop = true;

        while (loop) {
            loop = false;

            for (Triple triple : ut) {
                Node x = triple.getX();
                Node y = triple.getY();
                Node z = triple.getZ();

                if (G2.containsEdge(Edges.directedEdge(x, y))
                        && G2.getEndpoint(z, y) != Endpoint.ARROW
                        && !existsRecord(amb, y, z)
                        && !existsRecord(union(R0, R1), x, y, z)) {
                    G2.setEndpoint(y, z, Endpoint.ARROW);
                    addRecord(R1, x, y, z);
                    loop = true;
                }
            }

            for (List<Node> record : tri) {
                Node x = record.get(0);
                Node y = record.get(1);
                Node z = record.get(2);

                if (G2.containsEdge(Edges.directedEdge(x, y))
                        && G2.containsEdge(Edges.directedEdge(y, z))
                        && G2.getEndpoint(z, x) != Endpoint.ARROW
                        && !existsRecord(amb, y, z)
                        && !existsRecord(R2, y, x, z)) {
                    G2.setEndpoint(x, z, Endpoint.ARROW);
                    addRecord(R2, x, y, z);
                    loop = true;
                }
            }

            for (List<Node> record : kite) {
                Node x = record.get(0);
                Node y = record.get(1);
                Node z = record.get(2);
                Node w = record.get(3);

                if (G2.getEndpoint(y, z) != Endpoint.ARROW
                        && G2.containsEdge(Edges.undirectedEdge(x, y))
                        && G2.containsEdge(Edges.undirectedEdge(x, z))
                        && G2.containsEdge(Edges.directedEdge(y, w))
                        && G2.containsEdge(Edges.directedEdge(z, w))
                        && G2.getEndpoint(w, x) != Endpoint.ARROW
                        && G2.containsEdge(Edges.undirectedEdge(x, w))
                        && !existsRecord(amb, x, w)
                        && !existsRecord(R3, x, y, z, w)) {
                    G2.setEndpoint(y, z, Endpoint.ARROW);
                    addRecord(R3, x, y, z, w);
                    addRecord(R3, x, z, y, w);
                    loop = true;
                }
            }
        }

        for (List<Node> record : R3) {
            Node x = record.get(0);
            Node y = record.get(1);
            Node z = record.get(2);
            Node w = record.get(3);

            if (existsRecord(R2, x, y, w) || existsRecord(R2, x, z, w)) {
                R3.remove(record);
            }
        }

        // defining evidence of orientation

        Map<List<Node>, Set<Node>> e0 = new HashMap<>();
        Map<List<Node>, Set<Node>> e1 = new HashMap<>();
        Map<List<Node>, Set<Node>> e2 = new HashMap<>();
        Map<List<Node>, Set<Node>> e3 = new HashMap<>();
        Map<List<Node>, Set<Node>> eAll = new HashMap<>();

        return G2;
    }

    @Override
    public long getElapsedTime() {
        return 0;
    }

    private List<List<Node>> getC(Node x, Node y, Node z, Graph G) {
        List<List<Node>> c = new ArrayList<>();

        List<Node> adjx = G.getAdjacentNodes(x);

        DepthChoiceGenerator genx = new DepthChoiceGenerator(adjx.size(), adjx.size());
        int[] choicex;

        while ((choicex = genx.next()) != null) {
            List<Node> cond = GraphUtils.asList(choicex, adjx);
            if (cond.contains(z)) {
                c.add(cond);
            }
        }

        List<Node> adjy = G.getAdjacentNodes(y);

        DepthChoiceGenerator geny = new DepthChoiceGenerator(adjy.size(), adjy.size());
        int[] choicey;

        while ((choicey = geny.next()) != null) {
            List<Node> cond = GraphUtils.asList(choicey, adjy);
            if (cond.contains(z)) {
                c.add(cond);
            }
        }

        return c;
    }

    @SafeVarargs
    private final Set<List<Node>> union(Set<List<Node>>... r) {
        Set<List<Node>> union = new HashSet<>();

        for (Set<List<Node>> _r : r) {
            union.addAll(_r);
        }

        return union;
    }

    private List<Triple> getUT(Graph G1) {
        List<Triple> ut = new ArrayList<>();

        for (Node y : G1.getNodes()) {
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

        return ut;
    }

    private Set<List<Node>> getTri(Graph g1) {
        List<Node> nodes = g1.getNodes();
        Set<List<Node>> tri = new HashSet<>();

        for (Node x : nodes) {
            for (Node y : nodes) {
                for (Node z : nodes) {
                    if (g1.isAdjacentTo(y, z) && g1.isAdjacentTo(x, z) && g1.isAdjacentTo(y, z)) {
                        addRecord(tri, y, x, z);
                    }
                }
            }
        }

        return tri;
    }

    private Set<List<Node>> getKite(Graph g1) {
        List<Node> nodes = g1.getNodes();
        Set<List<Node>> tri = new HashSet<>();

        for (Node x : nodes) {
            for (Node y : nodes) {
                for (Node z : nodes) {
                    for (Node w : nodes) {
                        if (g1.isAdjacentTo(x, y)
                                && g1.isAdjacentTo(x, z)
                                && g1.isAdjacentTo(y, w)
                                && g1.isAdjacentTo(z, w)
                                && g1.isAdjacentTo(x, w)
                        ) {
                            addRecord(tri, x, y, z, w);
                        }
                    }
                }
            }
        }

        return tri;
    }

    private void addRecord(Set<List<Node>> R, Node... x) {
        List<Node> l = new ArrayList<>();
        addAll(l, x);
        R.add(l);
    }

    private boolean existsRecord(Set<List<Node>> R, Node... x) {
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
}





