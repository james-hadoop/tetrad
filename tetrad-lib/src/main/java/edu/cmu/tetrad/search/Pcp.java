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
 * Implements the PCP algorithm. The original idea for this was due to Eric Strobl and significantly revised by '
 * Wayne Lam and Peter Spirtes.
 *
 * @author Wayne Lam
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

        Map<List<Node>, List<Node>> S = new HashMap<>();
        Map<List<Node>, Set<Double>> V = new HashMap<>();
        Map<List<Node>, Double> P1 = new HashMap<>();

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

        for (List<Node> pair : new HashSet<>(V.keySet())) {
            if (!V.get(pair).isEmpty()) {
                setP(P1, pair, max(V.get(pair)));
            } else {
                P1.remove(pair);
            }
        }

        // algorithm 2

        Set<List<Node>> R0 = new HashSet<>();
        Set<Double> T = new HashSet<>();
        Map<List<Node>, Set<Double>> Tp = new HashMap<>();
        Map<List<Node>, Double> P2 = new HashMap<>();
        Set<List<Node>> amb = new HashSet<>();
        Set<List<Node>> ut = getUT(G1);

        for (List<Node> triple : ut) {
            Node x = triple.get(0);
            Node y = triple.get(1);
            Node z = triple.get(2);

            if (G1.isAdjacentTo(x, z)) continue;

            System.out.println("S = " + S);
            System.out.println("x = " + x + " z = " + z + " S.get(list(x, z)) = " + S.get(list(x, z)));

            if (!S.get(list(x, z)).contains(y)) {
                G1.setEndpoint(x, y, Endpoint.ARROW);
                G1.setEndpoint(z, y, Endpoint.ARROW);

                addRecord(R0, x, y, z);
                addRecord(R0, z, y, x);

                List<List<Node>> c = getC(x, z, y, G1);

                for (List<Node> cond : c) {
                    double p = pvalue(x, z, cond);
                    T.add(p);
                }

                addP(Tp, z, y, max(P1.get(list(x, y)), max(T)));
                addP(Tp, x, y, max(P1.get(list(z, y)), max(T)));
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
                    setP(P2, list(x, y), sum(Tp.get(list(x, y))));
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

            for (List<Node> triple : ut) {
                Node x = triple.get(0);
                Node y = triple.get(1);
                Node z = triple.get(2);

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
                Node y = record.get(0);
                Node z = record.get(1);
                Node x = record.get(2);

                if (G2.containsEdge(Edges.directedEdge(y, x))
                        && G2.containsEdge(Edges.directedEdge(x, z))
                        && G2.getEndpoint(z, x) != Endpoint.ARROW
                        && !existsRecord(amb, y, z)
                        && !existsRecord(R2, y, x, z)) {
                    G2.setEndpoint(y, z, Endpoint.ARROW);
                    addRecord(R2, y, x, z);
                    loop = true;
                }
            }

            for (List<Node> record : kite) {
                Node y = record.get(0);
                Node x = record.get(1);
                Node w = record.get(2);
                Node z = record.get(3);

                if (G2.getEndpoint(y, z) != Endpoint.ARROW
                        && G2.containsEdge(Edges.undirectedEdge(y, x))
                        && G2.containsEdge(Edges.undirectedEdge(y, w))
                        && G2.containsEdge(Edges.directedEdge(x, z))
                        && G2.containsEdge(Edges.directedEdge(w, z))
                        && G2.containsEdge(Edges.undirectedEdge(x, w))
                        && !existsRecord(amb, y, z)
                        && !existsRecord(R3, y, x, w, z)) {
                    G2.setEndpoint(y, z, Endpoint.ARROW);
                    addRecord(R3, y, x, w, z);
                    addRecord(R3, y, w, x, z);
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

        Map<List<Node>, Set<List<Node>>> e0 = new HashMap<>();
        Map<List<Node>, Set<List<Node>>> e1 = new HashMap<>();
        Map<List<Node>, Set<List<Node>>> e2 = new HashMap<>();
        Map<List<Node>, Set<List<Node>>> e3 = new HashMap<>();

        for (List<Node> record : R0) {
            Node y = record.get(0);
            Node z = record.get(1);
            Node x = record.get(2);

            addList(e0, y, z, list(x, z));
        }

        for (List<Node> record : R1) {
            Node x = record.get(0);
            Node y = record.get(1);
            Node z = record.get(2);

            addList(e1, y, z, list(x, y));
        }

        for (List<Node> record : R2) {
            Node y = record.get(0);
            Node x = record.get(1);
            Node z = record.get(2);

            addList(e2, y, z, list(y, x));
            addList(e2, y, z, list(x, z));
        }

        for (List<Node> record : R3) {
            Node y = record.get(0);
            Node x = record.get(1);
            Node w = record.get(2);
            Node z = record.get(3);

            addList(e3, y, z, list(y, x));
            addList(e3, y, z, list(y, w));
            addList(e3, y, z, list(x, z));
            addList(e3, y, z, list(w, z));
        }

        Map<List<Node>, Set<List<Node>>> eAll = union(e0, e1, e2, e3);
        Map<List<Node>, Set<List<Node>>> e123 = union(e1, e2, e3);

        Graph G3 = new EdgeListGraph(G2);

        for (Edge edge : G2.getEdges()) {
            if (!Edges.isBidirectedEdge(edge)) {
                continue;
            }

            Node y = edge.getNode1();
            Node z = edge.getNode2();

            G3.removeEdge(y, z);
            G3.addUndirectedEdge(y, z);

            addRecord(amb, y, z);
            addRecord(amb, z, y);

            for (List<Node> record : union(eAll.get(list(y, z)),
                    eAll.get(list(z, y)))) {
                Node x = record.get(0);
                Node w = record.get(1);

                G3.removeEdge(x, w);
                G3.addUndirectedEdge(x, w);

                addRecord(amb, x, w);
                addRecord(amb, w, x);
            }
        }

        // algorithm 4

        Set<List<Node>> directed = new HashSet<>();
        Set<List<Node>> undirected = new HashSet<>();

        for (Node y : nodes) {
            for (Node z : nodes) {
                if (G3.containsEdge(Edges.directedEdge(y, z))) directed.add(list(y, z));
                if (G3.containsEdge(Edges.undirectedEdge(y, z))) undirected.add(list(y, z));
            }
        }

        Map<List<Node>, Double> P3 = new HashMap<>();

        for (List<Node> pairyz : undirected) {
            Node y1 = pairyz.get(0);
            Node z1 = pairyz.get(1);

            if (!existsRecord(amb, y1, z1)) {
                P3.put(pairyz, P1.get(pairyz));
            }
        }

        Set<List<Node>> dup = new HashSet<>();

        for (List<Node> pairyz : directed) {
            Node y1 = pairyz.get(0);
            Node z1 = pairyz.get(1);

            if (e0.containsKey(pairyz) && !e123.containsKey(pairyz) && P1.containsKey(pairyz) && P2.containsKey(pairyz)) {
                P3.put(pairyz, max(P1.get(pairyz), P2.get(pairyz)));

                if (e0.get(pairyz).size() == 1 && !dup.contains(pairyz)) {
                    Node x = findR1X(R0, z1, y1);

                    List<Node> pairxz1 = list(x, z1);

                    if (e0.get(pairxz1).size() == 1) {
                        dup.add(pairxz1);
                    }
                }
            }
        }

        for (List<Node> pairyz : directed) {
            if (e0.get(pairyz).isEmpty()) {
                P2.put(pairyz, 0.0);
            }
        }

        Set<List<Node>> considered3 = new HashSet<>();
        Set<List<Node>> visited = new HashSet<>();

        while (!complement(directed, visited).isEmpty()) {
            for (List<Node> pairyz : complement(directed, visited)) {
                if (listComplement(eAll.get(pairyz), P3.keySet()).isEmpty()) {
                    Set<Double> U = new HashSet<>();

                    for (List<Node> R : R1) {
                        Node _x = R.get(0);
                        Node _y = R.get(1);

                        List<Node> pairxy = list(_x, _y);

                        U.add(P3.get(pairxy));
                    }

                    for (List<Node> R : R2) {
                        Node _y = R.get(0);
                        Node _x = R.get(1);
                        Node _z = R.get(2);

                        List<Node> pairyx = list(_y, _x);
                        List<Node> pairxz = list(_x, _z);

                        U.add(max(P3.get(pairyx), P3.get(pairxz)));
                    }

                    for (List<Node> R : R3) {
                        Node _y = R.get(0);
                        Node _x = R.get(1);
                        Node _w = R.get(2);
                        Node _z = R.get(3);

                        if (!existsRecord(considered3, _y, _x, _w, _z)) {
                            List<Node> pairyx = list(_y, _x);
                            List<Node> pairyw = list(_y, _w);
                            List<Node> pairxz = list(_x, _z);
                            List<Node> pairwz = list(_w, _z);

                            U.add(max(
                                    P3.get(pairyx),
                                    P3.get(pairyw),
                                    P3.get(pairxz),
                                    P3.get(pairwz)
                            ));

                            addRecord(considered3, _y, _x, _w, _z);
                        }
                    }

                    P3.put(pairyz, max(P1.get(pairyz), P2.get(pairyz), sum(U)));
                    visited.add(pairyz);
                }
            }
        }

        for (List<Node> pairyz : undirected) {
            Node y = pairyz.get(0);
            Node z = pairyz.get(1);

            if (!dup.contains(pairyz)) {
                List<Node> pairzy = list(z, y);
                dup.add(pairzy);
            }
        }

        // algorithm 5

        List<List<Node>> Pp = new ArrayList<>(P3.keySet());

        for (List<Node> pair : dup) {
            Pp.remove(pair);
        }

        int m = Pp.size();
        Pp.sort(Comparator.comparingDouble(P3::get));

        int R = Integer.MAX_VALUE;

        for (int i = m; i >= 1; i--) {
            if (P3.get(Pp.get(i - 1)) < alpha) {
                R = i;
                break;
            }
        }

        double sum = 0;

        for (int i = 1; i < m; i++) {
            sum += 1. / i;
        }

        double fdr = m * alpha * sum / max(R, 1);

        double[] q = new double[m + 1];
        double max = 0;
        int i = 0;

        for (int k = 1; k <= m; k++) {
            q[k] = m * P3.get(Pp.get(k - 1)) * sum / max(k, 1);

            if (q[k] > max) {
                max = q[k];
                i = k;
            }
        }

        Graph GStar = new EdgeListGraph(G3);

        if (i > 0) {
            double aStar = P3.get(Pp.get(i - 1));

            for (Node x : nodes) {
                for (Node y : nodes) {
                    if (P3.containsKey(list(x, y)) && P3.get(list(x, y)) > aStar) {
                        GStar.removeEdge(x, y);
                    }
                }
            }
        }

        return GStar;
    }

    private Node findR1X(Set<List<Node>> r0, Node y, Node z) {
        for (List<Node> R : r0) {
            Node _x = R.get(0);
            Node _y = R.get(1);
            Node _z = R.get(2);

            if (y == _y && z == _z) {
                return _x;
            }
        }

        throw new IllegalStateException();
    }

    private Set<List<Node>> complement(Set<List<Node>> A, Set<List<Node>> B) {
        Set<List<Node>> K = new HashSet<>(A);
        K.removeAll(B);
        return K;
    }

    private Set<List<Node>> listComplement(Set<List<Node>> A, Set<List<Node>> B) {
        Set<List<Node>> K = new HashSet<>(A);
        K.removeAll(B);
        return K;
    }

    private List<Node> list(Node... x) {
        List<Node> l = new ArrayList<>();
        addAll(l, x);
        return l;
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

    @SafeVarargs
    private final Map<List<Node>, Set<List<Node>>> union(Map<List<Node>, Set<List<Node>>>... r) {
        Map<List<Node>, Set<List<Node>>> union = new HashMap<>();

        for (Map<List<Node>, Set<List<Node>>> _r : r) {
            for (List<Node> pair : _r.keySet()) {
                for (List<Node> image : _r.get(pair)) {
                    Node x = pair.get(0);
                    Node y = pair.get(1);

                    addList(union, x, y, image);
                }
            }
        }

        return union;
    }

    private Set<List<Node>> getUT(Graph g) {
        List<Node> nodes = g.getNodes();
        Set<List<Node>> ut = new HashSet<>();

        for (Node x : nodes) {
            for (Node y : nodes) {
                for (Node z : nodes) {
                    if (x == z) continue;

                    if (g.isAdjacentTo(x, y) && g.isAdjacentTo(y, z) && !g.isAdjacentTo(x, z)) {
                        addRecord(ut, x, y, z);
                    }
                }
            }
        }

        return ut;
    }

    private Set<List<Node>> getTri(Graph g) {
        List<Node> nodes = g.getNodes();
        Set<List<Node>> tri = new HashSet<>();

        for (Node x : nodes) {
            for (Node y : nodes) {
                for (Node z : nodes) {
                    if (g.isAdjacentTo(y, z) && g.isAdjacentTo(x, z) && g.isAdjacentTo(y, z)) {
                        addRecord(tri, y, x, z);
                    }
                }
            }
        }

        return tri;
    }

    private Set<List<Node>> getKite(Graph g) {
        List<Node> nodes = g.getNodes();
        Set<List<Node>> kite = new HashSet<>();

        for (Node x : nodes) {
            for (Node y : nodes) {
                for (Node z : nodes) {
                    for (Node w : nodes) {
                        if (x == w) continue;
                        if (y == z) continue;

                        if (g.isAdjacentTo(x, y)
                                && g.isAdjacentTo(y, x)
                                && g.isAdjacentTo(y, w)
                                && g.isAdjacentTo(w, z)
                                && g.isAdjacentTo(y, z)
                                && !g.isAdjacentTo(x, w)
                        ) {
                            addRecord(kite, y, x, w, z);
                        }
                    }
                }
            }
        }

        return kite;
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

    private void clear(Map<List<Node>, Set<Double>> v, Node x, Node y) {
        v.computeIfAbsent(list(x, y), k -> new HashSet<>());
        v.get(list(x, y)).clear();
    }

    private void setList(Map<List<Node>, List<Node>> s, Node x, Node y, List<Node> SS) {
        s.put(list(x, y), SS);
    }

    private void addList(Map<List<Node>, Set<List<Node>>> s, Node x, Node y, List<Node> SS) {
        s.computeIfAbsent(list(x, y), k -> new HashSet<>());
        s.get(list(x, y)).add(SS);
    }

    private void addP(Map<List<Node>, Set<Double>> v, Node x, Node y, double p) {
        v.computeIfAbsent(list(x, y), k -> new HashSet<>());
        v.get(list(x, y)).add(p);
    }

    private void setP(Map<List<Node>, Double> V, List<Node> pair, double p) {
        V.put(pair, p);
    }

    private double max(Set<Double> p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double d : p) {
            if (d > max) max = d;
        }

        return max;
    }

    private double max(double... p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double _p : p) {
            if (_p > max) max = _p;
        }

        return max;
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





