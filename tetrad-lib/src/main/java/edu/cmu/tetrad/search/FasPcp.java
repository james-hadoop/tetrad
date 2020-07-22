///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014 by Peter Spirtes, Richard Scheines, Joseph   //
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

import java.util.*;

/**
 * Implements the "fast adjacency search" used in several causal algorithm in this package. In the fast adjacency
 * search, at a given stage of the search, an edge X*-*Y is removed from the graph if X _||_ Y | S, where S is a subset
 * of size d either of adj(X) or of adj(Y), where d is the depth of the search. The fast adjacency search performs this
 * procedure for each pair of adjacent edges in the graph and for each depth d = 0, 1, 2, ..., d1, where d1 is either
 * the maximum depth or else the first such depth at which no edges can be removed. The interpretation of this adjacency
 * search is different for different algorithm, depending on the assumptions of the algorithm. A mapping from {x, y} to
 * S({x, y}) is returned for edges x *-* y that have been removed.
 *
 * @author Joseph Ramsey.
 */
public class FasPcp {

    /**
     * The search nodes.
     */
    private final List<Node> nodes;

    /**
     * The independence test. This should be appropriate to the types
     */
    private final IndependenceTest test;

    /**
     * The maximum number of variables conditioned on in any conditional independence test. If the depth is -1, it will
     * be taken to be the maximum value, which is 1000. Otherwise, it should be set to a non-negative integer.
     */
    private int depth = 1000;

    /**
     * The sepsets found during the search.
     */
    private SepsetMap sepset = new SepsetMap();

    /**
     * True iff verbose output should be printed.
     */
    private boolean verbose = false;

    /**
     * List of p-values for each edge.
     */
    private final Map<NodePair, List<Double>> p1 = new HashMap<>();

    /**
     * List of unique identifiers for each edge.
     */
    private final Map<NodePair, Set<Object>> I = new HashMap<>();


    //==========================CONSTRUCTORS=============================//

    public FasPcp(IndependenceTest test) {
        this.test = test;
        this.nodes = test.getVariables();
    }

    //==========================PUBLIC METHODS===========================//

    /**
     * Discovers all adjacencies in data.  The procedure is to remove edges in the graph which connect pairs of
     * variables which are independent conditional on some other set of variables in the graph (the "sepset"). These are
     * removed in tiers.  First, edges which are independent conditional on zero other variables are removed, then edges
     * which are independent conditional on one other variable are removed, then two, then three, and so on, until no
     * more edges can be removed from the graph.  The edges which remain in the graph after this procedure are the
     * adjacencies in the data.
     *
     * @return a SepSet, which indicates which variables are independent conditional on which other variables
     */
    public Graph algorithm1() {
        sepset = new SepsetMap();

        int _depth = depth;

        if (_depth == -1) {
            _depth = 1000;
        }

        Graph graph = completeGraph(nodes);

        for (int d = 0; d <= _depth; d++) {
            boolean more = searchAtDepth(nodes, graph, d);

            if (!more) {
                break;
            }
        }

        for (NodePair key : p1.keySet()) {
            List<Double> value = p1.get(key);
            if (value != null && !value.isEmpty()) {
                double max = max(value);
                value.clear();
                value.add(max);
                p1.put(key, Collections.singletonList(max));
            }

            insertI(key, new Object());
        }

        return graph;
    }

    private void insertI(NodePair key, Object o) {
        I.computeIfAbsent(key, k -> new HashSet<>());
        I.get(key).add(o);
    }

    private double max(List<Double> p) {
        double max = Double.NEGATIVE_INFINITY;

        for (double d : p) {
            if (d > max) max = d;
        }

        return max;
    }

    public Map<NodePair, List<Double>> getP1() {
        return p1;
    }

    public Map<NodePair, Set<Object>> getI() {
        return I;
    }

    public int getDepth() {
        return depth;
    }

    public void setDepth(int depth) {
        if (depth < -1) {
            throw new IllegalArgumentException(
                    "Depth must be -1 (unlimited) or >= 0.");
        }

        this.depth = depth;
    }

    //==============================PRIVATE METHODS======================/

    private int freeDegree(List<Node> nodes, Graph graph) {
        int max = 0;

        for (Node x : nodes) {
            List<Node> opposites = graph.getAdjacentNodes(x);

            for (Node y : opposites) {
                Set<Node> adjx = new HashSet<>(opposites);
                adjx.remove(y);

                if (adjx.size() > max) {
                    max = adjx.size();
                }
            }
        }

        return max;
    }

    private boolean searchAtDepth(List<Node> nodes, Graph graph, int depth) {
        for (Node x : nodes) {
            List<Node> adjx = graph.getAdjacentNodes(x);

            EDGE:
            for (Node y : adjx) {
                List<Node> _adjx = graph.getAdjacentNodes(x);
                _adjx.remove(y);
                List<Node> ppx = possibleParents(graph, x, _adjx);

                if (ppx.size() >= depth) {
                    ChoiceGenerator cg = new ChoiceGenerator(ppx.size(), depth);
                    int[] choice;

                    while ((choice = cg.next()) != null) {
                        List<Node> condSet = GraphUtils.asList(choice, ppx);

                        double p;

                        synchronized (test) {
                            test.isIndependent(x, y, condSet);
                            p = test.getPValue();
                        }

                        if (p <= test.getAlpha()) {
                            p1.computeIfAbsent(new NodePair(x, y), k -> new ArrayList<>());
                            p1.get(new NodePair(x, y)).add(p);
                        } else {
                            graph.removeEdge(x, y);
                            getSepsets().set(x, y, condSet);
                            continue EDGE;
                        }
                    }
                }
            }
        }

        System.out.println("depth = " + depth + " " + graph);
        System.out.println("p1 = " + p1);
        System.out.println("I = " + I);

        return freeDegree(nodes, graph) > depth;
    }

    private List<Node> possibleParents(Graph graph, Node x, List<Node> adjx) {
        List<Node> pp = new ArrayList<>(adjx);
        pp.removeAll(graph.getChildren(x));
        return pp;
    }

    public SepsetMap getSepsets() {
        return sepset;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public static Graph completeGraph(List<Node> nodes) {
        Graph graph2 = new EdgeListGraph(nodes);

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node node1 = nodes.get(i);
                Node node2 = nodes.get(j);
                graph2.addNondirectedEdge(node1, node2);
            }
        }

        return graph2;
    }
}

