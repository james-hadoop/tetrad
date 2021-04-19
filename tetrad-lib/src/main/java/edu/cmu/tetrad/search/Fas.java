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

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.TetradLogger;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
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
public class Fas implements IFas {

    /**
     * The independence test. This should be appropriate to the types
     */
    private final IndependenceTest test;

    /**
     * Specification of which edges are forbidden or required.
     */
    private IKnowledge knowledge = new Knowledge2();

    /**
     * The maximum number of variables conditioned on in any conditional independence test. If the depth is -1, it will
     * be taken to be the maximum value, which is 1000. Otherwise, it should be set to a non-negative integer.
     */
    private int depth = 1000;

    /**
     * The number of independence tests.
     */
    private int numIndependenceTests;


    /**
     * The logger, by default the empty logger.
     */
    private final TetradLogger logger = TetradLogger.getInstance();

    /**
     * The number of dependence judgements. Temporary.
     */
    private int numDependenceJudgement;

    /**
     * The sepsets found during the search.
     */
    private SepsetMap sepset = new SepsetMap();

    private final NumberFormat nf = new DecimalFormat("0.00E0");

    /**
     * True iff verbose output should be printed.
     */
    private boolean verbose = false;

    private PrintStream out = System.out;
    private Graph initialGraph = null;
    private boolean h3 = true;

    //==========================CONSTRUCTORS=============================//

    /**
     * Constructs a new FastAdjacencySearch.
     */
    public Fas(Graph initialGraph, IndependenceTest test) {
        if (initialGraph != null) {
            this.initialGraph = new EdgeListGraph(initialGraph);
        }
        this.test = test;
    }

    public Fas(IndependenceTest test) {
        this.test = test;
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
    public Graph search() {
        this.logger.log("info", "Starting Fast Adjacency Search.");

        if (!h3) {
            return new FasOrig(initialGraph, test).search();
        } else {

            List<Edge> edges = new ArrayList<>();
            Map<Edge, Double> unconditionalEdgeScores = new HashMap<>();
            List<Node> nodes = test.getVariables();

            for (int i = 0; i < nodes.size(); i++) {
                for (int j = i + 1; j < nodes.size(); j++) {
                    edges.add(Edges.undirectedEdge(nodes.get(i), nodes.get(j)));
                }
            }

            sortedEdges(test, edges, unconditionalEdgeScores);

            sepset = new SepsetMap();

            int _depth = depth;

            if (_depth == -1) {
                _depth = 1000;
            }

            Map<Node, Set<Node>> adjacencies = new HashMap<>();

            for (Node node : nodes) {
                Set<Node> set = new LinkedHashSet<>();

                for (Node _node : nodes) {
                    if (_node == node) continue;
                    set.add(_node);
                }

                adjacencies.put(node, set);
            }

            for (Edge edge : new LinkedHashSet<>(unconditionalEdgeScores.keySet())) {
                if (unconditionalEdgeScores.get(edge) < 0) {
                    edges.remove(edge);
                    adjacencies.get(edge.getNode1()).remove(edge.getNode2());
                    adjacencies.get(edge.getNode2()).remove(edge.getNode1());
                    unconditionalEdgeScores.remove(edge);
                }
            }

            for (int d = 0; d <= _depth; d++) {
                boolean more;

                more = searchAtDepth(edges, nodes, test, adjacencies, d);

                if (!more) {
                    break;
                }
            }

            // The search graph. It is assumed going in that all of the true adjacencies of x are in this graph for every node
            // x. It is hoped (i.e. true in the large sample limit) that true adjacencies are never removed.
            Graph graph = new EdgeListGraph(nodes);

            for (int i = 0; i < nodes.size(); i++) {
                for (int j = i + 1; j < nodes.size(); j++) {
                    Node x = nodes.get(i);
                    Node y = nodes.get(j);

                    if (adjacencies.get(x).contains(y)) {
                        graph.addUndirectedEdge(x, y);
                    }
                }
            }

            this.logger.log("info", "Finishing Fast Adjacency Search.");

            return graph;
        }
    }

    private void sortedEdges(IndependenceTest test, List<Edge> edges, Map<Edge, Double> scores) {
        for (Edge edge : edges) {
            test.isIndependent(edge.getNode1(), edge.getNode2(), new ArrayList<>());
            scores.put(edge, test.getScore());
        }

        edges.sort(Comparator.comparing(scores::get));
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

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) {
            throw new NullPointerException("Cannot set knowledge to null");
        }
        this.knowledge = knowledge;
    }

    //==============================PRIVATE METHODS======================/

    private int freeDegree(List<Node> nodes, Map<Node, Set<Node>> adjacencies) {
        int max = 0;

        for (Node x : nodes) {
            Set<Node> opposites = adjacencies.get(x);

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

    private boolean searchAtDepth(List<Edge> edges, List<Node> nodes, final IndependenceTest test, Map<Node, Set<Node>> adjacencies, int depth) {
        for (Edge edge : edges) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            if (depth == 0 && initialGraph != null) {
                Node x2 = initialGraph.getNode(x.getName());
                Node y2 = initialGraph.getNode(y.getName());

                if (!initialGraph.isAdjacentTo(x2, y2)) {
                    continue;
                }
            }

            checkSide(test, adjacencies, depth, x, y);
            checkSide(test, adjacencies, depth, y, x);
        }

        return freeDegree(nodes, adjacencies) > depth;
    }

    private void checkSide(IndependenceTest test, Map<Node, Set<Node>> adjacencies, int depth, Node x, Node y) {
        if (!adjacencies.get(x).contains(y)) return;

        List<Node> _adjx = new ArrayList<>(adjacencies.get(x));
        _adjx.remove(y);
        List<Node> ppx = possibleParents(x, _adjx, knowledge, y);

        Map<Node, Double> scores = new HashMap<>();

        for (Node node : ppx) {
            test.isIndependent(x, node, new ArrayList<>());
            scores.put(node, test.getScore());
        }

        ppx.sort(Comparator.comparing(scores::get));
        Collections.reverse(ppx);

        if (ppx.size() >= depth) {
            ChoiceGenerator cg = new ChoiceGenerator(ppx.size(), depth);
            int[] choice;

            while ((choice = cg.next()) != null) {
                if (Thread.currentThread().isInterrupted()) {
                    break;
                }

                List<Node> Z = GraphUtils.asList(choice, ppx);

                boolean independent;

                try {
                    numIndependenceTests++;
                    independent = test.isIndependent(x, y, Z);
                } catch (Exception e) {
                    independent = false;
                }

                if (!independent) {
                    numDependenceJudgement++;
                }

                boolean noEdgeRequired =
                        knowledge.noEdgeRequired(x.getName(), y.getName());

                if (independent && noEdgeRequired) {
                    adjacencies.get(x).remove(y);
                    adjacencies.get(y).remove(x);

                    getSepsets().set(x, y, Z);

                    if (verbose) {
                        TetradLogger.getInstance().forceLogMessage(SearchLogUtils.independenceFact(x, y, Z) +
                                " score = " + nf.format(test.getScore()));
                        out.println(SearchLogUtils.independenceFactMsg(x, y, Z, test.getPValue()));
                    }
                }
            }
        }
    }

    private List<Node> possibleParents(Node x, List<Node> adjx,
                                       IKnowledge knowledge, Node y) {
        List<Node> possibleParents = new LinkedList<>();
        String _x = x.getName();

        for (Node z : adjx) {
            if (z == y) continue;
            String _z = z.getName();

            if (possibleParentOf(_z, _x, knowledge)) {
                possibleParents.add(z);
            }
        }

        return possibleParents;
    }

    private boolean possibleParentOf(String z, String x, IKnowledge knowledge) {
        return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
    }

    public int getNumIndependenceTests() {
        return numIndependenceTests;
    }

    public void setTrueGraph(Graph trueGraph) {
    }

    public int getNumDependenceJudgments() {
        return numDependenceJudgement;
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

    @Override
    public boolean isAggressivelyPreventCycles() {
        return false;
    }

    @Override
    public void setAggressivelyPreventCycles(boolean aggressivelyPreventCycles) {

    }

    @Override
    public IndependenceTest getIndependenceTest() {
        return null;
    }

    @Override
    public Graph search(List<Node> nodes) {
        return null;
    }

    @Override
    public long getElapsedTime() {
        return 0;
    }

    @Override
    public List<Node> getNodes() {
        return test.getVariables();
    }

    @Override
    public List<Triple> getAmbiguousTriples(Node node) {
        return null;
    }

    @Override
    public void setOut(PrintStream out) {
        this.out = out;
    }

    public void setH3(boolean fasOrig) {
        this.h3 = fasOrig;
    }
}

