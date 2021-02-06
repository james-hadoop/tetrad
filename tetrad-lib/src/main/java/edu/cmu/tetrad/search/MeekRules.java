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

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.TetradLogger;

import java.io.PrintStream;
import java.util.*;

/**
 * Implements Meek's complete orientation rule set for PC (Chris Meek (1995), "Causal inference and causal explanation
 * with background knowledge"), modified for Conservative PC to check noncolliders against recorded noncolliders before
 * orienting.
 * <p>
 * Rule R4 is only performed if knowledge is nonempty.
 *
 * @author Joseph Ramsey
 */
public class MeekRules implements ImpliedOrientation {

    private IKnowledge knowledge = new Knowledge2();

    //True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
    //graphs).
    private boolean aggressivelyPreventCycles = false;

    // If knowledge is available.
    boolean useRule4;

    //The logger to use.
    private final Map<Edge, Edge> changedEdges = new HashMap<>();

    // Whether verbose output should be generated.

    // Where verbose output should be sent.
    private PrintStream out;

    // The initial list of nodes to visit.
    private List<Node> nodes = new ArrayList<>();

    // The list of nodes actually visited.
    private final Set<Node> visited = new HashSet<>();

    // True if verbose output should be printed.
    private boolean verbose = false;

    /**
     * Constructs the <code>MeekRules</code> with no logging.
     */
    public MeekRules() {
        useRule4 = !knowledge.isEmpty();
    }

    //======================== Public Methods ========================//

    public void orientImplied(Graph graph) {
        orientImplied(graph, graph.getNodes());
    }

    public Set<Node> orientImplied(Graph graph, List<Node> nodes) {
        this.nodes = nodes;
        this.visited.addAll(nodes);

        TetradLogger.getInstance().log("impliedOrientations", "Starting Orientation Step D.");
        orientUsingMeekRulesLocally(knowledge, graph);
        TetradLogger.getInstance().log("impliedOrientations", "Finishing Orientation Step D.");

        return visited;
    }

    public void revertToColliders(List<Node> nodes, Graph graph) {
        boolean reverted = true;

        while (reverted) {
            reverted = false;

            for (Node node : nodes) {
                if (revertToColliders(node, graph)) {
                    reverted = true;
                }
            }
        }
    }

    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) throw new IllegalArgumentException();
        this.knowledge = knowledge;
    }


    public boolean isAggressivelyPreventCycles() {
        return aggressivelyPreventCycles;
    }

    public void setAggressivelyPreventCycles(boolean aggressivelyPreventCycles) {
        this.aggressivelyPreventCycles = aggressivelyPreventCycles;
    }

    public Map<Edge, Edge> getChangedEdges() {
        return changedEdges;
    }

    public void setOut(PrintStream out) {
        this.out = out;
    }

    public PrintStream getOut() {
        return out;
    }

    public Set<Node> getVisited() {
        return visited;
    }

    //============================== Private Methods ===================================//

    private void orientUsingMeekRulesLocally(IKnowledge knowledge, Graph graph) {

//        this.nodes = graph.getNodes();
//
//        boolean reverted = true;
//
//        while (reverted) {
//            reverted = false;
//
//            for (Node node : nodes) {
//                if (revertToColliders(node, graph)) {
//                    reverted = true;
//                }
//            }
//        }

        revertToColliders(this.nodes, graph);

        boolean oriented;

        do {
            oriented = false;

            for (Node y : graph.getNodes()) {
                if (meekR1(y, graph)) oriented = true;
                if (meekR2(y, graph)) oriented = true;
                if (meekR3(y, graph)) oriented = true;
//                if (meekR4(y, graph)) oriented = true;

            }
        } while (oriented);
    }

    /**
     * Meek's rule R1: if a-->node, node---c, and a not adj to c, then a-->c
     */
    private boolean meekR1(Node b, Graph graph) {
        List<Node> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            return false;
        }

        ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
        int[] choice;
        boolean oriented = false;

        while ((choice = cg.next()) != null) {
            List<Node> nodes = GraphUtils.asList(choice, adjacentNodes);
            Node a = nodes.get(0);
            Node c = nodes.get(1);

            if (r1Helper(a, b, c, graph)) {
                oriented = true;
            } else if (r1Helper(c, b, a, graph)) {
                oriented = true;
            }
        }

        return oriented;
    }

    private boolean r1Helper(Node a, Node b, Node c, Graph graph) {
        boolean directed = false;

        if (!graph.isAdjacentTo(a, c) && graph.isDirectedFromTo(a, b) && graph.isUndirectedFromTo(b, c)) {
            directed = direct(b, c, graph);
            String message = SearchLogUtils.edgeOrientedMsg(
                    "Meek R1 triangle (" + a + "-->" + b + "---" + c + ")", graph.getEdge(b, c));
            log(message);
        }

        return directed;
    }

    /**
     * If a-->b-->c, a--c, then b-->c.
     */
    private boolean meekR2(Node b, Graph graph) {
        List<Node> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            return false;
        }

        ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
        int[] choice;
        boolean oriented = false;

        while ((choice = cg.next()) != null) {
            List<Node> nodes = GraphUtils.asList(choice, adjacentNodes);
            Node a = nodes.get(0);
            Node c = nodes.get(1);

            if (r2Helper(a, b, c, graph)) {
                oriented = true;
            } else if (r2Helper(c, b, a, graph)) {
                oriented = true;
            } else if (r2Helper(a, c, b, graph)) {
                oriented = true;
            } else if (r2Helper(c, a, b, graph)) {
                oriented = true;
            }
        }

        return oriented;
    }

    private boolean r2Helper(Node a, Node b, Node c, Graph graph) {
        boolean directed = false;

        if (graph.isDirectedFromTo(a, b) &&
                graph.isDirectedFromTo(b, c) &&
                graph.isUndirectedFromTo(a, c)) {
            directed = direct(a, c, graph);
            log(SearchLogUtils.edgeOrientedMsg(
                    "Meek R2 triangle (" + a + "-->" + b + "-->" + c + ", " + a + "---" + c + ")", graph.getEdge(a, c)));
        }

        return directed;
    }

    /**
     * Meek's rule R3. If a--b, a--c, a--d, c-->b, d-->b, then orient a-->b.
     */
    private boolean meekR3(Node a, Graph graph) {
        List<Node> adjacentNodes = graph.getAdjacentNodes(a);

        if (adjacentNodes.size() < 3) {
            return false;
        }

        boolean oriented = false;

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

                    if (r3Helper(a, d, b, c, graph)) {
                        oriented = true;
                    }
                }
            }
        }

        return oriented;
    }

    private boolean r3Helper(Node a, Node d, Node b, Node c, Graph graph) {
        boolean directed = false;

        boolean b4 = graph.isUndirectedFromTo(d, c);
        boolean b5 = graph.isUndirectedFromTo(d, b);
        boolean b6 = graph.isDirectedFromTo(b, a);
        boolean b7 = graph.isDirectedFromTo(c, a);
        boolean b8 = graph.isUndirectedFromTo(d, a);

        if (b4 && b5 && b6 && b7 && b8) {
            directed = direct(d, a, graph);
            log(SearchLogUtils.edgeOrientedMsg("Meek R3 " + d + "--" + a + ", " + b + ", "
                    + c, graph.getEdge(d, a)));
        }

        return directed;
    }

    private boolean meekR4(Node node, Graph graph) {
        if (!useRule4) {
            return false;
        }

        List<Node> adjacentNodes = graph.getAdjacentNodes(node);

        if (adjacentNodes.size() < 3) {
            return false;
        }

        for (Node b : adjacentNodes) {
            List<Node> otherAdjacents = new LinkedList<>(adjacentNodes);
            otherAdjacents.remove(b);

            ChoiceGenerator cg = new ChoiceGenerator(otherAdjacents.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node a = otherAdjacents.get(combination[0]);
                Node c = otherAdjacents.get(combination[1]);

                if (graph.isAdjacentTo(a, b) && graph.isAdjacentTo(b, c)) {
                    boolean directed = false;

                    if (graph.isDirectedFromTo(c, node) && graph.isDirectedFromTo(node, a)) {
                        directed = direct(b, a, graph);
                        log(SearchLogUtils.edgeOrientedMsg("Meek R4", graph.getEdge(b, a)));
                    }

                    if (directed) return true;
                    directed = false;

                    if (graph.isDirectedFromTo(a, node) && graph.isDirectedFromTo(node, c)) {
                        directed = direct(b, c, graph);
                        log(SearchLogUtils.edgeOrientedMsg("Meek R4", graph.getEdge(b, c)));
                    }

                    if (directed) return true;
                }
            }
        }

        return false;
    }

    private boolean direct(Node a, Node c, Graph graph) {
        if (!isArrowpointAllowed(a, c, knowledge)) return false;

        Edge before = graph.getEdge(a, c);

        if (knowledge != null && knowledge.isForbidden(a.getName(), c.getName())) {
            return false;
        }

        Edge after = Edges.directedEdge(a, c);

        visited.add(a);
        visited.add(c);

        graph.removeEdge(before);
        graph.addEdge(after);

        return true;
    }

    private static boolean isArrowpointAllowed(Node from, Node to, IKnowledge knowledge) {
        if (knowledge.isEmpty()) return true;
        return !knowledge.isRequired(to.toString(), from.toString()) &&
                !knowledge.isForbidden(from.toString(), to.toString());
    }

    private boolean revertToColliders(Node y, Graph graph) {
//        List<Node> adj = graph.getAdjacentNodes(y);
//
//        if (adj.size() < 2) return false;
//
//        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
//        int[] choice;
//        Set<Node> keep = new HashSet<>();
//
//        if ((choice = gen.next()) != null) {
//            List<Node> _adj = GraphUtils.asList(choice, adj);
//            Node x = _adj.get(0);
//            Node z = _adj.get(1);
//
//            if (!graph.isAdjacentTo(x, z) && graph.isDefCollider(x, y, z)) {
//                keep.add(x);
//                keep.add(z);
//            }
//        }

//        boolean changed = false;
//
//        for (Node p : graph.getParents(y)) {
//            if (!keep.contains(p)) {
//                graph.removeEdge(p, y);
//                graph.addUndirectedEdge(p, y);
//                changed = true;
//            }
//        }
//
//        return changed;

//
        boolean did = false;
//
//        TRUE:
//        while (true) {
//            List<Node> parents = graph.getParents(y);
//
//            for (Node p : parents) {
//                Set<Node> _parents = new HashSet<>(parents);
//                _parents.remove(p);
//                if (graph.getAdjacentNodes(p).containsAll(_parents)) {
//                    graph.removeEdge(p, y);
//                    graph.addUndirectedEdge(p, y);
//                    did = true;
//                    continue TRUE;
//                }
//            }
//
//            break;
//        }

        List<Node> parents = graph.getParents(y);

        for (Node p : parents) {
            Set<Node> _parents = new HashSet<>(parents);
            _parents.remove(p);
            if (graph.getAdjacentNodes(p).containsAll(_parents)) {
//                String msg = "Unorienting " + graph.getEdge(p, y);

                graph.removeEdge(p, y);
                graph.addUndirectedEdge(p, y);

//                msg += " to " + graph.getEdge(p, y);
//
//                if (verbose) {
//                    TetradLogger.getInstance().forceLogMessage(msg);
//                }
//
                did = true;
            }
        }

        return did;
    }

    private void log(String message) {
        if (verbose) {
            TetradLogger.getInstance().forceLogMessage(message);
        }
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}



