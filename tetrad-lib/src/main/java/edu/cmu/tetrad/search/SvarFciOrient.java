///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015, 2022 by Peter Spirtes, Richard        //
// Scheines, Joseph Ramsey, and Clark Glymour.                               //
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
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;


/**
 * Extends Erin Korber's implementation of the Fast Causal Inference algorithm (found in FCI.java) with Jiji Zhang's
 * Augmented FCI rules (found in sec. 4.1 of Zhang's 2006 PhD dissertation, "Causal Inference and Reasoning in Causally
 * Insufficient Systems").
 * <p>
 * This class is based off a copy of FCI.java taken from the repository on 2008/12/16, revision 7306. The extension is
 * done by extending doFinalOrientation() with methods for Zhang's rules R5-R10 which implements the augmented search.
 * (By a remark of Zhang's, the rule applications can be staged in this way.)
 *
 * @author Erin Korber, June 2004
 * @author Alex Smith, December 2008
 * @author Joseph Ramsey
 * @author Choh-Man Teng
 * @author Daniel Malinsky
 * <p>
 * This is a copy of FciOrient.java for the SvarFCI algorithm. The main difference is that if an edge is orient, it will
 * also orient all homologous edges to preserve the time-repeating structure assumed by SvarFCI. Based on (but not
 * identicial to) code by Entner and Hoyer for their 2010 paper. Modified by DMalinsky 4/20/2016.
 */
public final class SvarFciOrient {

    /**
     * The SepsetMap being constructed.
     */
    private final SepsetProducer sepsets;

    private IKnowledge knowledge = new Knowledge2();

    private boolean changeFlag = true;

    /**
     * flag for complete rule set, true if should use complete rule set, false otherwise.
     */
    private boolean completeRuleSetUsed;

    /**
     * True iff the possible dsep search is done.
     */
    private boolean possibleDsepSearchDone = true;

    /**
     * The maximum length for any discriminating path. -1 if unlimited; otherwise, a positive integer.
     */
    private int maxPathLength = -1;

    /**
     * The logger to use.
     */
    private final TetradLogger logger = TetradLogger.getInstance();

    /**
     * True iff verbose output should be printed.
     */
    private boolean verbose;

    private Graph truePag;

    private final IndependenceTest independenceTest;

    //============================CONSTRUCTORS============================//

    /**
     * Constructs a new FCI search for the given independence test and background knowledge.
     */
    public SvarFciOrient(SepsetProducer sepsets, IndependenceTest independenceTest) {
        this.sepsets = sepsets;
        this.independenceTest = independenceTest;
    }

    //========================PUBLIC METHODS==========================//

    public Graph orient(Graph graph) {

        this.logger.log("info", "Starting FCI algorithm.");

        ruleR0(graph);

        if (this.verbose) {
            System.out.println("R0");
        }


        // Step CI D. (Zhang's step F4.)
        doFinalOrientation(graph);

//        graph.closeInducingPaths();   //to make sure it's a legal PAG
        if (this.verbose) {
            this.logger.log("graph", "Returning graph: " + graph);
        }

        return graph;
    }

    public SepsetProducer getSepsets() {
        return this.sepsets;
    }

    /**
     * The background knowledge.
     */
    public IKnowledge getKnowledge() {
        return this.knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) {
            throw new NullPointerException();
        }

        this.knowledge = knowledge;
    }

    /**
     * @return true if Zhang's complete rule set should be used, false if only R1-R4 (the rule set of the original FCI)
     * should be used. False by default.
     */
    public boolean isCompleteRuleSetUsed() {
        return this.completeRuleSetUsed;
    }

    /**
     * @param completeRuleSetUsed set to true if Zhang's complete rule set should be used, false if only R1-R4 (the rule
     *                            set of the original FCI) should be used. False by default.
     */
    public void setCompleteRuleSetUsed(boolean completeRuleSetUsed) {
        this.completeRuleSetUsed = completeRuleSetUsed;
    }

    //===========================PRIVATE METHODS=========================//

    private List<Node> getSepset(Node i, Node k) {
        return this.sepsets.getSepset(i, k);
    }

    /**
     * Orients colliders in the graph.  (FCI Step C)
     * <p>
     * Zhang's step F3, rule R0.
     */
    public void ruleR0(Graph graph) {
        graph.reorientAllWith(Endpoint.CIRCLE);
        fciOrientbk(this.knowledge, graph, graph.getNodes());

        List<Node> nodes = graph.getNodes();

        for (Node b : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(b);

            if (adjacentNodes.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node a = adjacentNodes.get(combination[0]);
                Node c = adjacentNodes.get(combination[1]);
                if (this.knowledge.isInWhichTier(a) == 0 && this.knowledge.isInWhichTier(b) == 0 && this.knowledge.isInWhichTier(c) == 0) {
                    System.out.println("Skipping triple a,b,c : " + a + " , " + b + " , " + c);
                    continue; // This is added as a temporary measure. Sepsets for lagged vars may be out of window, leading to incorrect collider orientations
                }
                // Skip triples that are shielded.
                if (graph.isAdjacentTo(a, c)) {
                    continue;
                }

                if (graph.isDefCollider(a, b, c)) {
                    continue;
                }

                if (this.sepsets.isCollider(a, b, c)) {
                    if (!isArrowpointAllowed(a, b, graph)) {
                        continue;
                    }

                    if (!isArrowpointAllowed(c, b, graph)) {
                        continue;
                    }

                    graph.setEndpoint(a, b, Endpoint.ARROW);
                    graph.setEndpoint(c, b, Endpoint.ARROW);
                    if (this.verbose) {
                        this.logger.log("colliderOrientations", SearchLogUtils.colliderOrientedMsg(a, b, c));
                        System.out.println(SearchLogUtils.colliderOrientedMsg(a, b, c));
                        final String location = "R0";

                        printWrongColliderMessage(a, b, c, graph);
                    }
                    this.orientSimilarPairs(graph, this.knowledge, a, b, Endpoint.ARROW);
                    this.orientSimilarPairs(graph, this.knowledge, c, b, Endpoint.ARROW);
                }
            }
        }
    }

    private void printWrongColliderMessage(Node a, Node b, Node c, Graph graph) {
        if (this.truePag != null && graph.isDefCollider(a, b, c) && !this.truePag.isDefCollider(a, b, c)) {
            System.out.println("R0" + ": Orienting collider by mistake: " + a + "*-&gt;" + b + "&lt;-*" + c);
        }
    }


    /**
     * Orients the graph according to rules in the graph (FCI step D).
     * <p>
     * Zhang's step F4, rules R1-R10.
     */
    public void doFinalOrientation(Graph graph) {
        if (this.completeRuleSetUsed) {
            zhangFinalOrientation(graph);
        } else {
            spirtesFinalOrientation(graph);
        }
    }

    private void spirtesFinalOrientation(Graph graph) {
        this.changeFlag = true;
        boolean firstTime = true;

        while (this.changeFlag) {
            this.changeFlag = false;
            rulesR1R2cycle(graph);
            ruleR3(graph);

            // R4 requires an arrow orientation.
            if (this.changeFlag || (firstTime && !this.knowledge.isEmpty())) {
                ruleR4B(graph);
                firstTime = false;
            }

            if (this.verbose) {
                System.out.println("Epoch");
            }
        }
    }

    private void zhangFinalOrientation(Graph graph) {
        this.changeFlag = true;
        boolean firstTime = true;

        while (this.changeFlag) {
            this.changeFlag = false;
            rulesR1R2cycle(graph);
            ruleR3(graph);

            // R4 requires an arrow orientation.
            if (this.changeFlag || (firstTime && !this.knowledge.isEmpty())) {
                ruleR4B(graph);
                firstTime = false;
            }

            if (this.verbose) {
                System.out.println("Epoch");
            }
        }

        if (isCompleteRuleSetUsed()) {
            // Now, by a remark on page 100 of Zhang's dissertation, we apply rule
            // R5 once.
            ruleR5(graph);

            // Now, by a further remark on page 102, we apply R6,R7 as many times
            // as possible.
            this.changeFlag = true;

            while (this.changeFlag) {
                this.changeFlag = false;
                ruleR6R7(graph);
            }

            // Finally, we apply R8-R10 as many times as possible.
            this.changeFlag = true;

            while (this.changeFlag) {
                this.changeFlag = false;
                rulesR8R9R10(graph);
            }

        }
    }

    //Does all 3 of these rules at once instead of going through all
    // triples multiple times per iteration of doFinalOrientation.

    public void rulesR1R2cycle(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node B : nodes) {
            List<Node> adj = graph.getAdjacentNodes(B);

            if (adj.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adj.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node A = adj.get(combination[0]);
                Node C = adj.get(combination[1]);

                //choice gen doesnt do diff orders, so must switch A & C around.
                ruleR1(A, B, C, graph);
                ruleR1(C, B, A, graph);
                ruleR2(A, B, C, graph);
                ruleR2(C, B, A, graph);
            }
        }
    }

    /// R1, away from collider
    // If a*-&gt;bo-*c and a, c not adjacent then a*-&gt;b->c
    private void ruleR1(Node a, Node b, Node c, Graph graph) {
        if (graph.isAdjacentTo(a, c)) {
            return;
        }

        if (graph.getEndpoint(a, b) == Endpoint.ARROW && graph.getEndpoint(c, b) == Endpoint.CIRCLE) {
            if (!isArrowpointAllowed(b, c, graph)) {
                return;
            }

            graph.setEndpoint(c, b, Endpoint.TAIL);
            graph.setEndpoint(b, c, Endpoint.ARROW);
            this.changeFlag = true;

            if (this.verbose) {
                this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Away from collider", graph.getEdge(b, c)));
                System.out.println(SearchLogUtils.edgeOrientedMsg("Away from collider", graph.getEdge(b, c)));
            }
            this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), b, c, Endpoint.ARROW);
        }
    }

    private boolean isNoncollider(Node a, Node b, Node c) {
        return this.sepsets.isNoncollider(a, b, c);
    }

    //if a*-oc and either a-->b*-&gt;c or a*-&gt;b-->c, then a*-&gt;c
    // This is Zhang's rule R2.
    private void ruleR2(Node a, Node b, Node c, Graph graph) {
        if ((graph.isAdjacentTo(a, c)) &&
                (graph.getEndpoint(a, c) == Endpoint.CIRCLE)) {

            if ((graph.getEndpoint(a, b) == Endpoint.ARROW) &&
                    (graph.getEndpoint(b, c) == Endpoint.ARROW) && (
                    (graph.getEndpoint(b, a) == Endpoint.TAIL) ||
                            (graph.getEndpoint(c, b) == Endpoint.TAIL))) {

                if (!isArrowpointAllowed(a, c, graph)) {
                    return;
                }

                graph.setEndpoint(a, c, Endpoint.ARROW);
                this.orientSimilarPairs(graph, this.getKnowledge(), a, c, Endpoint.ARROW);
                if (this.verbose) {
                    this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Away from ancestor", graph.getEdge(a, c)));
                    System.out.println(SearchLogUtils.edgeOrientedMsg("Away from ancestor", graph.getEdge(a, c)));
                }

                this.changeFlag = true;
            }
        }
    }

    /**
     * Implements the double-triangle orientation rule, which states that if D*-oB, A*-&gt;B&lt;-*C and A*-oDo-*C, then
     * D*-&gt;B.
     * <p>
     * This is Zhang's rule R3.
     */
    public void ruleR3(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node B : nodes) {

            List<Node> intoBArrows = graph.getNodesInTo(B, Endpoint.ARROW);
            List<Node> intoBCircles = graph.getNodesInTo(B, Endpoint.CIRCLE);

            for (Node D : intoBCircles) {
                if (intoBArrows.size() < 2) {
                    continue;
                }

                ChoiceGenerator gen = new ChoiceGenerator(intoBArrows.size(), 2);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    Node A = intoBArrows.get(choice[0]);
                    Node C = intoBArrows.get(choice[1]);

                    if (graph.isAdjacentTo(A, C)) {
                        continue;
                    }

                    if (!graph.isAdjacentTo(A, D) ||
                            !graph.isAdjacentTo(C, D)) {
                        continue;
                    }

                    if (graph.getEndpoint(A, D) != Endpoint.CIRCLE) {
                        continue;
                    }

                    if (graph.getEndpoint(C, D) != Endpoint.CIRCLE) {
                        continue;
                    }

                    if (!isArrowpointAllowed(D, B, graph)) {
                        continue;
                    }

                    graph.setEndpoint(D, B, Endpoint.ARROW);
                    this.orientSimilarPairs(graph, this.getKnowledge(), D, B, Endpoint.ARROW);
                    if (this.verbose) {
                        this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Double triangle", graph.getEdge(D, B)));
                        System.out.println(SearchLogUtils.edgeOrientedMsg("Double triangle", graph.getEdge(D, B)));
                    }

                    this.changeFlag = true;
                }
            }
        }
    }


    /**
     * The triangles that must be oriented this way (won't be done by another rule) all look like the ones below, where
     * the dots are a collider path from L to A with each node on the path (except L) a parent of C.
     * <pre>
     *          B
     *         xo           x is either an arrowhead or a circle
     *        /  \
     *       v    v
     * L....A --> C
     * </pre>
     * <p>
     * This is Zhang's rule R4, discriminating undirectedPaths.
     */
    private void ruleR4A(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node b : nodes) {

            //potential A and C candidate pairs are only those
            // that look like this:   A&lt;-*Bo-*C
            List<Node> possA = graph.getNodesOutTo(b, Endpoint.ARROW);
            List<Node> possC = graph.getNodesInTo(b, Endpoint.CIRCLE);

            for (Node a : possA) {
                for (Node c : possC) {
                    if (!graph.isParentOf(a, c)) {
                        continue;
                    }

                    if (graph.getEndpoint(b, c) != Endpoint.ARROW) {
                        continue;
                    }

                    LinkedList<Node> reachable = new LinkedList<>();
                    reachable.add(a);

                    if (this.verbose) {
                        System.out.println("Found CPDAG " + a + " " + b + " " + c);
                        reachablePathFind(a, b, c, reachable, graph);
                    }
                }
            }
        }
    }

    /**
     * a method to search "back from a" to find a DDP. It is called with a reachability list (first consisting only of
     * a). This is breadth-first, utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. The body of
     * a DDP consists of colliders that are parents of c.
     */
    private void reachablePathFind(Node a, Node b, Node c,
                                   LinkedList<Node> reachable, Graph graph) {

        Set<Node> cParents = new HashSet<>(graph.getParents(c));

        // Needed to avoid cycles in failure case.
        Set<Node> visited = new HashSet<>();
        visited.add(b);
        visited.add(c);

        Node e = reachable.getFirst();
        int distance = 0;

        // We don't want to include a,b,or c on the path, so they are added to
        // the "visited" set.  b and c are added explicitly here; a will be
        // added in the first while iteration.
        while (reachable.size() > 0) {
            Node x = reachable.removeFirst();
            visited.add(x);

            if (e == x) {
                e = x;
                distance++;

                int _maxPathLength = this.maxPathLength == -1 ? 1000 : this.maxPathLength;

                if (distance > 0 && distance > _maxPathLength) {
                    continue;
                }
            }

            // Possible DDP path endpoints.
            List<Node> pathExtensions = graph.getNodesInTo(x, Endpoint.ARROW);
            pathExtensions.removeAll(visited);

            for (Node d : pathExtensions) {
                // If d is reachable and not adjacent to c, its a DDP
                // endpoint, so do DDP orientation. Otherwise, if d &lt;-> c,
                // add d to the list of reachable nodes.
                if (!graph.isAdjacentTo(d, c)) {
                    // Check whether <a, b, c> should be reoriented given
                    // that d is not adjacent to c; if so, orient and stop.
                    doDdpOrientation(d, a, b, c, graph);
                    return;
                } else if (cParents.contains(d)) {
                    if (graph.getEndpoint(x, d) == Endpoint.ARROW) {
                        reachable.add(d);

                        // RFCI: only record the next node of the first (shortest) occurrence
                    }
                }
            }
        }
    }

    /**
     * Orients the edges inside the definte discriminating path triangle. Takes the left endpoint, and a,b,c as
     * arguments.
     */
    private void doDdpOrientation(Node d, Node a, Node b, Node c, Graph graph) {
        List<Node> sepset = getSepset(d, c);

        if (sepset == null) return;

        if (sepset.contains(b)) {
            graph.setEndpoint(c, b, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.TAIL);
            if (this.verbose) {
                this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Definite discriminating path d = " + d, graph.getEdge(b, c)));
                System.out.println(SearchLogUtils.edgeOrientedMsg("Definite discriminating path d = " + d, graph.getEdge(b, c)));
            }

        } else {
            if (!isArrowpointAllowed(a, b, graph)) {
                return;
            }

            if (!isArrowpointAllowed(c, b, graph)) {
                return;
            }

            graph.setEndpoint(a, b, Endpoint.ARROW);
            graph.setEndpoint(c, b, Endpoint.ARROW);
            this.orientSimilarPairs(graph, this.getKnowledge(), a, b, Endpoint.ARROW);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.ARROW);
            this.logger.log("colliderOrientations", SearchLogUtils.colliderOrientedMsg("Definite discriminating path.. d = " + d, a, b, c));
        }
        this.changeFlag = true;
    }

    /**
     * The triangles that must be oriented this way (won't be done by another rule) all look like
     * the ones below, where the dots are a collider path from L to A with each node on the path
     * (except L) a parent of C.
     * <pre>
     *          B
     *         xo           x is either an arrowhead or a circle
     *        /  \
     *       v    v
     * L....A --&gt; C
     * </pre>
     * <p>
     * This is Zhang's rule R4, discriminating undirectedPaths.
     */
    public void ruleR4B(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node b : nodes) {

            //potential A and C candidate pairs are only those
            // that look like this:   A&lt;-*Bo-*C
            List<Node> possA = graph.getNodesOutTo(b, Endpoint.ARROW);
            List<Node> possC = graph.getNodesInTo(b, Endpoint.CIRCLE);

            for (Node a : possA) {
                for (Node c : possC) {
                    if (!graph.isParentOf(a, c)) {
                        continue;
                    }

                    if (graph.getEndpoint(b, c) != Endpoint.ARROW) {
                        continue;
                    }

                    ddpOrient(a, b, c, graph);
                }
            }
        }
    }

    /**
     * a method to search "back from a" to find a DDP. It is called with a reachability list (first consisting only of
     * a). This is breadth-first, utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. The body of
     * a DDP consists of colliders that are parents of c.
     */
    public void ddpOrient(Node a, Node b, Node c, Graph graph) {
        Queue<Node> Q = new ArrayDeque<>();
        Set<Node> V = new HashSet<>();

        Node e = null;
        int distance = 0;

        Map<Node, Node> previous = new HashMap<>();

        List<Node> cParents = graph.getParents(c);

        Q.offer(a);
        V.add(a);
        V.add(b);
        previous.put(a, b);

        while (!Q.isEmpty()) {
            Node t = Q.poll();

            if (e == null || e == t) {
                e = t;
                distance++;
                if (distance > 0 && distance > (this.maxPathLength == -1 ? 1000 : this.maxPathLength)) return;
            }

            List<Node> nodesInTo = graph.getNodesInTo(t, Endpoint.ARROW);

            for (Node d : nodesInTo) {
                if (V.contains(d)) continue;

                previous.put(d, t);
                Node p = previous.get(t);

                if (!graph.isDefCollider(d, t, p)) {
                    continue;
                }

                previous.put(d, t);

                if (!graph.isAdjacentTo(d, c)) {
                    if (doDdpOrientation(d, a, b, c, previous, graph)) {
                        return;
                    }
                }

                if (cParents.contains(d)) {
                    Q.offer(d);
                    V.add(d);
                }
            }
        }
    }

    /**
     * Orients the edges inside the definte discriminating path triangle. Takes the left endpoint, and a,b,c as
     * arguments.
     */
    private boolean doDdpOrientation(Node d, Node a, Node b, Node c, Map<Node, Node> previous, Graph graph) {
        if (graph.isAdjacentTo(d, c)) {
            throw new IllegalArgumentException();
        }

        List<Node> path = getPath(d, previous);

        boolean ind = getSepsets().isIndependent(d, c, path);

        List<Node> path2 = new ArrayList<>(path);

        path2.remove(b);

        boolean ind2 = getSepsets().isIndependent(d, c, path2);

        if (!ind && !ind2) {
            List<Node> sepset = getSepsets().getSepset(d, c);

            if (this.verbose) {
                System.out.println("Sepset for d = " + d + " and c = " + c + " = " + sepset);
            }

            if (sepset == null) {
                if (this.verbose) {
                    System.out.println("Must be a sepset: " + d + " and " + c + "; they're non-adjacent.");
                }
                return false;
            }

            ind = sepset.contains(b);
        }

//        printDdp(d, path, a, b, c, graph);

        if (ind) {
//            if (sepset.contains(b)) {
            graph.setEndpoint(c, b, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.TAIL);
            if (this.verbose) {
                this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Definite discriminating path d = " + d, graph.getEdge(b, c)));
                System.out.println(SearchLogUtils.edgeOrientedMsg("Definite discriminating path d = " + d, graph.getEdge(b, c)));
            }

        } else {
            if (!isArrowpointAllowed(a, b, graph)) {
                return false;
            }

            if (!isArrowpointAllowed(c, b, graph)) {
                return false;
            }

            graph.setEndpoint(a, b, Endpoint.ARROW);
            graph.setEndpoint(c, b, Endpoint.ARROW);
            this.orientSimilarPairs(graph, this.getKnowledge(), a, b, Endpoint.ARROW);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.ARROW);
            if (this.verbose) {
                this.logger.log("impliedOrientations", SearchLogUtils.colliderOrientedMsg("Definite discriminating path.. d = " + d, a, b, c));
                System.out.println(SearchLogUtils.colliderOrientedMsg("Definite discriminating path.. d = " + d, a, b, c));
            }

        }
        this.changeFlag = true;
        return true;
    }

    private void printDdp(Node d, List<Node> path, Node a, Node b, Node c, Graph graph) {
        List<Node> nodes = new ArrayList<>();
        nodes.add(d);
        nodes.addAll(path);
        nodes.add(a);
        nodes.add(b);
        nodes.add(c);

        if (this.verbose) {
            System.out.println("DDP subgraph = " + graph.subgraph(nodes));
        }
    }

    private List<Node> getPath(Node c, Map<Node, Node> previous) {
        List<Node> l = new ArrayList<>();

        Node p = c;

        do {
            p = previous.get(p);

            if (p != null) {
                l.add(p);
            }
        } while (p != null);

        return l;
    }

    /**
     * Implements Zhang's rule R5, orient circle undirectedPaths: for any Ao-oB, if there is an uncovered circle path u =
     * [A,C,...,D,B] such that A,D nonadjacent and B,C nonadjacent, then A---B and orient every edge on u undirected.
     */
    public void ruleR5(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node a : nodes) {
            List<Node> adjacents = graph.getNodesInTo(a, Endpoint.CIRCLE);

            for (Node b : adjacents) {
                if (!(graph.getEndpoint(a, b) == Endpoint.CIRCLE)) continue;
                // We know Ao-oB.

                List<List<Node>> ucCirclePaths = getUcCirclePaths(a, b, graph);

                for (List<Node> u : ucCirclePaths) {
                    if (u.size() < 3) continue;

                    Node c = u.get(1);
                    Node d = u.get(u.size() - 2);

                    if (graph.isAdjacentTo(a, d)) continue;
                    if (graph.isAdjacentTo(b, c)) continue;
                    // We know u is as required: R5 applies!

                    this.logger.log("colliderOrientations", SearchLogUtils.edgeOrientedMsg("Orient circle path", graph.getEdge(a, b)));

                    graph.setEndpoint(a, b, Endpoint.TAIL);
                    this.orientSimilarPairs(graph, this.getKnowledge(), a, b, Endpoint.TAIL);
                    graph.setEndpoint(b, a, Endpoint.TAIL);
                    this.orientSimilarPairs(graph, this.getKnowledge(), b, a, Endpoint.TAIL);
                    orientTailPath(u, graph);
                    this.changeFlag = true;
                }
            }
        }
    }

    /**
     * Implements Zhang's rules R6 and R7, applies them over the graph once. Orient single tails. R6: If A---Bo-*C then
     * A---B--*C. R7: If A--oBo-*C and A,C nonadjacent, then A--oB--*C
     */
    public void ruleR6R7(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node b : nodes) {
            List<Node> adjacents = graph.getAdjacentNodes(b);

            if (adjacents.size() < 2) continue;

            ChoiceGenerator cg = new ChoiceGenerator(adjacents.size(), 2);

            for (int[] choice = cg.next(); choice != null; choice = cg.next()) {
                Node a = adjacents.get(choice[0]);
                Node c = adjacents.get(choice[1]);

                if (graph.isAdjacentTo(a, c)) continue;

                if (!(graph.getEndpoint(b, a) == Endpoint.TAIL)) continue;
                if (!(graph.getEndpoint(c, b) == Endpoint.CIRCLE)) continue;
                // We know A--*Bo-*C.

                if (graph.getEndpoint(a, b) == Endpoint.TAIL) {

                    // We know A---Bo-*C: R6 applies!
                    graph.setEndpoint(c, b, Endpoint.TAIL);
                    this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.TAIL);
                    this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Single tails (tail)", graph.getEdge(c, b)));

                    this.changeFlag = true;
                }

                if (graph.getEndpoint(a, b) == Endpoint.CIRCLE) {
//                    if (graph.isAdjacentTo(a, c)) continue;

                    this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Single tails (tail)", graph.getEdge(c, b)));

                    // We know A--oBo-*C and A,C nonadjacent: R7 applies!
                    graph.setEndpoint(c, b, Endpoint.TAIL);
                    this.orientSimilarPairs(graph, this.getKnowledge(), c, b, Endpoint.TAIL);
                    this.changeFlag = true;
                }

            }
        }
    }

    /**
     * Implements Zhang's rules R8, R9, R10, applies them over the graph once. Orient arrow tails. I.e., tries R8, R9,
     * and R10 in that sequence on each Ao-&gt;C in the graph.
     */
    public void rulesR8R9R10(Graph graph) {
        List<Node> nodes = graph.getNodes();

        for (Node c : nodes) {
            List<Node> intoCArrows = graph.getNodesInTo(c, Endpoint.ARROW);

            for (Node a : intoCArrows) {
                if (!(graph.getEndpoint(c, a) == Endpoint.CIRCLE)) continue;
                // We know Ao->C.

                // Try each of R8, R9, R10 in that order, stopping ASAP.
                if (!ruleR8(a, c, graph)) {
                    boolean b = ruleR9(a, c, graph);

                    if (!b) {
                        ruleR10(a, c, graph);
                    }
                }
            }
        }

    }

    /**
     * Orients every edge on a path as undirected (i.e. A---B).
     * <p>
     * DOES NOT CHECK IF SUCH EDGES ACTUALLY EXIST: MAY DO WEIRD THINGS IF PASSED AN ARBITRARY LIST OF NODES THAT IS NOT
     * A PATH.
     *
     * @param path The path to orient as all tails.
     */
    private void orientTailPath(List<Node> path, Graph graph) {
        for (int i = 0; i < path.size() - 1; i++) {
            Node n1 = path.get(i);
            Node n2 = path.get(i + 1);

            graph.setEndpoint(n1, n2, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), n1, n2, Endpoint.TAIL);
            graph.setEndpoint(n2, n1, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), n2, n1, Endpoint.TAIL);
            this.changeFlag = true;

            this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("Orient circle undirectedPaths", graph.getEdge(n1, n2)));
        }
    }

    /**
     * Gets a list of every uncovered partially directed path between two nodes in the graph.
     * <p>
     * Probably extremely slow.
     *
     * @param n1 The beginning node of the undirectedPaths.
     * @param n2 The ending node of the undirectedPaths.
     * @return A list of uncovered partially directed undirectedPaths from n1 to n2.
     */
    private List<List<Node>> getUcPdPaths(Node n1, Node n2, Graph graph) {
        List<List<Node>> ucPdPaths = new LinkedList<>();

        LinkedList<Node> soFar = new LinkedList<>();
        soFar.add(n1);

        List<Node> adjacencies = graph.getAdjacentNodes(n1);
        for (Node curr : adjacencies) {
            getUcPdPsHelper(curr, soFar, n2, ucPdPaths, graph);
        }

        return ucPdPaths;
    }

    /**
     * Used in getUcPdPaths(n1,n2) to perform a breadth-first search on the graph.
     * <p>
     * ASSUMES soFar CONTAINS AT LEAST ONE NODE!
     * <p>
     * Probably extremely slow.
     *
     * @param curr      The getModel node to test for addition.
     * @param soFar     The getModel partially built-up path.
     * @param end       The node to finish the undirectedPaths at.
     * @param ucPdPaths The getModel list of uncovered p.d. undirectedPaths.
     */
    private void getUcPdPsHelper(Node curr, List<Node> soFar, Node end,
                                 List<List<Node>> ucPdPaths, Graph graph) {

        if (soFar.contains(curr)) return;

        Node prev = soFar.get(soFar.size() - 1);
        if (graph.getEndpoint(prev, curr) == Endpoint.TAIL ||
                graph.getEndpoint(curr, prev) == Endpoint.ARROW) {
            return; // Adding curr would make soFar not p.d.
        } else if (soFar.size() >= 2) {
            Node prev2 = soFar.get(soFar.size() - 2);
            if (graph.isAdjacentTo(prev2, curr)) {
                return; // Adding curr would make soFar not uncovered.
            }
        }

        soFar.add(curr); // Adding curr is OK, so let's do it.

        if (curr.equals(end)) {
            // We've reached the goal! Save soFar as a path.
            ucPdPaths.add(new LinkedList<>(soFar));
        } else {
            // Otherwise, try each node adjacent to the getModel one.
            List<Node> adjacents = graph.getAdjacentNodes(curr);
            for (Node next : adjacents) {
                getUcPdPsHelper(next, soFar, end, ucPdPaths, graph);
            }
        }

        soFar.remove(soFar.get(soFar.size() - 1)); // For other recursive calls.
    }

    /**
     * Gets a list of every uncovered circle path between two nodes in the graph by iterating through the uncovered
     * partially directed undirectedPaths and only keeping the circle undirectedPaths.
     * <p>
     * Probably extremely slow.
     *
     * @param n1 The beginning node of the undirectedPaths.
     * @param n2 The ending node of the undirectedPaths.
     * @return A list of uncovered circle undirectedPaths between n1 and n2.
     */
    private List<List<Node>> getUcCirclePaths(Node n1, Node n2, Graph graph) {
        List<List<Node>> ucCirclePaths = new LinkedList<>();
        List<List<Node>> ucPdPaths = getUcPdPaths(n1, n2, graph);

        for (List<Node> path : ucPdPaths) {
            for (int i = 0; i < path.size() - 1; i++) {
                Node j = path.get(i);
                Node sj = path.get(i + 1);

                if (!(graph.getEndpoint(j, sj) == Endpoint.CIRCLE)) break;
                if (!(graph.getEndpoint(sj, j) == Endpoint.CIRCLE)) break;
                // This edge is OK, it's all circles.

                if (i == path.size() - 2) {
                    // We're at the last edge, so this is a circle path.
                    ucCirclePaths.add(path);
                }
            }
        }

        return ucCirclePaths;
    }

    /**
     * Tries to apply Zhang's rule R8 to a pair of nodes A and C which are assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R8: If Ao->C and A-->B-->C or A--oB-->C, then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R8 was successfully applied.
     */
    private boolean ruleR8(Node a, Node c, Graph graph) {
        List<Node> intoCArrows = graph.getNodesInTo(c, Endpoint.ARROW);

        for (Node b : intoCArrows) {
            // We have B*-&gt;C.
            if (!graph.isAdjacentTo(a, b)) continue;
            if (!graph.isAdjacentTo(b, c)) continue;

            // We have A*-*B*-&gt;C.
            if (!(graph.getEndpoint(b, a) == Endpoint.TAIL)) continue;
            if (!(graph.getEndpoint(c, b) == Endpoint.TAIL)) continue;
            // We have A--*B-->C.

            if (graph.getEndpoint(a, b) == Endpoint.TAIL) continue;
            // We have A-->B-->C or A--oB-->C: R8 applies!

            this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("R8", graph.getEdge(c, a)));

            graph.setEndpoint(c, a, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, a, Endpoint.TAIL);
            this.changeFlag = true;
            return true;
        }

        return false;
    }

    /**
     * Tries to apply Zhang's rule R9 to a pair of nodes A and C which are assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R9: If Ao->C and there is an uncovered p.d. path u=<A,B,..,C> such that C,B nonadjacent, then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R9 was succesfully applied.
     */
    private boolean ruleR9(Node a, Node c, Graph graph) {
        List<List<Node>> ucPdPsToC = getUcPdPaths(a, c, graph);

        for (List<Node> u : ucPdPsToC) {
            Node b = u.get(1);
            if (graph.isAdjacentTo(b, c)) continue;
            if (b == c) continue;
            // We know u is as required: R9 applies!

            this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("R9", graph.getEdge(c, a)));

            graph.setEndpoint(c, a, Endpoint.TAIL);
            this.orientSimilarPairs(graph, this.getKnowledge(), c, a, Endpoint.TAIL);
            this.changeFlag = true;
            return true;
        }

        return false;
    }

    /**
     * Tries to apply Zhang's rule R10 to a pair of nodes A and C which are assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R10: If Ao->C, B-->C&lt;--D, there is an uncovered p.d. path u1=<A,M,...,B> and an uncovered p.d. path
     * u2=<A,N,...,D> with M != N and M,N nonadjacent then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R10 was successfully applied.
     */
    private boolean ruleR10(Node a, Node c, Graph graph) {
        List<Node> intoCArrows = graph.getNodesInTo(c, Endpoint.ARROW);

        for (Node b : intoCArrows) {
            if (b == a) continue;

            if (!(graph.getEndpoint(c, b) == Endpoint.TAIL)) continue;
            // We know Ao->C and B-->C.

            for (Node d : intoCArrows) {
                if (d == a || d == b) continue;

                if (!(graph.getEndpoint(d, c) == Endpoint.TAIL)) continue;
                // We know Ao->C and B-->C<--D.

                List<List<Node>> ucPdPsToB = getUcPdPaths(a, b, graph);
                List<List<Node>> ucPdPsToD = getUcPdPaths(a, d, graph);
                for (List<Node> u1 : ucPdPsToB) {
                    Node m = u1.get(1);
                    for (List<Node> u2 : ucPdPsToD) {
                        Node n = u2.get(1);

                        if (m.equals(n)) continue;
                        if (graph.isAdjacentTo(m, n)) continue;
                        // We know B,D,u1,u2 as required: R10 applies!

                        this.logger.log("impliedOrientations", SearchLogUtils.edgeOrientedMsg("R10", graph.getEdge(c, a)));

                        graph.setEndpoint(c, a, Endpoint.TAIL);
                        this.changeFlag = true;
                        this.orientSimilarPairs(graph, this.getKnowledge(), c, a, Endpoint.TAIL);
                        return true;
                    }
                }
            }
        }

        return false;
    }

    /**
     * Orients according to background knowledge
     */
    private void fciOrientbk(IKnowledge bk, Graph graph, List<Node> variables) {
        this.logger.log("info", "Starting BK Orientation.");

        for (Iterator<KnowledgeEdge> it =
             bk.forbiddenEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in the graph.
            Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
            Node to = SearchGraphUtils.translate(edge.getTo(), variables);


            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            // Orient to*-&gt;from
            graph.setEndpoint(to, from, Endpoint.ARROW);
            graph.setEndpoint(from, to, Endpoint.CIRCLE);
            this.changeFlag = true;
            this.logger.log("knowledgeOrientation", SearchLogUtils.edgeOrientedMsg("Knowledge", graph.getEdge(from, to)));
        }

        for (Iterator<KnowledgeEdge> it =
             bk.requiredEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in this graph
            Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
            Node to = SearchGraphUtils.translate(edge.getTo(), variables);

            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            graph.setEndpoint(to, from, Endpoint.TAIL);
            graph.setEndpoint(from, to, Endpoint.ARROW);
            this.changeFlag = true;
            this.logger.log("knowledgeOrientation", SearchLogUtils.edgeOrientedMsg("Knowledge", graph.getEdge(from, to)));
        }

        this.logger.log("info", "Finishing BK Orientation.");
    }


    /**
     * Helper method. Appears to check if an arrowpoint is permitted by background knowledge.
     *
     * @param x The possible other node.
     * @param y The possible point node.
     * @return Whether the arrowpoint is allowed.
     */
    private boolean isArrowpointAllowed(Node x, Node y, Graph graph) {
        if (graph.getEndpoint(x, y) == Endpoint.ARROW) {
            return true;
        }

        if (graph.getEndpoint(x, y) == Endpoint.TAIL) {
            return false;
        }

        if (graph.getEndpoint(y, x) == Endpoint.ARROW) {
            return true; // changed by DMalinsky 02/11/2019, if already an arrowhead at x, arrowpoint is allowed at y
//            if (!knowledge.isForbidden(x.getName(), y.getName())) return true;
        }

        if (graph.getEndpoint(y, x) == Endpoint.TAIL) {
            if (!this.knowledge.isForbidden(x.getName(), y.getName())) return true;
        }

        return graph.getEndpoint(y, x) == Endpoint.CIRCLE;
    }

    public boolean isPossibleDsepSearchDone() {
        return this.possibleDsepSearchDone;
    }

    public void setPossibleDsepSearchDone(boolean possibleDsepSearchDone) {
        this.possibleDsepSearchDone = possibleDsepSearchDone;
    }

    /**
     * @return the maximum length of any discriminating path, or -1 of unlimited.
     */
    public int getMaxPathLength() {
        return this.maxPathLength;
    }

    /**
     * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
     */
    public void setMaxPathLength(int maxPathLength) {
        if (maxPathLength < -1) {
            throw new IllegalArgumentException("Max path length must be -1 (unlimited) or >= 0: " + maxPathLength);
        }

        this.maxPathLength = maxPathLength;
    }

    /**
     * True iff verbose output should be printed.
     */
    public boolean isVerbose() {
        return this.verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void setTruePag(Graph truePag) {
        this.truePag = truePag;
    }

    /**
     * The true PAG if available. Can be null.
     */
    public Graph getTruePag() {
        return this.truePag;
    }

    public void setChangeFlag(boolean changeFlag) {
        this.changeFlag = changeFlag;
    }

    /**
     * change flag for repeat rules
     */
    public boolean isChangeFlag() {
        return this.changeFlag;
    }

    private void orientSimilarPairs(Graph graph, IKnowledge knowledge, Node x, Node y, Endpoint mark) {
        if (x.getName().equals("time") || y.getName().equals("time")) {
            return;
        }
        System.out.println("Entering orient similar pairs method for x and y: " + x + ", " + y);
        int ntiers = knowledge.getNumTiers();
        int indx_tier = knowledge.isInWhichTier(x);
        int indy_tier = knowledge.isInWhichTier(y);
        int tier_diff = Math.max(indx_tier, indy_tier) - Math.min(indx_tier, indy_tier);
        int indx_comp = -1;
        int indy_comp = -1;
        List tier_x = knowledge.getTier(indx_tier);
//        Collections.sort(tier_x);
        List tier_y = knowledge.getTier(indy_tier);
//        Collections.sort(tier_y);

        int i;
        for (i = 0; i < tier_x.size(); ++i) {
            if (getNameNoLag(x.getName()).equals(getNameNoLag(tier_x.get(i)))) {
                indx_comp = i;
                break;
            }
        }

        for (i = 0; i < tier_y.size(); ++i) {
            if (getNameNoLag(y.getName()).equals(getNameNoLag(tier_y.get(i)))) {
                indy_comp = i;
                break;
            }
        }

        if (indx_comp == -1) System.out.println("WARNING: indx_comp = -1!!!! ");
        if (indy_comp == -1) System.out.println("WARNING: indy_comp = -1!!!! ");

        for (i = 0; i < ntiers - tier_diff; ++i) {
            if (knowledge.getTier(i).size() == 1) continue;
            String A;
            Node x1;
            String B;
            Node y1;
            if (indx_tier >= indy_tier) {
                List tmp_tier1 = knowledge.getTier(i + tier_diff);
//                Collections.sort(tmp_tier1);
                List tmp_tier2 = knowledge.getTier(i);
//                Collections.sort(tmp_tier2);
                A = (String) tmp_tier1.get(indx_comp);
                B = (String) tmp_tier2.get(indy_comp);
                if (A.equals(B)) continue;
                if (A.equals(tier_x.get(indx_comp)) && B.equals(tier_y.get(indy_comp))) continue;
                if (B.equals(tier_x.get(indx_comp)) && A.equals(tier_y.get(indy_comp))) continue;
                x1 = this.independenceTest.getVariable(A);
                y1 = this.independenceTest.getVariable(B);

                if (graph.isAdjacentTo(x1, y1) && graph.getEndpoint(x1, y1) == Endpoint.CIRCLE) {
                    System.out.print("Orient edge " + graph.getEdge(x1, y1).toString());
                    graph.setEndpoint(x1, y1, mark);
                    System.out.println(" by structure knowledge as: " + graph.getEdge(x1, y1).toString());
                }
            } else {
//                System.out.println("############## WARNING (orientSimilarPairs): did not catch x,y pair " + x + ", " + y);
            }
        }

    }


    public String getNameNoLag(Object obj) {
        String tempS = obj.toString();
        if (tempS.indexOf(':') == -1) {
            return tempS;
        } else return tempS.substring(0, tempS.indexOf(':'));
    }


}

