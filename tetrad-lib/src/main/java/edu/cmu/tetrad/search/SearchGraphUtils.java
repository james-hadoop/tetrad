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
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;

import java.util.*;

/**
 * Graph utilities for search algorithm. Lots of orientation method, for
 * instance.
 *
 * @author Joseph Ramsey
 */
public final class SearchGraphUtils {

    /**
     * Orients according to background knowledge.
     */
    public static void pcOrientbk(IKnowledge bk, Graph graph, List<Node> nodes) {
        log("Staring BK Orientation.");
        for (Iterator<KnowledgeEdge> it = bk.forbiddenEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in the graph.
            Node from = translate(edge.getFrom(), nodes);
            Node to = translate(edge.getTo(), nodes);

            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            // Orient to-->from
            graph.removeEdge(from, to);
            graph.addDirectedEdge(to, from);
//            graph.setEndpoint(from, to, Endpoint.TAIL);
//            graph.setEndpoint(to, from, Endpoint.ARROW);

            log(SearchLogUtils.edgeOrientedMsg("IKnowledge", graph.getEdge(to, from)));
        }

        for (Iterator<KnowledgeEdge> it = bk.requiredEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in this graph
            Node from = translate(edge.getFrom(), nodes);
            Node to = translate(edge.getTo(), nodes);

            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            // Orient from-->to
            graph.removeEdges(from, to);
            graph.addDirectedEdge(from, to);

//            graph.setEndpoint(to, from, Endpoint.TAIL);
//            graph.setEndpoint(from, to, Endpoint.ARROW);
            log(SearchLogUtils.edgeOrientedMsg("IKnowledge", graph.getEdge(from, to)));
        }

        log("Finishing BK Orientation.");
    }

    /**
     * Performs step C of the algorithm, as indicated on page xxx of CPS, with
     * the modification that X--W--Y is oriented as X-->W<--Y if W is
     * *determined by* the sepset of (X, Y), rather than W just being *in* the
     * sepset of (X, Y).
     */
    public static void pcdOrientC(SepsetMap set, IndependenceTest test,
                                  IKnowledge knowledge, Graph graph) {
        log("Starting Collider Orientation:");

        List<Node> nodes = graph.getNodes();

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

                List<Node> sepset = sepset(graph, x, z, new HashSet<Node>(), new HashSet<Node>(),
                        -1, test);
                //set.get(x, z);

                if (sepset == null) {
                    continue;
                }

                if (sepset.contains(y)) {
                    continue;
                }

                List<Node> augmentedSet = new LinkedList<>(sepset);

                if (!augmentedSet.contains(y)) {
                    augmentedSet.add(y);
                }

                if (test.determines(sepset, x)) {
//                    System.out.println(SearchLogUtils.determinismDetected(sepset, x));
                    continue;
                }

                if (test.determines(sepset, z)) {
//                    System.out.println(SearchLogUtils.determinismDetected(sepset, z));
                    continue;
                }

                if (test.determines(augmentedSet, x)) {
//                    System.out.println(SearchLogUtils.determinismDetected(augmentedSet, x));
                    continue;
                }

                if (test.determines(augmentedSet, z)) {
//                    System.out.println(SearchLogUtils.determinismDetected(augmentedSet, z));
                    continue;
                }

                if (!isArrowpointAllowed(x, y, knowledge)
                        || !isArrowpointAllowed(z, y, knowledge)) {
                    continue;
                }

                graph.setEndpoint(x, y, Endpoint.ARROW);
                graph.setEndpoint(z, y, Endpoint.ARROW);

                System.out.println(SearchLogUtils.colliderOrientedMsg(x, y, z) + " sepset = " + sepset);
                log(SearchLogUtils.colliderOrientedMsg(x, y, z));
            }
        }

        log("Finishing Collider Orientation.");
    }

    private static void log(String s) {
        System.out.println(s);
    }

    private static List<Node> sepset(Graph graph, Node a, Node c, Set<Node> containing, Set<Node> notContaining, int depth,
                                     IndependenceTest independenceTest) {
        List<Node> adj = graph.getAdjacentNodes(a);
        adj.addAll(graph.getAdjacentNodes(c));
        adj.remove(c);
        adj.remove(a);

        for (int d = 0; d <= Math.min((depth == -1 ? 1000 : depth), Math.max(adj.size(), adj.size())); d++) {
            if (d <= adj.size()) {
                ChoiceGenerator gen = new ChoiceGenerator(adj.size(), d);
                int[] choice;

                WHILE:
                while ((choice = gen.next()) != null) {
                    Set<Node> v2 = GraphUtils.asSet(choice, adj);
                    v2.addAll(containing);
                    v2.removeAll(notContaining);
                    v2.remove(a);
                    v2.remove(c);

//                    if (isForbidden(a, c, new ArrayList<>(v2)))
                    independenceTest.isIndependent(a, c, new ArrayList<>(v2));
                    double p2 = independenceTest.getScore();

                    if (p2 < 0) {
                        return new ArrayList<>(v2);
                    }
                }
            }
        }

        return null;
    }

    //    /**
//     * Performs step D of the algorithm, as indicated on page xxx of CPS. This
//     * method should be called again if it returns true.
//     *
//     * <pre>
//     * Meek-Orient(G, t)
//     * 1.while orientations can be made, for arbitrary a, b, c, and d:
//     * 2.    If a --> b, b --> c, a not in adj(c), and Is-Noncollider(a, b, c) then orient b --> c.
//     * 3.    If a --> b, b --> c, a --- c, then orient a --> c.
//     * 4.    If a --- b, a --- c, a --- d, c --> b, d --> b, then orient a --> b.
//     * 5.    If a --> b, b in adj(d) a in adj(c), a --- d, b --> c, c --> d, then orient a --> d.
//     * </pre>
//     */
//    public static void orientUsingMeekRules(IKnowledge knowledge, Graph graph) {
//        LogUtils.getInstance().info("Starting Orientation Step D.");
//        boolean changed;
//
//        do {
//            changed = meekR1(graph, knowledge) || meekR2(graph, knowledge) ||
//                    meekR3(graph, knowledge) || meekR4(graph, knowledge);
//        } while (changed);
//
//        LogUtils.getInstance().info("Finishing Orientation Step D.");
//    }

    /**
     * Orients using Meek rules, double checking noncolliders locally.
     */
    public static void orientUsingMeekRulesLocally(IKnowledge knowledge,
                                                   Graph graph, IndependenceTest test, int depth) {
        log("Starting Orientation Step D.");
        boolean changed;

        do {
            changed = meekR1Locally(graph, knowledge, test, depth)
                    || meekR2(graph, knowledge) || meekR3(graph, knowledge)
                    || meekR4(graph, knowledge);
        } while (changed);

        log("Finishing Orientation Step D.");
    }

    public static void orientUsingMeekRulesLocally2(IKnowledge knowledge,
                                                    Graph graph, IndependenceTest test, int depth) {
        log("Starting Orientation Step D.");
        boolean changed;

        do {
            changed = meekR1Locally2(graph, knowledge, test, depth)
                    || meekR2(graph, knowledge) || meekR3(graph, knowledge)
                    || meekR4(graph, knowledge);
        } while (changed);

        log("Finishing Orientation Step D.");
    }

    /**
     * Step C of PC; orients colliders using specified sepset. That is, orients
     * x *-* y *-* z as x *-> y <-* z just in case y is in Sepset({x, z}).
     */
    public static List<Triple> orientCollidersUsingSepsets(SepsetMap set, IKnowledge knowledge, Graph graph, boolean verbose,
                                                           boolean enforcePattern) {
        log("Starting Collider Orientation:");
        List<Triple> colliders = new ArrayList<>();

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

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(a, c)) {
                    continue;
                }

                List<Node> sepset = set.get(a, c);

                //I think the null check needs to be here --AJ
                if (sepset != null && !sepset.contains(b)
                        && isArrowpointAllowed(a, b, knowledge)
                        && isArrowpointAllowed(c, b, knowledge)) {
                    if (verbose) {
                        System.out.println("Collider orientation <" + a + ", " + b + ", " + c + "> sepset = " + sepset);
                    }

                    if (enforcePattern) {
                        if (graph.getEndpoint(b, a) == Endpoint.ARROW || graph.getEndpoint(b, c) == Endpoint.ARROW) {
                            continue;
                        }
                    }

//                    graph.setEndpoint(a, b, Endpoint.ARROW);
//                    graph.setEndpoint(c, b, Endpoint.ARROW);
                    graph.removeEdge(a, b);
                    graph.removeEdge(c, b);

                    graph.addDirectedEdge(a, b);
                    graph.addDirectedEdge(c, b);

                    colliders.add(new Triple(a, b, c));
                    log(SearchLogUtils.colliderOrientedMsg(a, b, c, sepset));
                }
            }
        }

        log("Finishing Collider Orientation.");

        return colliders;
    }

    private static List<Node> union(List<Node> nodes, Node a) {
        List<Node> union = new ArrayList<>(nodes);
        union.add(a);
        return union;
    }

    public static void orientCollidersUsingSepsets(SepsetProducer sepset, IKnowledge knowledge, Graph graph, boolean verbose) {
        log("Starting Collider Orientation:");

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

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(a, c)) {
                    continue;
                }

                // Skip triple already oriented.
                if (graph.getEdge(a, b).pointsTowards(b) && graph.getEdge(b, c).pointsTowards(b)) {
                    continue;
                }

                final List<Node> sepset1 = sepset.getSepset(a, c);

                if (!sepset1.contains(b) && isArrowpointAllowed(a, b, knowledge)
                        && isArrowpointAllowed(c, b, knowledge)) {
                    graph.setEndpoint(a, b, Endpoint.ARROW);
                    graph.setEndpoint(c, b, Endpoint.ARROW);

                    log(SearchLogUtils.colliderOrientedMsg(a, b, c, sepset1));

                    if (verbose) {
                        System.out.println(SearchLogUtils.colliderOrientedMsg(a, b, c, sepset1));
                    }
                }
            }
        }

        log("Finishing Collider Orientation.");
    }

    //use this for oritentation with an initial graph if using null trick for unconditional independence
    //AJ
    public static List<Triple> orientCollidersUsingSepsets(SepsetMap set, IKnowledge knowledge, Graph graph, Graph initialGraph, boolean verbose) {
        log("Starting Collider Orientation:");
        List<Triple> colliders = new ArrayList<>();

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

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(a, c)) {
                    continue;
                }

                //for vanilla pc, I think a check if already oriented might need to be here -AJ
                // Skip triples with parents not adjacent in initialGraph
                // may need a similar check for knowledge... -AJ
                if (initialGraph != null && !initialGraph.isAdjacentTo(a, c)) {
                    continue;
                }

                List<Node> sepset = set.get(a, c);

                //Null check needs to be here if sepsets.setReturnEmptyIfNotSet(false)--AJ
                if (sepset != null && !sepset.contains(b)
                        && isArrowpointAllowed(a, b, knowledge)
                        && isArrowpointAllowed(c, b, knowledge)) {
                    if (verbose) {
                        System.out.println("Collider orientation <" + a + ", " + b + ", " + c + "> sepset = " + sepset);
                    }

                    graph.setEndpoint(a, b, Endpoint.ARROW);
                    graph.setEndpoint(c, b, Endpoint.ARROW);
                    colliders.add(new Triple(a, b, c));
                    log(SearchLogUtils.colliderOrientedMsg(a, b, c, sepset));
                }
            }
        }

        log("Finishing Collider Orientation.");

        return colliders;
    }

    private static boolean createsBidirectedEdge(Node a, Node b, Node c, Graph graph) {
        if (graph.getEdge(b, a).getDistalEndpoint(b) == Endpoint.ARROW) {
            return true;
        }
        if (graph.getEdge(b, c).getDistalEndpoint(b) == Endpoint.ARROW) {
            return true;
        }
        return false;
    }

    // Tests whether adding a for b--a--c to the sepset (if it's not there) yields independence. Poor man's CPC.
    public static void orientCollidersUsingSepsets(SepsetMap set,
                                                   IKnowledge knowledge, Graph graph,
                                                   IndependenceTest test) {

        log("Starting Collider Orientation:");

//        verifySepsetIntegrity(set, graph);
        List<Node> nodes = graph.getNodes();

        for (Node b : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(b);
            Collections.sort(adjacentNodes);

            if (adjacentNodes.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node a = adjacentNodes.get(combination[0]);
                Node c = adjacentNodes.get(combination[1]);

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(a, c)) {
                    continue;
                }

                List<Node> sepset = set.get(a, c);

                List<Node> sepset2 = new ArrayList<>(sepset);

                if (!sepset2.contains(b)) {
                    System.out.println("\nADDING " + b);

                    sepset2.add(b);
                    double alpha = test.getAlpha();
                    test.setAlpha(alpha);

                    if (test.isIndependent(a, c, sepset2)) {
                        sepset = sepset2;
                    }
                }

                if (!sepset.contains(b)
                        && isArrowpointAllowed(a, b, knowledge)
                        && isArrowpointAllowed(c, b, knowledge)) {
                    System.out.println("Collider orientation <" + a + ", " + b + ", " + c + "> sepset = " + sepset);

                    graph.setEndpoint(a, b, Endpoint.ARROW);
                    graph.setEndpoint(c, b, Endpoint.ARROW);
                    log(SearchLogUtils.colliderOrientedMsg(a, b, c, sepset));
                }
            }
        }

        log("Finishing Collider Orientation.");
    }

    public static void orientCollidersLocally(IKnowledge knowledge, Graph graph,
                                              IndependenceTest test,
                                              int depth) {
        orientCollidersLocally(knowledge, graph, test, depth, null);
    }

    public static void orientCollidersLocally(IKnowledge knowledge, Graph graph,
                                              IndependenceTest test,
                                              int depth, Set<Node> nodesToVisit) {
        log("Starting Collider Orientation:");

        if (nodesToVisit == null) {
            nodesToVisit = new HashSet<>(graph.getNodes());
        }

        for (Node a : nodesToVisit) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(a);

            if (adjacentNodes.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node b = adjacentNodes.get(combination[0]);
                Node c = adjacentNodes.get(combination[1]);

                // Skip triples that are shielded.
                if (graph.isAdjacentTo(b, c)) {
                    continue;
                }

                if (isArrowpointAllowed1(b, a, knowledge)
                        && isArrowpointAllowed1(c, a, knowledge)) {
                    if (!existsLocalSepsetWith(b, a, c, test, graph, depth)) {
                        graph.setEndpoint(b, a, Endpoint.ARROW);
                        graph.setEndpoint(c, a, Endpoint.ARROW);
                        log(SearchLogUtils.colliderOrientedMsg(b, a, c));
                    }
                }
            }
        }

        log("Finishing Collider Orientation.");
    }

    public static boolean existsLocalSepsetWith(Node x, Node y, Node z,
                                                IndependenceTest test, Graph graph, int depth) {
        Set<Node> __nodes = new HashSet<>(graph.getAdjacentNodes(x));
        __nodes.addAll(graph.getAdjacentNodes(z));
        __nodes.remove(x);
        __nodes.remove(z);

        List<Node> _nodes = new LinkedList<>(__nodes);
        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        int _depth = depth;
        if (_depth == -1) {
            _depth = 1000;
        }
        _depth = Math.min(_depth, _nodes.size());

        for (int d = 1; d <= _depth; d++) {
            if (_nodes.size() >= d) {
                ChoiceGenerator cg2 = new ChoiceGenerator(_nodes.size(), d);
                int[] choice;

                while ((choice = cg2.next()) != null) {
                    List<Node> condSet = GraphUtils.asList(choice, _nodes);

                    if (!condSet.contains(y)) {
                        continue;
                    }

//                    LogUtils.getInstance().finest("Trying " + condSet);
                    if (test.isIndependent(x, z, condSet)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    public static boolean existsLocalSepsetWithout(Node x, Node y, Node z,
                                                   IndependenceTest test, Graph graph, int depth) {
        Set<Node> __nodes = new HashSet<>(graph.getAdjacentNodes(x));
        __nodes.addAll(graph.getAdjacentNodes(z));
        __nodes.remove(x);
        __nodes.remove(z);
        List<Node> _nodes = new LinkedList<>(__nodes);
        log(
                "Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        int _depth = depth;
        if (_depth == -1) {
            _depth = 1000;
        }
        _depth = Math.min(_depth, _nodes.size());

        for (int d = 0; d <= _depth; d++) {
            if (_nodes.size() >= d) {
                ChoiceGenerator cg2 = new ChoiceGenerator(_nodes.size(), d);
                int[] choice;

                while ((choice = cg2.next()) != null) {
                    List<Node> condSet = GraphUtils.asList(choice, _nodes);

                    if (condSet.contains(y)) {
                        continue;
                    }

                    //            LogUtils.getInstance().finest("Trying " + condSet);
                    if (test.isIndependent(x, z, condSet)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    public static boolean existsLocalSepsetWithoutDet(Node x, Node y, Node z,
                                                      IndependenceTest test, Graph graph, int depth) {
        Set<Node> __nodes = new HashSet<>(graph.getAdjacentNodes(x));
        __nodes.addAll(graph.getAdjacentNodes(z));
        __nodes.remove(x);
        __nodes.remove(z);
        List<Node> _nodes = new LinkedList<>(__nodes);
        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        int _depth = depth;
        if (_depth == -1) {
            _depth = 1000;
        }
        _depth = Math.min(_depth, _nodes.size());

        for (int d = 0; d <= _depth; d++) {
            if (_nodes.size() >= d) {
                ChoiceGenerator cg2 = new ChoiceGenerator(_nodes.size(), d);
                int[] choice;

                while ((choice = cg2.next()) != null) {
                    List<Node> condSet = GraphUtils.asList(choice, _nodes);

                    if (condSet.contains(y)) {
                        continue;
                    }

                    if (test.determines(condSet, y)) {
                        continue;
                    }

                    //        LogUtils.getInstance().finest("Trying " + condSet);
                    if (test.isIndependent(x, z, condSet)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    /**
     * Orient away from collider.
     */
    public static boolean meekR1Locally(Graph graph, IKnowledge knowledge,
                                        IndependenceTest test, int depth) {
        List<Node> nodes = graph.getNodes();
        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node a : nodes) {
                List<Node> adjacentNodes = graph.getAdjacentNodes(a);

                if (adjacentNodes.size() < 2) {
                    continue;
                }

                ChoiceGenerator cg
                        = new ChoiceGenerator(adjacentNodes.size(), 2);
                int[] combination;

                while ((combination = cg.next()) != null) {
                    Node b = adjacentNodes.get(combination[0]);
                    Node c = adjacentNodes.get(combination[1]);

                    // Skip triples that are shielded.
                    if (graph.isAdjacentTo(b, c)) {
                        continue;
                    }

                    if (graph.getEndpoint(b, a) == Endpoint.ARROW
                            && graph.isUndirectedFromTo(a, c)) {
                        if (existsLocalSepsetWithout(b, a, c, test, graph,
                                depth)) {
                            continue;
                        }

                        if (isArrowpointAllowed(a, c, knowledge)) {
                            graph.setEndpoint(a, c, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek R1", graph.getEdge(a, c)));
                            changed = true;
                        }
                    } else if (graph.getEndpoint(c, a) == Endpoint.ARROW
                            && graph.isUndirectedFromTo(a, b)) {
                        if (existsLocalSepsetWithout(b, a, c, test, graph,
                                depth)) {
                            continue;
                        }

                        if (isArrowpointAllowed(a, b, knowledge)) {
                            graph.setEndpoint(a, b, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek R1", graph.getEdge(a, b)));
                            changed = true;
                        }
                    }
                }
            }
        }

        return changed;
    }

    public static boolean meekR1Locally2(Graph graph, IKnowledge knowledge,
                                         IndependenceTest test, int depth) {
        List<Node> nodes = graph.getNodes();
        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node a : nodes) {
                List<Node> adjacentNodes = graph.getAdjacentNodes(a);

                if (adjacentNodes.size() < 2) {
                    continue;
                }

                ChoiceGenerator cg
                        = new ChoiceGenerator(adjacentNodes.size(), 2);
                int[] combination;

                while ((combination = cg.next()) != null) {
                    Node b = adjacentNodes.get(combination[0]);
                    Node c = adjacentNodes.get(combination[1]);

                    // Skip triples that are shielded.
                    if (graph.isAdjacentTo(b, c)) {
                        continue;
                    }

                    if (graph.getEndpoint(b, a) == Endpoint.ARROW
                            && graph.isUndirectedFromTo(a, c)) {
                        if (existsLocalSepsetWithoutDet(b, a, c, test, graph,
                                depth)) {
                            continue;
                        }

                        if (isArrowpointAllowed(a, c, knowledge)) {
                            graph.setEndpoint(a, c, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek R1", graph.getEdge(a, c)));
                            changed = true;
                        }
                    } else if (graph.getEndpoint(c, a) == Endpoint.ARROW
                            && graph.isUndirectedFromTo(a, b)) {
                        if (existsLocalSepsetWithoutDet(b, a, c, test, graph,
                                depth)) {
                            continue;
                        }

                        if (isArrowpointAllowed(a, b, knowledge)) {
                            graph.setEndpoint(a, b, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek R1", graph.getEdge(a, b)));
                            changed = true;
                        }
                    }
                }
            }
        }

        return changed;
    }

    /**
     * If
     */
    public static boolean meekR2(Graph graph, IKnowledge knowledge) {
        List<Node> nodes = graph.getNodes();
        boolean changed = false;

        for (Node a : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(a);

            if (adjacentNodes.size() < 2) {
                continue;
            }

            ChoiceGenerator cg = new ChoiceGenerator(adjacentNodes.size(), 2);
            int[] combination;

            while ((combination = cg.next()) != null) {
                Node b = adjacentNodes.get(combination[0]);
                Node c = adjacentNodes.get(combination[1]);

                if (graph.isDirectedFromTo(b, a)
                        && graph.isDirectedFromTo(a, c)
                        && graph.isUndirectedFromTo(b, c)) {
                    if (isArrowpointAllowed(b, c, knowledge)) {
                        graph.setEndpoint(b, c, Endpoint.ARROW);
                        log(SearchLogUtils.edgeOrientedMsg("Meek R2", graph.getEdge(b, c)));
                    }
                } else if (graph.isDirectedFromTo(c, a)
                        && graph.isDirectedFromTo(a, b)
                        && graph.isUndirectedFromTo(c, b)) {
                    if (isArrowpointAllowed(c, b, knowledge)) {
                        graph.setEndpoint(c, b, Endpoint.ARROW);
                        log(SearchLogUtils.edgeOrientedMsg("Meek R2", graph.getEdge(c, b)));
                    }
                }
            }
        }

        return changed;
    }

    /**
     * Meek's rule R3. If a--b, a--c, a--d, c-->b, c-->b, then orient a-->b.
     */
    public static boolean meekR3(Graph graph, IKnowledge knowledge) {

        List<Node> nodes = graph.getNodes();
        boolean changed = false;

        for (Node a : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(a);

            if (adjacentNodes.size() < 3) {
                continue;
            }

            for (Node b : adjacentNodes) {
                List<Node> otherAdjacents = new LinkedList<>(adjacentNodes);
                otherAdjacents.remove(b);

                if (!graph.isUndirectedFromTo(a, b)) {
                    continue;
                }

                ChoiceGenerator cg
                        = new ChoiceGenerator(otherAdjacents.size(), 2);
                int[] combination;

                while ((combination = cg.next()) != null) {
                    Node c = otherAdjacents.get(combination[0]);
                    Node d = otherAdjacents.get(combination[1]);

                    if (graph.isAdjacentTo(c, d)) {
                        continue;
                    }

                    if (!graph.isUndirectedFromTo(a, c)) {
                        continue;
                    }

                    if (!graph.isUndirectedFromTo(a, d)) {
                        continue;
                    }

                    if (graph.isDirectedFromTo(c, b)
                            && graph.isDirectedFromTo(d, b)) {
                        if (isArrowpointAllowed(a, b, knowledge)) {
                            graph.setEndpoint(a, b, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek R3", graph.getEdge(a, b)));
                            changed = true;
                            break;
                        }
                    }
                }
            }
        }

        return changed;
    }

    public static boolean meekR4(Graph graph, IKnowledge knowledge) {
        if (knowledge == null) {
            return false;
        }

        List<Node> nodes = graph.getNodes();
        boolean changed = false;

        for (Node a : nodes) {
            List<Node> adjacentNodes = graph.getAdjacentNodes(a);

            if (adjacentNodes.size() < 3) {
                continue;
            }

            for (Node d : adjacentNodes) {
                if (!graph.isAdjacentTo(a, d)) {
                    continue;
                }

                List<Node> otherAdjacents = new LinkedList<>(adjacentNodes);
                otherAdjacents.remove(d);

                ChoiceGenerator cg
                        = new ChoiceGenerator(otherAdjacents.size(), 2);
                int[] combination;

                while ((combination = cg.next()) != null) {
                    Node b = otherAdjacents.get(combination[0]);
                    Node c = otherAdjacents.get(combination[1]);

                    if (!graph.isUndirectedFromTo(a, b)) {
                        continue;
                    }

                    if (!graph.isUndirectedFromTo(a, c)) {
                        continue;
                    }

//                    if (!isUnshieldedNoncollider(c, a, b, graph)) {
//                        continue;
//                    }
                    if (graph.isDirectedFromTo(b, c)
                            && graph.isDirectedFromTo(d, c)) {
                        if (isArrowpointAllowed(a, c, knowledge)) {
                            graph.setEndpoint(a, c, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek T1", graph.getEdge(a, c)));
                            changed = true;
                            break;
                        }
                    } else if (graph.isDirectedFromTo(c, d)
                            && graph.isDirectedFromTo(d, b)) {
                        if (isArrowpointAllowed(a, b, knowledge)) {
                            graph.setEndpoint(a, b, Endpoint.ARROW);
                            log(SearchLogUtils.edgeOrientedMsg("Meek T1", graph.getEdge(a, b)));
                            changed = true;
                            break;
                        }
                    }
                }
            }
        }

        return changed;
    }

    /**
     * Checks if an arrowpoint is allowed by background knowledge.
     */
    public static boolean isArrowpointAllowed(Object from, Object to,
                                              IKnowledge knowledge) {
        if (knowledge == null) {
            return true;
        }
        return !knowledge.isRequired(to.toString(), from.toString())
                && !knowledge.isForbidden(from.toString(), to.toString());
    }

    /**
     * Transforms a maximally directed pattern (PDAG) represented in graph
     * <code>g</code> into an arbitrary DAG by modifying <code>g</code> itself.
     * Based on the algorithm described in </p> Chickering (2002) "Optimal
     * structure identification with greedy search" Journal of Machine Learning
     * Research. </p> R. Silva, June 2004
     */
    public static void pdagToDag(Graph g) {
        Graph p = new EdgeListGraph(g);
        List<Edge> undirectedEdges = new ArrayList<>();

        for (Edge edge : g.getEdges()) {
            if (edge.getEndpoint1() == Endpoint.TAIL
                    && edge.getEndpoint2() == Endpoint.TAIL
                    && !undirectedEdges.contains(edge)) {
                undirectedEdges.add(edge);
            }
        }
        g.removeEdges(undirectedEdges);
        List<Node> pNodes = p.getNodes();

        do {
            Node x = null;

            for (Node pNode : pNodes) {
                x = pNode;

                if (p.getChildren(x).size() > 0) {
                    continue;
                }

                Set<Node> neighbors = new HashSet<>();

                for (Edge edge : p.getEdges()) {
                    if (edge.getNode1() == x || edge.getNode2() == x) {
                        if (edge.getEndpoint1() == Endpoint.TAIL
                                && edge.getEndpoint2() == Endpoint.TAIL) {
                            if (edge.getNode1() == x) {
                                neighbors.add(edge.getNode2());
                            } else {
                                neighbors.add(edge.getNode1());
                            }
                        }
                    }
                }
                if (neighbors.size() > 0) {
                    Collection<Node> parents = p.getParents(x);
                    Set<Node> all = new HashSet<>(neighbors);
                    all.addAll(parents);
                    if (!GraphUtils.isClique(all, p)) {
                        continue;
                    }
                }

                for (Node neighbor : neighbors) {
                    Node node1 = g.getNode(neighbor.getName());
                    Node node2 = g.getNode(x.getName());

                    g.addDirectedEdge(node1, node2);
                }
                p.removeNode(x);
                break;
            }
            pNodes.remove(x);
        } while (pNodes.size() > 0);
    }

    /**
     * Get a graph and direct only the unshielded colliders.
     */
    public static void basicPattern(Graph graph, boolean orientInPlace) {
        Set<Edge> undirectedEdges = new HashSet<>();

        NEXT_EDGE:
        for (Edge edge : graph.getEdges()) {
            if (!edge.isDirected()) {
                continue;
            }

            Node x = Edges.getDirectedEdgeTail(edge);
            Node y = Edges.getDirectedEdgeHead(edge);

            for (Node parent : graph.getParents(y)) {
                if (parent != x) {
                    if (!graph.isAdjacentTo(parent, x)) {
                        continue NEXT_EDGE;
                    }
                }
            }

            undirectedEdges.add(edge);
        }

        for (Edge nextUndirected : undirectedEdges) {
            if (orientInPlace) {
                nextUndirected.setEndpoint1(Endpoint.TAIL);
                nextUndirected.setEndpoint2(Endpoint.TAIL);
            } else {
                Node node1 = nextUndirected.getNode1();
                Node node2 = nextUndirected.getNode2();

                graph.removeEdges(node1, node2);
                graph.addUndirectedEdge(node1, node2);
            }
        }
    }

    public static void basicPatternRestricted(Graph graph, Set<Edge> edges) {
        Set<Edge> undirectedEdges = new HashSet<>();

        NEXT_EDGE:
        for (Edge edge : edges) {
            if (!edge.isDirected()) {
                continue;
            }

            Node _x = Edges.getDirectedEdgeTail(edge);
            Node _y = Edges.getDirectedEdgeHead(edge);

            for (Node parent : graph.getParents(_y)) {
                if (parent != _x) {
                    if (!graph.isAdjacentTo(parent, _x)) {
                        continue NEXT_EDGE;
                    }
                }
            }

            undirectedEdges.add(edge);
        }

        for (Edge nextUndirected : undirectedEdges) {
            Node node1 = nextUndirected.getNode1();
            Node node2 = nextUndirected.getNode2();

            graph.removeEdge(nextUndirected);
            graph.addUndirectedEdge(node1, node2);
        }
    }

    public static void basicPatternRestricted2(Graph graph, Node node) {
        Set<Edge> undirectedEdges = new HashSet<>();

        NEXT_EDGE:
        for (Edge edge : graph.getEdges(node)) {
            if (!edge.isDirected()) {
                continue;
            }

            Node _x = Edges.getDirectedEdgeTail(edge);
            Node _y = Edges.getDirectedEdgeHead(edge);

            for (Node parent : graph.getParents(_y)) {
                if (parent != _x) {
                    if (!graph.isAdjacentTo(parent, _x)) {
                        continue NEXT_EDGE;
                    }
                }
            }

            undirectedEdges.add(edge);
        }

        for (Edge nextUndirected : undirectedEdges) {
            Node node1 = nextUndirected.getNode1();
            Node node2 = nextUndirected.getNode2();

            graph.removeEdge(nextUndirected);
            graph.addUndirectedEdge(node1, node2);
        }
    }

    /**
     * Transforms a DAG represented in graph <code>graph</code> into a maximally
     * directed pattern (PDAG) by modifying <code>g</code> itself. Based on the
     * algorithm described in </p> Chickering (2002) "Optimal structure
     * identification with greedy search" Journal of Machine Learning Research.
     * It works for both BayesNets and SEMs.
     * </p> R. Silva, June 2004
     */
    public static void dagToPdag(Graph graph) {
        //do topological sort on the nodes
        Graph graphCopy = new EdgeListGraph(graph);
        Node orderedNodes[] = new Node[graphCopy.getNodes().size()];
        int count = 0;
        while (graphCopy.getNodes().size() > 0) {
            Set<Node> exogenousNodes = new HashSet<>();

            for (Node next : graphCopy.getNodes()) {
                if (graphCopy.isExogenous(next)) {
                    exogenousNodes.add(next);
                    orderedNodes[count++] = graph.getNode(next.getName());
                }
            }

            graphCopy.removeNodes(new ArrayList<>(exogenousNodes));
        }
        //ordered edges - improvised, inefficient implementation
        count = 0;
        Edge edges[] = new Edge[graph.getNumEdges()];
        boolean edgeOrdered[] = new boolean[graph.getNumEdges()];
        Edge orderedEdges[] = new Edge[graph.getNumEdges()];

        for (Edge edge : graph.getEdges()) {
            edges[count++] = edge;
        }

        for (int i = 0; i < edges.length; i++) {
            edgeOrdered[i] = false;
        }

        while (count > 0) {
            for (Node orderedNode : orderedNodes) {
                for (int k = orderedNodes.length - 1; k >= 0; k--) {
                    for (int q = 0; q < edges.length; q++) {
                        if (!edgeOrdered[q]
                                && edges[q].getNode1() == orderedNodes[k]
                                && edges[q].getNode2() == orderedNode) {
                            edgeOrdered[q] = true;
                            orderedEdges[orderedEdges.length - count]
                                    = edges[q];
                            count--;
                        }
                    }
                }
            }
        }

        //label edges
        boolean compelledEdges[] = new boolean[graph.getNumEdges()];
        boolean reversibleEdges[] = new boolean[graph.getNumEdges()];
        for (int i = 0; i < graph.getNumEdges(); i++) {
            compelledEdges[i] = false;
            reversibleEdges[i] = false;
        }
        for (int i = 0; i < graph.getNumEdges(); i++) {
            if (compelledEdges[i] || reversibleEdges[i]) {
                continue;
            }
            Node x = orderedEdges[i].getNode1();
            Node y = orderedEdges[i].getNode2();
            for (int j = 0; j < orderedEdges.length; j++) {
                if (orderedEdges[j].getNode2() == x && compelledEdges[j]) {
                    Node w = orderedEdges[j].getNode1();
                    if (!graph.isParentOf(w, y)) {
                        for (int k = 0; k < orderedEdges.length; k++) {
                            if (orderedEdges[k].getNode2() == y) {
                                compelledEdges[k] = true;
                                break;
                            }
                        }
                    } else {
                        for (int k = 0; k < orderedEdges.length; k++) {
                            if (orderedEdges[k].getNode1() == w
                                    && orderedEdges[k].getNode2() == y) {
                                compelledEdges[k] = true;
                                break;
                            }
                        }
                    }
                }
                if (compelledEdges[i]) {
                    break;
                }
            }
            if (compelledEdges[i]) {
                continue;
            }
            boolean foundZ = false;

            for (Edge orderedEdge : orderedEdges) {
                Node z = orderedEdge.getNode1();
                if (z != x && orderedEdge.getNode2() == y
                        && !graph.isParentOf(z, x)) {
                    compelledEdges[i] = true;
                    for (int k = i + 1; k < graph.getNumEdges(); k++) {
                        if (orderedEdges[k].getNode2() == y
                                && !reversibleEdges[k]) {
                            compelledEdges[k] = true;
                        }
                    }
                    foundZ = true;
                    break;
                }
            }

            if (!foundZ) {
                reversibleEdges[i] = true;

                for (int j = i + 1; j < orderedEdges.length; j++) {
                    if (!compelledEdges[j] && orderedEdges[j].getNode2() == y) {
                        reversibleEdges[j] = true;
                    }
                }
            }
        }

        //undirect edges that are reversible
        for (int i = 0; i < reversibleEdges.length; i++) {
            if (reversibleEdges[i]) {
                graph.setEndpoint(orderedEdges[i].getNode1(),
                        orderedEdges[i].getNode2(), Endpoint.TAIL);
                graph.setEndpoint(orderedEdges[i].getNode2(),
                        orderedEdges[i].getNode1(), Endpoint.TAIL);
            }
        }
    }

    public static Graph patternFromEPattern(Graph ePattern) {
        ePattern = new EdgeListGraph(ePattern);

        MeekRules rules = new MeekRules();
        rules.orientImplied(ePattern);

        List<Triple> ambiguousTriples = new ArrayList<>(ePattern.getAmbiguousTriples());
        removeExtraAmbiguousTriples(ePattern, new ArrayList<>(ambiguousTriples));

        while (!ambiguousTriples.isEmpty()) {
            Triple triple = ambiguousTriples.get(0);

            Node x = triple.getX();
            Node y = triple.getY();
            Node z = triple.getZ();

            ePattern.removeAmbiguousTriple(x, y, z);
            ambiguousTriples.remove(triple);

//            if (!ePattern.isDefCollider(x, y, z)) {
//                ePattern.removeEdge(x, y);
//                ePattern.removeEdge(z, y);
//                ePattern.addDirectedEdge(x, y);
//                ePattern.addDirectedEdge(z, y);
//            }

            rules.orientImplied(ePattern);
            removeExtraAmbiguousTriples(ePattern, ambiguousTriples);
        }

        return ePattern;
    }

    private static void removeExtraAmbiguousTriples(Graph graph, List<Triple> ambiguousTriples) {
        Set<Triple> ambiguities = graph.getAmbiguousTriples();

        for (Triple triple : new HashSet<>(ambiguities)) {
            final Node x = triple.getX();
            final Node y = triple.getY();
            final Node z = triple.getZ();

            if (!graph.isAdjacentTo(x, y) || !graph.isAdjacentTo(y, x)) {
                graph.removeAmbiguousTriple(x, y, z);
                ambiguousTriples.remove(triple);
            }

            if (graph.isDefCollider(x, y, z)) {
                graph.removeAmbiguousTriple(x, y, z);
                ambiguousTriples.remove(triple);
            }

            if (graph.getEdge(x, y).pointsTowards(x) || graph.getEdge(y, z).pointsTowards(z)) {
                graph.removeAmbiguousTriple(x, y, z);
                ambiguousTriples.remove(triple);
            }
        }
    }

    private static List<Triple> asList(int[] indices, List<Triple> nodes) {
        List<Triple> list = new LinkedList<>();

        for (int i : indices) {
            list.add(nodes.get(i));
        }

        return list;
    }

    private static void direct(Node a, Node c, Graph graph) {
        Edge before = graph.getEdge(a, c);
        Edge after = Edges.directedEdge(a, c);
        graph.removeEdge(before);
        graph.addEdge(after);
    }

    private static boolean orientOneCircle(Graph graph) {
        for (Edge edge : graph.getEdges()) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            if (graph.getEndpoint(x, y) == Endpoint.CIRCLE) {
                graph.setEndpoint(x, y, Endpoint.ARROW);
                return true;
            }

            if (graph.getEndpoint(y, x) == Endpoint.CIRCLE) {
                graph.setEndpoint(y, x, Endpoint.ARROW);
                return true;
            }
        }

        return false;
    }

    /**
     * Double checks a sepset map against a pattern to make sure that X is
     * adjacent to Y in the pattern iff {X, Y} is not in the domain of the
     * sepset map.
     *
     * @param sepset  a sepset map, over variables v.
     * @param pattern a pattern over variables W, v subset of W.
     * @return true if the sepset map is consistent with the pattern.
     */
    public static boolean verifySepsetIntegrity(SepsetMap sepset, Graph pattern) {
        for (Node x : pattern.getNodes()) {
            for (Node y : pattern.getNodes()) {
                if (x == y) {
                    continue;
                }

                if ((pattern.isAdjacentTo(y, x)) != (sepset.get(x, y) == null)) {
                    System.out.println("Sepset not consistent with graph for {" + x + ", " + y + "}");
                    return false;
                }
            }
        }

        return true;
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

    public static List<Set<Node>> powerSet(List<Node> nodes) {
        List<Set<Node>> subsets = new ArrayList<>();
        int total = (int) Math.pow(2, nodes.size());
        for (int i = 0; i < total; i++) {
            Set<Node> newSet = new HashSet<>();
            String selection = Integer.toBinaryString(i);
            for (int j = selection.length() - 1; j >= 0; j--) {
                if (selection.charAt(j) == '1') {
                    newSet.add(nodes.get(selection.length() - j - 1));
                }
            }
            subsets.add(newSet);
        }
        return subsets;
    }

    /**
     * Checks if an arrowpoint is allowed by background knowledge.
     */
    public static boolean isArrowpointAllowed1(Node from, Node to,
                                               IKnowledge knowledge) {
        if (knowledge == null) {
            return true;
        }

        return !knowledge.isRequired(to.toString(), from.toString())
                && !knowledge.isForbidden(from.toString(), to.toString());
    }

    public static boolean isArrowpointAllowed2(Node from, Node to,
                                               IKnowledge knowledge, Graph graph) {
        if (knowledge == null) {
            return true;
        }

        if (!graph.getNodesInTo(to, Endpoint.ARROW).isEmpty()) {
            return false;
        }

        return !knowledge.isRequired(to.toString(), from.toString())
                && !knowledge.isForbidden(from.toString(), to.toString());
    }

    // The published version.
    public static CpcTripleType getCpcTripleType(Node x, Node y, Node z,
                                                 IndependenceTest test, int depth,
                                                 Graph graph, boolean verbose) {
        int numSepsetsContainingY = 0;
        int numSepsetsNotContainingY = 0;

        List<Node> _nodes = graph.getAdjacentNodes(x);
        _nodes.remove(z);
        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        int _depth = depth;
        if (_depth == -1) {
            _depth = 1000;
        }
        _depth = Math.min(_depth, _nodes.size());

        for (int d = 0; d <= _depth; d++) {
            ChoiceGenerator cg = new ChoiceGenerator(_nodes.size(), d);
            int[] choice;

            while ((choice = cg.next()) != null) {
                List<Node> cond = GraphUtils.asList(choice, _nodes);

                if (test.isIndependent(x, z, cond)) {
                    if (verbose) {
//                        System.out.println("Indep: " + x + " _||_ " + z + " | " + cond);
                    }

                    if (cond.contains(y)) {
                        numSepsetsContainingY++;
                    } else {
                        numSepsetsNotContainingY++;
                    }
                }

                if (numSepsetsContainingY > 0 && numSepsetsNotContainingY > 0) {
                    return CpcTripleType.AMBIGUOUS;
                }
            }
        }

        _nodes = graph.getAdjacentNodes(z);
        _nodes.remove(x);
        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        if (_depth == -1) {
            _depth = 1000;
        }

        _depth = Math.min(_depth, _nodes.size());

        for (int d = 0; d <= _depth; d++) {
            ChoiceGenerator cg = new ChoiceGenerator(_nodes.size(), d);
            int[] choice;

            while ((choice = cg.next()) != null) {
                List<Node> cond = GraphUtils.asList(choice, _nodes);

                if (test.isIndependent(x, z, cond)) {
                    if (cond.contains(y)) {
                        numSepsetsContainingY++;
                    } else {
                        numSepsetsNotContainingY++;
                    }
                }

                if (numSepsetsContainingY > 0 && numSepsetsNotContainingY > 0) {
                    return CpcTripleType.AMBIGUOUS;
                }
            }
        }

        if (numSepsetsContainingY > 0) {
            return CpcTripleType.NONCOLLIDER;
        } else {
            return CpcTripleType.COLLIDER;
        }
    }

    public static CpcTripleType getCpcTripleType3(Node x, Node y, Node z,
                                                  IndependenceTest test, int depth,
                                                  Graph graph) {
        int numSepsetsContainingY = 0;
        int numSepsetsNotContainingY = 0;

        Set<Set<Node>> withY = new HashSet<>();
        Set<Set<Node>> withoutY = new HashSet<>();

        Set<Node> n = new HashSet<>(graph.getAdjacentNodes(x));
        n.addAll(graph.getAdjacentNodes(z));
        List<Node> _nodes = new ArrayList<>(n);

//        List<Node> _nodes = graph.getAdjacentNodes(x);
//        _nodes.addAll(graph.getAdjacentNodes(z));
        _nodes.remove(x);
        _nodes.remove(z);
        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);

        int _depth = depth;
        if (_depth == -1) {
            _depth = 1000;
        }
        _depth = Math.min(_depth, _nodes.size());

        for (int d = 0; d <= _depth; d++) {
            ChoiceGenerator cg = new ChoiceGenerator(_nodes.size(), d);
            int[] choice;

            while ((choice = cg.next()) != null) {
                List<Node> cond = GraphUtils.asList(choice, _nodes);

                if (test.isIndependent(x, z, cond)) {
                    if (cond.contains(y)) {
                        numSepsetsContainingY++;
                        withY.add(new HashSet<>(cond));
                    } else {
                        numSepsetsNotContainingY++;
                        withoutY.add(new HashSet<>(cond));
                    }
                }
            }
        }

//        _nodes = graph.getAdjacentNodes(z);
//        _nodes.remove(x);
//        log("Adjacents for " + x + "--" + y + "--" + z + " = " + _nodes);
//
//        _depth = depth;
//        if (_depth == -1) {
//            _depth = 1000;
//        }
//        _depth = Math.min(_depth, _nodes.size());
//        for (int d = 0; d <= _depth; d++) {
//            ChoiceGenerator cg = new ChoiceGenerator(_nodes.size(), d);
//            int[] choice;
//
//            while ((choice = cg.next()) != null) {
//                List<Node> cond = DataGraphUtils.asList(choice, _nodes);
//
//                if (test.isIndependent(x, z, cond)) {
//                    if (cond.contains(y)) {
//                        numSepsetsContainingY++;
//                        withY.add(new HashSet<Node>(cond));
//                    } else {
//                        numSepsetsNotContainingY++;
//                        withoutY.add(new HashSet<Node>(cond));
//                    }
//                }
//            }
//        }
//        int factor = 3;
        int factor = 1;

        numSepsetsContainingY = withY.size();
        numSepsetsNotContainingY = withoutY.size();

        if (numSepsetsContainingY > factor * numSepsetsNotContainingY) {
            return CpcTripleType.NONCOLLIDER;
        } else if (numSepsetsNotContainingY > factor * numSepsetsContainingY) {
            return CpcTripleType.COLLIDER;
        } else {
            return CpcTripleType.AMBIGUOUS;
        }
    }

    public static void orientRequired(IKnowledge bk, Graph graph, List<Node> nodes) {
        log("Staring BK Orientation.");
        for (Iterator<KnowledgeEdge> it = bk.forbiddenEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in the graph.
            Node from = translate(edge.getFrom(), nodes);
            Node to = translate(edge.getTo(), nodes);

            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            // Orient to-->from
            graph.removeEdge(from, to);
            graph.addDirectedEdge(to, from);
//            graph.setEndpoint(from, to, Endpoint.TAIL);
//            graph.setEndpoint(to, from, Endpoint.ARROW);

            log(SearchLogUtils.edgeOrientedMsg("IKnowledge", graph.getEdge(to, from)));
        }

        for (Iterator<KnowledgeEdge> it = bk.requiredEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge edge = it.next();

            //match strings to variables in this graph
            Node from = translate(edge.getFrom(), nodes);
            Node to = translate(edge.getTo(), nodes);

            if (from == null || to == null) {
                continue;
            }

            if (graph.getEdge(from, to) == null) {
                continue;
            }

            // Orient from-->to
            graph.removeEdge(from, to);
            graph.addDirectedEdge(from, to);

//            graph.setEndpoint(to, from, Endpoint.TAIL);
//            graph.setEndpoint(from, to, Endpoint.ARROW);
            log(SearchLogUtils.edgeOrientedMsg("IKnowledge", graph.getEdge(from, to)));
        }

        log("Finishing BK Orientation.");
    }

    public static int structuralHammingDistance(Graph trueGraph, Graph estGraph) {
        int error = 0;

        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());

        Set<Node> _allNodes = new HashSet<>();

        List<Node> trueLatents = trueGraph.getNodes();
        List<Node> estLatents = estGraph.getNodes();

//        List<Node> trueLatents = GraphUtils.getLatents(trueGraph);
//        List<Node> estLatents = GraphUtils.getLatents(graph);
        Graph u = trueGraph.subgraph(trueLatents);
        Graph t = estGraph.subgraph(estLatents);

        Graph G = u; //patternForDag(u);
        Graph H = t; //patternForDag(t);

//        System.out.println("Pattern of true graph over latents = " + G);
        _allNodes.addAll(trueLatents);
        _allNodes.addAll(estLatents);

        List<Node> allNodes = new ArrayList<>(_allNodes);

        for (int i1 = 0; i1 < allNodes.size(); i1++) {
            for (int i2 = i1 + 1; i2 < allNodes.size(); i2++) {
                Node l1 = allNodes.get(i1);
                Node l2 = allNodes.get(i2);

                Edge e1 = G.getEdge(l1, l2);
                Edge e2 = H.getEdge(l1, l2);

                int shd = structuralHammingDistanceOneEdge(e1, e2);

                error += shd;
            }
        }
        return error;
    }

    private static int structuralHammingDistanceOneEdge(Edge e1, Edge e2) {
        if (noEdge(e1) && undirected(e2)) {
            return 1;
        } else if (noEdge(e2) && undirected(e1)) {
            return 1;
        } else if (noEdge(e1) && directed(e2)) {
            return 2;
        } else if (noEdge(e2) && directed(e1)) {
            return 2;
        } else if (undirected(e1) && directed(e2)) {
            return 1;
        } else if (undirected(e2) && directed(e1)) {
            return 1;
        } else if (directed(e1) && directed(e2)) {
            if (Edges.getDirectedEdgeHead(e1) == Edges.getDirectedEdgeTail(e2)) {
                return 1;
            }
        } else if (bidirected(e1) || bidirected(e2)) {
            return 2;
        }

        return 0;
    }

    private static boolean directed(Edge e2) {
        return e2 != null && Edges.isDirectedEdge(e2);
    }

    private static boolean bidirected(Edge e2) {
        return e2 != null && Edges.isBidirectedEdge(e2);
    }

    private static boolean undirected(Edge e2) {
        return e2 != null && Edges.isUndirectedEdge(e2);
    }

    private static boolean noEdge(Edge e1) {
        return e1 == null;
    }

    public static int structuralHammingDistance3a(Graph trueGraph, Graph estGraph) {
        int error = 0;

        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());

        List<Node> _allNodes = estGraph.getNodes();

        List<Node> allNodes = new ArrayList<>(_allNodes);

        for (int i1 = 0; i1 < allNodes.size(); i1++) {
            for (int i2 = i1 + 1; i2 < allNodes.size(); i2++) {
                Node l1 = allNodes.get(i1);
                Node l2 = allNodes.get(i2);

                Edge e1 = trueGraph.getEdge(l1, l2);
                Edge e2 = estGraph.getEdge(l1, l2);

                int shd = structuralHammingDistanceOneEdge(e1, e2);
                error += shd;
            }
        }
        return error;
    }

    public static int structuralHammingDistance3(Graph trueGraph, Graph estGraph) {
        int error = 0;

        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());

        List<Node> _allNodes = estGraph.getNodes();

        List<Node> allNodes = new ArrayList<>(_allNodes);

        for (int i1 = 0; i1 < allNodes.size(); i1++) {
            for (int i2 = i1 + 1; i2 < allNodes.size(); i2++) {
                Node l1 = allNodes.get(i1);
                Node l2 = allNodes.get(i2);

                Edge e1 = trueGraph.getEdge(l1, l2);
                Edge e2 = estGraph.getEdge(l1, l2);

                int shd = structuralHammingDistanceOneEdge3(e1, e2);
                error += shd;
            }
        }
        return error;
    }

    private static int structuralHammingDistanceOneEdge3(Edge e1, Edge e2) {
        if (noEdge3(e1) && nondirected3(e2)) {
            return 1;
        } else if (noEdge3(e2) && nondirected3(e1)) {
            return 1;
        } else if (noEdge3(e1) && directed3(e2)) {
            return 2;
        } else if (noEdge3(e2) && directed3(e1)) {
            return 2;
        } else if (nondirected3(e1) && directed3(e2)) {
            return 1;
        } else if (nondirected3(e2) && directed3(e1)) {
            return 1;
        } else if (directed3(e1) && directed3(e2)) {
            if (Edges.getDirectedEdgeHead(e1) == Edges.getDirectedEdgeTail(e2)) {
                return 1;
            }
        }

        return 0;
    }

    private static boolean directed3(Edge e2) {
        return e2 != null && Edges.isDirectedEdge(e2);
    }

    private static boolean nondirected3(Edge e2) {
        return e2 != null && Edges.isNondirectedEdge(e2);
    }

    private static boolean noEdge3(Edge e1) {
        return e1 == null;
    }

    /**
     * Simple class to store edges for the reachability search.
     *
     * @author Joseph Ramsey
     */
    private static class ReachabilityEdge {

        private Node from;
        private Node to;

        public ReachabilityEdge(Node from, Node to) {
            this.from = from;
            this.to = to;
        }

        public int hashCode() {
            int hash = 17;
            hash += 63 * getFrom().hashCode();
            hash += 71 * getTo().hashCode();
            return hash;
        }

        public boolean equals(Object obj) {
            ReachabilityEdge edge = (ReachabilityEdge) obj;

            if (!(edge.getFrom().equals(this.getFrom()))) {
                return false;
            }

            return edge.getTo().equals(this.getTo());
        }

        public Node getFrom() {
            return from;
        }

        public Node getTo() {
            return to;
        }
    }

    public enum CpcTripleType {
        COLLIDER, NONCOLLIDER, AMBIGUOUS
    }
}
