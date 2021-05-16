package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.ChoiceGenerator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Attempts to improve an initial graph using a global score the adding, removing or reversing
 * edges.
 *
 * @author josephramsey
 */
public class GlobalScoreSearch {
    private final Scorer scorer;
    private IKnowledge knowledge = new Knowledge2();

    public GlobalScoreSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    public Graph improve(Graph graph) {
        EdgeListGraph _graph = new EdgeListGraph(graph);
        addRequiredEdges(_graph);
        graph = SearchGraphUtils.dagFromPattern(_graph);
        return decreaseScoreLoop(graph);
    }

    public Score scoreGraph(Graph graph, Scorer scorer) {
        Graph dag = new EdgeListGraph(graph);
        scorer.score(dag);
        return new Score(scorer);
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    private Graph decreaseScoreLoop(Graph graph) {
        double score0 = scoreGraph(graph, scorer).getScore();

        MOVES:
        while (true) {
            List<Move> moves = new ArrayList<>();
//            moves.addAll(getRemoveColliderMoves(graph));
//            moves.addAll(getDoubleRemoveMoves(graph));
//            moves.addAll(getRemoveTriangleMoves(graph));
//            moves.addAll(getSwapMoves(graph));
//            moves.addAll(getRedirectMoves(graph));
//            moves.addAll(getAddColliderMoves(graph));
            moves.addAll(getRemoveMoves(graph));
            moves.addAll(getAddMoves(graph));



            for (Move move : moves) {
                Graph g2 = makeMove(graph, move);

                if (getKnowledge().isViolatedBy(g2)) continue;

                Score _score = scoreGraph(g2, scorer);

                if (_score.getScore() < score0) {
                    score0 = _score.getScore();
                    graph = g2;

                    System.out.println("Decreases score (" + move.getType() + "): score = " + _score.getScore());

                    continue MOVES;
                }
            }

            break;
        }

        return graph;
    }

    private Graph makeMove(Graph graph, Move move) {
        graph = new EdgeListGraph(graph);
        Edge firstEdge = move.getFirstEdge();
        Edge secondEdge = move.getSecondEdge();

        if (firstEdge != null && move.getType() == Move.Type.ADD) {
            graph.removeEdge(firstEdge.getNode1(), firstEdge.getNode2());
            graph.addEdge(firstEdge);
        } else if (firstEdge != null && move.getType() == Move.Type.REMOVE) {
            graph.removeEdge(firstEdge);
        } else if (firstEdge != null && move.getType() == Move.Type.DOUBLE_REMOVE) {
            graph.removeEdge(firstEdge);
            graph.removeEdge(secondEdge);
        } else if (firstEdge != null && move.getType() == Move.Type.REDIRECT) {
            graph.removeEdge(graph.getEdge(firstEdge.getNode1(), firstEdge.getNode2()));
            graph.addEdge(firstEdge);
        } else if (firstEdge != null && secondEdge != null && move.getType() == Move.Type.ADD_COLLIDER) {
            Edge existingEdge1 = graph.getEdge(firstEdge.getNode1(), firstEdge.getNode2());
            Edge existingEdge2 = graph.getEdge(secondEdge.getNode1(), secondEdge.getNode2());

            if (existingEdge1 != null) {
                graph.removeEdge(existingEdge1);
            }

            if (existingEdge2 != null) {
                graph.removeEdge(existingEdge2);
            }

            graph.addEdge(firstEdge);
            graph.addEdge(secondEdge);
        } else if (firstEdge != null && secondEdge != null && move.getType() == Move.Type.REMOVE_COLLIDER) {
            graph.removeEdge(firstEdge);
            graph.removeEdge(secondEdge);
        } else if (firstEdge != null && secondEdge != null && move.getType() == Move.Type.SWAP) {
            graph.removeEdge(firstEdge);
            Edge secondEdgeStar = graph.getEdge(secondEdge.getNode1(), secondEdge.getNode2());

            if (secondEdgeStar != null) {
                graph.removeEdge(secondEdgeStar);
            }

            graph.addEdge(secondEdge);
        }

        return graph;
    }

    private List<Move> getRemoveMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        // Remove moves:
        List<Edge> edges = new ArrayList<>(graph.getEdges());

        for (Edge edge : edges) {
            Node i = edge.getNode1();
            Node j = edge.getNode2();

            if (getKnowledge().isRequired(i.getName(), j.getName())) {
                continue;
            }

            moves.add(new Move(edge, Move.Type.REMOVE));
        }

        return moves;
    }

    private List<Move> getRedirectMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        // Reverse moves:
        List<Edge> edges = new ArrayList<>(graph.getEdges());

        for (Edge edge : edges) {
            Node i = edge.getNode1();
            Node j = edge.getNode2();
            if (knowledge.isForbidden(j.getName(), i.getName())) {
                continue;
            }

            if (getKnowledge().isRequired(i.getName(), j.getName())) {
                continue;
            }

            if (graph.isAncestorOf(j, i)) {
                continue;
            }

            moves.add(new Move(Edges.directedEdge(j, i), Move.Type.REDIRECT));
        }

        return moves;
    }

    private List<Move> getAddMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        // Add moves:
        List<Node> nodes = graph.getNodes();
        Collections.sort(nodes);

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    continue;
                }

                if (graph.isAdjacentTo(nodes.get(i), nodes.get(j))) {
                    continue;
                }

                if (getKnowledge().isForbidden(nodes.get(i).getName(), nodes.get(j).getName())) {
                    continue;
                }

                if (getKnowledge().isRequired(nodes.get(j).getName(), nodes.get(i).getName())) {
                    continue;
                }

//                if (graph.isAncestorOf(nodes.get(j), nodes.get(i))) continue;

//                if (graph.getDegree(nodes.get(i)) > 0 || graph.getDegree(nodes.get(j)) > 0) {
                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
                moves.add(new Move(edge, Move.Type.ADD));
//                }

//                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
//                moves.add(new Move(edge, Move.Type.ADD));
//                }
            }
        }

        return moves;
    }

    private List<Move> getAddColliderMoves(Graph graph) {
//         Make collider moves:
        List<Move> moves = new ArrayList<>();

        for (Node b : graph.getNodes()) {
//            if (graph.getAdjacentNodes(b).isEmpty()) {
                List<Node> nodes = graph.getNodes();
                nodes.remove(b);
                nodes.removeAll(graph.getAdjacentNodes(b));

                if (nodes.size() >= 2) {
                    ChoiceGenerator gen = new ChoiceGenerator(nodes.size(), 2);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        List<Node> _nodes = GraphUtils.asList(choice, nodes);
                        Node a = _nodes.get(0);
                        Node c = _nodes.get(1);

                        if (a == b || c == b) continue;

                        Edge edge1 = Edges.directedEdge(a, b);
                        Edge edge2 = Edges.directedEdge(c, b);

                        if (getKnowledge().isForbidden(edge1.getNode1().getName(), edge1.getNode2().getName())) {
                            continue;
                        }

                        if (getKnowledge().isForbidden(edge2.getNode1().getName(), edge2.getNode2().getName())) {
                            continue;
                        }

                        moves.add(new Move(edge1, edge2, Move.Type.ADD_COLLIDER));
                    }
                }
//            }
        }

        return moves;
    }

    private List<Move> getSwapMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        for (Node b : graph.getNodes()) {
            List<Node> adj = graph.getAdjacentNodes(b);

            if (adj.size() < 2) continue;

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> set = GraphUtils.asList(choice, adj);

                Node a = set.get(0);
                Node c = set.get(1);

                if (graph.getEdge(a, b) != null && graph.getEdge(b, c) != null &&
                        graph.getEdge(a, b).pointsTowards(b) && graph.getEdge(b, c).pointsTowards(c)) {
                    moves.add(new Move(Edges.directedEdge(a, b), Edges.directedEdge(b, c), Move.Type.SWAP));
                } else if (graph.getEdge(b, a) != null && graph.getEdge(a, c) != null &&
                        graph.getEdge(b, a).pointsTowards(a) && graph.getEdge(a, c).pointsTowards(c)) {
                    moves.add(new Move(Edges.directedEdge(b, a), Edges.directedEdge(a, c), Move.Type.SWAP));
                }
            }
        }

        return moves;
    }

    private List<Move> getRemoveTriangleMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        for (Node b : graph.getNodes()) {
            List<Node> adj = graph.getAdjacentNodes(b);

            if (adj.size() < 2) continue;

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> set = GraphUtils.asList(choice, adj);

                Node a = set.get(0);
                Node c = set.get(1);

                Edge edge1 = graph.getEdge(a, b);
                Edge edge2 = graph.getEdge(b, c);
                Edge edge3 = graph.getEdge(a, c);

                if (edge1 != null && edge2 != null && edge3 != null &&
                        edge1.pointsTowards(a) && edge3.pointsTowards(c) &&
                        edge2.pointsTowards(c)) {
                    moves.add(new Move(Edges.directedEdge(b, c), Edges.directedEdge(c, a), Move.Type.SWAP));
                } else if (edge1 != null && edge2 != null && edge3 != null &&
                        edge3.pointsTowards(a) && edge1.pointsTowards(b) &&
                        edge2.pointsTowards(b)) {
                    moves.add(new Move(Edges.directedEdge(b, c), Edges.directedEdge(b, a), Move.Type.SWAP));
                }
            }
        }

        return moves;
    }

    private List<Move> getRemoveColliderMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        for (Node b : graph.getNodes()) {
            List<Node> adj = graph.getAdjacentNodes(b);

            if (adj.size() < 2) continue;

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> set = GraphUtils.asList(choice, adj);

                Node a = set.get(0);
                Node c = set.get(1);

                if (graph.isDefCollider(a, b, c)) {
                    Edge edge1 = Edges.directedEdge(a, b);
                    Edge edge2 = Edges.directedEdge(c, b);

                    moves.add(new Move(edge1, edge2, Move.Type.REMOVE_COLLIDER));
                }
            }
        }

        return moves;
    }

    private List<Move> getDoubleRemoveMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();
        List<Edge> edges = new ArrayList<>(graph.getEdges());

        // Remove moves:
        for (int i = 0; i < edges.size(); i++) {
            for (int j = i + 1; j < edges.size(); j++) {
                moves.add(new Move(edges.get(i), edges.get(j), Move.Type.DOUBLE_REMOVE));
            }
        }

        return moves;
    }







    private void addRequiredEdges(Graph graph) {
        for (Iterator<KnowledgeEdge> it =
             this.getKnowledge().requiredEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge next = it.next();
            String a = next.getFrom();
            String b = next.getTo();
            Node nodeA = null, nodeB = null;
            Iterator<Node> itn = graph.getNodes().iterator();
            while (itn.hasNext() && (nodeA == null || nodeB == null)) {
                Node nextNode = itn.next();
                if (nextNode.getName().equals(a)) {
                    nodeA = nextNode;
                }
                if (nextNode.getName().equals(b)) {
                    nodeB = nextNode;
                }
            }
//            if (!graph.isAncestorOf(nodeB, nodeA)) {
            graph.removeEdge(nodeA, nodeB);
            graph.addDirectedEdge(nodeA, nodeB);
//                TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeA, nodeB));
//            }
        }
        for (Iterator<KnowledgeEdge> it =
             getKnowledge().forbiddenEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge next = it.next();
            String a = next.getFrom();
            String b = next.getTo();
            Node nodeA = null, nodeB = null;
            Iterator<Node> itn = graph.getNodes().iterator();
            while (itn.hasNext() && (nodeA == null || nodeB == null)) {
                Node nextNode = itn.next();
                if (nextNode.getName().equals(a)) {
                    nodeA = nextNode;
                }
                if (nextNode.getName().equals(b)) {
                    nodeB = nextNode;
                }
            }
            if (nodeA != null && nodeB != null && graph.isAdjacentTo(nodeA, nodeB) &&
                    !graph.isChildOf(nodeA, nodeB)) {
                if (!graph.isAncestorOf(nodeA, nodeB)) {
                    graph.removeEdges(nodeA, nodeB);
                    graph.addDirectedEdge(nodeB, nodeA);
//                    TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeB, nodeA));
                }
            }
        }
    }

    public static class Score {
        private final Scorer scorer;

        public Score(Scorer scorer) {
            this.scorer = scorer;
        }

        public double getScore() {
            return scorer.getScore();
        }
    }

    private static class Move {
        private final Edge edge;
        private final Type type;
        private Edge secondEdge;

        public Move(Edge edge, Type type) {
            this.edge = edge;
            this.type = type;
        }

        public Move(Edge edge, Edge secondEdge, Type type) {
            this.edge = edge;
            this.secondEdge = secondEdge;
            this.type = type;
        }

        public Edge getFirstEdge() {
            return this.edge;
        }

        public Edge getSecondEdge() {
            return secondEdge;
        }

        public Type getType() {
            return this.type;
        }

        public String toString() {
            String s = (secondEdge != null) ? (secondEdge + ", ") : "";
            return "<" + edge + ", " + s + type + ">";

        }

        public boolean equals(Object o) {
            if ((!(o instanceof Move))) return false;
            if (o == this) return true;
            Move m = (Move) o;
            return m.type.equals(type) && m.edge.equals(edge) && (secondEdge == null || m.secondEdge.equals(secondEdge));
        }

        public enum Type {
            ADD, REMOVE, REDIRECT, ADD_COLLIDER, REMOVE_COLLIDER, SWAP, DOUBLE_REMOVE
        }
    }
}
