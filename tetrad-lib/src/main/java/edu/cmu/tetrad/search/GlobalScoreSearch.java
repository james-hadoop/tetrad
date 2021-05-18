package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.Scorer;

import java.util.*;

/**
 * Searches for a DAG in a pattern by adding or removing directed edges starting
 * with an empty graph.
 *
 * @author josephramsey
 */
public class GlobalScoreSearch {
    private final Scorer scorer;
    private IKnowledge knowledge = new Knowledge2();

    public GlobalScoreSearch(Scorer scorer) {
        this.scorer = scorer;
    }

    public Graph search() {
        EdgeListGraph _graph = new EdgeListGraph(scorer.getVariables());
        addRequiredEdges(_graph);
        return decreaseScoreLoop(_graph);
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    private Graph decreaseScoreLoop(Graph graph) {
        double score0 = scoreGraph(graph, scorer).getScore();

        List<Move> addMoves = getAddMoves(graph);
        List<Move> addMovesRecall = new ArrayList<>(addMoves);

        MOVES:
        while (true) {
            for (Move move : getRemoveMoves(graph)) {
                Edge edge = makeMove(graph, move);

                if (knowledge.isViolatedBy(graph)) continue;

                Score _score = scoreGraph(graph, scorer);

                if (_score.getScore() < score0) {
                    score0 = _score.getScore();

                    System.out.println("Decreases score (" + move.getType() + "): score = " + _score.getScore());

                    addMoves.remove(move);
                    continue MOVES;
                }

                graph.addEdge(edge);
            }

            if (addMoves.isEmpty()) {
                addMoves = new ArrayList<>(addMovesRecall);
            }

            for (Move move : new LinkedList<>(addMoves)) {
                if (graph.isAdjacentTo(move.edge.getNode1(), move.edge.getNode2())) continue;
                if (graph.existsDirectedPathFromTo(move.edge.getNode2(), move.edge.getNode1())) continue;

                Edge edge = makeMove(graph, move);

                if (knowledge.isViolatedBy(graph)) continue;

                Score _score = scoreGraph(graph, scorer);

                if (_score.getScore() < score0) {
                    score0 = _score.getScore();

                    System.out.println("Decreases score (" + move.getType() + "): score = " + _score.getScore());

                    addMoves.remove(move);
                    continue MOVES;
                }

                graph.removeEdge(edge);
            }

            break;
        }

        return graph;
    }

    private Score scoreGraph(Graph graph, Scorer scorer) {
        scorer.score(graph);
        return new Score(scorer);
    }

    private Edge makeMove(Graph graph, Move move) {
        Edge firstEdge = move.getFirstEdge();

        if (firstEdge != null && move.getType() == Move.Type.ADD) {
            graph.addEdge(firstEdge);
            return firstEdge;
        } else if (firstEdge != null && move.getType() == Move.Type.REMOVE) {
            graph.removeEdge(firstEdge);
            return firstEdge;
        }

        throw new IllegalArgumentException();
    }

    private List<Move> getRemoveMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

        List<Edge> edges = new ArrayList<>(graph.getEdges());

        for (Edge edge : edges) {
            Node i = edge.getNode1();
            Node j = edge.getNode2();

            if (knowledge.isRequired(i.getName(), j.getName())) {
                continue;
            }

            moves.add(new Move(edge, Move.Type.REMOVE));
        }

        return moves;
    }

    private List<Move> getAddMoves(Graph graph) {
        List<Move> moves = new ArrayList<>();

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

                if (knowledge.isForbidden(nodes.get(i).getName(), nodes.get(j).getName())) {
                    continue;
                }

                if (knowledge.isRequired(nodes.get(j).getName(), nodes.get(i).getName())) {
                    continue;
                }

                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
                moves.add(new Move(edge, Move.Type.ADD));
            }
        }

        return moves;
    }

    private void addRequiredEdges(Graph graph) {
        for (Iterator<KnowledgeEdge> it =
            knowledge.requiredEdgesIterator(); it.hasNext(); ) {
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
            graph.removeEdge(nodeA, nodeB);
            graph.addDirectedEdge(nodeA, nodeB);
        }
        for (Iterator<KnowledgeEdge> it =
            knowledge.forbiddenEdgesIterator(); it.hasNext(); ) {
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
                }
            }
        }
    }

    private static class Move {
        private final Edge edge;
        private final Type type;

        public Move(Edge edge, Type type) {
            this.edge = edge;
            this.type = type;
        }

        public Edge getFirstEdge() {
            return this.edge;
        }

        public Type getType() {
            return this.type;
        }

        public String toString() {
            return "<" + edge + ", " + type + ">";
        }

        public boolean equals(Object o) {
            if ((!(o instanceof Move))) return false;
            if (o == this) return true;
            Move m = (Move) o;
            return m.type.equals(type) && m.edge.equals(edge);
        }

        public enum Type {
            ADD, REMOVE
        }
    }

    private static class Score {
        private final Scorer scorer;
        public Score(Scorer scorer) {
            this.scorer = scorer;
        }
        public double getScore() {
            return scorer.getScore();
        }
    }
}
