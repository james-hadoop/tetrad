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

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.*;

import java.util.*;

/**
 * Heuristic Best Significant Model Search using a beam search.
 *
 * @author Joseph Ramsey
 */
public final class HbsmsBeam2 implements Hbsms {
    private final Graph initialGraph;
    private final Scorer scorer;
    private final CovarianceMatrix cov;
    private IKnowledge knowledge;
    private Graph graph;
    private double alpha = 0.05;
    private SemIm originalSemIm;
    private SemIm newSemIm;
    private int beamWidth = 1;

    public HbsmsBeam2(Graph graph, DataSet data, IKnowledge knowledge) {
        if (graph == null) graph = new EdgeListGraph(data.getVariables());

        this.knowledge = knowledge;
        this.graph = graph;
        this.initialGraph = new EdgeListGraph(graph);
        CovarianceMatrix cov = new CovarianceMatrix(data);
        this.cov = cov;
        this.scorer = new FmlBicScorer(cov);
    }

    public Graph search() {
        EdgeListGraph _graph = new EdgeListGraph(initialGraph);
        addRequiredEdges(_graph);
        Graph best = SearchGraphUtils.dagFromPattern(_graph);

        Score score0 = scoreGraph(best, scorer);
        MeekRules meekRules = new MeekRules();
        meekRules.setKnowledge(getKnowledge());

        best = decreaseScoreLoop2(best);

        Score score = scoreGraph(best, scorer);

        SemPm pm = new SemPm(best);
        SemEstimator est = new SemEstimator(cov, pm);
        this.newSemIm = est.getEstimatedSem();

        return best;
    }

    private Graph decreaseScoreLoop(Graph best) {
        Set<Graph> visited = new HashSet<>();

        double initialScore = scoreGraph(best, scorer).getScore();

        Map<Graph, Double> S = new HashMap<>();
        S.put(best, initialScore);

        LinkedList<Graph> bracket = new LinkedList<>();
        bracket.add(best);

        MOVES:
        while (true) {
            for (Graph g : bracket) {
                List<Move> moves = new ArrayList<>();
                moves.addAll(getRemoveMoves(g));
                moves.addAll(getRedirectMoves(g));
                moves.addAll(getAddMoves(g));

                for (Move move : moves) {
                    Graph g2 = makeMove(g, move);

                    if (visited.contains(g2)) continue;
                    visited.add(g2);

                    if (getKnowledge().isViolatedBy(g2)) continue;

                    Score _score = scoreGraph(g2, scorer);

                    if (bracket.size() != 1) throw new IllegalArgumentException();

                    double _s = _score.getScore();

                    if (decreasesScore(bracket, S, _s)) {
                        System.out.println("Inserting new score in top bracket (" + move.getType() + "): score = " + _s);

                        LinkedList<Graph> bracket2 = new LinkedList<>(bracket);

                        bracket.add(g2);
                        S.put(g2, _s);

                        if (bracket.size() > this.beamWidth) {
                            removeMaximalScore(bracket, S);
                        }

                        bracket = bracket2;

                        continue MOVES;
                    }
                }
            }

            break;
        }

        this.graph = maximumScore(S);
        return maximumScore(S);
    }

    private Graph decreaseScoreLoop2(Graph best) {
        double score0 = scoreGraph(best, scorer).getScore();

        MOVES:
        while (true) {
            List<Move> moves = new ArrayList<>();
            moves.addAll(getRemoveMoves(best));
            moves.addAll(getRedirectMoves(best));
            moves.addAll(getAddMoves(best));

            for (Move move : moves) {
                Graph g2 = makeMove(best, move);

                if (getKnowledge().isViolatedBy(g2)) continue;

                Score _score = scoreGraph(g2, scorer);

                if (_score.getScore() < score0) {
                    score0 = _score.getScore();
                    best = g2;

                    System.out.println("Decreases score (" + move.getType() + "): score = " + _score.getScore());

                    continue MOVES;
                }
            }

            break;
        }

        this.graph = best;
        return best;
    }

    private Graph decreaseScoreLoop3(Graph best) {
        double score0 = scoreGraph(best, scorer).getScore();

        LinkedList<Graph> bracket = new LinkedList<>();
        Map<Graph, Double> scores = new HashMap<>();
        bracket.add(best);
        scores.put(best, score0);

        MOVES:
        while (true) {
            Set<Graph> graphs = new HashSet<>();

            for (Graph graph : bracket) {
                List<Move> moves = new ArrayList<>();
                moves.addAll(getRemoveMoves(best));
                moves.addAll(getRedirectMoves(best));
                moves.addAll(getAddMoves(best));

                for (Move move : moves) {
                    graphs.add(makeMove(best, move));
                }
            }

            for (Graph g2 : graphs) {
                if (getKnowledge().isViolatedBy(g2)) continue;

                Score _score = scoreGraph(g2, scorer);

                if (decreasesScore(bracket, scores, _score.getScore())) {

//                if (_score.getScore() < score0) {
                    double score1 = _score.getScore();
                    scores.put(g2, score1);
                    bracket.add(g2);

                    removeMaximalScore(bracket, scores);

                    System.out.println("Decreases score " + _score.getScore());

                    continue MOVES;
                }
            }

            break;
        }

        this.graph = best;
        return best;
    }

    private boolean decreasesScore(LinkedList<Graph> bracket, Map<Graph, Double> s, double score) {
        double max = Double.NEGATIVE_INFINITY;

        for (Graph graph : bracket) {
            Double aDouble = s.get(graph);

            if (aDouble > max) {
                max = aDouble;
            }
        }

        return score < max;
    }

    private Graph maximumScore(Map<Graph, Double> s) {
        double maxScore = Double.NEGATIVE_INFINITY;
        Graph maxGraph = null;

        for (Graph graph : s.keySet()) {
            if (graph == null) {
                throw new NullPointerException();
            }

            double score = s.get(graph);

            if (score > maxScore) {
                maxScore = score;
                maxGraph = graph;
            }
        }

        return maxGraph;
    }

    private void removeMaximalScore(LinkedList<Graph> bracket, Map<Graph, Double> S) {
        double max = Double.NEGATIVE_INFINITY;
        Graph maxGraph = null;

        for (Graph graph : bracket) {
            if (S.get(graph) > max) {
                max = S.get(graph);
                maxGraph = graph;
            }
        }

        if (maxGraph != null && bracket.size() > 1) {
            bracket.remove(maxGraph);
        }
    }

    private Graph makeMove(Graph graph, Move move) {
        graph = new EdgeListGraph(graph);
        Edge firstEdge = move.getFirstEdge();

        if (firstEdge != null && move.getType() == Move.Type.ADD) {
            graph.removeEdge(firstEdge.getNode1(), firstEdge.getNode2());
            graph.addEdge(firstEdge);
        } else if (firstEdge != null && move.getType() == Move.Type.REMOVE) {
            graph.removeEdge(firstEdge);
        } else if (firstEdge != null && move.getType() == Move.Type.REDIRECT) {
            graph.removeEdge(graph.getEdge(firstEdge.getNode1(), firstEdge.getNode2()));
            graph.addEdge(firstEdge);
        }

        return graph;
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

                if (graph.getDegree(nodes.get(i)) > 0 || graph.getDegree(nodes.get(j)) > 0) {
                    Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
                    moves.add(new Move(edge, Move.Type.ADD));
                }

//                Edge edge = Edges.directedEdge(nodes.get(i), nodes.get(j));
//                moves.add(new Move(edge, Move.Type.ADD));
//                }
            }
        }

        return moves;
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

//            if (graph.isAncestorOf(j, i)) {
//                continue;
//            }

            moves.add(new Move(Edges.directedEdge(j, i), Move.Type.REDIRECT));
        }

        return moves;
    }

    public Graph getGraph() {
        return graph;
    }

    public SemIm getOriginalSemIm() {
        return originalSemIm;
    }

    public SemIm getNewSemIm() {
        return newSemIm;
    }

    public Score scoreGraph(Graph graph, Scorer scorer) {
        Graph dag = new EdgeListGraph(graph);// SearchGraphUtils.dagFromPattern(graph, getKnowledge());
        scorer.score(dag);
        return new Score(scorer);
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setBeamWidth(int beamWidth) {
        if (beamWidth < 1) throw new IllegalArgumentException();
        this.beamWidth = beamWidth;
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;

        if (knowledge.isViolatedBy(graph)) {
            throw new IllegalArgumentException("Graph violates knowledge.");
        }
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
            if (!graph.isAncestorOf(nodeB, nodeA)) {
                graph.removeEdge(nodeA, nodeB);
                graph.addDirectedEdge(nodeA, nodeB);
//                TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeA, nodeB));
            }
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
            ADD, REMOVE, REDIRECT, ADD_COLLIDER, REMOVE_COLLIDER, SWAP, DOUBLE_REMOVE
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
}


