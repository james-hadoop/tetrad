package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.PermutationGenerator;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Double.POSITIVE_INFINITY;


/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class Boss2 {
    private final List<Node> variables;
    private Score score;
    private IndependenceTest test;
    private boolean cachingScores = true;
    private int numStarts = 1;
    private Method method = Method.BOSS;
    private boolean verbose = false;
    private boolean breakTies = false;
    private TeyssierScorer.ScoreType scoreType = TeyssierScorer.ScoreType.Edge;
    private IKnowledge knowledge = new Knowledge2();

    public Boss2(Score score) {
        this.score = score;
        this.variables = new ArrayList<>(score.getVariables());
    }

    public Boss2(IndependenceTest test) {
        this.test = test;
        this.variables = new ArrayList<>(test.getVariables());
    }

    public List<Node> bestOrder(List<Node> order) {

        long start = System.currentTimeMillis();
        order = new ArrayList<>(order);


        TeyssierScorer scorer;

        if (score != null && !(score instanceof GraphScore)) {
            scorer = new TeyssierScorer(score);
        } else if (test != null) {
            scorer = new TeyssierScorer(test);
        } else {
            throw new IllegalArgumentException("Need a score (not GraphScore) or a test.");
        }


        scorer.setKnowledge(knowledge);
        scorer.setScoreType(scoreType);

        scorer.setCachingScores(cachingScores);

        double best = POSITIVE_INFINITY;
        List<Node> bestPerm = null;

//        if (knowledge != null) {
        makeValidKnowledgeOrder(order);
//        }

        System.out.println("valid knowledge order = " + order);
        scorer.score(order);

        for (int r = 0; r < numStarts; r++) {
            if (r > 0) scorer.shuffleVariables();

            List<Node> perm;

            if (method == Method.BOSS) {
                perm = boss(scorer);
            } else if (method == Method.SP) {
                perm = sp(scorer);
            } else if (method == Method.GSP) {
                perm = gsp(scorer);
            } else {
                throw new IllegalArgumentException("Unrecognized method: " + method);
            }

            if (scorer.score() < best) {
                best = scorer.score();
                bestPerm = perm;
            }
        }

        long stop = System.currentTimeMillis();

        if (verbose) {
            System.out.println("Final order = " + scorer.getOrder());
            System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
        }

        return bestPerm;
    }

    private void makeValidKnowledgeOrder(List<Node> order) {
        order.sort((o1, o2) -> {
            if (o1.getName().equals(o2.getName())) {
                return 0;
            } else if (knowledge.isRequired(o1.getName(), o2.getName())) {
                return 1;
            } else if (knowledge.isRequired(o2.getName(), o1.getName())) {
                return -1;
            } else if (knowledge.isForbidden(o2.getName(), o1.getName())) {
                return -1;
            } else if (knowledge.isForbidden(o1.getName(), o2.getName())) {
                return 1;
            } else {
                return 1;
            }
        });
    }

    private boolean satisfiesKnowledge(List<Node> order) {
        if (!knowledge.isEmpty()) {
            for (int i = 0; i < order.size(); i++) {
                for (int j = i + 1; j < order.size(); j++) {
                    if (knowledge.isForbidden(order.get(i).getName(), order.get(j).getName())) {
                        return false;
                    }
                }
            }

        }

        return true;
    }

    public List<Node> boss(TeyssierScorer scorer) {

        // Take each variable in turn and try moving it to each position to the left (i.e.,
        // try promoting it in the causal order). Score each causal order by building
        // a DAG using something like the K2 method. (We actually use Grow-Shrink, a
        // Markov blanket algorithm, to find putative parents for each node.) Finally,
        // place the node in whichever position yielded the highest score. Do this
        // for each variable in turn. Once you're done, do it all again, until no more
        // variables can be relocated.

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Initial)");
        }

        if (breakTies) {
            do {
                while (relocateLoop(scorer)) ;
            } while (twoStepLoop(scorer));
        } else {
            while (relocateLoop(scorer)) ;
        }

        return scorer.getOrder();
    }

    // For ties choose the most senior.
    private boolean relocateLoop(TeyssierScorer scorer) {
        scorer.bookmark();

        double score;

        do {
            score = scorer.score();

            for (Node v : scorer.getOrder()) {
                scorer.moveToLast(v);

                double score2 = POSITIVE_INFINITY;

                if (scorer.score() < score2 && satisfiesKnowledge(scorer.getOrder())) {
                    scorer.bookmark();
                    score2 = scorer.score();
                }

                while (scorer.promote(v)) {
                    if (scorer.score() < score2 && satisfiesKnowledge(scorer.getOrder())) {
                        scorer.bookmark();
                        score2 = scorer.score();
                    }
                }

                scorer.goToBookmark();
            }
        } while (scorer.score() < score);

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (Single Moves)");
        }

        return scorer.score() < score;
    }

    private boolean twoStepLoop(TeyssierScorer scorer) {
        List<Node> b = scorer.getOrder();

        for (Node v : b) {
            for (int _r1 = 0; _r1 < b.size(); _r1++) {
                for (int _r2 = _r1 + 1; _r2 < b.size(); _r2++) {
                    Node r1 = scorer.get(_r1);
                    Node r2 = scorer.get(_r2);

                    if (((scorer.getParents(v).contains(r1) && scorer.getParents(v).contains(r2))
                            || (scorer.getParents(r1).contains(v) && scorer.getParents(r2).contains(v)))
                            && triangle(v, r1, r2, scorer) && satisfiesKnowledge(scorer.getOrder())) {
                        double score = scorer.score();
                        scorer.bookmark();

                        scorer.swap(v, r1);
                        scorer.swap(v, r2);

                        if (scorer.score() < score && satisfiesKnowledge(scorer.getOrder())) {
                            return true;
                        }

                        scorer.goToBookmark();
                    }
                }
            }
        }

        if (verbose) {
            System.out.println("# Edges = " + scorer.getNumEdges() + " Score = " + scorer.score() + " (" + "Simultaneous Moves" + ")");
        }

        return false;
    }

    private boolean triangle(Node v, Node r1, Node r2, TeyssierScorer scorer) {
        return adjacent(v, r1, scorer) && adjacent(v, r2, scorer) && adjacent(r1, r2, scorer);
    }

    private boolean adjacent(Node v, Node r1, TeyssierScorer scorer) {
        return scorer.getParents(v).contains(r1) || scorer.getParents(r1).contains(v);
    }

    private List<Node> gsp(TeyssierScorer scorer) {
        int maxDepth = Integer.MAX_VALUE;

        double score = scorer.score();

        while (true) {
            gspVisit(scorer, 0, maxDepth, new HashSet<>(), null);

            if (scorer.score() == score) {
                break;
            }

            score = scorer.score();
        }

        return scorer.getOrder();
    }

    private boolean gspVisit(TeyssierScorer scorer, int depth, int maxDepth, Set<Node> path, Node ww) {
        if (depth > maxDepth) return false;
        if (path.contains(ww)) return false;
        System.out.println("GSP visit depth = " + depth);
        path.add(ww);

        double s = scorer.score();

        for (Node w : scorer.getOrder()) {
            Set<Node> parentsw = scorer.getParents(scorer.indexOf(w));

            for (Node v : parentsw) {
                Set<Node> parentsv = scorer.getParents(scorer.indexOf(v));
                Set<Node> _parentsw = new HashSet<>(parentsw);
                _parentsw.remove(v);

                if (parentsv.equals(_parentsw)) {
                    scorer.swap(v, w);

                    if (scorer.score() < s) {
                        return gspVisit(scorer, depth + 1, maxDepth, path, w);
                    }

                    scorer.swap(v, w);
                }
            }
        }

        path.remove(ww);
        return false;
    }

    public List<Node> sp(TeyssierScorer scorer) {
        double minScore = POSITIVE_INFINITY;
        List<Node> minP = null;

        List<Node> variables = scorer.getOrder();
        PermutationGenerator gen = new PermutationGenerator(variables.size());
        int[] perm;

        while ((perm = gen.next()) != null) {
            List<Node> p = GraphUtils.asList(perm, variables);
            double score = scorer.score(p);

            if (score < minScore) {
                minScore = score;
                minP = p;
            }
        }

        return minP;
    }

    @NotNull
    public Graph getGraph(List<Node> order, boolean cpdag) {
        TeyssierScorer scorer;

        if (score != null) {
            scorer = new TeyssierScorer(score);
        } else {
            scorer = new TeyssierScorer(test);
        }

        scorer.setScoreType(scoreType);

        scorer.score(order);
        return scorer.getGraph(cpdag);
    }

    public void setCacheScores(boolean cachingScores) {
        this.cachingScores = cachingScores;
    }

    public void setNumStarts(int numStarts) {
        this.numStarts = numStarts;
    }

    public Method getMethod() {
        return method;
    }

    public void setMethod(Method method) {
        this.method = method;
    }

    public List<Node> getVariables() {
        return this.variables;
    }

    public void setBreakTies(boolean breakTies) {
        this.breakTies = breakTies;
    }

    public void setScoreType(TeyssierScorer.ScoreType scoreType) {
        this.scoreType = scoreType;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public enum Method {BOSS, SP, GSP}
}
