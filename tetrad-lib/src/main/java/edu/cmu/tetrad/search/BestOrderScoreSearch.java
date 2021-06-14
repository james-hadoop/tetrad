package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
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
public class BestOrderScoreSearch {
    private final Score score;
    private final List<Node> variables;
    private final Map<ScoreProblem, Double> scoresHash = new HashMap<>();
    private IKnowledge knowledge = new Knowledge2();

    /**
     * Constructs a GSS search
     */
    public BestOrderScoreSearch(Score score) {
        this.score = score;
        this.variables = score.getVariables();
    }

    @NotNull
    public Graph search(List<Node> b0) {
        long start = System.currentTimeMillis();

        LinkedList<Node> br = new LinkedList<>(b0);

        // Modeled on an insertion sort. Each insertions a BIC comparison of models over the same
        // variables.
        for (int i = 0; i < br.size(); i++) {

            // Pick a various nodes beyond i and move it up into the i.
            for (int j = i; j < br.size(); j++) {
                Node v = br.get(j);
                move(br, v, 0);

                double[] scores = new double[i];

                for (int k = 0; k < i; k++) {
                    scores[k] = getScore(br, k);
                }

                double min = sumScores(scores);
                int b = -0;

                for (int l = 1; l < i; l++) {
                    swap(l, scores, br, i);

                    if (sumScores(scores) < min) {
                        b = l;
                        min = sumScores(scores);
                    }
                }

                move(br, v, b);
            }
        }

        long stop = System.currentTimeMillis();
        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");
        System.out.println("Final order = " + br);
        return SearchGraphUtils.patternForDag(new K3(score).search(br));
    }

    private void move(LinkedList<Node> br, Node v, int b) {
        br.remove(v);
        br.add(b, v);
    }

    private double sumScores(double[] scores) {
        double s2 = 0;
        for (double _s : scores) s2 += _s;
        return s2;
    }

    private double getScore(LinkedList<Node> br, int p) {
        Node n = br.get(p);
        Set<Node> pi = new HashSet<>();
        boolean changed = true;
        double s_new = POSITIVE_INFINITY;
        double s_node = score(n, new HashSet<>());

        while (changed) {
            changed = false;

            // Let z be the node that maximizes the score...
            Node z = null;

            for (Node z0 : new HashSet<>(getPrefix(br, p))) {
                if (pi.contains(z0)) continue;
                if (knowledge.isForbidden(z0.getName(), n.getName())) continue;
                pi.add(z0);

                double s2 = score(n, pi);

                if (s2 < s_new) {
                    s_new = s2;
                    z = z0;
                }

                pi.remove(z0);
            }

            if (s_new < s_node) {
                pi.add(z);
                s_node = s_new;
                changed = true;

                boolean changed2 = true;

                while (changed2) {
                    changed2 = false;

                    for (Node z0 : new HashSet<>(pi)) {
                        if (z0 == z) continue;
                        pi.remove(z0);

                        double s2 = score(n, pi);

                        if (s2 < s_new) {
                            s_new = s2;
                            z = z0;
                        }

                        pi.add(z0);
                    }

                    if (s_new < s_node) {
                        pi.remove(z);
                        s_node = s_new;
                        changed2 = true;
                    }
                }
            }
        }

        return s_node;
    }

    private void swap(int ww, double[] scores, LinkedList<Node> br, int prefixSize) {
        if (ww < 0) throw new IllegalArgumentException("ww < 1");
        if (ww >= prefixSize) throw new IllegalArgumentException("ww >= prefixSize");

        move(br, br.get(ww - 1), ww);
        recalculate(ww - 1, scores, br);
        recalculate(ww, scores, br);
    }

    private void recalculate(int x, double[] scores, LinkedList<Node> br) {
        double score = getScore(br, x);
        scores[x] = score;
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    private double score(Node n, Set<Node> pi) {
        if (scoresHash.containsKey(new ScoreProblem(n, pi))) {
            return scoresHash.get(new ScoreProblem(n, pi));
        }

        int[] parentIndices = new int[pi.size()];

        int k = 0;

        for (Node p : pi) {
            parentIndices[k++] = variables.indexOf(p);
        }

        double score = -this.score.localScore(variables.indexOf(n), parentIndices);

        scoresHash.put(new ScoreProblem(n, pi), score);

        return score;
    }

    public static class ScoreProblem {
        private final Node y;
        private final Set<Node> prefix;

        public ScoreProblem(Node y, Set<Node> prefix) {
            this.y = y;
            this.prefix = new HashSet<>(prefix);
        }

        public int hashCode() {
            return 17 * y.hashCode() + 3 * prefix.hashCode();
        }

        public boolean equals(Object o) {
            if (!(o instanceof ScoreProblem)) {
                return false;
            }

            ScoreProblem spec = (ScoreProblem) o;
            return spec.y.equals(this.y) && spec.getPrefix().equals(this.prefix);
        }

        public Node getY() {
            return y;
        }

        public Set<Node> getPrefix() {
            return prefix;
        }
    }
}
