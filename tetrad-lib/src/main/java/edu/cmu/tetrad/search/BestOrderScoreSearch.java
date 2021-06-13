package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.util.Collections.shuffle;


/**
 * Searches for a DAG by adding or removing directed edges starting
 * with an empty graph, using a score (default FML BIC). Implements
 * the Global Score Search (GSS) algorithm.
 *
 * @author josephramsey
 */
public class BestOrderScoreSearch {
    // The score used (default FML BIC). Lower is better.
    private final ForwardScore forwardScore;
    long stop = System.currentTimeMillis();
    private double scoreOriginalOrder = Double.NaN;
    private double scoreLearnedOrder = Double.NaN;
    private Method method = Method.NONRECURSIVE;
    private int numRestarts = 1;
    /**
     * Constructs a GSS search
     */
    public BestOrderScoreSearch(Score score) {
        this.forwardScore = new K3(score);
    }

    @NotNull
    public Graph search(List<Node> b0) {
        long start = System.currentTimeMillis();

        List<Node> br = new ArrayList<>(b0);

        for (int k = 0; k < numRestarts; k++) {

            for (int i = 0; i < b0.size(); i++) {
                List<Node> best = new ArrayList<>(br);
                double sss = forwardScore.score(getPrefix(br, i));

                // ...and move v to an optimal position.
                for (int j = i; j < br.size(); j++) {
                    List<Node> b1 = new ArrayList<>(br);
                    Node v = b1.get(j);
                    b1.remove(v);
                    b1.add(i, v);

                    for (int w = 0; w <= i; w++) {
                        List<Node> b12 = new ArrayList<>(b1);
                        Node v2 = b1.get(j);
                        b12.remove(v2);
                        b12.add(w, v2);

                        double ss = forwardScore.score(getPrefix(b12, i));

                        if (ss < sss) {
                            best = b12;
                            sss = ss;
                        }
                    }
                }

                br = best;
            }

            long stop = System.currentTimeMillis();
            System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

            b0 = br;
        }

        System.out.println("Final b0 = " + br);
        return forwardScore.search(b0);
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j <= i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }

    public Graph search2(List<Node> variables) {
        long start = System.currentTimeMillis();

        List<Node> b0 = new ArrayList<>(variables);
        double s0 = forwardScore.score(b0);
        int m = variables.size();

        scoreOriginalOrder = s0;

        boolean changed = true;

        List<Node> br = new ArrayList<>(b0);
        double sr = forwardScore.score(br);

        for (int r = 0; r < numRestarts; r++) {
            if (method == Method.NONRECURSIVE) {

                // Until you can't do it anymore...
                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(br);
                    double s2 = sr;

                    // ...pick a variable v...
                    for (int i = 0; i < m; i++) {
                        Node v = b0.get(i);

                        List<Node> bt = new ArrayList<>(b2);
                        double st = forwardScore.score(bt);

                        // ...and move v to an optimal position.
                        for (int j = 0; j < m; j++) {
                            List<Node> b1 = new ArrayList<>(b2);
                            b1.remove(v);
                            b1.add(j, v);

                            double s1 = forwardScore.score(b1);

                            if (s1 < st) {
                                st = s1;
                                bt = b1;
                            }
                        }

                        if (st < s2) {
                            s2 = st;
                            b2 = bt;
                        }
                    }

                    // Output the best order you find.
                    if (s2 < sr) {
                        br = b2;
                        sr = s2;
                        System.out.println("br = " + br);
                        changed = true;
                    }
                }
            } else if (method == Method.RECURSIVE) {
                br = bossRecursive(br, sr, variables);
                sr = forwardScore.score(br);
            } else {
                throw new IllegalStateException("Unrecognized method: " + method);
            }

            if (sr < s0) {
                s0 = sr;
                b0 = new ArrayList<>(br);
            }

            System.out.println("Completed round " + (r + 1) + " elapsed = " + (System.currentTimeMillis() - start) / 1000.0 + " s");
            System.out.println("br = " + br);

            shuffle(br);
            System.out.println("shufflng br");
            sr = forwardScore.score(br);
        }

        scoreLearnedOrder = s0;

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return forwardScore.search(b0);
    }

    public Graph search3(List<Node> variables) {
        long start = System.currentTimeMillis();

        List<Node> b0 = new ArrayList<>(variables);
        double s0 = forwardScore.score(b0);
        int m = variables.size();

        scoreOriginalOrder = s0;

        boolean changed = true;

        List<Node> br = new ArrayList<>(b0);
        double sr = forwardScore.score(br);

        for (int r = 0; r < numRestarts; r++) {
            if (method == Method.NONRECURSIVE) {

                // Until you can't do it anymore...
                while (changed) {
                    changed = false;

                    List<Node> b2 = new ArrayList<>(br);
                    double s2 = sr;

                    // ...pick a variable v...
                    for (int i = 0; i < m; i++) {
                        Node v = b0.get(i);

                        List<Node> bt = new ArrayList<>(b2);
                        double st = forwardScore.score(bt);

                        // ...and move v to an optimal position.
                        for (int j = 0; j < m; j++) {
                            List<Node> b1 = new ArrayList<>(b2);
                            b1.remove(v);
                            b1.add(j, v);

                            double s1 = forwardScore.score(b1);

                            if (s1 < st) {
                                st = s1;
                                bt = b1;
                            }
                        }

                        if (st < s2) {
                            s2 = st;
                            b2 = bt;
                        }
                    }

                    // Output the best order you find.
                    if (s2 < sr) {
                        br = b2;
                        sr = s2;
                        System.out.println("br = " + br);
                        changed = true;
                    }
                }
            } else if (method == Method.RECURSIVE) {
                br = bossRecursive(br, sr, variables);
                sr = forwardScore.score(br);
            } else {
                throw new IllegalStateException("Unrecognized method: " + method);
            }

            if (sr < s0) {
                s0 = sr;
                b0 = new ArrayList<>(br);
            }

            System.out.println("Completed round " + (r + 1) + " elapsed = " + (System.currentTimeMillis() - start) / 1000.0 + " s");
            System.out.println("br = " + br);

            shuffle(br);
            System.out.println("shufflng br");
            sr = forwardScore.score(br);
        }

        scoreLearnedOrder = s0;

        long stop = System.currentTimeMillis();

        System.out.println("BOSS Elapsed time = " + (stop - start) / 1000.0 + " s");

        return forwardScore.search(b0);
    }

    private List<Node> bossRecursive(List<Node> b0, double s0, List<Node> variables) {
        int m = variables.size();
        List<Node> b2 = new ArrayList<>(b0);

        // ...pick a variable v...
        for (Node v : variables) {

            List<Node> bt = new ArrayList<>(b2);
            double st = forwardScore.score(bt);

            // ...and move v to an optimal position.
            for (int j = 0; j < m; j++) {
                List<Node> b1 = new ArrayList<>(b2);
                b1.remove(v);
                b1.add(j, v);

                double s1 = forwardScore.score(b1);

                if (s1 < st) {
                    st = s1;
                    bt = b1;
                }
            }

            if (st < s0) {
                return bossRecursive(bt, st, variables);
            }
        }

        return b0;
    }

    public double getScoreOriginalOrder() {
        return scoreOriginalOrder;
    }

    public double getScoreLearnedOrder() {
        return scoreLearnedOrder;
    }

    public void setMethod(Method method) {
        this.method = method;
    }

    public void setNumRestarts(int numRestarts) {
        this.numRestarts = numRestarts;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.forwardScore.setKnowledge(knowledge);
    }

    public enum Method {NONRECURSIVE, RECURSIVE}
}
