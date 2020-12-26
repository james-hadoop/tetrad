package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.RandomUtil;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.io.File;
import java.util.*;

/**

 */
public class TestEIModel {

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());
        Graph g = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        Map<Edge, Double> times = getRandomTimes(g);
        Map<Edge, Integer> excitements = getRandomExcitements(g);

        Node x5 = node(g, "X5");
        Node x15 = node(g, "X15");
        List<Node> cond = new ArrayList<>();

        EIModel wellen = new EIModel(g, times, excitements, 300);

        int wellenPrediction = wellen.getWellenPrediction(x5, x15, cond);

        switch (wellenPrediction) {
            case 1:
                System.out.println("Excitatory");
                break;
            case 2:
                System.out.println("Inhibitory");
                break;
            case 3:
                System.out.println("Mixed");
                break;
            default:
                throw new IllegalStateException("Prediction should be 1, 2, 3: " + wellenPrediction);
        }
    }

    @NotNull
    private Map<Edge, Double> getRandomTimes(Graph graph) {
        Map<Edge, Double> times = new HashMap<>();

        for (Edge edge : graph.getEdges()) {
            times.put(edge, RandomUtil.getInstance().nextInt(60) + 20.);
        }
        return times;
    }

    @NotNull
    private Map<Edge, Integer> getRandomExcitements(Graph graph) {
        Map<Edge, Integer> excitements = new HashMap<>();

        for (Edge edge : graph.getEdges()) {
            if (RandomUtil.getInstance().nextDouble() < 0.2) {
                excitements.put(edge, 2);
            } else {
                excitements.put(edge, 1);
            }
        }

        return excitements;
    }

    public Node node(Graph graph, String x1) {
        return graph.getNode(x1);
    }

}
