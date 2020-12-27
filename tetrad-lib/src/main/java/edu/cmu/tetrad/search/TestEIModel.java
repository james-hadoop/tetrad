package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 *
 */
public class TestEIModel {

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(38284848383L);

        System.out.println("seed = " + RandomUtil.getInstance().getSeed());
        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        EIModel.Records records = new EIModel.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    getRandomTime(), getRandomExcitement(.9));
        }

        try {
            records.toFile(new File("/Users/user/Downloads/records-ei-model.txt").toPath());
            records = EIModel.Records.fromFile(new File("/Users/user/Downloads/records-ei-model.txt").toPath());
            System.out.println(records);
        } catch (IOException e) {
            e.printStackTrace();
        }

        graph = records.getGraph();
        Map<Edge, Double> times = records.getTimes();
        Map<Edge, Integer> excitements = records.getExcitements();

        Node x5 = node(graph, "X5");
        Node x15 = node(graph, "X15");
        List<Node> cond = new ArrayList<>();
//        cond.add(node(graph, "X6"));

        EIModel wellen = new EIModel(graph, times, excitements, 100);

        int wellenPrediction = wellen.getWellenPrediction(x5, x15, cond);

        switch (wellenPrediction) {
            case 0:
                System.out.println("Disconnected");
                break;
            case 1:
                System.out.println("Excitatory");
                break;
            case 2:
                System.out.println("Inhibitory");
                break;
            case 3:
                System.out.println("Uncertain");
                break;
            default:
                throw new IllegalStateException("Prediction should be 0, 1, 2, or 3: " + wellenPrediction);
        }
    }

    private double getRandomTime() {
        return RandomUtil.getInstance().nextUniform(20, 80);
    }

    private int getRandomExcitement(double cutoff) {
        if (cutoff < 0 || cutoff > 1) throw new IllegalStateException("Cutoff should be in [0, 1]: " + cutoff);

        if (RandomUtil.getInstance().nextDouble() < cutoff) {
            return 1;
        } else {
            return 2;
        }
    }

    public Node node(Graph graph, String x1) {
        return graph.getNode(x1);
    }


}
