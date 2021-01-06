package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * Tests the EI model.
 *
 * @author jdramsey@andrew.cmu.edu 2021.1.4
 */
public class TestEIModel {

    private double cutoff;

    @Test
    public void test1() {
        setCutoff(0.8);

        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph7.txt"));

        EIModel.Records records = new EIModel.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    getRandomTime(), getRandomExcitement(getCutoff()));
        }

        try {
            records.toFile(new File("/Users/user/Downloads/records-acyclic.txt").toPath());
        } catch (IOException e) {
            e.printStackTrace();
        }

        EIModel.main("/Users/user/Downloads/records-acyclic.txt", "60", "X3", "X8", "X12");
    }

    @Test
    public void test2() {
        setCutoff(0.9);

        Graph graph = GraphUtils.loadGraphTxt(new File("src/test/resources/graph8.txt"));

        EIModel.Records records = new EIModel.Records();

        for (Edge edge : graph.getEdges()) {
            records.addRecord(edge.getNode1().toString(), edge.getNode2().toString(),
                    getRandomTime(), getRandomExcitement(getCutoff()));
        }

        try {
            records.toFile(new File("/Users/user/Downloads/records-cyclic.txt").toPath());
        } catch (IOException e) {
            e.printStackTrace();
        }

        EIModel.main("/Users/user/Downloads/records-cyclic.txt", "100", "X2", "X4");
    }

    private double getCutoff() {
        return cutoff;
    }

    private void setCutoff(double cutoff) {
        if (cutoff < 0 || cutoff > 1) throw new IllegalStateException("Cutoff should be in [0, 1]: " + cutoff);
        this.cutoff = cutoff;
    }

    private double getRandomTime() {
        return RandomUtil.getInstance().nextUniform(3, 20);
    }

    private int getRandomExcitement(double cutoff) {
        if (RandomUtil.getInstance().nextDouble() < cutoff) {
            return 1;
        } else {
            return RandomUtil.getInstance().nextDouble() < 0.2 ? 3 : 2;
        }
    }

}
