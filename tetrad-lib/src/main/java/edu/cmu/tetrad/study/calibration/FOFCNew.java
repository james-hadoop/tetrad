package edu.cmu.tetrad.study.calibration;

import edu.cmu.tetrad.data.CorrelationMatrix;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphNode;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.RandomUtil;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;

import java.io.File;
import java.util.*;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class FOFCNew {

    // From the constructor.
    private DataSet data;
    private IndependenceTest test;
    private ContinuousTetradTest tetradTest;
    private double alpha;

    // Internal, don't expose.
    private CorrelationMatrix corr;
    private DeltaTetradTest deltaTest;
    private Set<CommonCause> C = new HashSet<>();
    private int latentIndex = 0;

    public FOFCNew(DataSet data, IndependenceTest test, ContinuousTetradTest tetradTest, double alpha) {
        this.data = data;
        this.test = test;
        this.tetradTest = tetradTest;
        this.alpha = alpha;

        // Internal.
        this.corr = new CorrelationMatrix(data);
        this.deltaTest = new DeltaTetradTest(data);
    }

    public Graph search() {
        PcStable pc = new PcStable(test);
        Graph pcGraph = pc.search();

        List<Node> variables = pcGraph.getNodes();

        ChoiceGenerator gen = new ChoiceGenerator(variables.size(), 3);
        int[] choice;

        C1:
        while ((choice = gen.next()) != null) {

            int i = choice[0];
            int j = choice[1];
            int k = choice[2];

            Node vi = variables.get(i);
            Node vj = variables.get(j);
            Node vk = variables.get(k);

            List<Node> other = new ArrayList<>(variables);
            other.remove(vi);
            other.remove(vj);
            other.remove(vk);

            // Usual purity check.
            for (int l = 0; l < variables.size(); l++) {
                if (l == i) continue;
                if (l == j) continue;
                if (l == k) continue;

                if (!threeTetrads(i, j, k, l)) {
                    continue C1;
                }
            }

            // step 1a, genius second part.
            ChoiceGenerator gen2 = new ChoiceGenerator(other.size(), 2);
            int[] choice2;

            while ((choice2 = gen2.next()) != null) {

                int l = variables.indexOf(other.get(choice2[0]));
                int m = variables.indexOf(other.get(choice2[1]));

                if (threeTetrads(i, j, l, m)) continue C1;
                if (threeTetrads(i, k, l, m)) continue C1;
                if (threeTetrads(j, k, l, m)) continue C1;

                C.add(new CommonCause(new GraphNode("L" + ++latentIndex), vi, vj, vk));
                break;
            }
        }

        // Step 1c, beautiful.
        gen = new ChoiceGenerator(variables.size(), 2);

        C2:
        while ((choice = gen.next()) != null) {
            int i = choice[0];
            int j = choice[1];

            Node vi = variables.get(i);
            Node vj = variables.get(j);

            List<Node> other = new ArrayList<>(variables);
            other.remove(vi);
            other.remove(vj);

            if (correlation(vi, vj)) continue;

            ChoiceGenerator gen2 = new ChoiceGenerator(other.size(), 2);
            int[] choice2;

            while ((choice2 = gen2.next()) != null) {
                int k = variables.indexOf(other.get(choice2[0]));
                int l = variables.indexOf(other.get(choice2[1]));

                if (oneTetradTest(i, j, k, l)) continue C2;

                C.add(new CommonCause(new GraphNode("L" + ++latentIndex), vi, vj));
                break;
            }
        }

        // Step 1d seems wrong; mixed clusters will satisfy this. Is OneTetradTest wrong?

        // step 2, doesn't work. (Generates a bunch of mixed clusters.)
        ChoiceGenerator gen2 = new ChoiceGenerator(variables.size(), 4);
        int[] choice2;

        while ((choice2 = gen2.next()) != null) {

            int i = choice2[0];
            int j = choice2[1];
            int k = choice2[2];
            int l = choice2[3];

            Node vi = variables.get(i);
            Node vj = variables.get(j);
            Node vk = variables.get(k);
            Node vl = variables.get(l);

            if (tetradTest.tetradPValue(i, k, j, l) > alpha) {
                CommonCause c1 = new CommonCause(new GraphNode("L" + ++latentIndex), vi, vl);
                CommonCause c2 = new CommonCause(new GraphNode("L" + ++latentIndex), vj, vk);

                if (!C.contains(c1) || !C.contains(c2)) {
//                    CI.add(c1);
//                    CI.add(c2);
                }
            }
        }

        boolean changed = true;

        // step 3, combine clusters that share two measures in common.
        while (changed) {
            changed = false;

            List<CommonCause> CIL = new ArrayList<>(C);

            C1:
            for (int s = 0; s < CIL.size(); s++) {
                for (int t = s + 1; t < CIL.size(); t++) {
                    CommonCause c1 = CIL.get(s);
                    CommonCause c2 = CIL.get(t);

                    Set<Node> intersection = c1.getNodes();
                    intersection.retainAll(c2.getNodes());

                    if (intersection.size() >= 2) {
                        C.remove(c1);
                        C.remove(c2);

                        Set<Node> union = new HashSet<>();
                        union.addAll(c1.getNodes());
                        union.addAll(c2.getNodes());

                        C.add(new CommonCause(new GraphNode("L" + ++latentIndex), union));
                        changed = true;
                        continue C1;
                    }
                }
            }
        }

        latentIndex = 0;

        for (CommonCause c : C) {
            c.rename("L" + ++latentIndex);
        }

        // debugging check
        for (CommonCause c : C) {
            System.out.println(c);
        }


        Graph graph = GraphUtils.emptyGraph(data.getNumColumns());
        graph = GraphUtils.replaceNodes(graph, data.getVariables());

        for (CommonCause c : C) {
            graph.addNode(c.getLatent());

            for (Node n : c.getNodes()) {
                graph.addDirectedEdge(c.getLatent(), n);
            }
        }

        List<CommonCause> cl = new ArrayList<>(C);

        for (int s = 0; s < cl.size(); s++) {
            for (int t = s + 1; t < cl.size(); t++) {
                graph.addNondirectedEdge(cl.get(s).getLatent(), cl.get(t).getLatent());
            }
        }

        boolean correlated;

        for (int s = 0; s < cl.size(); s++) {
            CommonCause c1 = cl.get(s);
            Set<Node> n1Indices = c1.getNodes();

            correlated = false;

            for (int t = s + 1; t < cl.size(); t++) {
                CommonCause c2 = cl.get(t);
                Set<Node> n2Indices = c2.getNodes();

                for (Node n1Index : n1Indices) {
                    for (Node n2Index : n2Indices) {
                        if (correlation(n1Index, n2Index)) {
                            correlated = true;
                            break;
                        }

                    }

                    if (correlated) {
                        break;
                    }
                }

                if (!correlated) {
                    graph.removeEdge(c1.getLatent(), c2.getLatent());
                } else {
                    correlated = false;
                }

            }
        }

        for (int s = 0; s < cl.size(); s++) {
            CommonCause c1 = cl.get(s);
            Set<Node> n1Indices = c1.getNodes();

            for (int t = s + 1; t < cl.size() - 1; t++) {
                CommonCause c2 = cl.get(t);
                Set<Node> n2Indices = c2.getNodes();

                for (int u = t + 1; u < cl.size(); u++) {
                    CommonCause c3 = cl.get(u);
                    Set<Node> n3Indices = c3.getNodes();

                    if (graph.isAdjacentTo(c1.getLatent(), c2.getLatent()) && graph.isAdjacentTo(c2.getLatent(), c3.getLatent()) && !graph.isAdjacentTo(c1.getLatent(), c3.getLatent())) {

                        if (collider(n1Indices, n2Indices, n3Indices)) {

                            graph.removeEdge(c1.getLatent(), c2.getLatent());
                            graph.addDirectedEdge(c1.getLatent(), c2.getLatent());

                            graph.removeEdge(c3.getLatent(), c2.getLatent());
                            graph.addDirectedEdge(c3.getLatent(), c2.getLatent());
                        }

                    }

                    if (graph.isAdjacentTo(c1.getLatent(), c2.getLatent()) && graph.isAdjacentTo(c1.getLatent(), c3.getLatent()) && !graph.isAdjacentTo(c2.getLatent(), c3.getLatent())) {

                        if (collider(n2Indices, n1Indices, n3Indices)) {

                            graph.removeEdge(c1.getLatent(), c2.getLatent());
                            graph.addDirectedEdge(c2.getLatent(), c1.getLatent());

                            graph.removeEdge(c3.getLatent(), c1.getLatent());
                            graph.addDirectedEdge(c3.getLatent(), c1.getLatent());
                        }

                    }

                    if (graph.isAdjacentTo(c1.getLatent(), c3.getLatent()) && graph.isAdjacentTo(c2.getLatent(), c3.getLatent()) && !graph.isAdjacentTo(c1.getLatent(), c2.getLatent())) {

                        if (collider(n1Indices, n3Indices, n2Indices)) {

                            graph.removeEdge(c1.getLatent(), c3.getLatent());
                            graph.addDirectedEdge(c1.getLatent(), c3.getLatent());

                            graph.removeEdge(c3.getLatent(), c2.getLatent());
                            graph.addDirectedEdge(c2.getLatent(), c3.getLatent());
                        }

                    }
                }
            }

        }

        return graph;
    }

    private boolean collider(Set<Node> nodes1, Set<Node> nodes2, Set<Node> nodes3) {

        boolean correlated1 = false;
        boolean correlated2 = false;

        for (Node n1 : nodes1) {
            for (Node n2 : nodes2) {
                if (correlation(n1, n2)) {
                    correlated1 = true;
                    break;
                }

            }

            if (correlated1) {
                break;
            }
        }

        for (Node n3 : nodes3) {
            for (Node n2 : nodes2) {
                if (correlation(n3, n2)) {
                    correlated2 = true;
                    break;
                }

            }

            if (correlated2) {
                break;
            }
        }

        return (correlated1 && correlated2);
    }

    private boolean threeTetrads(int i, int j, int k, int l) {
        if (zeroCorr(i, j, k, l)) {
            return false;
        }

        List<Node> variables = data.getVariables();
        Tetrad t1 = new Tetrad(variables.get(i), variables.get(j), variables.get(k), variables.get(l));
        return deltaTest.getPValue(t1) > alpha;
    }

    private boolean oneTetradTest(int i, int j, int k, int l) {
        return tetradTest.tetradPValue(i, j, k, l) > alpha
                || tetradTest.tetradPValue(i, k, j, l) > alpha
                || tetradTest.tetradPValue(i, j, l, k) > alpha;
    }

    public boolean correlation(Node node1, Node node2) {

        double r = this.corr.getValue(data.getColumn(node1), data.getColumn(node2));
        int N = this.corr.getSampleSize();
        double f = sqrt(N) * Math.log((1. + r) / (1. - r));
        double p = 2.0 * (1.0 - RandomUtil.getInstance().normalCdf(0, 1, abs(f)));
        return (p > alpha);

    }

    private boolean zeroCorr(int i, int j, int k, int l) {
        int count = 0;

        ArrayList<Integer> cluster = new ArrayList<>();
        cluster.add(i);
        cluster.add(j);
        cluster.add(k);
        cluster.add(l);

        for (int m = 0; m < cluster.size(); m++) {
            for (int n = m + 1; n < cluster.size(); n++) {
                double r = this.corr.getValue(cluster.get(m), cluster.get(n));
                int N = this.corr.getSampleSize();
                double f = sqrt(N) * Math.log((1. + r) / (1. - r));
                double p = 2.0 * (1.0 - RandomUtil.getInstance().normalCdf(0, 1, abs(f)));
                if (p > alpha) count++;
            }
        }

        return count >= 1;
    }

    private static class CommonCause {
        private Node latent;
        private Set<Node> nodes;

        public CommonCause(Node latent, Node... nodes) {
            this.latent = latent;
            Set<Node> _nodes = new HashSet<>();
            Collections.addAll(_nodes, nodes);
            if (_nodes.size() < 2)
                throw new IllegalArgumentException("It's not a common cause if it's not a parent of at least two measures.");
            this.nodes = _nodes;
        }

        public CommonCause(Node latent, Set<Node> nodes) {
            if (nodes.size() < 2)
                throw new IllegalArgumentException("It's not a common cause if it's not a parent of at least two measures.");
            this.latent = latent;
            this.nodes = new HashSet<>(nodes);
        }

        public String toString() {
            StringBuilder s = new StringBuilder();

            for (Node n : this.nodes) {
                s.append(n).append("\t");
            }

            return s.toString();
        }

        public boolean equals(Object o) {
            if (!(o instanceof CommonCause)) return false;
            return getNodes().equals(((CommonCause) o).getNodes());
        }

        public int hashCode() {
            return getNodes().hashCode();
        }

        public Set<Node> getNodes() {
            return new HashSet<>(nodes);
        }

        public Node getLatent() {
            return latent;
        }

        public void rename(String s) {
            latent.setName(s);
        }
    }

    public static void main(String... args) {
//        File dataFile = new File("/Users/user/downloads/latentmodel1.txt");
        File dataFile = new File("/Users/user/Downloads/11.txt");
        ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(dataFile.toPath(), Delimiter.TAB);
        try {
            Data dataSet = reader.readInData();
            DataModel model = DataConvertUtils.toDataModel(dataSet);
            DataSet data = (DataSet) model;

            double alpha = .1;

            IndTestFisherZ fisherZ = new IndTestFisherZ(data, alpha);
            ContinuousTetradTest test = new ContinuousTetradTest(data, TestType.TETRAD_WISHART, alpha);

            FOFCNew fofc = new FOFCNew(data, fisherZ, test, alpha);
            Graph graph = fofc.search();

            System.out.println(graph);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
