package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.*;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;


public class HillStats {
    private Graph trueGraph;
    private Graph estGraph;
    private List<Node> tp;
    private List<Node> tn;
    private List<Node> fp;
    private List<Node> fn;
    private List<Node> all;
    private List<Node> positives;
    private List<Node> negatives;

    public static DataSet mTorGold = null;

    public HillStats(Graph trueGraph, Graph estGraph) {

        this.trueGraph = trueGraph;
        this.estGraph = estGraph;

        if (mTorGold == null) {
            try {
                mTorGold = getDiscreteData(
                        "/Users/user/Box/data/4cellLineData/ground.truth/gold_descendants_AZD8055_mTOR2.txt",
                        Delimiter.TAB);
            } catch (IOException e) {
                e.printStackTrace();
            }

            System.out.println(mTorGold);

            for (int i = 0; i < 35; i++) {
                System.out.println(i + ". " + mTorGold.getInt(i, 27));
            }
        }

        Graph _trueGraph = GraphUtils.replaceNodes(this.getTrueGraph(), this.getEstGraph().getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(this.getEstGraph(), _trueGraph.getNodes());

        List<Node> allVars = new ArrayList<>();

        DiscreteVariable vars = (DiscreteVariable) mTorGold.getVariable(0);

        for (int i = 0; i < mTorGold.getNumRows(); i++) {
            final String category = vars.getCategory(mTorGold.getInt(i, 0));
            allVars.add(estGraph.getNode(category));
        }

        List<Node> trues = new ArrayList<>();

        for (int i = 0; i < allVars.size(); i++) {
            int count = 0;

            for (int j = 25; j <= 32; j++) {
                if (mTorGold.getInt(i, j) == 1) {
                    count++;

                    if (count == 1) {
                        final String category = vars.getCategory(mTorGold.getInt(i, 0));
                        trues.add(estGraph.getNode(category));
                        break;
                    }
                }
            }
        }

        System.out.println("Trues = " + trues);


        Node mTor = _estGraph.getNode("mTOR_pS2448");
        this.all = allVars;


        List<Node> desc = getDescendants(_estGraph, Collections.singletonList(mTor));

        this.positives = new ArrayList<>();
        this.negatives = new ArrayList<>();

        tp = new ArrayList<>();
        tn = new ArrayList<>();
        fp = new ArrayList<>();
        fn = new ArrayList<>();

        for (Node d : allVars) {
            if (trues.contains(d) && desc.contains(d)) {
                tp.add(d);
                this.positives.add(d);
            } else if (trues.contains(d) && !desc.contains(d)) {
                tn.add(d);
                this.negatives.add(d);
            } else if (!trues.contains(d) && desc.contains(d)) {
                fp.add(d);
                this.positives.add(d);
            } else if (!trues.contains(d) && !desc.contains(d)) {
                fn.add(d);
                this.negatives.add(d);
            }
        }

    }

    public List<Node> getTp() {
        return tp;
    }

    public List<Node> getTn() {
        return tn;
    }

    public List<Node> getFp() {
        return fp;
    }

    public List<Node> getFn() {
        return fn;
    }

    public List<Node> getPositives() {
        return positives;
    }

    public List<Node> getAll() {
        return all;
    }

    public List<Node> getNegatives() {
        return negatives;
    }


    public static DataSet getDiscreteData(String path, Delimiter delimiter) throws IOException {
        Path dataFile = Paths.get(path);
        VerticalDiscreteTabularDatasetReader dataReader = new VerticalDiscreteTabularDatasetFileReader(dataFile, delimiter);
        return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
    }

    public Graph getEstGraph() {
        return estGraph;
    }

    public Graph getTrueGraph() {
        return trueGraph;
    }

    private void collectDescendantsVisit(Graph graph, Node node, Set<Node> descendants) {
        descendants.add(node);
        List<Node> children = getChildren(graph, node);

        if (!children.isEmpty()) {
            for (Object aChildren : children) {
                Node child = (Node) aChildren;
                doChildClosureVisit(graph, child, descendants, 1);
            }
        }
    }

    /**
     * closure under the child relation
     */
    private void doChildClosureVisit(Graph graph, Node node, Set<Node> closure, int depth) {
        if (depth > 3) return;

        if (!closure.contains(node)) {
            closure.add(node);

            for (Edge edge1 : graph.getEdges(node)) {
                Node sub = Edges.traverseDirected(node, edge1);

                if (sub == null) {
                    continue;
                }

                doChildClosureVisit(graph, sub, closure, depth + 1);
            }
        }
    }

    /**
     * @return the list of children for a node.
     */
    public List<Node> getChildren(Graph graph, Node node) {
        List<Node> children = new ArrayList<>();

        for (Object o : graph.getEdges(node)) {
            Edge edge = (Edge) (o);
            Node sub = Edges.traverseDirected(node, edge);

            if (sub != null) {
                children.add(sub);
            }
        }

        return children;
    }

    public List<Node> getDescendants(Graph graph, List<Node> nodes) {
        Set<Node> descendants = new HashSet<>();

        for (Object node1 : nodes) {
            Node node = (Node) node1;
            collectDescendantsVisit(graph, node, descendants);
        }

        return new LinkedList<>(descendants);
    }
}
