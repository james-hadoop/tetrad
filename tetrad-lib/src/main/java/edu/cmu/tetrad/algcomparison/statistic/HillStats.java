package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.DataColumn;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.*;

import java.io.File;
import java.io.IOException;
import java.util.*;


public class HillStats {
//    public static String vars = "4EBP1_pS65,ACC_pS79,Akt_pS473,Akt_pT308,AMPK_pT172,Bad_pS112,c_Met_pY1235," +
//            "c_Raf_pS338,Chk1_pS345,Chk2_pT68,EGFR_pY1068,EGFR_pY1173,ER_alpha_pS118,GSK3_alpha_beta_pS21_S9," +
//            "HER2_pY1248,JNK_pT183_pT185,MAPK_pT202_Y204,MEK1_pS217_S221,mTOR_pS2448,NF_kB_p65_pS536,p27_pT157," +
//            "p27_pT198,p38_pT180_Y182,p70S6K_pT389,p90RSK_pT359_S363,PDK1_pS241,PKC_alpha_pS657,PRAS40_pT246," +
//            "Rb_pS807_S811,S6_pS235_S236,Src_pY416,Src_pY527,STAT3_pY705,YAP_pS127,YB_1_PS102";

//    public static String[] bt549 = new String[]{"YAP_pS127", "PDK1_pS241", "PKC_alpha_pS657", "GSK3_alpha_beta_pS21_S9",
//            "PRAS40_pT246", "4EBP1_pS65", "p70S6K_pT389", "S6_pS235_S236", "mTOR_pS2448", "Akt_pT308", "Akt_pS473",
//            "Bad_pS112", "AMPK_pT172", "ACC_pS79", "Src_pY416", "c_Met_pY1235", "ER_alpha_pS118", "Chk1_pS345",
//            "Chk2_pT68", "p27_pT198", "MEK1_pS217_S221", "p90RSK_pT359_S363", "YB_1_PS102", "p38_pT180_Y182",
//            "JNK_pT183_pT185", "STAT3_pY705", "Src_pY416", "EGFR_pY1068", "HER2_pY1248", "EGFR_pY1173"
////            , "Rb_pS807_S811", "MAPK_pT202_Y204"
//
//    };
//
//    public static String[] canonical = new String[]{
//            "YAP_pS127", "PDK1_pS241", "PKC_alpha_pS657", "GSK3_alpha_beta_pS21_S9",
//            "PRAS40_pT246", "4EBP1_pS65", "p70S6K_pT389", "S6_pS235_S236", "mTOR_pS2448", "Akt_pT308", "Akt_pS473",
//            "Bad_pS112"
//
//    };

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
                        Delimiter.TAB, '\"');
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

//        String[] tokens = vars.split(",");

        List<Node> allVars = new ArrayList<>();

        DiscreteVariable vars = (DiscreteVariable) mTorGold.getVariable(0);

        for (int i = 0; i < mTorGold.getNumRows(); i++) {
            final String category = vars.getCategory(mTorGold.getInt(i, 0));
            allVars.add(estGraph.getNode(category));
        }

        List<Node> trues = new ArrayList<>();

//        final String[] trueGuys = HillStats.bt549;
//        for (int i = 0; i < trueGuys.length; i++) trues.add(_estGraph.getNode(trueGuys[i]));

//        List<String> trueGuys = new ArrayList<>();

//        for (int i = 0; i < allVars.size(); i++) {
//            for (int j = 1; j <= 8; j++) {
//                if (mTorGold.getInt(i, j) == 1) {
//                    final String category = vars.getCategory(i);
//                    trues.add(estGraph.getNode(category));
//                    break;
//                }
//            }
//        }
//
//        for (int i = 0; i < allVars.size(); i++) {
//            for (int j = 9; j <= 16; j++) {
//                if (mTorGold.getInt(i, j) == 1) {
//                    final String category = vars.getCategory(i);
//                    trues.add(estGraph.getNode(category));
//                }
//            }
//        }
//
//        for (int i = 0; i < allVars.size(); i++) {
//            for (int j = 17; j <= 24; j++) {
//                if (mTorGold.getInt(i, j) == 1) {
//                    final String category = vars.getCategory(i);
//                    trues.add(estGraph.getNode(category));
//                }
//            }
//        }

        for (int i = 0; i < allVars.size(); i++) {
//            if (mTorGold.getInt(i, 27) == 1) {
//                final String category = vars.getCategory(mTorGold.getInt(i, 0));
////                System.out.println("True " + i + ". " + category);
//                trues.add(estGraph.getNode(category));
//            }


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


    public static DataSet getDiscreteData(String path, Delimiter delimiter, char quoteChar) throws IOException {
        TabularColumnFileReader colReader = new TabularColumnFileReader(
                new File(path).toPath(),
                delimiter
        );

        colReader.setQuoteCharacter(quoteChar);

        DataColumn[] cols = colReader.generateColumns(new int[0], true);

        TabularDataReader reader = new TabularDataFileReader(
                new File(path).toPath(),
                delimiter
        );

//        MixedTabularDatasetReader reader = new MixedTabularDatasetFileReader(
//                new File(path).toPath(),
//                delimiter,
//                100
//        );


        Data data = reader.read(cols, true);

        final DataSet dataSet = (DataSet) DataConvertUtils.toDataModel(data);

        return dataSet;
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
