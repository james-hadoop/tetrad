package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class HillStats {
    public static String vars = "4EBP1_pS65,ACC_pS79,Akt_pS473,Akt_pT308,AMPK_pT172,Bad_pS112,c_Met_pY1235," +
            "c_Raf_pS338,Chk1_pS345,Chk2_pT68,EGFR_pY1068,EGFR_pY1173,ER_alpha_pS118,GSK3_alpha_beta_pS21_S9," +
            "HER2_pY1248,JNK_pT183_pT185,MAPK_pT202_Y204,MEK1_pS217_S221,mTOR_pS2448,NF_kB_p65_pS536,p27_pT157," +
            "p27_pT198,p38_pT180_Y182,p70S6K_pT389,p90RSK_pT359_S363,PDK1_pS241,PKC_alpha_pS657,PRAS40_pT246," +
            "Rb_pS807_S811,S6_pS235_S236,Src_pY416,Src_pY527,STAT3_pY705,YAP_pS127,YB_1_PS102";

    public static String[] bt549s = new String[]{"YAP_pS127", "PDK1_pS241", "PKC_alpha_pS657", "GSK3_alpha_beta_pS21_S9",
            "PRAS40_pT246", "4EBP1_pS65", "p70S6K_pT389", "S6_pS235_S236", "mTOR_pS2448", "Akt_pS473", "Akt_pS473",
            "Bad_pS112", "AMPK_pT172", "ACC_pS79", "Src_pY416", "c_Met_pY1235", "ER_alpha_pS118", "Chk1_pS345",
            "Chk2_pT68", "p27_pT198", "MEK1_pS217_S221", "p90RSK_pT359_S363", "YB_1_PS102", "p38_pT180_Y182",
            "JNK_pT183_pT185", "STAT3_pY705", "Src_pY416", "EGFR_pY1068", "HER2_pY1248", "EGFR_pY1173"
//            , "Rb_pS807_S811", "MAPK_pT202_Y204"

    };

    private Graph trueGraph;
    private Graph estGraph;
    private List<Node> tp;
    private List<Node> tn;
    private List<Node> fp;
    private List<Node> fn;
    private List<Node> all;
    private List<Node> positives;
    private List<Node> negatives;

    public HillStats(Graph trueGraph, Graph estGraph) {
        this.trueGraph = trueGraph;
        this.estGraph = estGraph;
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


    public HillStats invoke() {
        Graph _trueGraph = GraphUtils.replaceNodes(trueGraph, estGraph.getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(estGraph, _trueGraph.getNodes());

        String[] tokens = vars.split(",");

        List<Node> allVars = new ArrayList<>();

        for (int i = 0; i < tokens.length; i++) allVars.add(_estGraph.getNode(tokens[i]));

        List<Node> trues = new ArrayList<>();

        for (int i = 0; i < bt549s.length; i++) trues.add(_estGraph.getNode(bt549s[i]));

        Node mTor = _estGraph.getNode("mTOR_pS2448");
        this.all = allVars;
        this.positives = _estGraph.getDescendants(Collections.singletonList(mTor));
        List<Node> temp = new ArrayList<>(getAll());
        temp.removeAll(positives);
        this.negatives = temp;

        tp = new ArrayList<>();
        tn = new ArrayList<>();
        fp = new ArrayList<>();
        fn = new ArrayList<>();

        for (Node d : allVars) {
            if (trues.contains(d) && getPositives().contains(d)) {
                tp.add(d);
            } else if (trues.contains(d) && !getPositives().contains(d)) {
                tn.add(d);
            } else if (!trues.contains(d) && getPositives().contains(d)) {
                fp.add(d);
            } else if (!trues.contains(d) && !getPositives().contains(d)) {
                fn.add(d);
            }
        }

        return this;
    }

}
