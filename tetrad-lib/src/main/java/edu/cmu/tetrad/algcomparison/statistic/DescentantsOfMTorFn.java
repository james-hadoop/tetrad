package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The arrow precision. This counts arrowheads maniacally, wherever they occur in the graphs.
 * The true positives are the number of arrowheads in both the true and estimated graphs.
 * Thus, if the true contains X*->Y and estimated graph either does not contain an edge from
 * X to Y or else does not contain an arrowhead at X for an edge from X to Y, one false
 * positive is counted. Similarly for false negatives.
 *
 * @author jdramsey
 */
public class DescentantsOfMTorFn implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "mtorDescFn";
    }

    @Override
    public String getDescription() {
        return "Descendants of mTor FN";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        Graph _trueGraph = GraphUtils.replaceNodes(trueGraph, estGraph.getNodes());
        Graph _estGraph = GraphUtils.replaceNodes(estGraph, _trueGraph.getNodes());

        String vars = "4EBP1_pS65,ACC_pS79,Akt_pS473,Akt_pT308,AMPK_pT172,Bad_pS112,c_Met_pY1235," +
                "c_Raf_pS338,Chk1_pS345,Chk2_pT68,EGFR_pY1068,EGFR_pY1173,ER_alpha_pS118,GSK3_alpha_beta_pS21_S9," +
                "HER2_pY1248,JNK_pT183_pT185,MAPK_pT202_Y204,MEK1_pS217_S221,mTOR_pS2448,NF_kB_p65_pS536,p27_pT157," +
                "p27_pT198,p38_pT180_Y182,p70S6K_pT389,p90RSK_pT359_S363,PDK1_pS241,PKC_alpha_pS657,PRAS40_pT246," +
                "Rb_pS807_S811,S6_pS235_S236,Src_pY416,Src_pY527,STAT3_pY705,YAP_pS127,YB_1_PS102";

        String[] tokens = vars.split(",");

        List<Node> allVars = new ArrayList<>();

        for (int i = 0; i < tokens.length; i++) allVars.add(_estGraph.getNode(tokens[i]));

        String[] bt549s = new String[]{"YAP_pS127", "PDK1_pS241", "PKC_alpha_pS657", "GSK3_alpha_beta_pS21_S9",
                "PRAS40_pT246", "4EBP1_pS65", "p70S6K_pT389", "S6_pS235_S236", "mTOR_pS2448", "Akt_pS473", "Akt_pS473",
                "Bad_pS112", "AMPK_pT172", "ACC_pS79", "Src_pY416", "c_Met_pY1235", "ER_alpha_pS118", "Chk1_pS345",
                "Chk2_pT68", "p27_pT198", "MEK1_pS217_S221", "p90RSK_pT359_S363", "YB_1_PS102", "p38_pT180_Y182",
                "JNK_pT183_pT185", "STAT3_pY705", "Src_pY416", "EGFR_pY1068", "HER2_pY1248", "EGFR_pY1173"
        };

        List<Node> truthbt549 = new ArrayList<>();

        for (int i = 0; i < bt549s.length; i++) truthbt549.add(_estGraph.getNode(bt549s[i]));

        Node mTor = _estGraph.getNode("mTOR_pS2448");
        List<Node> est = _estGraph.getDescendants(Collections.singletonList(mTor));

        List<Node> fn = new ArrayList<>();

        for (Node d : truthbt549) {
            if (!est.contains(d)) {
                fn.add(d);
                System.out.println("FN: " + d);
            }
        }

        return fn.size();
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
