package edu.cmu.tetrad.algcomparison.statistic.utils;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

import java.util.List;

/**
 * A confusion matrix for arrows--i.e. TP, FP, TN, FN for counts of arrow endpoints.
 * A true positive arrow is counted for X*-&gt;Y in the estimated graph if X is not adjacent
 * to Y or X--Y or X&lt;--Y.
 *
 * @author jdramsey, rubens (November, 2016)
 */
public class ArrowConfusion {

    // For arrowhead FP's, don't count an error unless the variables are adj in the true graph.
    private final boolean truthAdj;

    private int arrowsTp;
    private int arrowsTpc;
    private int arrowsFp;
    private int arrowsFpc;
    private int arrowsFn;
    private int arrowsFnc;
    private int arrowsTn;
    private int arrowsTnc;
    private int TCtp;
    private int TCfn;
    private int TCfp;

    public ArrowConfusion(Graph truth, Graph est) {
        this(truth, est, false);
    }

    public ArrowConfusion(Graph truth, Graph est, boolean truthAdj) {
        Graph truth1 = truth;
        Graph est1 = est;
        this.arrowsTp = 0;
        this.arrowsTpc = 0;
        this.arrowsFp = 0;
        this.arrowsFpc = 0;
        this.arrowsFn = 0;
        this.arrowsFnc = 0;
        this.TCtp = 0; //for the two-cycle accuracy
        this.TCfn = 0;
        this.TCfp = 0;
        this.truthAdj = truthAdj;


        est1 = GraphUtils.replaceNodes(est, truth.getNodes());
        truth1 = GraphUtils.replaceNodes(truth, est.getNodes());


        // Get edges from the true Graph to compute TruePositives, TrueNegatives and FalseNeagtives
        //    System.out.println(this.truth.getEdges());

        for (Edge edge : truth1.getEdges()) {

            List<Edge> edges1 = est1.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est1.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            //      System.out.println(edge1 + "(est)");

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(edge.getNode1());
                e2Est = edge1.getProximalEndpoint(edge.getNode2());
            }

            List<Edge> edges2 = truth1.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
//                if (Edges.isUndirectedEdge(edge2)) continue;
            } else {
                edge2 = truth1.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            //       System.out.println(edge2 + "(truth)");

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(edge.getNode1());
                e2True = edge2.getProximalEndpoint(edge.getNode2());
            }


            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                this.arrowsFn++;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                this.arrowsFn++;
            }

            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsFnc = getArrowsFnc() + 1;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsFnc = getArrowsFnc() + 1;
            }


            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW) {
                this.arrowsTp++;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW) {
                this.arrowsTp++;
            }

            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsTpc = getArrowsTpc() + 1;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsTpc = getArrowsTpc() + 1;
            }

            if (e1True != Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                this.arrowsTn++;
            }

            if (e2True != Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                this.arrowsTn++;
            }

            if (e1True != Endpoint.ARROW && e1Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsTnc = getArrowsTnc() + 1;
            }

            if (e2True != Endpoint.ARROW && e2Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2()) && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                this.arrowsTnc = getArrowsTnc() + 1;
            }
        }
// Get edges from the estimated graph to compute only FalsePositives
        // System.out.println(this.est.getEdges());

        for (Edge edge : est1.getEdges()) {

            List<Edge> edges1 = est1.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est1.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }
            //      System.out.println(edge1 + "(est)");

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(edge.getNode1());
                e2Est = edge1.getProximalEndpoint(edge.getNode2());
            }


            List<Edge> edges2 = truth1.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
//                if (Edges.isUndirectedEdge(edge2)) continue;
            } else {
                edge2 = truth1.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            //          System.out.println(edge2 + "(truth)");

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(edge.getNode1());
                e2True = edge2.getProximalEndpoint(edge.getNode2());
            }

            if (isTruthAdj()) {
                if (truth.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                    if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                        this.arrowsFp++;
                    }

                    if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                        this.arrowsFp++;
                    }
                }
            } else {
                if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                    this.arrowsFp++;
                }

                if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                    this.arrowsFp++;
                }
            }

            if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW && edge1 != null && edge2 != null) {
                this.arrowsFpc = getArrowsFpc() + 1;
            }

            if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW && edge1 != null && edge2 != null) {
                this.arrowsFpc = getArrowsFpc() + 1;
            }

        }


        // test for 2-cycle

        for (Edge edge : truth1.getEdges()) {


            List<Edge> TwoCycle1 = truth1.getEdges(edge.getNode1(), edge.getNode2());
            List<Edge> TwoCycle2 = est1.getEdges(edge.getNode1(), edge.getNode2());

            if (TwoCycle1.size() == 2 && TwoCycle2.size() == 2) {
                //              System.out.println("2-cycle correctly inferred " + TwoCycle1);
                this.TCtp++;
            }

            if (TwoCycle1.size() == 2 && TwoCycle2.size() != 2) {
                //             System.out.println("2-cycle not inferred " + TwoCycle1);
                this.TCfn++;
            }
        }

        for (Edge edge : est1.getEdges()) {

            List<Edge> TwoCycle1 = truth1.getEdges(edge.getNode1(), edge.getNode2());
            List<Edge> TwoCycle2 = est1.getEdges(edge.getNode1(), edge.getNode2());

            if (TwoCycle1.size() != 2 && TwoCycle2.size() == 2) {
                //              System.out.println("2-cycle falsely inferred" + TwoCycle2);
                this.TCfp++;
            }
        }

        //divide by 2, the 2cycle accuracy is duplicated due to how getEdges is used
        this.TCtp = this.TCtp / 2;
        this.TCfn = this.TCfn / 2;
        this.TCfp = this.TCfp / 2;

    }


    public int getArrowsTp() {
        return this.arrowsTp;
    }

    public int getArrowsFp() {
        return this.arrowsFp;
    }

    public int getArrowsFn() {
        return this.arrowsFn;
    }

    public int getArrowsTn() {
        return this.arrowsTn;
    }

    public int getTwoCycleTp() {
        return this.TCtp;
    }

    public int getTwoCycleFp() {
        return this.TCfp;
    }

    public int getTwoCycleFn() {
        return this.TCfn;
    }

    /**
     * Two positives for common edges.
     */
    public int getArrowsTpc() {
        return this.arrowsTpc;
    }

    /**
     * False positives for common edges.
     */
    public int getArrowsFpc() {
        return this.arrowsFpc;
    }

    /**
     * False negatives for common edges.
     */
    public int getArrowsFnc() {
        return this.arrowsFnc;
    }

    /**
     * True Negatives for common edges.
     */
    public int getArrowsTnc() {
        return this.arrowsTnc;
    }

    public boolean isTruthAdj() {
        return this.truthAdj;
    }
}
