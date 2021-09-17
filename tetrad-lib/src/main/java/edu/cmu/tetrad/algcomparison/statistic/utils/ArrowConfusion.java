package edu.cmu.tetrad.algcomparison.statistic.utils;

import edu.cmu.tetrad.graph.*;

import java.util.List;

/**
 * A confusion matrix for arrows--i.e. TP, FP, TN, FN for counts of arrow endpoints.
 * A true positive arrow is counted for X*->Y in the estimated graph if X is not adjacent
 * to Y or X--Y or X<--Y.
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
        arrowsTp = 0;
        arrowsTpc = 0;
        arrowsFp = 0;
        arrowsFpc = 0;
        arrowsFn = 0;
        arrowsFnc = 0;
        TCtp = 0; //for the two-cycle accuracy
        TCfn = 0;
        TCfp = 0;
        this.truthAdj = truthAdj;

        est = GraphUtils.replaceNodes(est, truth.getNodes());
        assert est != null;
        truth = GraphUtils.replaceNodes(truth, est.getNodes());
        assert truth != null;

        // Get edges from the true Graph to compute TruePositives, TrueNegatives and FalseNeagtives
        //    System.out.println(this.truth2.getEdges());

        for (Edge edge : truth.getEdges()) {

            List<Edge> edges1 = est.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            //      System.out.println(edge1 + "(est2)");

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(edge.getNode1());
                e2Est = edge1.getProximalEndpoint(edge.getNode2());
            }

            List<Edge> edges2 = truth.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
            } else {
                edge2 = truth.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(edge.getNode1());
                e2True = edge2.getProximalEndpoint(edge.getNode2());
            }

            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsFnc = getArrowsFnc() + 1;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsFnc = getArrowsFnc() + 1;
            }


            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW) {
                arrowsTp++;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW) {
                arrowsTp++;
            }

            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsTpc = getArrowsTpc() + 1;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsTpc = getArrowsTpc() + 1;
            }

            if (e1True != Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                arrowsTn++;
            }

            if (e2True != Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                arrowsTn++;
            }

            if (e1True != Endpoint.ARROW && e1Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsTnc = getArrowsTnc() + 1;
            }

            if (e2True != Endpoint.ARROW && e2Est != Endpoint.ARROW && truth.isAdjacentTo(edge.getNode1(), edge.getNode2())
                    && est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                arrowsTnc = getArrowsTnc() + 1;
            }
        }

        // Get edges from the estimated graph to compute only FalsePositives
        // System.out.println(this.est2.getEdges());

        for (Edge edge : est.getEdges()) {

            List<Edge> edges1 = est.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(edge.getNode1());
                e2Est = edge1.getProximalEndpoint(edge.getNode2());
            }

            List<Edge> edges2 = truth.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
            } else {
                edge2 = truth.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(edge.getNode1());
                e2True = edge2.getProximalEndpoint(edge.getNode2());
            }

            if (isTruthAdj()) {
                if (truth.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                    if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                        arrowsFp++;
                    }

                    if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                        arrowsFp++;
                    }
                }
            } else {
                if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                    arrowsFp++;
                }

                if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                    arrowsFp++;
                }
            }

            if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW && edge2 != null) {
                arrowsFpc = getArrowsFpc() + 1;
            }

            if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW && edge2 != null) {
                arrowsFpc = getArrowsFpc() + 1;
            }

        }

        // test for 2-cycle
        for (Edge edge : truth.getEdges()) {
            List<Edge> TwoCycle1 = truth.getEdges(edge.getNode1(), edge.getNode2());
            List<Edge> TwoCycle2 = est.getEdges(edge.getNode1(), edge.getNode2());

            if (TwoCycle1.size() == 2 && TwoCycle2.size() == 2) {
                TCtp++;
            }

            if (TwoCycle1.size() == 2 && TwoCycle2.size() != 2) {
                TCfn++;
            }
        }

        for (Edge edge : est.getEdges()) {

            List<Edge> TwoCycle1 = truth.getEdges(edge.getNode1(), edge.getNode2());
            List<Edge> TwoCycle2 = est.getEdges(edge.getNode1(), edge.getNode2());

            if (TwoCycle1.size() != 2 && TwoCycle2.size() == 2) {
                TCfp++;
            }
        }

        //divide by 2, the 2cycle accuracy is duplicated due to how getEdges is used
        TCtp = TCtp / 2;
        TCfn = TCfn / 2;
        TCfp = TCfp / 2;
    }


    public int getArrowsTp() {
        return arrowsTp;
    }

    public int getArrowsFp() {
        return arrowsFp;
    }

    public int getArrowsFn() {
        return arrowsFn;
    }

    public int getArrowsTn() {
        return arrowsTn;
    }

    public int getTwoCycleTp() {
        return TCtp;
    }

    public int getTwoCycleFp() {
        return TCfp;
    }

    public int getTwoCycleFn() {
        return TCfn;
    }

    /**
     * Two positives for common edges.
     */
    public int getArrowsTpc() {
        return arrowsTpc;
    }

    /**
     * False positives for common edges.
     */
    public int getArrowsFpc() {
        return arrowsFpc;
    }

    /**
     * False negatives for common edges.
     */
    public int getArrowsFnc() {
        return arrowsFnc;
    }

    /**
     * True Negatives for common edges.
     */
    public int getArrowsTnc() {
        return arrowsTnc;
    }

    public boolean isTruthAdj() {
        return truthAdj;
    }
}
