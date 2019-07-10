package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 *
 *
 * @author jdramsey
 */
public class MySuperDooperStatistic implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "TrueRate";
    }

    @Override
    public String getDescription() {
        return "Super Dooper Statistic";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
//        System.out.println("Still Working!");

//        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(trueGraph, estGraph);
//        int adjTp = adjConfusion.getAdjTp();
//        int adjFp = adjConfusion.getAdjFp();
////        int adjFn = adjConfusion.getAdjFn();
////        int adjTn = adjConfusion.getAdjTn();
//        return adjTp / (double) (adjTp + adjFp);
//        if (estGraph.getNumEdges()==trueGraph.getNumEdges()){
////        if (estGraph==trueGraph){
//            System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
//            return 1;
//
//        }else {
//            return 0;
//        }
        return estGraph.getNumNodes();

    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
