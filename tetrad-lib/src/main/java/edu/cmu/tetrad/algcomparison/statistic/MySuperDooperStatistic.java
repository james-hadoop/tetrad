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
//        System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);


//        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(trueGraph, estGraph);
//        int adjTp = adjConfusion.getAdjTp();
//        int adjFp = adjConfusion.getAdjFp();
////        int adjFn = adjConfusion.getAdjFn();
////        int adjTn = adjConfusion.getAdjTn();
//        return adjTp / (double) (adjTp + adjFp);

//            if (estGraph.getChildren(estGraph.getNode("_L1")).size()==trueGraph.getNumEdges()){
//    //////        if (estGraph==trueGraph){
//    ////            System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
//                return 1;
//    ////
//            }else {
//                if(estGraph.getChildren(estGraph.getNode("_L1")).size()>0){
//
//                    System.out.println(" est graph cluster Number = " + estGraph.getChildren(estGraph.getNode("_L1")).size());
//                }
//                return 0;
//            }
          return estGraph.getNumEdges();
//        return estGraph.getChildren(estGraph.getNode("_L1")).size();

    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
