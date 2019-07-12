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
public class ClusterCounter implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "Cluster";
    }

    @Override
    public String getDescription() {
        return "Counter Clusters";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
//        System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);

        int CusterNumber = 0;


//        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(trueGraph, estGraph);
//        int adjTp = adjConfusion.getAdjTp();
//        int adjFp = adjConfusion.getAdjFp();
////        int adjFn = adjConfusion.getAdjFn();
////        int adjTn = adjConfusion.getAdjTn();
//        return adjTp / (double) (adjTp + adjFp);
        if(estGraph.getNumNodes() <= trueGraph.getNumNodes()) {

            return estGraph.getNumEdges();


        }else{

            if ( estGraph.getParents(estGraph.getNode("X5")).size()>0 && estGraph.getParents(estGraph.getNode("X1")).size()>0){

//                System.out.println( estGraph.getParents(estGraph.getNode("X1")).size());
//
//                System.out.println( estGraph.getParents(estGraph.getNode("X5")).size());

                int estCluster_one = estGraph.getChildren(estGraph.getParents(estGraph.getNode("X1")).get(0)).size();

                System.out.println( "estCluster_one: "+ estCluster_one);

                int estCluster_two = estGraph.getChildren(estGraph.getNode(estGraph.getParents(estGraph.getNode("X5")).get(0).getName())).size();

                System.out.println( "estCluster_two: "+ estCluster_two);

                int trueCluster_one = trueGraph.getChildren(trueGraph.getNode(trueGraph.getParents(trueGraph.getNode("X1")).get(0).getName())).size();

                System.out.println( "trueCluster_one: "+ trueCluster_one);

                int trueCluster_two = trueGraph.getChildren(trueGraph.getNode(trueGraph.getParents(trueGraph.getNode("X5")).get(0).getName())).size();

                System.out.println( "trueCluster_two: "+ trueCluster_two);

                if(trueGraph.isAdjacentTo(trueGraph.getParents(trueGraph.getNode("X1")).get(0), trueGraph.getParents(trueGraph.getNode("X5")).get(0))){

                    estCluster_one ++;

                }

                if (estCluster_one == trueCluster_one && estCluster_two == trueCluster_two){
                    return 1;
                }else {
                    System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
                    return 0;
                }
            }else{
                System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
                return  0;
            }



//            return estGraph.getChildren(estGraph.getNode("_L1")).size()*estGraph.getChildren(estGraph.getNode("_L2")).size();

        }
//            if (estGraph.getChildren(estGraph.getNode("_L1")).size()==trueGraph.getNumEdges()){
//////////        if (estGraph==trueGraph){
////////            System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);
////            return 1;
////////
////        }else {
////            if(estGraph.getChildren(estGraph.getNode("_L1")).size()>0){
////
////                System.out.println(" est graph cluster Number = " + estGraph.getChildren(estGraph.getNode("_L1")).size());
////            }
////            return 0;
////        }
//                return estGraph.getNumEdges();
//        return estGraph.getChildren(estGraph.getNode("_L1")).size();
//        }


    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
