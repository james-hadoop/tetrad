package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.List;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 *
 *
 * @author jdramsey
 */
public class ClusterRecall implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "Recall";
    }

    @Override
    public String getDescription() {
        return "Maximum of true edges over all true edges";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {

//      Recall = tp/(tp+fn)

        if(estGraph.getNumEdges()>trueGraph.getNumEdges()){//MGMgraph

            List<List<Node>> trueClusters = GrabCluster(trueGraph);

//            List<List<Node>> estClusters = GrabCluster(estGraph);



            double recall = 0;

            for (List<Node> truecluster : trueClusters){


                if (isComplete(estGraph.subgraph(truecluster))){

                    recall ++;
                }



            }

            return recall/trueClusters.size();

        }



        List<List<Node>> trueClusters = GrabCluster(trueGraph);
        List<List<Node>> estClusters = GrabCluster(estGraph);



        double recall = 0;
        double optRecall =1;

        for (List<Node> truecluster : trueClusters){

            double sRecall = 0;

            for (List<Node> estcluster : estClusters){


                double nodeCounter = 0;

                for (Node node : estcluster){

                    if (truecluster.contains(node)){
                        nodeCounter++;
                    }


                }

                double trueRatio = nodeCounter/truecluster.size();

                if (sRecall<trueRatio){  sRecall = trueRatio; }//pick the highest accuracy of recovering

            }

            if (optRecall>sRecall) { optRecall = sRecall; }

            recall += sRecall;//for each estimated cluster, pick the recall as the highest accuracy of recovering a true cluster

        }


        return recall/trueClusters.size();//for all true clusters, how good they are reflected


    }

    private List<List<Node>> GrabCluster(Graph graph){

        List<List<Node>> Clusters = new ArrayList<List<Node>>();

        for (Node d : graph.getNodes()){

            if ( d.getName().contains("L") ){//if it is a latent variable

                List<Node> legalChildren = new ArrayList<Node>();


                for (Node c : graph.getChildren(d)){

                    if(!c.getName().contains("L")){//checking if the child is a measured variable

                         if(graph.getNumEdges(c) ==1 ){//if the child is forms a pure model

                             if(graph.getChildren(c).size() ==0){
                                 legalChildren.add(c);
                             }

                         }
// else{
//                             System.out.println(c.getName());
//                         }

                    }

                }


                if(legalChildren.size()>0) Clusters.add(legalChildren);

            }


        }

        return Clusters;

    }

    private boolean isComplete(Graph graph){

        //System.out.println("For recall: "+graph);

        int numEdges = graph.getNumEdges();
        int numNodes = graph.getNumNodes();
        int fullyNumEdges = numNodes*(numNodes-1)/2;

        //System.out.println("cluster edges " +numEdges+" nodes "+numNodes+" connectivity "+ fullyNumEdges);

        if (numEdges == fullyNumEdges){
            return true;
        }else {
            return false;
        }

    }

    private boolean isMeaningful(Node node, Graph graph){

        if(graph.getAdjacentNodes(node).size()>0){
            return true;
        }else {
            return false;
        }

    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
