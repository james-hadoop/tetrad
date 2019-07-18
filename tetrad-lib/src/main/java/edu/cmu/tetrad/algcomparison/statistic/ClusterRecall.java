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

        //Recall = tp/(tp+fn)
         if(trueGraph.getNode("L.X2")!= null){

                    if(trueGraph.isAdjacentTo(trueGraph.getNode("L.X1"),trueGraph.getNode("L.X2"))){

                        trueGraph.removeEdge(trueGraph.getNode("L.X1"),trueGraph.getNode("L.X2"));

                    }
                }

                if(trueGraph.getNode("L.X3")!= null){

                    if(trueGraph.isAdjacentTo(trueGraph.getNode("L.X3"),trueGraph.getNode("L.X2"))){

                        trueGraph.removeEdge(trueGraph.getNode("L.X2"),trueGraph.getNode("L.X2"));

                    }
                }

                if(trueGraph.getNode("L.X4")!= null){

                    if(trueGraph.isAdjacentTo(trueGraph.getNode("L.X3"),trueGraph.getNode("L.X4"))){

                        trueGraph.removeEdge(trueGraph.getNode("L.X3"),trueGraph.getNode("L.X4"));

                    }
                }

                if(trueGraph.getNode("L.X5")!= null){

                    if(trueGraph.isAdjacentTo(trueGraph.getNode("L.X4"),trueGraph.getNode("L.X5"))){

                        trueGraph.removeEdge(trueGraph.getNode("L.X4"),trueGraph.getNode("L.X5"));

                    }
                }

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

        for (List<Node> estcluster : estClusters){

            double sRecall = 0;

            for (List<Node> truecluster : trueClusters){

                double nodeCounter = 0;

                for (Node node : estcluster){

                    if (truecluster.contains(node)){
                        nodeCounter++;
                    }


                }

                double trueRatio = nodeCounter/truecluster.size();

                if (sRecall<trueRatio){  sRecall = trueRatio; }

            }

            recall += sRecall;


        }

        return recall/estClusters.size();

    }

    private List<List<Node>> GrabCluster(Graph graph){

        List<List<Node>> Clusters = new ArrayList<List<Node>>();

        for (Node d : graph.getNodes()){

            if (graph.getChildren(d).size() > 0 ){
                Clusters.add(graph.getChildren(d));

            }

        }

        return Clusters;

    }

    private boolean isComplete(Graph graph){

        System.out.println("For recall: "+graph);

        int numEdges = graph.getNumEdges();
        int numNodes = graph.getNumNodes();
        int fullyNumEdges = numNodes*(numNodes-1)/2;

        System.out.println("cluster edges " +numEdges+" nodes "+numNodes+" connectivity "+ fullyNumEdges);

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
