package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.dbmi.data.Dataset;

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
public class ClusterPrecision implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "Precision";
    }

    @Override
    public String getDescription() {
        return "Maximum of true edges over all edges";
    }


    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {

        //Precision = tp/(tp+fp)

        if(trueGraph.getNode("L.X2")!= null){

            if(trueGraph.isAdjacentTo(trueGraph.getNode("L.X1"),trueGraph.getNode("L.X2"))){

                trueGraph.removeEdge(trueGraph.getNode("L.X1"),trueGraph.getNode("L.X2"));

            }
        }

        List<List<Node>> trueClusters = GrabCluster(trueGraph);
        List<List<Node>> estClusters = GrabCluster(estGraph);

        double precision = 0;


        for (List<Node> truecluster : trueClusters){

            double sprecision = 0;

            for (List<Node> estcluster : estClusters){

                double nodeCounter = 0;

                for (Node node : truecluster){

                    if (estcluster.contains(node)){
                        nodeCounter++;
                    }


                }

                double trueRatio = nodeCounter/estcluster.size();

                if (sprecision<trueRatio){  sprecision = trueRatio; }

            }

            precision += sprecision;


        }

        return precision/trueClusters.size();
//

    }

    @Override
    public double getNormValue(double value) {
        return value;
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
}
