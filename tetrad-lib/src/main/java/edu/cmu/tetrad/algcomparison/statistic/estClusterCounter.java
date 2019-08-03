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
public class estClusterCounter implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "eCC";
    }

    @Override
    public String getDescription() {
        return "Estimated Counter Clusters";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {


        double estClusters = GrabCluster(estGraph);

        return estClusters;

    }

    private double GrabCluster(Graph graph){

        double Clusters = 0;

        for (Node d : graph.getNodes()){

            if (graph.getChildren(d).size()>0){

                Clusters ++;

            }

        }

        return Clusters;


    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
