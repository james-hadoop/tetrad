package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.AdjacencyConfusion;
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
        return "SDS";
    }

    @Override
    public String getDescription() {
        return "Super Dooper Statistic";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        System.out.println("True graph = " + trueGraph + " est graph = " + estGraph);

//        AdjacencyConfusion adjConfusion = new AdjacencyConfusion(trueGraph, estGraph);
//        int adjTp = adjConfusion.getAdjTp();
//        int adjFp = adjConfusion.getAdjFp();
////        int adjFn = adjConfusion.getAdjFn();
////        int adjTn = adjConfusion.getAdjTn();
//        return adjTp / (double) (adjTp + adjFp);
        if (estGraph.getNumEdges()==trueGraph.getNumEdges()){
            return 1;
        }else {
            return 0;
        }

    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
