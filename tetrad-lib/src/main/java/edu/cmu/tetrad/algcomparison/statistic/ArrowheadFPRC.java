package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.algcomparison.statistic.utils.ArrowConfusion;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Graph;

/**
 * The adjacency true positive rate. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class ArrowheadFPRC implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "AHFPRC";
    }

    @Override
    public String getDescription() {
        return "Arrowhead False Positive Rate Common Edges";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        ArrowConfusion adjConfusion = new ArrowConfusion(trueGraph, estGraph);
        int adjFp = adjConfusion.getArrowsFpc();
        int adjTn = adjConfusion.getArrowsTnc();
        return adjFp / (double) (adjFp + adjTn);
    }

    @Override
    public double getNormValue(double value) {
        return value;
    }
}
