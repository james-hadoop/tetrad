package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;

import static java.lang.Math.abs;
import static java.lang.Math.tanh;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class DeltaBic implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "DeltaBic";
    }

    @Override
    public String getDescription() {
        return "True BIC - Estimated BIC";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph g = SearchGraphUtils.dagFromPattern(estGraph);
        double bicEst = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search
                .SemBicScore((DataSet) dataModel)).scoreDag(g);

        double bicTrue = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search
                .SemBicScore((DataSet) dataModel)).scoreDag(trueGraph);

        return bicTrue - bicEst;
    }

    @Override
    public double getNormValue(double value) {
        return 1.0 - tanh(abs(value) / 1000000);
    }
}
