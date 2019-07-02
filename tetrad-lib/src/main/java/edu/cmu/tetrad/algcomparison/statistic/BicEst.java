package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;

import static java.lang.Math.tanh;

/**
 * The adjacency precision. The true positives are the number of adjacencies in both
 * the true and estimated graphs.
 *
 * @author jdramsey
 */
public class BicEst implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "BicEst";
    }

    @Override
    public String getDescription() {
        return "BIC of the estimated pattern";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph g = SearchGraphUtils.dagFromPattern(estGraph);
        return new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search
                .SemBicScore((DataSet) dataModel, false, false)).scoreDag(g);
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value);
    }
}
