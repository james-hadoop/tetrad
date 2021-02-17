package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.SemBicScore;
import edu.cmu.tetrad.search.SemBicScorer;
import edu.cmu.tetrad.search.SearchGraphUtils;

import java.util.HashSet;

import static java.lang.Math.tanh;

/**
 * Difference between the true and estiamted BIC scores.
 *
 * @author jdramsey
 */
public class BicDiff implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "BicDiff";
    }

    @Override
    public String getDescription() {
        return "Difference between the true and estimated BIC scores";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());

        assert estGraph != null;
        if (!new HashSet<>(estGraph.getNodes()).equals(new HashSet<>(trueGraph.getNodes()))) {
            throw new IllegalArgumentException("Variables must be the same between true and estimate to calculate a BIC score.");
        }

        SemBicScore score = new SemBicScore((DataSet) dataModel);
        score.setRuleType(SemBicScore.RuleType.CHICKERING);
        score.setPenaltyDiscount(1);
        Fges fges = new Fges(score);
        double _true = fges.scoreDag(trueGraph);
        double est = fges.scoreDag(SearchGraphUtils.dagFromPattern(estGraph));
        return (est - _true);
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value / 1e6);
    }
}

