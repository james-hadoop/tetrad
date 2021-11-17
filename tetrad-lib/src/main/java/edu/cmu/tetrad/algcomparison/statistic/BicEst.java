package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;

import static java.lang.Math.tanh;

/**
 * Estimated BIC score.
 *
 * @author jdramsey
 */
public class BicEst implements Statistic {
    static final long serialVersionUID = 23L;
    private Score score;
    private IndependenceTest test;

    public BicEst() {

    }

    public BicEst(Score score) {
        this.score = score;
    }

    public BicEst(IndependenceTest test) {
        this.test = test;
    }

    @Override
    public String getAbbreviation() {
        return "BicEst";
    }

    @Override
    public String getDescription() {
        return "BIC of the estimated cpdag";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph, DataModel dataModel) {
        Graph dag = SearchGraphUtils.dagFromCpdag(estGraph);

        {
            DataModel data = null;

            if (dataModel != null && test != null) {
                throw new IllegalArgumentException("Should only specify one of test, score, or data model.");
            }

            if (dataModel != null) {
                data = dataModel;
            } else if (test != null) {
                data = test.getData();
            }

            if (data != null) {
                Score score;

                if (data instanceof CovarianceMatrix) {
                    return new Fges(new SemBicScore((CovarianceMatrix) data)).scoreDag(dag);
                } else if (data instanceof DataSet) {
                    DataSet dataSet = (DataSet) data;

                    if (dataSet.isContinuous()) {
                        score = new SemBicScore(dataSet);
                    } else if (dataSet.isDiscrete()) {
                        score = new BicScore(dataSet);
                    } else {
                        score = new ConditionalGaussianScore(dataSet, 1, 0, true);
                    }

                    return new Fges(score).scoreDag(dag);
                }
            }
        }

        return Double.NaN;
    }

    @Override
    public double getNormValue(double value) {
        return tanh(value / 1e6);
    }
}

