package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.ConditionalCorrelationIndependence;
import edu.cmu.tetrad.search.IndTestConditionalCorrelation;
import edu.cmu.tetrad.search.IndTestUncorrelatedGaussian;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

// Can't change the name of this yet.

/**
 * Wrapper for Daudin Conditional Independence test.
 *
 * @author jdramsey
 */
@TestOfIndependence(
        name = "Uncorrelated Gaussian Test",
        command = "uncor-gauss-test",
        dataType = DataType.Continuous
)
public class IndTestUncorrGaussian implements IndependenceWrapper {

    static final long serialVersionUID = 23L;
    private Graph initialGraph = null;

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        final IndTestUncorrelatedGaussian cci = new IndTestUncorrelatedGaussian(DataUtils.getContinuousDataSet(dataSet),
                parameters.getDouble("alpha"));
        cci.setKernelRegressionSampleSize(parameters.getInt("kernelRegressionSampleSize"));
        return cci;
    }

    @Override
    public String getDescription() {
        return "UCG Test";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add("alpha");
        params.add("kernelRegressionSampleSize");
        return params;
    }
}