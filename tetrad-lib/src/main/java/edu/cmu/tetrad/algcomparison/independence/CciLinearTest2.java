package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.ConditionalCorrelationLinear;
import edu.cmu.tetrad.search.IndTestConditionalCorrelationLinear;
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
        name = "Conditional Correlation (CCI) Test, Linear",
        command = "cci-test-linear",
        dataType = DataType.Continuous
)
public class CciLinearTest2 implements IndependenceWrapper {

    static final long serialVersionUID = 23L;
    private Graph initialGraph = null;

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        final IndTestConditionalCorrelationLinear cci = new IndTestConditionalCorrelationLinear(DataUtils.getContinuousDataSet(dataSet),
                parameters.getDouble("alpha"));
        if (parameters.getInt("kernelType") == 1) {
            cci.setKernel(ConditionalCorrelationLinear.Kernel.Gaussian);

        } else if (parameters.getInt("kernelType") == 2) {
            cci.setKernel(ConditionalCorrelationLinear.Kernel.Epinechnikov);
        } else {
            throw new IllegalStateException("Kernel not configured.");
        }

        if (parameters.getInt("basisType") == 1) {
            cci.setBasis(ConditionalCorrelationLinear.Basis.Polynomial);
        } else if (parameters.getInt("basisType") == 2) {
            cci.setBasis(ConditionalCorrelationLinear.Basis.Cosine);
        } else {
            throw new IllegalStateException("Basis not configured.");
        }

        cci.setNumFunctions(parameters.getInt("numBasisFunctions"));
        cci.setKernelMultiplier(parameters.getDouble("kernelMultiplier"));
        cci.setKernelRegressionSampleSize(parameters.getInt("kernelRegressionSampleSize"));
        cci.setNumDependenceSpotChecks(parameters.getInt("numDependenceSpotChecks"));
        cci.setEarlyReturn(true);

        return cci;
    }

    @Override
    public String getDescription() {
        return "CCI Test Linear";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add("alpha");
        params.add("numBasisFunctions");
        params.add("kernelType");
        params.add("kernelMultiplier");
        params.add("basisType");
        return params;
    }
}