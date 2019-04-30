package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.DagToPag2;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.List;

import static edu.cmu.tetrad.util.StatUtils.correlation;

/**
 * FCI.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FCI",
        command = "fci",
        algoType = AlgType.allow_latent_common_causes
)
public class Fci implements Algorithm, TakesInitialGraph, HasKnowledge, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public Fci() {
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    public Fci(IndependenceWrapper test) {
        this.test = test;
    }

    public Fci(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
    	if (parameters.getInt("numberResampling") < 1) {
            if (algorithm != null) {
                initialGraph = algorithm.search(dataSet, parameters);
            }

            knowledge = getKnowledge((DataSet) dataSet);

            edu.cmu.tetrad.search.Fci search = new edu.cmu.tetrad.search.Fci(test.getTest(dataSet, parameters));
            search.setDepth(parameters.getInt("depth"));
            search.setKnowledge(knowledge);
            search.setMaxPathLength(parameters.getInt("maxPathLength"));
            search.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed"));
            search.setVerbose(parameters.getBoolean("verbose"));

//            if (initialGraph != null) {
//                search.setInitialGraph(initialGraph);
//            }
            return search.search();
        } else {
            Fci algorithm = new Fci(test);
            //algorithm.setKnowledge(knowledge);
//          if (initialGraph != null) {
//      		algorithm.setInitialGraph(initialGraph);
//  		}

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, algorithm, parameters.getInt("numberResampling"));
            search.setKnowledge(knowledge);
            
            search.setPercentResampleSize(parameters.getDouble("percentResampleSize"));
            search.setResamplingWithReplacement(parameters.getBoolean("resamplingWithReplacement"));
            
            ResamplingEdgeEnsemble edgeEnsemble = ResamplingEdgeEnsemble.Highest;
            switch (parameters.getInt("resamplingEnsemble", 1)) {
                case 0:
                    edgeEnsemble = ResamplingEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = ResamplingEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = ResamplingEdgeEnsemble.Majority;
            }
            search.setEdgeEnsemble(edgeEnsemble);
            search.setAddOriginalDataset(parameters.getBoolean("addOriginalDataset"));
            
            search.setParameters(parameters);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        }
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return new DagToPag2(new EdgeListGraph(graph)).convert();
    }

    public String getDescription() {
        return "FCI (Fast Causal Inference) using " + test.getDescription()
                + (algorithm != null ? " with initial graph from "
                        + algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        parameters.add("depth");
        parameters.add("maxPathLength");
        parameters.add("completeRuleSetUsed");
        // Resampling
        parameters.add("numberResampling");
        parameters.add("percentResampleSize");
        parameters.add("resamplingWithReplacement");
        parameters.add("resamplingEnsemble");
        parameters.add("addOriginalDataset");
        parameters.add("verbose");
        return parameters;
    }

    @Override
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    @Override
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    @Override
    public Graph getInitialGraph() {
        return initialGraph;
    }

    @Override
    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    @Override
    public void setInitialGraph(Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    @Override
    public void setIndependenceWrapper(IndependenceWrapper test) {
        this.test = test;
    }

    @Override
    public IndependenceWrapper getIndependenceWrapper() {
        return test;
    }

    public static IKnowledge getKnowledge(DataSet dataSet) {
        dataSet = DataUtils.standardizeData(dataSet);
        double[][] colData = dataSet.getDoubleData().transpose().toArray();

        List<Node> variables = dataSet.getVariables();

//        Collections.sort(variables, new Comparator<Node>() {
//            @Override
//            public int compare(Node o1, Node o2) {
//                if (o1 == o2) return 0;
//                int i = variables.indexOf(o1);
//                int j = variables.indexOf(o2);
//                final double[] x = colData[i];
//                final double[] y = colData[j];
//                return leftRightMinnesota(x, y) ? +1 : -1;
//            }
//        });

        Graph graph = new EdgeListGraph(variables);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                final double[] x = colData[i];
                final double[] y = colData[j];

                if (leftRightMinnesota(x, y)) {
                    graph.addDirectedEdge(variables.get(i), variables.get(j));
                } else {
                    graph.addDirectedEdge(variables.get(j), variables.get(i));
                }
            }
        }

        IKnowledge knowledge = new Knowledge2();
//
//        for (int i = 0; i < variables.size(); i++) {
//            knowledge.addToTier(i + 1, variables.get(i).getName());
//        }

        int numOfNodes = variables.size();
        for (int i = 0; i < numOfNodes; i++) {
            for (int j = i + 1; j < numOfNodes; j++) {
                Node n1 = variables.get(i);
                Node n2 = variables.get(j);

                if (n1.getName().startsWith("E_") || n2.getName().startsWith("E_")) {
                    continue;
                }

                Edge edge = graph.getEdge(n1, n2);
                if (edge != null && edge.isDirected()) {
                    knowledge.setForbidden(edge.getNode2().getName(), edge.getNode1().getName());
                }
            }
        }

        return knowledge;
    }



    private static double[] correctSkewness(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }

    private static boolean leftRightMinnesota(double[] x, double[] y) {
        x = correctSkewness(x);
        y = correctSkewness(y);

        final double cxyx = cov(x, y, x);
        final double cxyy = cov(x, y, y);
        final double cxxx = cov(x, x, x);
        final double cyyx = cov(y, y, x);
        final double cxxy = cov(x, x, y);
        final double cyyy = cov(y, y, y);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        if (sk_ey < 0) {
            lr *= -1;
        }

        final double a = correlation(x, y);

        if (a < 0) {// && sk_ey > 0) {
            lr *= -1;
        }

        return lr > 0;
    }

    private static double cov(double[] x, double[] y, double[] condition) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }

    private static double[] residuals(double[] _y, double[][] _x) {
        TetradMatrix y = new TetradMatrix(new double[][]{_y}).transpose();
        TetradMatrix x = new TetradMatrix(_x).transpose();

        TetradMatrix xT = x.transpose();
        TetradMatrix xTx = xT.times(x);
        TetradMatrix xTxInv = xTx.inverse();
        TetradMatrix xTy = xT.times(y);
        TetradMatrix b = xTxInv.times(xTy);

        TetradMatrix yHat = x.times(b);
        if (yHat.columns() == 0) yHat = y.copy();

        return y.minus(yHat).getColumn(0).toArray();
    }
}
