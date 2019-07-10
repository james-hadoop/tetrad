package edu.cmu.tetrad.algcomparison.algorithm.mixed;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.dbmi.algo.bootstrap.BootstrapEdgeEnsemble;
import edu.pitt.dbmi.algo.bootstrap.GeneralBootstrapTest;

import java.util.*;

/**
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "MGM",
        command = "mgm",
        algoType = AlgType.produce_undirected_graphs
)
public class Mgm implements Algorithm {

    static final long serialVersionUID = 23L;

    public Mgm() {
    }

    @Override
    public Graph search(DataModel ds, Parameters parameters) {
    	// Notify the user that you need at least one continuous and one discrete variable to run MGM
        List<Node> variables = ds.getVariables();
        boolean hasContinuous = false;
        boolean hasDiscrete = false;

        for (Node node : variables) {
            if (node instanceof ContinuousVariable) {
                hasContinuous = true;
            }

            if (node instanceof DiscreteVariable) {
                hasDiscrete = true;
            }
        }

        if (!hasContinuous || !hasDiscrete) {
        //    throw new IllegalArgumentException("You need at least one continuous and one discrete variable to run MGM.");
            if(!hasDiscrete){
                ds = ds.copy();
                DataModel dataSetC = ds.copy();


                Map<Node, DiscretizationSpec> map = new HashMap<>();
                List<String> categories = new ArrayList<>();
                categories.add("0");
                categories.add("1");
                //categories.add("0");
//                if(square){
//                    //System.out.println(dataSet);
//                    dataSet = DataUtils.convertToSquareValue((DataSet)dataSet);
//                    //System.out.println("Square: "+ dataSet);
//                }

                double[] cutoffList = new double[]{};
//                if (randomCutoff) {
//                    cutoffList = new double[]{scale*Math.random(), scale*Math.random(), scale*Math.random(), scale*Math.random()};
//                    System.out.println("scale: "+scale+" "+Arrays.toString(cutoffList));
//                }else{
//                    if(square){
//
//                        cutoffList = new double[]{0.5, 0.8, 1.0, 1.40};
//
//                    }
//                    else{
                        cutoffList = new double[]{0, 0, 0, 0};
//                    }
//                }

                for (int i = 0; i < dataSetC.getVariables().size(); i++) {
                    Node node = dataSetC.getVariables().get(i);
                    map.put(node, new ContinuousDiscretizationSpec(new double[]{cutoffList[i]}, categories));
                }

                Discretizer discretizer = new Discretizer((DataSet) dataSetC, map);
                dataSetC = discretizer.discretize();

//                if(mixed) {
                    //System.out.println(dataSet);
                    int[] column = new int[]{0, 2};
                    ds = DataUtils.ColumnMontage((DataSet) ds, (DataSet) dataSetC, column);
                    //System.out.println("Mixed: "+ ds);
//                }

      //          ds = DataUtils.convertNumericalDiscreteToContinuous((DataSet) ds);
            }
        }
        
        if (parameters.getInt("bootstrapSampleSize") < 1) {
            DataSet _ds = DataUtils.getMixedDataSet(ds);

            double mgmParam1 = parameters.getDouble("mgmParam1");
            double mgmParam2 = parameters.getDouble("mgmParam2");
            double mgmParam3 = parameters.getDouble("mgmParam3");

            double[] lambda = {
                mgmParam1,
                mgmParam2,
                mgmParam3
            };

            MGM m = new MGM(_ds, lambda);

            return m.search();
        } else {
            Mgm algorithm = new Mgm();

            DataSet data = (DataSet) ds;

            GeneralBootstrapTest search = new GeneralBootstrapTest(data, algorithm, parameters.getInt("bootstrapSampleSize"));

            BootstrapEdgeEnsemble edgeEnsemble = BootstrapEdgeEnsemble.Highest;
            switch (parameters.getInt("bootstrapEnsemble", 1)) {
                case 0:
                    edgeEnsemble = BootstrapEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = BootstrapEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = BootstrapEdgeEnsemble.Majority;
            }
            search.setEdgeEnsemble(edgeEnsemble);
            search.setParameters(parameters);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        }
    }

    // Need to marry the parents on this.
    @Override
    public Graph getComparisonGraph(Graph graph) {
        return GraphUtils.undirectedGraph(graph);
    }

    @Override
    public String getDescription() {
        return "Returns the output of the MGM (Mixed Graphical Model) algorithm (a Markov random field)";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("mgmParam1");
        parameters.add("mgmParam2");
        parameters.add("mgmParam3");
        // Bootstrapping
        parameters.add("bootstrapSampleSize");
        parameters.add("bootstrapEnsemble");
        parameters.add("verbose");
        return parameters;
    }
}
