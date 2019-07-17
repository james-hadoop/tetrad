package edu.cmu.tetrad.algcomparison.algorithm.cluster;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FindOneFactorClusters;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.search.TestType;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.bootstrap.BootstrapEdgeEnsemble;
import edu.pitt.dbmi.algo.bootstrap.GeneralBootstrapTest;
import edu.pitt.dbmi.data.Dataset;

import java.util.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.lang.Double;


/**
 * FOFC.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FOFC",
        command = "fofc",
        algoType = AlgType.search_for_structure_over_latents
)
public class Fofc implements Algorithm, TakesInitialGraph, HasKnowledge, ClusterAlgorithm {

    static final long serialVersionUID = 23L;
    private boolean discretize = false;
    private boolean randomCutoff = false;
    private double scale = 1;
    private boolean square = false;
    private boolean mixed = false;
    private Graph initialGraph = null;
    private Algorithm algorithm = null;
    private IKnowledge knowledge = new Knowledge2();

    public Fofc() {

        this(false,false,1,false,false);
    }

    public Fofc(boolean disretize, boolean randomCutoff, double scale, boolean square, boolean mixed) {

        this.discretize = disretize;
        this.randomCutoff = randomCutoff;
        this.scale = scale;
        this.square = square;
        this.mixed = mixed;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (discretize) {
            System.out.println("discretize!");
            int measuredNum = dataSet.getVariables().size();
            dataSet = dataSet.copy();
            DataModel dataSetC = dataSet.copy();


            Map<Node, DiscretizationSpec> map = new HashMap<>();
//            List<String> categories = new ArrayList<>();
//            categories.add("0");
//            categories.add("1");
//            categories.add("2");
//            categories.add("3");
//            categories.add("4");
//            categories.add("5");
//            categories.add("6");
            //categories.add("0");
            if (square) {
                //System.out.println(dataSet);
                dataSet = DataUtils.convertToSquareValue((DataSet) dataSet);
                //System.out.println("Square: "+ dataSet);
            }

//            double[] cutoffList = new double[measuredNum];
            if (scale < 2) { //binary data

                List<String> categories = new ArrayList<>();
                categories.add("0");
                categories.add("1");


                double[] cutoffList = new double[measuredNum];

                if (randomCutoff) {
                    for (int i = 0; i < measuredNum; i++) {
                        cutoffList[i] = scale * Math.random();
                    }
                    ;
                    //                System.out.println("scale: "+scale+" "+Arrays.toString(cutoffList));
                } else {
                    //                if(square){
                    //
                    //                    cutoffList = new double[]{0.5, 0.8, 1.0, 1.40};
                    //
                    //                }
                    //                else{
                    for (int i = 0; i < measuredNum; i++) {
                        cutoffList[i] = 0;
                    }
                    // cutoffList = new double[]{0, 0, 0, 0};
                    //}
                }

                for (int i = 0; i < dataSet.getVariables().size(); i++) {
                    Node node = dataSet.getVariables().get(i);
                    map.put(node, new ContinuousDiscretizationSpec(new double[]{cutoffList[i]}, categories));
                }

                Discretizer discretizer = new Discretizer((DataSet) dataSet, map);
                dataSet = discretizer.discretize();

            } else {

                List<String> categories = new ArrayList<>();

                for (int i =0; i < scale +1;i++){

                    categories.add(Integer.toString(i));

                }

//                System.out.println("scale: "+scale +"categories"+ categories.size());


                double[][] cutoffList = new double[measuredNum][categories.size()-1];

                for (int i = 0; i < measuredNum; i++) {


                    for (int j = 0; j < categories.size() - 1; j++) {

                        if (j == 0) {
                            cutoffList[i][j] = -Math.random();

 //

                        } else {
                            cutoffList[i][j] = cutoffList[i][j - 1] + Math.random();
 //

                        }


                    }

                    for (int w = 0; w < categories.size() - 1; w++) {

                        if (cutoffList[i][w] > 0){

                            //System.out.println("Cutoff[" + w + "]" + cutoffList[i][w]);
                            cutoffList[i][w] = 1*cutoffList[i][w]/cutoffList[i][categories.size() - 2];

//                            System.out.println("Cutoff["+"largest"+"]"+ cutoffList[i][categories.size() - 2]);
//                            System.out.println("Cutoff[" + w + "]" + cutoffList[i][w]);

                        }


                }



//                    System.out.println("category:"+ categories.size());

//                    cutoffList[i] = scale * Math.random();
                }


                for (int i = 0; i < dataSet.getVariables().size(); i++) {
                    Node node = dataSet.getVariables().get(i);
                    map.put(node, new ContinuousDiscretizationSpec(cutoffList[i], categories));
                }

                Discretizer discretizer = new Discretizer((DataSet) dataSet, map);
                dataSet = discretizer.discretize();
//                System.out.println("Tichotomied: "+ dataSet);

            }

            if(mixed) {
//                System.out.println("binary:"+ dataSet);
//                System.out.println("Continuous:"+ dataSetC);
                //int[] column = new int[]{0, 2, 4, 5};
                int[] column = new int[measuredNum/2];
                for(int i =0; i<measuredNum/2;i++){
                    column[i] = i+i;
                }
                dataSet = DataUtils.ColumnMontage((DataSet) dataSetC, (DataSet) dataSet, column);
//
            }

            dataSet = DataUtils.convertNumericalDiscreteToContinuous((DataSet) dataSet);

//            System.out.println(((DataSet) dataSet).getCorrelationMatrix());


        }else{
            System.out.println("Continue Matrix!");

//            System.out.println(((DataSet) dataSet).getCorrelationMatrix());


        }

    	if (parameters.getInt("bootstrapSampleSize") < 1) {

//            System.out.println("Mixed: "+ dataSet);

            ICovarianceMatrix cov = DataUtils.getCovMatrix(dataSet);
            double alpha = parameters.getDouble("alpha");

            boolean wishart = parameters.getBoolean("useWishart", true);
            TestType testType;

            if (wishart) {
                testType = TestType.TETRAD_WISHART;
            } else {
                testType = TestType.TETRAD_DELTA;
            }

            boolean gap = parameters.getBoolean("useGap", true);
            FindOneFactorClusters.Algorithm algorithm;

            if (gap) {
                algorithm = FindOneFactorClusters.Algorithm.GAP;
            } else {
                algorithm = FindOneFactorClusters.Algorithm.SAG;
            }

            edu.cmu.tetrad.search.FindOneFactorClusters search
                    = new edu.cmu.tetrad.search.FindOneFactorClusters(cov, testType, algorithm, alpha);
            search.setVerbose(parameters.getBoolean("verbose"));

            return search.search();
        } else {
            Fofc algorithm = new Fofc();

            //algorithm.setKnowledge(knowledge);
//          if (initialGraph != null) {
//      		algorithm.setInitialGraph(initialGraph);
//  		}



            DataSet data = (DataSet) dataSet;


            GeneralBootstrapTest search = new GeneralBootstrapTest(data, algorithm, parameters.getInt("bootstrapSampleSize"));
            search.setKnowledge(knowledge);

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

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return SearchGraphUtils.patternForDag(new EdgeListGraph(graph));
    }

    @Override
    public String getDescription() {
        return "FOFC (Find One Factor Clusters)";
    }


    public double getScale(){return scale;}

    public boolean ifDiscrete(){
        return discretize;
    }

    @Override
    public DataType getDataType() {
        if (discretize) {return DataType.Discrete;}
        else {return DataType.Continuous;}
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("alpha");
        parameters.add("useWishart");
        parameters.add("useGap");
        parameters.add("verbose");
        // Bootstrapping
        parameters.add("bootstrapSampleSize");
        parameters.add("bootstrapEnsemble");
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

}
