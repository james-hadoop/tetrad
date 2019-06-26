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
    private Graph initialGraph = null;
    private Algorithm algorithm = null;
    private IKnowledge knowledge = new Knowledge2();

    public Fofc() {

        this(false,false,1,false);
    }

    public Fofc(boolean disretize, boolean randomCutoff, double scale, boolean square) {

        this.discretize = disretize;
        this.randomCutoff = randomCutoff;
        this.scale = scale;
        this.square = square;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (discretize) {
            dataSet = dataSet.copy();

            Map<Node, DiscretizationSpec> map = new HashMap<>();
            List<String> categories = new ArrayList<>();
            categories.add("0");
            categories.add("1");
            //categories.add("0");
            if(square){
                //System.out.println(dataSet);
                dataSet = DataUtils.convertToSquareValue((DataSet)dataSet);
                //System.out.println("Square: "+ dataSet);
            }
            double[] cutoffList = new double[]{};
            if (randomCutoff) {
                cutoffList = new double[]{scale*Math.random(), scale*Math.random(), scale*Math.random(), scale*Math.random()};
                System.out.println("scale: "+scale+" "+Arrays.toString(cutoffList));
            }else{
                if(square){

                    cutoffList = new double[]{0.5, 0.8, 1.0, 1.40};

                }
                else {cutoffList = new double[]{0, 0, 0, 0};}
            }

            for (int i = 0; i < dataSet.getVariables().size(); i++) {
                Node node = dataSet.getVariables().get(i);
                map.put(node, new ContinuousDiscretizationSpec(new double[]{cutoffList[i]}, categories));
            }

            Discretizer discretizer = new Discretizer((DataSet) dataSet, map);
            dataSet = discretizer.discretize();

            dataSet = DataUtils.convertNumericalDiscreteToContinuous((DataSet) dataSet);

            //System.out.println(dataSet);
        }

    	if (parameters.getInt("bootstrapSampleSize") < 1) {
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
