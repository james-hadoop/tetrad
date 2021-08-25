package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.graph.SingleGraph;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.SemGraph;
import edu.cmu.tetrad.sem.LinearSemIm;
import edu.cmu.tetrad.sem.LinearSemPm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
public class LinearSemSimulation implements Simulation {

    static final long serialVersionUID = 23L;
    private RandomGraph randomGraph;
    private LinearSemPm pm;
    private LinearSemIm im;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();
    private List<LinearSemIm> ims = new ArrayList<>();

    public LinearSemSimulation(RandomGraph graph) {
        this.randomGraph = graph;
    }

    public LinearSemSimulation(LinearSemPm pm) {
        SemGraph graph = pm.getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.pm = pm;
    }

    public LinearSemSimulation(LinearSemIm im) {
        SemGraph graph = im.getSemPm().getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.im = im;
        this.pm = im.getSemPm();
        this.ims = new ArrayList<>();
        ims.add(im);
    }

    @Override
    public void createData(Parameters parameters, boolean newModel) {
        if (!newModel && !dataSets.isEmpty()) return;

        Graph graph = randomGraph.createGraph(parameters);

        dataSets = new ArrayList<>();
        graphs = new ArrayList<>();
        ims = new ArrayList<>();

        for (int i = 0; i < parameters.getInt(Params.NUM_RUNS); i++) {
            System.out.println("Simulating dataset #" + (i + 1));

            if (parameters.getBoolean(Params.DIFFERENT_GRAPHS) && i > 0) {
                graph = randomGraph.createGraph(parameters);
            }

            graphs.add(graph);

            DataSet dataSet = simulate(graph, parameters);

            if (parameters.getBoolean(Params.STANDARDIZE)) {
                dataSet = DataUtils.standardizeData(dataSet);
            }

            double variance = parameters.getDouble(Params.MEASUREMENT_VARIANCE);

            if (variance > 0) {
                for (int k = 0; k < dataSet.getNumRows(); k++) {
                    for (int j = 0; j < dataSet.getNumColumns(); j++) {
                        double d = dataSet.getDouble(k, j);
                        double norm = RandomUtil.getInstance().nextNormal(0, Math.sqrt(variance));
                        dataSet.setDouble(k, j, d + norm);
                    }
                }
            }

            if (parameters.getBoolean(Params.RANDOMIZE_COLUMNS)) {
                dataSet = DataUtils.reorderColumns(dataSet);
            }

            dataSet.setName("" + (i + 1));
            dataSets.add(dataSet);
        }
    }

    @Override
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }

    @Override
    public Graph getTrueGraph(int index) {
        return graphs.get(index);
    }

    @Override
    public String getDescription() {
        return "Linear, Gaussian SEM simulation using " + randomGraph.getDescription();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        if (!(randomGraph instanceof SingleGraph)) {
            parameters.addAll(randomGraph.getParameters());
        }

//        if (pm == null) {
//            parameters.addAll(SemPm.getParameterNames());
//        }
        if (im == null) {
            parameters.addAll(LinearSemIm.getParameterNames());
        }

        parameters.add(Params.MEASUREMENT_VARIANCE);
        parameters.add(Params.NUM_RUNS);
        parameters.add(Params.DIFFERENT_GRAPHS);
        parameters.add(Params.RANDOMIZE_COLUMNS);
        parameters.add(Params.SAMPLE_SIZE);
        parameters.add(Params.SAVE_LATENT_VARS);
        parameters.add(Params.STANDARDIZE);

        return parameters;
    }

    @Override
    public int getNumDataModels() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    private DataSet simulate(Graph graph, Parameters parameters) {
        boolean saveLatentVars = parameters.getBoolean(Params.SAVE_LATENT_VARS);

        LinearSemIm im = this.im;

        if (im == null) {
            LinearSemPm pm = this.pm;

            if (pm == null) {
                pm = new LinearSemPm(graph);
                im = new LinearSemIm(pm, parameters);
                ims.add(im);
                return im.simulateData(parameters.getInt(Params.SAMPLE_SIZE), saveLatentVars);
            } else {
                im = new LinearSemIm(pm, parameters);
                ims.add(im);
                return im.simulateData(parameters.getInt(Params.SAMPLE_SIZE), saveLatentVars);
            }
        } else {
            ims.add(im);
            return im.simulateData(parameters.getInt(Params.SAMPLE_SIZE), saveLatentVars);
        }
    }

    public List<LinearSemIm> getSemIms() {
        return ims;
    }
}
