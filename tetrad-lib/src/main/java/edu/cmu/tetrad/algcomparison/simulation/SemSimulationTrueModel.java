package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.graph.SingleGraph;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.SemGraph;
import edu.cmu.tetrad.sem.LinearSemIm;
import edu.cmu.tetrad.sem.LinearSemPm;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
public class SemSimulationTrueModel implements Simulation {

    static final long serialVersionUID = 23L;
    private RandomGraph randomGraph;
    private LinearSemPm pm;
    private LinearSemIm im;
    private List<DataModel> dataSets = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();
    private List<LinearSemIm> ims = new ArrayList<>();

    public SemSimulationTrueModel(RandomGraph graph) {
        this.randomGraph = graph;
    }

    public SemSimulationTrueModel(LinearSemPm pm) {
        SemGraph graph = pm.getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.pm = pm;
    }

    public SemSimulationTrueModel(LinearSemIm im) {
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

            makeModels(graph, parameters);

            Matrix impl = getSemIms().get(0).getImplCovar(getSemIms().get(0).getMeasuredNodes());
            ICovarianceMatrix cov = new CovarianceMatrix(getSemIms().get(0).getMeasuredNodes(), impl,
                    Integer.MAX_VALUE);
            if (parameters.getBoolean(Params.RANDOMIZE_COLUMNS)) {
                cov = DataUtils.reorderColumns(cov);
            }

            cov.setName("" + (i + 1));

            dataSets.add(cov);
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

    private void makeModels(Graph graph, Parameters parameters) {
        boolean saveLatentVars = parameters.getBoolean(Params.SAVE_LATENT_VARS);

        LinearSemIm im = this.im;

        if (im == null) {
            LinearSemPm pm = this.pm;

            if (pm == null) {
                pm = new LinearSemPm(graph);
                im = new LinearSemIm(pm, parameters);
                ims.add(im);
            } else {
                im = new LinearSemIm(pm, parameters);
                ims.add(im);
            }
        } else {
            ims.add(im);
        }
    }

    public List<LinearSemIm> getSemIms() {
        return ims;
    }
}
