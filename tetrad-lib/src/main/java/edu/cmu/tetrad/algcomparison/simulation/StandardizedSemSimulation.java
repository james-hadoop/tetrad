package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.graph.SingleGraph;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.SemGraph;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.sem.StandardizedSemIm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
public class StandardizedSemSimulation implements Simulation {
    static final long serialVersionUID = 23L;
    private final RandomGraph randomGraph;
    private SemPm pm;
    private StandardizedSemIm standardizedIm;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();

    public StandardizedSemSimulation(RandomGraph graph) {
        this.randomGraph = graph;
    }

    public StandardizedSemSimulation(SemPm pm) {
        SemGraph graph = pm.getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.pm = pm;
    }

    public StandardizedSemSimulation(StandardizedSemIm im) {
        this.randomGraph = new SingleGraph(im.getSemPm().getGraph());
        this.standardizedIm = im;
        this.pm = im.getSemPm();
    }

    @Override
    public void createData(Parameters parameters, boolean newModel) {
//        if (!newModel && !dataSets.isEmpty()) return;

        Graph graph = this.randomGraph.createGraph(parameters);

        this.dataSets = new ArrayList<>();
        this.graphs = new ArrayList<>();

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            System.out.println("Simulating dataset #" + (i + 1));

            if (parameters.getBoolean("differentGraphs") && i > 0) {
                graph = this.randomGraph.createGraph(parameters);
            }

            this.graphs.add(graph);

            DataSet dataSet = simulate(graph, parameters);
            dataSet.setName("" + (i + 1));
            this.dataSets.add(dataSet);
        }
    }

    @Override
    public DataModel getDataModel(int index) {
        return this.dataSets.get(index);
    }

    @Override
    public Graph getTrueGraph(int index) {
        return this.graphs.get(index);
    }

    @Override
    public String getDescription() {
        return "Linear, Gaussian SEM simulation using " + this.randomGraph.getDescription();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        if (!(this.randomGraph instanceof SingleGraph)) {
            parameters.addAll(this.randomGraph.getParameters());
        }

        parameters.add(Params.NUM_RUNS);
        parameters.add(Params.DIFFERENT_GRAPHS);
        parameters.add(Params.SAMPLE_SIZE);
        return parameters;
    }

    @Override
    public int getNumDataModels() {
        return this.dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    private DataSet simulate(Graph graph, Parameters parameters) {
        if (this.standardizedIm == null) {
            SemPm pm = this.pm;

            if (pm == null) {
                pm = new SemPm(graph);
            }

            this.standardizedIm = new StandardizedSemIm(new SemIm(pm), parameters);
        }

        return this.standardizedIm.simulateData(parameters.getInt(Params.SAMPLE_SIZE), false);
    }
}
