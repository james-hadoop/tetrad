package edu.cmu.tetrad.algcomparison.algorithm;

import edu.cmu.tetrad.algcomparison.utils.TakesExternalGraph;
import edu.cmu.tetrad.data.BootstrapSampler;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Stability selection.
 *
 * @author jdramsey
 */
public class StabilitySelection implements Algorithm, TakesExternalGraph {

    static final long serialVersionUID = 23L;
    private final Algorithm algorithm;
    private Graph externalGraph;

    public StabilitySelection(Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        DataSet _dataSet = (DataSet) dataSet;

        double percentageB = parameters.getDouble("percentSubsampleSize");
        int numSubsamples = parameters.getInt("numSubsamples");

        Map<Edge, Integer> counts = new HashMap<>();

        List<Graph> graphs = new ArrayList<>();

        ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {

            private final int chunk;
            private final int from;
            private final int to;

            private StabilityAction(int chunk, int from, int to) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            @Override
            protected void compute() {
                if (this.to - this.from <= this.chunk) {
                    for (int s = this.from; s < this.to; s++) {
                        BootstrapSampler sampler = new BootstrapSampler();
                        sampler.setWithoutReplacements(true);
                        DataSet sample = sampler.sample(_dataSet, (int) (percentageB * _dataSet.getNumRows()));
                        Graph graph = StabilitySelection.this.algorithm.search(sample, parameters);
                        graphs.add(graph);
                    }
                } else {
                    int mid = (this.to + this.from) / 2;

                    StabilityAction left = new StabilityAction(this.chunk, this.from, mid);
                    StabilityAction right = new StabilityAction(this.chunk, mid, this.to);

                    left.fork();
                    right.compute();
                    left.join();
                }
            }
        }

        final int chunk = 2;

        pool.invoke(new StabilityAction(chunk, 0, numSubsamples));

        for (Graph graph : graphs) {
            for (Edge edge : graph.getEdges()) {
                increment(edge, counts);
            }
        }

        this.externalGraph = new EdgeListGraph(dataSet.getVariables());
        double percentStability = parameters.getDouble("percentStability");

        for (Edge edge : counts.keySet()) {
            if (counts.get(edge) > percentStability * numSubsamples) {
                this.externalGraph.addEdge(edge);
            }
        }

        return this.externalGraph;
    }

    private void increment(Edge edge, Map<Edge, Integer> counts) {
        counts.putIfAbsent(edge, 0);
        counts.put(edge, counts.get(edge) + 1);
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return this.algorithm.getComparisonGraph(graph);
    }

    @Override
    public String getDescription() {
        return "Stability selection for " + this.algorithm.getDescription();
    }

    @Override
    public DataType getDataType() {
        return this.algorithm.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = this.algorithm.getParameters();
        parameters.add("depth");
        parameters.add("verbose");
        parameters.add("numSubsamples");
        parameters.add("percentSubsampleSize");
        parameters.add("percentStability");

        return parameters;
    }

    @Override
    public Graph getExternalGraph() {
        return this.externalGraph;
    }

    @Override
    public void setExternalGraph(Graph externalGraph) {
        this.externalGraph = externalGraph;
    }

    @Override
    public void setExternalGraph(Algorithm algorithm) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
