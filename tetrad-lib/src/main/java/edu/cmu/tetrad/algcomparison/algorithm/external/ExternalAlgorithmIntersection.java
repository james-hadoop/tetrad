package edu.cmu.tetrad.algcomparison.algorithm.external;

import edu.cmu.tetrad.algcomparison.algorithm.ExternalAlgorithm;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.Parameters;

import java.util.Set;

/**
 * An API to allow results from external algorithms to be included in a report through the algrorithm
 * comparison tool. This one is for matrix generated by PC in pcalg. See below. This script can generate
 * the files in R.
 * <p>
 * library("MASS");
 * library("pcalg");
 * <p>
 * path<-"/Users/user/tetrad/comparison-final";
 * simulation<-1;
 * <p>
 * subdir<-"pc.solve.confl.TRUE";
 * dir.create(paste(path, "/save/", simulation, "/", subdir, sep=""));
 * <p>
 * for (i in 1:10) {
 * data<-read.table(paste(path, "/save/", simulation, "/data/data.", i, ".txt", sep=""), header=TRUE)
 * n<-nrow(data)
 * C<-cor(data)
 * v<-names(data)
 * suffStat<-list(C = C, n=n)
 * pc.fit<-pc(suffStat=suffStat, indepTest=gaussCItest, alpha=0.001, labels=v,
 * solve.conf=TRUE)
 * A<-as(pc.fit, "amat")
 * name<-paste(path, "/save/", simulation, "/", subdir, "/graph.", i, ".txt", sep="")
 * print(name)
 * write.matrix(A, file=name, sep="\t")
 * }
 *
 * @author jdramsey
 */
public class ExternalAlgorithmIntersection extends ExternalAlgorithm {
    static final long serialVersionUID = 23L;
    private final ExternalAlgorithm[] algorithms;
    private String shortDescription = null;
    private long elapsed = -99;

    public ExternalAlgorithmIntersection(String shortDescription, ExternalAlgorithm... algorithms) {
        this.algorithms = algorithms;
        this.shortDescription = shortDescription;
    }

    /**
     * Reads in the relevant graph from the file (see above) and returns it.
     */
    public Graph search(DataModel dataSet, Parameters parameters, Graph trueGraph) {
        this.elapsed = 0;

        for (ExternalAlgorithm algorithm : algorithms) {
            algorithm.setPath(this.path);
            algorithm.setSimIndex(this.simIndex);
            algorithm.setSimulation(this.simulation);
            elapsed += algorithm.getElapsedTime((DataSet) dataSet, parameters);
        }

        Graph graph0 = algorithms[0].search(dataSet, parameters, trueGraph);
        Set<Edge> edges = graph0.getEdges();

        for (int i = 1; i < algorithms.length; i++) {
            edges.retainAll(algorithms[i].search(dataSet, parameters, trueGraph).getEdges());
        }

        EdgeListGraph intersection = new EdgeListGraph(graph0.getNodes());

        for (Edge edge : edges) {
            intersection.addEdge(edge);
        }

        return intersection;
    }

    /**
     * Returns the cpdag of the supplied DAG.
     */
    public Graph getComparisonGraph(Graph graph) {
        return algorithms[0].getComparisonGraph(graph);
    }

    public String getDescription() {
        return shortDescription;
    }

    public DataType getDataType() {
        return DataType.Continuous;
    }

    public long getElapsedTime(DataModel dataSet, Parameters parameters) {
        return this.elapsed;
    }

}
