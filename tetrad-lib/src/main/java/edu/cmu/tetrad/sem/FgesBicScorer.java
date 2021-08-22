package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.Score;

import java.util.List;

public class FgesBicScorer implements Scorer {
    private final DataSet dataSet;
    private final Fges fges;
    private Graph dag = null;

    public FgesBicScorer(Score score, DataSet dataSet) {
        this.fges = new Fges(score);
        this.dataSet = dataSet;
    }

    @Override
    public double score(Graph dag) {
        this.dag = dag;
        return -fges.scoreDag(dag);
    }

    @Override
    public double score(Edge edge) {
        return -fges.scoreDag(dag);
    }

    @Override
    public DataSet getDataSet() {
        return null;
    }

    @Override
    public int getSampleSize() {
        return dataSet.getNumRows();
    }

    @Override
    public List<Node> getMeasuredNodes() {
        return dataSet.getVariables();
    }

    @Override
    public List<Node> getVariables() {
        return dataSet.getVariables();
    }

    @Override
    public DataType getDataType() {
        if (dataSet.isDiscrete()) return DataType.Discrete;
        else return DataType.Continuous;
    }

    @Override
    public void resetParameters(Edge edge) {
    }
}
