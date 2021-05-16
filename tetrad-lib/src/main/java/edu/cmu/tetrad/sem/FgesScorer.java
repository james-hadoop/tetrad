package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.SemBicScore;

import java.util.List;

public class FgesScorer implements Scorer {
    private final DataSet dataSet;
    private final Fges fges;
    private double score = Double.NaN;

    public FgesScorer(DataSet data) {
        this.dataSet = data;
        this.fges = new Fges(new SemBicScore(dataSet));
    }

    @Override
    public double score(Graph dag) {
        double score = -fges.scoreDag(dag);
        this.score = score;
        return score;
    }

    @Override
    public double getScore() {
        return score;
    }

    @Override
    public DataSet getDataSet() {
        return dataSet;
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
}
