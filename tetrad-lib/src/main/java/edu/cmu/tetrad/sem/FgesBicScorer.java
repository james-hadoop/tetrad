package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.SemBicScore;

import java.util.List;

public class FgesBicScorer implements Scorer {
    private DataSet dataSet;
    private Fges fges;
    private double score = Double.NaN;

    public FgesBicScorer(DataSet dataSet, double penaltyDiscount) {

        if (dataSet == null) throw new NullPointerException();

        this.dataSet = dataSet;
        SemBicScore score = new SemBicScore(new CovarianceMatrix(dataSet));
        score.setPenaltyDiscount(penaltyDiscount);
        this.fges = new Fges(score);
    }

    public FgesBicScorer(Score score, DataSet dataSet) {
        this.fges = new Fges(score);
        this.dataSet = dataSet;
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
}
