package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.SemBicScore;

import java.util.List;

public class FgesBicScorer implements Scorer {
    private final CovarianceMatrix cov;
    private final Fges fges;
    private double score = Double.NaN;

    public FgesBicScorer(CovarianceMatrix cov) {
        this.cov = cov;
        this.fges = new Fges(new SemBicScore(cov));
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
        return cov.getSampleSize();
    }

    @Override
    public List<Node> getMeasuredNodes() {
        return cov.getVariables();
    }

    @Override
    public List<Node> getVariables() {
        return cov.getVariables();
    }
}
