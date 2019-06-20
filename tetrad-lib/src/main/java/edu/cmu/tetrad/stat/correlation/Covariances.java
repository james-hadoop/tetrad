package edu.cmu.tetrad.stat.correlation;

import edu.cmu.tetrad.util.TetradMatrix;

public interface Covariances {
    double covariance(int i, int j);

    int size();

    void setCovariance(int i, int j, double v);

    double[][] getMatrix();

    double[][] getSubMatrix(int[] rows, int[] cols);
}
