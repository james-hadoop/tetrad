package edu.cmu.tetrad.util;

public interface Matrix {
    Matrix eye(int rows);

    Matrix sqrt();

    int rows();

    int columns();

    Matrix getSelection(int[] rows, int[] cols);

    Matrix copy();

    Vector getColumn(int j);

    Matrix times(Matrix m);

    Vector times(Vector v);

    double[][] toArray();

    double get(int i, int j);

    Matrix like();

    void set(int i, int j, double v);

    Vector getRow(int i);

    Matrix getPart(int i, int j, int k, int l);

    Matrix inverse();

    Matrix symmetricInverse();

    Matrix ginverse();

    void assignRow(int row, Vector doubles);

    void assignColumn(int col, Vector doubles);

    double trace();

    double det();

    Matrix transpose();

    Matrix transposeWithoutCopy();

    boolean zeroDimension();

    boolean equals(Matrix m, double tolerance);

    boolean isSquare();

    boolean isSymmetric(double tolerance);

    double zSum();

    Matrix minus(Matrix mb);

    Matrix plus(Matrix mb);

    Matrix scalarMult(double scalar);

    int rank();

    double norm1();

    Vector diag();

    String toString();

    void assign(Matrix matrix);

    Vector sum(int direction);
}
