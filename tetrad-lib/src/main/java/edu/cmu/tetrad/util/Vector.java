package edu.cmu.tetrad.util;

public interface Vector {
    double[] toArray();

    int size();

    void set(int j, double v);

    double dotProduct(Vector v2);

    Vector like();

    double get(int i);

    Vector copy();

    Vector viewSelection(int[] selection);

    Vector minus(Vector mb);

    Vector plus(Vector mb);

    Vector scalarMult(double scalar);

    Matrix diag();

    String toString();

    boolean equals(Object o);

    void assign(double value);
}
