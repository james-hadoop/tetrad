///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines,    //
// Joseph Ramsey, and Clark Glymour.                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.Vector;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static edu.cmu.tetrad.data.Discretizer.discretize;
import static edu.cmu.tetrad.data.Discretizer.getEqualFrequencyBreakPoints;


/**
 * Calculates Mixed Variables Polynomial likelihood.
 *
 * @author Bryan Andrews
 */
public class MVPLikelihood {

    private final DataSet dataSet;

    // The variables of the dataset.
    private final List<Node> variables;

    // The variables of the discrete dataset.
    private List<Node> discreteVariables;

    // Indices of variables.
    private final Map<Node, Integer> nodesHash;

    // Continuous data only.
    private final double[][] continuousData;

    // Discrete data only.
    private final int[][] discreteData;

    // Partitions
    private final AdLeafTree adTree;

    // Fix degree
    private final int fDegree;

    // Structure Prior
    private final double structurePrior;

    // Discretize
    private final boolean discretize;

    public MVPLikelihood(DataSet dataSet, double structurePrior, int fDegree, boolean discretize) {

        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.dataSet = dataSet;
        this.variables = dataSet.getVariables();
        this.structurePrior = structurePrior;
        this.fDegree = fDegree;
        this.discretize = discretize;

        this.continuousData = new double[dataSet.getNumColumns()][];
        this.discreteData = new int[dataSet.getNumColumns()][];
        for (int j = 0; j < dataSet.getNumColumns(); j++) {
            Node v = dataSet.getVariable(j);
            if (v instanceof ContinuousVariable) {
                double[] col = new double[dataSet.getNumRows()];
                for (int i = 0; i < dataSet.getNumRows(); i++) {
                    col[i] = dataSet.getDouble(i, j);
                }
                this.continuousData[j] = col;
            } else if (v instanceof DiscreteVariable) {
                int[] col = new int[dataSet.getNumRows()];
                for (int i = 0; i < dataSet.getNumRows(); i++) {
                    col[i] = dataSet.getInt(i, j);
                }
                this.discreteData[j] = col;
            }
        }

        this.nodesHash = new HashMap<>();
        for (int j = 0; j < dataSet.getNumColumns(); j++) {
            Node v = dataSet.getVariable(j);
            this.nodesHash.put(v, j);
        }

        if (discretize) {
            DataSet discreteDataSet = useErsatzVariables();
            this.discreteVariables = discreteDataSet.getVariables();
            this.adTree = new AdLeafTree(discreteDataSet);
            // All discrete data
            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                int[] col = new int[discreteDataSet.getNumRows()];
                for (int i = 0; i < discreteDataSet.getNumRows(); i++) {
                    col[i] = discreteDataSet.getInt(i, j);
                }
                this.discreteData[j] = col;
            }
        } else {
            this.adTree = new AdLeafTree(dataSet);
        }

    }

    private double multipleRegression(Vector Y, Matrix X) {

        int n = X.rows();
        Vector r;
        if (X.columns() >= n) {
            Vector ones = new Vector(n);
            for (int i = 0; i < n; i++) ones.set(i, 1);
            r = ones.scalarMult(ones.dotProduct(Y) / (double) n).minus(Y);
        } else {
            try {
                Matrix Xt = X.transpose();
                Matrix XtX = Xt.times(X);
                r = X.times(XtX.inverse().times(Xt.times(Y))).minus(Y);
            } catch (Exception e) {
                Vector ones = new Vector(n);
                for (int i = 0; i < n; i++) ones.set(i, 1);
                r = ones.scalarMult(ones.dotProduct(Y) / (double) n).minus(Y);
            }
        }

        double sigma2 = r.dotProduct(r) / n;
        double lik;

        if (sigma2 < 0) {
            Vector ones = new Vector(n);
            for (int i = 0; i < n; i++) ones.set(i, 1);
            r = ones.scalarMult(ones.dotProduct(Y) / (double) Math.max(n, 2)).minus(Y);
            sigma2 = r.dotProduct(r) / n;
            lik = -(n / 2.) * (Math.log(2 * Math.PI) + Math.log(sigma2) + 1);
        } else if (sigma2 == 0) {
            lik = 0;
        } else {
            lik = -(n / 2.) * (Math.log(2 * Math.PI) + Math.log(sigma2) + 1);
        }


        if (Double.isInfinite(lik) || Double.isNaN(lik)) {
            System.out.println(lik);
        }

        return lik;
    }

    private double approxMultinomialRegression(Matrix Y, Matrix X) {

        int n = X.rows();
        int d = Y.columns();
        double lik = 0.0;
        Matrix P;


        if (d >= n || X.columns() >= n) {
            Matrix ones = new Matrix(n, 1);
            for (int i = 0; i < n; i++) ones.set(i, 0, 1);
            P = ones.times(ones.transpose().times(Y).scalarMult(1 / (double) n));
        } else {
            try {
                Matrix Xt = X.transpose();
                Matrix XtX = Xt.times(X);
                P = X.times(XtX.inverse().times(Xt.times(Y)));
            } catch (Exception e) {
                Matrix ones = new Matrix(n, 1);
                for (int i = 0; i < n; i++) ones.set(i, 0, 1);
                P = ones.times(ones.transpose().times(Y).scalarMult(1 / (double) n));
            }

            for (int i = 0; i < n; i++) {
                double min = 1;
                double center = 1 / (double) d;
                double bound = 1 / (double) n;
                for (int j = 0; j < d; j++) {
                    min = Math.min(min, P.get(i, j));
                }
                if (X.columns() > 1 && min < bound) {
                    min = (bound - center) / (min - center);
                    for (int j = 0; j < d; j++) {
                        P.set(i, j, min * P.get(i, j) + center * (1 - min));
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            lik += Math.log(P.getRow(i).dotProduct(Y.getRow(i)));
        }

        if (Double.isInfinite(lik) || Double.isNaN(lik)) {
            System.out.println(lik);
        }

        return lik;
    }


    public double getLik(int child_index, int[] parents) {

        double lik = 0;
        Node c = this.variables.get(child_index);
        List<ContinuousVariable> continuous_parents = new ArrayList<>();
        List<DiscreteVariable> discrete_parents = new ArrayList<>();

        if (c instanceof DiscreteVariable && this.discretize) {
            for (int p : parents) {
                Node parent = this.discreteVariables.get(p);
                discrete_parents.add((DiscreteVariable) parent);
            }
        } else {
            for (int p : parents) {
                Node parent = this.variables.get(p);
                if (parent instanceof ContinuousVariable) {
                    continuous_parents.add((ContinuousVariable) parent);
                } else {
                    discrete_parents.add((DiscreteVariable) parent);
                }
            }
        }

        int p = continuous_parents.size();

        List<List<Integer>> cells = this.adTree.getCellLeaves(discrete_parents);

        int[] continuousCols = new int[p];
        for (int j = 0; j < p; j++) continuousCols[j] = this.nodesHash.get(continuous_parents.get(j));

        for (List<Integer> cell : cells) {
//            for (int[] cell : cells) {
            int r = cell.size();
//                int r = cell.length;
            if (r > 1) {

                double[] mean = new double[p];
                double[] var = new double[p];
                for (int i = 0; i < p; i++) {
                    for (Integer integer : cell) {
                        mean[i] += this.continuousData[continuousCols[i]][integer];
                        var[i] += Math.pow(this.continuousData[continuousCols[i]][integer], 2);
                    }
                    mean[i] /= r;
                    var[i] /= r;
                    var[i] -= Math.pow(mean[i], 2);
                    var[i] = Math.sqrt(var[i]);

                    if (Double.isNaN(var[i])) {
                        System.out.println(var[i]);
                    }
                }

                int degree = this.fDegree;
                if (this.fDegree < 1) {
                    degree = (int) Math.floor(Math.log(r));
                }
                Matrix subset = new Matrix(r, p * degree + 1);
                for (int i = 0; i < r; i++) {
                    subset.set(i, p * degree, 1);
                    for (int j = 0; j < p; j++) {
                        for (int d = 0; d < degree; d++) {
                            subset.set(i, p * d + j, Math.pow((this.continuousData[continuousCols[j]][cell.get(i)] - mean[j]) / var[j], d + 1));
                        }
                    }
                }

                if (c instanceof ContinuousVariable) {
                    Vector target = new Vector(r);
                    for (int i = 0; i < r; i++) {
                        target.set(i, this.continuousData[child_index][cell.get(i)]);
//                        target.set(i, continuousData[child_index][cell[i]]);
                    }
                    lik += multipleRegression(target, subset);
                } else {
                    assert c instanceof DiscreteVariable;
                    Matrix target = new Matrix(r, ((DiscreteVariable) c).getNumCategories());
                    for (int i = 0; i < r; i++) {
                        target.set(i, this.discreteData[child_index][cell.get(i)], 1);
                    }
                    lik += approxMultinomialRegression(target, subset);
                }
            }
        }

        return lik;
    }

    public double getDoF(int child_index, int[] parents) {

        double dof = 0;
        Node c = this.variables.get(child_index);
        List<ContinuousVariable> continuous_parents = new ArrayList<>();
        List<DiscreteVariable> discrete_parents = new ArrayList<>();

        if (c instanceof DiscreteVariable && this.discretize) {
            for (int p : parents) {
                Node parent = this.discreteVariables.get(p);
                discrete_parents.add((DiscreteVariable) parent);
            }
        } else {
            for (int p : parents) {
                Node parent = this.variables.get(p);
                if (parent instanceof ContinuousVariable) {
                    continuous_parents.add((ContinuousVariable) parent);
                } else {
                    discrete_parents.add((DiscreteVariable) parent);
                }
            }
        }

        List<List<Integer>> cells = this.adTree.getCellLeaves(discrete_parents);

        for (List<Integer> cell : cells) {
            int r = cell.size();
            if (r > 0) {

                int degree = this.fDegree;
                if (this.fDegree < 1) {
                    degree = (int) Math.floor(Math.log(r));
                }
                if (c instanceof ContinuousVariable) {
                    dof += degree * continuous_parents.size() + 1;
                } else {
                    assert c instanceof DiscreteVariable;
                    dof += ((degree * continuous_parents.size()) + 1) * (((DiscreteVariable) c).getNumCategories() - 1);
                }
            }
        }

        return dof;
    }

    public double getStructurePrior(int k) {

        if (this.structurePrior < 0) {
            return getEBICprior();
        }

        double n = this.dataSet.getNumColumns() - 1;
        double p = this.structurePrior / n;

        if (this.structurePrior == 0) {
            return 0;
        }
        return k * Math.log(p) + (n - k) * Math.log(1 - p);

    }

    public double getEBICprior() {

        double n = this.dataSet.getNumColumns();
        double gamma = -this.structurePrior;
        return gamma * Math.log(n);

    }

    private DataSet useErsatzVariables() {
        List<Node> nodes = new ArrayList<>();
        // Number of categories to use to discretize continuous mixedVariables.
        int numCategories = 3;

        for (Node x : this.variables) {
            if (x instanceof ContinuousVariable) {
                nodes.add(new DiscreteVariable(x.getName(), numCategories));
            } else {
                nodes.add(x);
            }
        }

        DataSet replaced = new BoxDataSet(new VerticalIntDataBox(this.dataSet.getNumRows(), this.dataSet.getNumColumns()), nodes);

        for (int j = 0; j < this.variables.size(); j++) {
            if (this.variables.get(j) instanceof DiscreteVariable) {
                for (int i = 0; i < this.dataSet.getNumRows(); i++) {
                    replaced.setInt(i, j, this.dataSet.getInt(i, j));
                }
            } else {
                double[] column = this.continuousData[j];

                double[] breakpoints = getEqualFrequencyBreakPoints(column, numCategories);

                List<String> categoryNames = new ArrayList<>();

                for (int i = 0; i < numCategories; i++) {
                    categoryNames.add("" + i);
                }

                Discretizer.Discretization d = discretize(column, breakpoints, this.variables.get(j).getName(), categoryNames);

                for (int i = 0; i < this.dataSet.getNumRows(); i++) {
                    replaced.setInt(i, j, d.getData()[i]);
                }
            }
        }

        return replaced;
    }

}
