///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015, 2022 by Peter Spirtes, Richard        //
// Scheines, Joseph Ramsey, and Clark Glymour.                               //
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

package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.LogisticRegression;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.search.IndependenceResult;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.SearchLogUtils;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.ProbUtils;
import edu.cmu.tetrad.util.TetradLogger;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * Performs a test of conditional independence X _||_ Y | Z1...Zn where all searchVariables are either continuous or discrete.
 * This test is valid for both ordinal and non-ordinal discrete searchVariables.
 *
 * This logisticRegression makes multiple assumptions: 1. IIA 2. Large sample size (multiple regressions needed on subsets of
 * sample)
 *
 * @author Joseph Ramsey
 * @author Augustus Mayo.
 */
public class IndTestMultinomialLogisticRegressionWald implements IndependenceTest {
    private final DataSet originalData;
    private final List<Node> searchVariables;
    private final DataSet internalData;
    private double alpha;
    private double lastP;
    private final Map<Node, List<Node>> variablesPerNode = new HashMap<>();
    private final LogisticRegression logisticRegression;
    private final RegressionDataset regression;
    private final boolean preferLinear;
    private boolean verbose;

    public IndTestMultinomialLogisticRegressionWald(DataSet data, double alpha, boolean preferLinear) {
        if (!(alpha >= 0 && alpha <= 1)) {
            throw new IllegalArgumentException("Alpha mut be in [0, 1]");
        }

        this.searchVariables = data.getVariables();
        this.originalData = data.copy();
        DataSet internalData = data.copy();
        this.alpha = alpha;
        this.preferLinear = preferLinear;

        List<Node> variables = internalData.getVariables();

        for (Node node : variables) {
            List<Node> nodes = expandVariable(internalData, node);
            this.variablesPerNode.put(node, nodes);
        }

        this.internalData = internalData;
        this.logisticRegression = new LogisticRegression(internalData);
        this.regression = new RegressionDataset(internalData);
    }

    /**
     * @return an Independence test for a subset of the searchVariables.
     */
    public IndependenceTest indTestSubset(List<Node> vars) {
        throw new UnsupportedOperationException();
    }

    /**
     * @return true if the given independence question is judged true, false if not. The independence question is of the
     * form x _||_ y | z, z = &lt;z1,...,zn&gt;, where x, y, z1,...,zn are searchVariables in the list returned by
     * getVariableNames().
     */
    public IndependenceResult checkIndependence(Node x, Node y, List<Node> z) {
        if (x instanceof DiscreteVariable && y instanceof DiscreteVariable) {
            return isIndependentMultinomialLogisticRegression(x, y, z);
        } else if (!this.preferLinear) {
            if (x instanceof DiscreteVariable)
                return isIndependentMultinomialLogisticRegression(x, y, z);
            else if (y instanceof DiscreteVariable)
                return isIndependentMultinomialLogisticRegression(y, x, z);
            else
                return isIndependentRegression(x, y, z);

        } else {
            if (x instanceof DiscreteVariable)
                return isIndependentRegression(y, x, z);
            else
                return isIndependentRegression(x, y, z);
        }
    }

    private List<Node> expandVariable(DataSet dataSet, Node node) {
        if (node instanceof ContinuousVariable) {
            return Collections.singletonList(node);
        }

        if (node instanceof DiscreteVariable && ((DiscreteVariable) node).getNumCategories() < 3) {
            return Collections.singletonList(node);
        }

        if (!(node instanceof DiscreteVariable)) {
            throw new IllegalArgumentException();
        }

        List<String> varCats = new ArrayList<>(((DiscreteVariable) node).getCategories());
        varCats.remove(0);
        List<Node> variables = new ArrayList<>();

        for (String cat : varCats) {

            Node newVar;

            do {
                String newVarName = node.getName() + "MULTINOM" + "." + cat;
                newVar = new DiscreteVariable(newVarName, 2);
            } while (dataSet.getVariable(newVar.getName()) != null);

            variables.add(newVar);

            dataSet.addVariable(newVar);
            int newVarIndex = dataSet.getColumn(newVar);
            int numCases = dataSet.getNumRows();

            for (int l = 0; l < numCases; l++) {
                Object dataCell = dataSet.getObject(l, dataSet.getColumn(node));
                int dataCellIndex = ((DiscreteVariable) node).getIndex(dataCell.toString());

                if (dataCellIndex == ((DiscreteVariable) node).getIndex(cat))
                    dataSet.setInt(l, newVarIndex, 1);
                else
                    dataSet.setInt(l, newVarIndex, 0);
            }
        }

        return variables;
    }

    private IndependenceResult isIndependentMultinomialLogisticRegression(Node x, Node y, List<Node> z) {
        if (!this.variablesPerNode.containsKey(x)) {
            throw new IllegalArgumentException("Unrecogized node: " + x);
        }

        if (!this.variablesPerNode.containsKey(y)) {
            throw new IllegalArgumentException("Unrecogized node: " + y);
        }

        for (Node node : z) {
            if (!this.variablesPerNode.containsKey(x)) {
                throw new IllegalArgumentException("Unrecogized node: " + node);
            }
        }

        List<Double> pValues = new ArrayList<>();

        int[] _rows = getNonMissingRows(x, y, z);
        this.logisticRegression.setRows(_rows);

        boolean independent;

        double p = 1.0;
        for (Node _x : this.variablesPerNode.get(x)) {

            // Without y
//            List<Node> regressors0 = new ArrayList<Node>();
//
//            for (Node _z : z) {
//                regressors0.addAll(variablesPerNode.get(_z));
//            }
//
//            LogisticRegression.Result result0 = logisticRegression.regress((DiscreteVariable) _x, regressors0);

            // With y.
            List<Node> regressors1 = new ArrayList<>(this.variablesPerNode.get(y));

            for (Node _z : z) {
                regressors1.addAll(this.variablesPerNode.get(_z));
            }

            LogisticRegression.Result result1 = this.logisticRegression.regress((DiscreteVariable) _x, regressors1);

            // Returns -2 LL
            int n = this.originalData.getNumRows();
            int k = regressors1.size() + 1;

            for (int i = 0; i < this.variablesPerNode.get(y).size(); i++) {
                double wald = Math.abs(result1.getCoefs()[i + 1] / result1.getStdErrs()[i + 1]);
                //double val = (1.0 - new NormalDistribution(0,1).cumulativeProbability(wald))*2;//two-tailed test
                //double val = 1-result1.getProbs()[i+1];

                //this is exactly the same test as the linear case
                double val = (1.0 - ProbUtils.tCdf(wald, n - k)) * 2;

                //System.out.println("My p: " + val + " Their p: " + otherVal + "1-their p:" + (1-otherVal));
                if (val < p) {
                    p = val;
                }

                //faster but won't find min p
                if (p <= this.alpha) {
                    independent = false;
                    this.lastP = p;

                    if (independent) {
                        TetradLogger.getInstance().log("independencies", SearchLogUtils.independenceFactMsg(x, y, z, p));
                    } else {
                        TetradLogger.getInstance().log("dependencies", SearchLogUtils.dependenceFactMsg(x, y, z, p));
                    }

                    return new IndependenceResult(new IndependenceFact(x, y, z), independent, p);
                }
            }

        }

        // Choose the minimum of the p-values
        // This is only one method that can be used, this requires every coefficient to be significant

        independent = p > this.alpha;

        this.lastP = p;

        if (this.verbose) {
            if (independent) {
                TetradLogger.getInstance().forceLogMessage(
                        SearchLogUtils.independenceFactMsg(x, y, z, p));
            }
        }

        return new IndependenceResult(new IndependenceFact(x, y, z), independent, p);
    }

    // This takes an inordinate amount of time. -jdramsey 20150929
    private int[] getNonMissingRows(Node x, Node y, List<Node> z) {

        int[] _rows = new int[this.internalData.getNumRows()];
        for (int k = 0; k < _rows.length; k++) _rows[k] = k;

        return _rows;
    }

    private boolean isMissing(Node x, int i) {
        int j = this.internalData.getColumn(x);

        if (x instanceof DiscreteVariable) {
            int v = this.internalData.getInt(i, j);

            if (v == -99) {
                return true;
            }
        }

        if (x instanceof ContinuousVariable) {
            double v = this.internalData.getDouble(i, j);

            return Double.isNaN(v);
        }

        return false;
    }

    private IndependenceResult isIndependentRegression(Node x, Node y, List<Node> z) {
        if (!this.variablesPerNode.containsKey(x)) {
            throw new IllegalArgumentException("Unrecogized node: " + x);
        }

        if (!this.variablesPerNode.containsKey(y)) {
            throw new IllegalArgumentException("Unrecogized node: " + y);
        }

        for (Node node : z) {
            if (!this.variablesPerNode.containsKey(node)) {
                throw new IllegalArgumentException("Unrecogized node: " + node);
            }
        }

        List<Node> regressors = new ArrayList<>();
        if (y instanceof ContinuousVariable) {
            regressors.add(this.internalData.getVariable(y.getName()));
        } else {
            regressors.addAll(this.variablesPerNode.get(y));
        }

        for (Node _z : z) {
            regressors.addAll(this.variablesPerNode.get(_z));
        }

        int[] _rows = getNonMissingRows(x, y, z);
        this.regression.setRows(_rows);

        RegressionResult result;

        try {
            result = this.regression.regress(x, regressors);
        } catch (Exception e) {
            return new IndependenceResult(new IndependenceFact(x, y, z), false, Double.NaN);
        }

        double p = 1;
        if (y instanceof ContinuousVariable) {
            p = result.getP()[1];
        } else {
            for (int i = 0; i < this.variablesPerNode.get(y).size(); i++) {
                double val = result.getP()[1 + i];
                if (val < p)
                    p = val;
            }
        }
        this.lastP = p;

        boolean independent = p > this.alpha;

        if (this.verbose) {
            if (independent) {
                TetradLogger.getInstance().forceLogMessage(
                        SearchLogUtils.independenceFactMsg(x, y, z, p));
            }
        }

        return new IndependenceResult(new IndependenceFact(x, y, z), independent, p);
    }

    /**
     * @return the probability associated with the most recently executed independence test, of Double.NaN if p value is
     * not meaningful for tis test.
     */
    public double getPValue() {
        return this.lastP; //STUB
    }

    /**
     * @return the list of searchVariables over which this independence checker is capable of determinining independence
     * relations.
     */
    public List<Node> getVariables() {
        return this.searchVariables; // Make sure the variables from the ORIGINAL data set are returned, not the modified dataset!
    }

    /**
     * @return the list of variable varNames.
     */
    public List<String> getVariableNames() {
        List<Node> variables = getVariables();
        List<String> variableNames = new ArrayList<>();
        for (Node variable1 : variables) {
            variableNames.add(variable1.getName());
        }
        return variableNames;
    }

    public Node getVariable(String name) {
        for (int i = 0; i < getVariables().size(); i++) {
            Node variable = getVariables().get(i);
            if (variable.getName().equals(name)) {
                return variable;
            }
        }

        return null;
    }

    /**
     * @return true if y is determined the variable in z.
     */
    public boolean determines(List<Node> z, Node y) {
        return false; //stub
    }

    /**
     * @return the significance level of the independence test.
     * @throws UnsupportedOperationException if there is no significance level.
     */
    public double getAlpha() {
        return this.alpha; //STUB
    }

    /**
     * Sets the significance level.
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public DataSet getData() {
        return this.originalData;
    }

    @Override
    public ICovarianceMatrix getCov() {
        return null;
    }

    @Override
    public List<DataSet> getDataSets() {
        return null;
    }

    @Override
    public int getSampleSize() {
        return 0;
    }

    @Override
    public List<Matrix> getCovMatrices() {
        return null;
    }

    @Override
    public double getScore() {
        return getPValue();
    }

    /**
     * @return a string representation of this test.
     */
    public String toString() {
        NumberFormat nf = new DecimalFormat("0.0000");
        return "Multinomial Logistic Regression, alpha = " + nf.format(getAlpha());
    }

    @Override
    public boolean isVerbose() {
        return this.verbose;
    }

    @Override
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}

