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

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.LogisticRegression;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.TetradLogger;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * Performs a test of conditional independence X _||_ Y | Z1...Zn where all searchVariables are either continuous or discrete.
 * This test is valid for both ordinal and non-ordinal discrete searchVariables.
 * <p>
 * This logisticRegression makes multiple assumptions: 1. IIA 2. Large sample size (multiple regressions needed on subsets of
 * sample)
 *
 * @author Joseph Ramsey
 * @author Augustus Mayo.
 */
public class IndTestMultinomialLogisticRegression implements IndependenceTest {
    private final DataSet originalData;
    private final List<Node> searchVariables;
    private final DataSet internalData;
    private double alpha;
    private double lastP;
    private final Map<Node, List<Node>> variablesPerNode = new HashMap<>();
    private final LogisticRegression logisticRegression;
    private final RegressionDataset regression;
    private boolean verbose;

    public IndTestMultinomialLogisticRegression(DataSet data, double alpha) {
        this.searchVariables = data.getVariables();
        this.originalData = data.copy();
        DataSet internalData = data.copy();
        this.alpha = alpha;

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
     * form x _||_ y | z, z = [z1,...,zn], where x, y, z1,...,zn are searchVariables in the list returned by
     * getVariableNames().
     */
    public IndependenceResult checkIndependence(Node x, Node y, List<Node> z) {
        if (x instanceof DiscreteVariable) {
            return isIndependentMultinomialLogisticRegression(x, y, z);
        } else if (y instanceof DiscreteVariable) {
            return isIndependentMultinomialLogisticRegression(y, x, z);
        } else {
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

        int[] _rows = getNonMissingRows();
        this.logisticRegression.setRows(_rows);

        for (Node _x : this.variablesPerNode.get(x)) {

            // Without y
            List<Node> regressors0 = new ArrayList<>();

            for (Node _z : z) {
                regressors0.addAll(this.variablesPerNode.get(_z));
            }

            LogisticRegression.Result result0 = logisticRegression.regress((DiscreteVariable) _x, regressors0);

            // With y.
            List<Node> regressors1 = new ArrayList<>(variablesPerNode.get(y));

            for (Node _z : z) {
                regressors1.addAll(variablesPerNode.get(_z));
            }

            LogisticRegression.Result result1 = this.logisticRegression.regress((DiscreteVariable) _x, regressors1);

            // Returns -2 LL
            double ll0 = result0.getLogLikelihood();
            double ll1 = result1.getLogLikelihood();

            double chisq = (ll0 - ll1);
            int df = this.variablesPerNode.get(y).size();
            double p = 1.0 - new ChiSquaredDistribution(df).cumulativeProbability(chisq);
            pValues.add(p);
        }

        double p = 1.0;

        // Choose the minimum of the p-values
        // This is only one method that can be used, this requires every coefficient to be significant
        for (double val : pValues) {
            if (val < p) p = val;
        }

        boolean independent = p > this.alpha;

        this.lastP = p;

        if (this.verbose) {
            if (independent) {
                TetradLogger.getInstance().forceLogMessage(
                        SearchLogUtils.independenceFactMsg(x, y, z, p));
            }
        }

        return new IndependenceResult(new IndependenceFact(x, y, z), independent, p);
    }

    int[] _rows;

    private int[] getNonMissingRows() {
        if (this._rows == null) {
            this._rows = new int[this.internalData.getNumRows()];
            for (int k = 0; k < this._rows.length; k++) this._rows[k] = k;
        }

        return this._rows;
    }

    private IndependenceResult isIndependentRegression(Node x, Node y, List<Node> z) {
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

        List<Node> regressors = new ArrayList<>();
        regressors.add(this.internalData.getVariable(y.getName()));

        for (Node _z : z) {
            regressors.addAll(this.variablesPerNode.get(_z));
        }

        int[] _rows = getNonMissingRows();
        this.regression.setRows(_rows);

        RegressionResult result;

        try {
            result = this.regression.regress(x, regressors);
        } catch (Exception e) {
            return new IndependenceResult(new IndependenceFact(x, y, z), false, Double.NaN);
        }

        double p = result.getP()[1];
        this.lastP = p;

        boolean indep = p > this.alpha;

        if (this.verbose) {
            if (indep) {
                TetradLogger.getInstance().log("independencies", SearchLogUtils.independenceFactMsg(x, y, z, p));
            } else {
                TetradLogger.getInstance().log("dependencies", SearchLogUtils.dependenceFactMsg(x, y, z, p));
            }
        }

        return new IndependenceResult(new IndependenceFact(x, y, z), indep, p);
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
        return 0;
    }

    /**
     * @return a string representation of this test.
     */
    public String toString() {
        NumberFormat nf = new DecimalFormat("0.0000");
        return "Multinomial Logistic Regression, alpha = " + nf.format(getAlpha());
    }

    public boolean isVerbose() {
        return this.verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}


