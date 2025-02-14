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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.NumberFormatUtil;
import edu.cmu.tetrad.util.TetradLogger;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Checks the conditional independence X _||_ Y | S, where S is a set of discrete variable, and X and Y are discrete
 * variable not in S, by applying a conditional G Square test. A description of such a test is given in Fienberg, "The
 * Analysis of Cross-Classified Categorical Data," 2nd edition. The formula for degrees of freedom used in this test are
 * equivalent to the formulation on page 142 of Fienberg.
 *
 * @author Joseph Ramsey
 * @see GSquareTest
 */
public final class IndTestGSquare implements IndependenceTest {

    /**
     * The G Square tester.
     */
    private final GSquareTest gSquareTest;

    /**
     * The variables in the discrete data sets or which conditional independence judgements are desired.
     */
    private final List<Node> variables;

    /**
     * The dataset of discrete variables.
     */
    private final DataSet dataSet;

    /**
     * The significance level for the test.
     */
    private final double alpha;

    /**
     * The p value associated with the most recent call of isIndependent.
     */
    private double pValue;

    /**
     * The lower bound of percentages of observation of some category in the data, given some particular combination of
     * values of conditioning variables, that coefs as 'determining."
     */
    private double determinationP = 0.99;

    /**
     * The standard number formatter for Tetrad.
     */
    private static final NumberFormat nf = NumberFormatUtil.getInstance().getNumberFormat();


    private boolean verbose;

    /**
     * Constructs a new independence checker to check conditional independence facts for discrete data using a g square
     * test.
     *
     * @param dataSet the discrete data set.
     * @param alpha   the significance level of the tests.
     */
    public IndTestGSquare(DataSet dataSet, double alpha) {

        // The g square test requires as parameters: (a) the data set
        // itself, (b) an array containing the number of values for
        // each variable in order, and (c) the significance level of
        // the test. Also, in order to perform specific conditional
        // independence tests, it is necessary to construct an array
        // containing the variables of the requested test, in
        // order. Specifically, to test whether X _||_ Y | Z1, ...,
        // Zn, an array is constructed with the indices, in order of
        // X, Y, Z1, ..., Zn. Therefore, the indices of these
        // variables must be stored. We do this by storing the
        // variables themselves in a List.
        this.dataSet = dataSet;
        this.alpha = alpha;

        this.variables = new ArrayList<>(dataSet.getVariables());
        this.gSquareTest = new GSquareTest(dataSet, alpha);
    }

    /**
     * Creates a new IndTestGSquare for a subset of the variables.
     */
    public IndependenceTest indTestSubset(List<Node> vars) {
        if (vars.isEmpty()) {
            throw new IllegalArgumentException("Subset may not be empty.");
        }

        int[] indices = new int[vars.size()];
        int j = -1;

        for (int i = 0; i < this.variables.size(); i++) {
            if (!vars.contains(this.variables.get(i))) {
                continue;
            }

            indices[++j] = i;
        }

        DataSet newDataSet = this.dataSet.subsetColumns(indices);
        return new IndTestGSquare(newDataSet, this.alpha);
    }

    /**
     * @return the p value associated with the most recent call of isIndependent.
     */
    public double getPValue() {
        return this.pValue;
    }

    /**
     * Determines whether variable x is independent of variable y given a list of conditioning varNames z.
     *
     * @param x the one variable being compared.
     * @param y the second variable being compared.
     * @param z the list of conditioning varNames.
     * @return true iff x _||_ y | z.
     */
    public IndependenceResult checkIndependence(Node x, Node y, List<Node> z) {
        if (x == null) {
            throw new NullPointerException();
        }

        if (y == null) {
            throw new NullPointerException();
        }

        if (z == null) {
            throw new NullPointerException();
        }

        for (Node node : z) {
            if (node == null) {
                throw new NullPointerException();
            }
        }

        // For testing x, y given z1,...,zn, set up an array of length
        // n + 2 containing the indices of these variables in order.
        int[] testIndices = new int[2 + z.size()];

        testIndices[0] = this.variables.indexOf(x);
        testIndices[1] = this.variables.indexOf(y);

        for (int i = 0; i < z.size(); i++) {
            testIndices[i + 2] = this.variables.indexOf(z.get(i));
        }

        // the following is lame code--need a better test
        for (int i = 0; i < testIndices.length; i++) {
            if (testIndices[i] < 0) {
                throw new IllegalArgumentException(
                        "Variable " + i + " was not used in the constructor.");
            }
        }

        //        System.out.println("Testing " + x + " _||_ " + y + " | " + z);

        GSquareTest.Result result = this.gSquareTest.calcGSquare(testIndices);
        this.pValue = result.getPValue();

        if (this.verbose) {
            if (result.isIndep()) {
                TetradLogger.getInstance().forceLogMessage(
                        SearchLogUtils.independenceFactMsg(x, y, z, getPValue()));
            }
        }

        return new IndependenceResult(new IndependenceFact(x, y, z),
                result.isIndep(), result.getPValue());
    }

    /**
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     *
     * @param alpha the new significance level.
     */
    public void setAlpha(double alpha) {
        this.gSquareTest.setAlpha(alpha);
    }

    /**
     * Gets the getModel significance level.
     *
     * @return this number.
     */
    public double getAlpha() {
        return this.gSquareTest.getAlpha();
    }

    /**
     * @return the list of variables over which this independence checker is capable of determinine independence
     * relations-- that is, all the variables in the given graph or the given data set.
     */
    public List<Node> getVariables() {
        return Collections.unmodifiableList(this.variables);
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

    public String toString() {
        return "G Square, alpha = " + IndTestGSquare.nf.format(getAlpha());
    }

    /**
     * @param z  The list of variables z1,...,zn with respect to which we want to know whether z determines x oir z.
     * @param x1 The one variable whose determination by z we want to know.
     * @return true if it is estimated that z determines x or z determines y.
     */
    public boolean determines(List<Node> z, Node x1) {
        if (z == null) {
            throw new NullPointerException();
        }

        for (Node node : z) {
            if (node == null) {
                throw new NullPointerException();
            }
        }

        // For testing x, y given z1,...,zn, set up an array of length
        // n + 2 containing the indices of these variables in order.
        int[] testIndices = new int[1 + z.size()];
        testIndices[0] = this.variables.indexOf(x1);

        for (int i = 0; i < z.size(); i++) {
            testIndices[i + 1] = this.variables.indexOf(z.get(i));
        }

        // the following is lame code--need a better test
        for (int i = 0; i < testIndices.length; i++) {
            if (testIndices[i] < 0) {
                throw new IllegalArgumentException(
                        "Variable " + i + "was not used in the constructor.");
            }
        }

        //        System.out.println("Testing " + x + " _||_ " + y + " | " + z);

        boolean determined =
                this.gSquareTest.isDetermined(testIndices, getDeterminationP());

        if (determined) {
            StringBuilder sb = new StringBuilder();
            sb.append("Determination found: ").append(x1).append(
                    " is determined by {");

            for (int i = 0; i < z.size(); i++) {
                sb.append(z.get(i));

                if (i < z.size() - 1) {
                    sb.append(", ");
                }
            }

            sb.append("}");

            TetradLogger.getInstance().log("independencies", sb.toString());
        }

        return determined;
    }

    public double getDeterminationP() {
        return this.determinationP;
    }

    public void setDeterminationP(double determinationP) {
        this.determinationP = determinationP;
    }

    public DataSet getData() {
        return this.dataSet;
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

    @Override
    public boolean isVerbose() {
        return this.verbose;
    }

    @Override
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}




