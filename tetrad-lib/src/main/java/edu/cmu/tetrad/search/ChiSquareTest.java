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

import edu.cmu.tetrad.data.CellTable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.util.CombinationIterator;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.util.Arrays;

/**
 * Calculates marginal chi square test results for a discrete dataset.
 *
 * @author Frank Wimberly original version
 * @author Joseph Ramsey revision 10/01, believe it or now also 2020.
 */
public class ChiSquareTest {

    /**
     * The data set this test uses.
     */
    private final DataSet dataSet;

    /**
     * The number of values for each variable in the data.
     */
    private final int[] dims;

    /**
     * Stores the data in the form of a cell table.
     */
    private final CellTable cellTable;

    /**
     * The significance level of the test.
     */
    private double alpha;

    //==============================CONSTRUCTORS=========================//

    /**
     * Constructs a test using the given data set and significance level.
     *
     * @param dataSet A data set consisting entirely of discrete variables.
     * @param alpha   The significance level, usually 0.05.
     */
    public ChiSquareTest(DataSet dataSet, double alpha) {
        if (alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("Significance level must be in " +
                    "[0, 1]: " + alpha);
        }

        this.dims = new int[dataSet.getNumColumns()];

        for (int i = 0; i < getDims().length; i++) {
            DiscreteVariable variable =
                    (DiscreteVariable) dataSet.getVariable(i);
            this.getDims()[i] = variable.getNumCategories();
        }

        this.dataSet = dataSet;
        this.alpha = alpha;
        this.cellTable = new CellTable(null);
        this.getCellTable().setMissingValue(DiscreteVariable.MISSING_VALUE);
    }

    /**
     * Calculates chi square for a conditional crosstabulation table for independence question 0 _||_ 1 | 2, 3, ...max
     * by summing up chi square and degrees of freedom for each conditional table in turn, where rows or columns that
     * consist entirely of zeros have been removed.
     */
    public synchronized Result calcChiSquare(int[] testIndices) {

        // Reset the cell table for the columns referred to in
        // 'testIndices.' Do cell coefs for those columns.
        this.getCellTable().addToTable(getDataSet(), testIndices);

        // Indicator arrays to tell the cell table which margins
        // to calculate. For x _||_ y | z1, z2, ..., we want to
        // calculate the margin for x, the margin for y, and the
        // margin for x and y. (These will be used later.)
        int[] firstVar = {0};
        int[] secondVar = {1};
        int[] bothVars = {0, 1};

        double xSquare = 0.0;
        int df = 0;

        int[] condDims = new int[testIndices.length - 2];
        System.arraycopy(selectFromArray(getDims(), testIndices), 2, condDims, 0,
                condDims.length);

        int[] coords = new int[testIndices.length];
        int numRows = this.getCellTable().getNumValues(0);
        int numCols = this.getCellTable().getNumValues(1);

        CombinationIterator combinationIterator =
                new CombinationIterator(condDims);

        // Make a chi square table for each condition combination, strike zero rows and columns and calculate
        // chi square and degrees of freedom for the remaining rows and columns in the table. See Friedman.
        while (combinationIterator.hasNext()) {
            boolean[] attestedRows = new boolean[numRows];
            boolean[] attestedCols = new boolean[numCols];

            Arrays.fill(attestedRows, true);
            Arrays.fill(attestedCols, true);

            int[] combination = combinationIterator.next();

            System.arraycopy(combination, 0, coords, 2, combination.length);

            long total = getCellTable().calcMargin(coords, bothVars);

            if (total == 0) continue;

            double _xSquare = 0.0;

            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    coords[0] = i;
                    coords[1] = j;

                    long sumRow = getCellTable().calcMargin(coords, firstVar);
                    long sumCol = getCellTable().calcMargin(coords, secondVar);
                    long observed = getCellTable().getValue(coords);

                    if (sumRow == 0L) {
                        attestedRows[i] = false;
                        continue;
                    }

                    if (sumCol == 0L) {
                        attestedCols[j] = false;
                        continue;
                    }

                    double expected = (sumRow * sumCol) / (double) total;
                    _xSquare += Math.pow(observed - expected, 2.0) / expected;
                }
            }

            int numAttestedRows = 0;
            int numAttestedCols = 0;

            for (boolean attestedRow : attestedRows) {
                if (attestedRow) {
                    numAttestedRows++;
                }
            }

            for (boolean attestedCol : attestedCols) {
                if (attestedCol) {
                    numAttestedCols++;
                }
            }

            if (numAttestedRows > 0 && numAttestedCols > 0) {
                df += (numAttestedRows - 1) * (numAttestedCols - 1);
                xSquare += _xSquare;
            }
        }

        // If df == 0, this is definitely an indepedent table.
        if (df == 0) {
            final double pValue = 1.0;
            return new Result(xSquare, pValue, 0, true);
        }

        double pValue = 1.0 - new ChiSquaredDistribution(df).cumulativeProbability(xSquare);
        boolean indep = (pValue > getAlpha());
        return new Result(xSquare, pValue, df, indep);
    }

    /**
     * @param testIndices An array of indices for variables in the dataset supplied in the constructor.
     * @param p           The probability that some marginal for some table dominates. A good value is 0.99.
     * @return True if the variable at index 0 is determined by the variables at the other indices.
     */
    public boolean isDetermined(int[] testIndices, double p) {

        // Reset the cell table for the columns referred to in
        // 'testIndices.' Do cell coefs for those columns.
        this.getCellTable().addToTable(getDataSet(), testIndices);

        // Indicator arrays to tell the cell table which margins
        // to calculate. For x _||_ y | z1, z2, ..., we want to
        // calculate the margin for x, the margin for y, and the
        // margin for x and y. (These will be used later.)
        int[] firstVar = {0};

        int[] condDims = new int[testIndices.length - 1];
        System.arraycopy(selectFromArray(getDims(), testIndices), 1, condDims, 0,
                condDims.length);

        int[] coords = new int[testIndices.length];
        int numValues = this.getCellTable().getNumValues(0);

        CombinationIterator combinationIterator =
                new CombinationIterator(condDims);

        while (combinationIterator.hasNext()) {
            int[] combination = combinationIterator.next();
            System.arraycopy(combination, 0, coords, 1, combination.length);

            long total = this.getCellTable().calcMargin(coords, firstVar);

            if (total == 0) {
                continue;
            }

            boolean dominates = false;

            for (int i = 0; i < numValues; i++) {
                coords[0] = i;

                long numi = this.getCellTable().getValue(coords);

                if ((double) numi / total >= p) {
                    dominates = true;
                }
            }

            if (!dominates) {
                return false;
            }
        }

        return true;
    }

    /**
     * @return the getModel significance level being used for tests.
     */
    public double getAlpha() {
        return this.alpha;
    }

    /**
     * Sets the significance level to be used for tests.
     */
    public void setAlpha(double alpha) {
        if (alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("Significance level must be in " +
                    "[0, 1]: " + alpha);
        }

        this.alpha = alpha;
    }

    //================================PRIVATE==============================//

    public int[] selectFromArray(int[] arr, int[] indices) {
        int[] retArr = new int[indices.length];

        for (int i = 0; i < indices.length; i++) {
            retArr[i] = arr[indices[i]];
        }

        return retArr;
    }

    public DataSet getDataSet() {
        return this.dataSet;
    }

    public int[] getDims() {
        return this.dims;
    }

    public CellTable getCellTable() {
        return this.cellTable;
    }

    //===============================CLASSES==============================//

    /**
     * Simple class to store the parameters of the result returned by the G Square test.
     *
     * @author Frank Wimberly
     */
    public static class Result {

        /**
         * The chi square value.
         */
        private final double chiSquare;

        /**
         * The pValue of the result.
         */
        private final double pValue;

        /**
         * The adjusted degrees of freedom.
         */
        private final int df;

        /**
         * Whether the conditional independence holds or not. (True if it does, false if it doesn't.
         */
        private final boolean isIndep;

        /**
         * Constructs a new g square result using the given parameters.
         */
        public Result(double chiSquare, double pValue, int df, boolean isIndep) {
            this.chiSquare = chiSquare;
            this.pValue = pValue;
            this.df = df;
            this.isIndep = isIndep;
        }

        public double getXSquare() {
            return this.chiSquare;
        }

        public double getPValue() {
            return this.pValue;
        }

        public int getDf() {
            return this.df;
        }

        public boolean isIndep() {
            return this.isIndep;
        }
    }
}





