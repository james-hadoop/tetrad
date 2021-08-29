///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
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

import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * Implements the continuous BIC score for FGES.
 *
 * @author Joseph Ramsey
 */
public class KimEtAlScores implements Score {

    // The dataset.
    private DataSet dataSet;

    // The correlation matrix.
    private ICovarianceMatrix covariances;

    // The variables of the covariance matrix.
    private List<Node> variables;

    // The sample size of the covariance matrix.
    private final int sampleSize;

    // True if verbose output should be sent to out.
    private boolean verbose = false;

    // The true error variance
    private double trueErrorVariance = 1.0;

    // The rule type to use.
    private RuleType ruleType = RuleType.MANUAL;

    // Sample size or equivalent sample size.
    private double N;

    // Manually set lambda, by default log(n);
    private double lambda;
    private boolean calculateRowSubsets = false;
    Matrix data;
    //    private double penaltyDiscount;
    private double correlationThreshold = 1.0;
    private boolean takeLog = true;
    private boolean calculateSquareEuclideanNorms = false;
    private double penaltyDiscount = 1;

    /**
     * Constructs the score using a covariance matrix.
     */
    public KimEtAlScores(ICovarianceMatrix covariances, double correlationThreshold) {
        if (covariances == null) {
            throw new NullPointerException();
        }

        this.correlationThreshold = correlationThreshold;

        setCovariances(covariances);
        this.variables = covariances.getVariables();
        this.sampleSize = covariances.getSampleSize();
        this.setLambda(log(this.sampleSize));
    }

    /**
     * Constructs the score using a covariance matrix.
     */
    public KimEtAlScores(DataSet dataSet, double correlationThreshold) {
        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.correlationThreshold = correlationThreshold;

        dataSet = DataUtils.center(dataSet);

        double[][] cov = new double[dataSet.getNumColumns()][dataSet.getNumColumns()];

        for (int i = 0; i < dataSet.getNumColumns(); i++) {
            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                double sum = 0.0;

                for (int k = 0; k < dataSet.getNumRows(); k++) {
                    sum += dataSet.getDouble(k, i) * dataSet.getDouble(k, j);
                }

                cov[i][j] = sum / dataSet.getNumRows();
            }
        }

        CovarianceMatrix covarianceMatrix = new CovarianceMatrix(dataSet.getVariables(), cov, dataSet.getNumRows());

        this.data = dataSet.getDoubleData();
        this.dataSet = dataSet;

        if (!dataSet.existsMissingValue()) {
            setCovariances(covarianceMatrix);// new CovarianceMatrix(dataSet, false));
            this.variables = covariances.getVariables();
            this.sampleSize = covariances.getSampleSize();
            calculateRowSubsets = false;
            return;
        }

        this.variables = dataSet.getVariables();
        this.sampleSize = dataSet.getNumRows();
        calculateRowSubsets = true;
    }

    @Override
    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScoreDiff(x, y, new int[0]);
    }

    public double localScore(int i, int... parents) {
        double sn = 12;

        if (parents.length > sn) return Double.NEGATIVE_INFINITY;
        final int k = parents.length;

        // Only do this once.
        double pn = variables.size();
        pn = min(pn, sn);
        double n = N;

        double varry = SemBicScore.getVarRy(i, parents, data, covariances, calculateRowSubsets);

        double lambda;

        // Defaults to the manually set lambda.
        if (ruleType == RuleType.MANUAL) {
            lambda = this.lambda;
        } else if (ruleType == RuleType.BIC) {
            lambda = log(n);
        } else if (ruleType == RuleType.GIC2) {

            // Following Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high dimensions.
            // The Journal 0of Machine Learning Research, 13(1), 1037-1057.
            lambda = pow(n, .33);
        } else if (ruleType == RuleType.RIC) {

            // Following Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high dimensions.
            // The Journal 0of Machine Learning Research, 13(1), 1037-1057.
            lambda = 2.2 * (log(pn));
        } else if (ruleType == RuleType.RICc) {

            // Following Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high dimensions.
            // The Journal of Machine Learning Research, 13(1), 1037-1057.
            lambda = 2 * (log(pn) + log(log(pn)));
        } else if (ruleType == RuleType.GIC5) {

            // Following Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high dimensions.
            // The Journal of Machine Learning Research, 13(1), 1037-1057.
            lambda = log(log(n)) * (log(pn));
        } else if (ruleType == RuleType.GIC6) {

            // Following Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high dimensions.
            // The Journal of Machine Learning Resjearch, 13(1), 1037-1057.
            lambda = log(n) * log(pn);
        } else {
            throw new IllegalStateException("That lambda rule is not configured: " + ruleType);
        }

//        double c = penaltyDiscount;

        if (takeLog) {
            return -n * log(varry) - lambda * getPenaltyDiscount() * k;
        } else {
            return -n * (varry) - lambda * getPenaltyDiscount() * k * trueErrorVariance;
        }

    }

    /**
     * Specialized scoring method for a single parent. Used to speed up the effect edges search.
     */
    public double localScore(int i, int parent) {
        return localScore(i, new int[]{parent});
    }

    /**
     * Specialized scoring method for no parents. Used to speed up the effect edges search.
     */
    public double localScore(int i) {
        return localScore(i, new int[0]);
    }

    public double getTrueErrorVariance() {
        return trueErrorVariance;
    }

    public ICovarianceMatrix getCovariances() {
        return covariances;
    }

    public int getSampleSize() {
        return sampleSize;
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    public DataSet getDataSet() {
        return dataSet;
    }

    public void setTrueErrorVariance(double trueErrorVariance) {
        this.trueErrorVariance = trueErrorVariance;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    @Override
    public List<Node> getVariables() {
        return variables;
    }

    public void setVariables(List<Node> variables) {
        if (covariances != null) {
            covariances.setVariables(variables);
        }

        this.variables = variables;
    }

    @Override
    public Node getVariable(String targetName) {
        for (Node node : variables) {
            if (node.getName().equals(targetName)) {
                return node;
            }
        }

        return null;
    }

    @Override
    public int getMaxDegree() {
        return (int) Math.ceil(log(sampleSize));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        int i = variables.indexOf(y);

        int[] k = new int[z.size()];

        for (int t = 0; t < z.size(); t++) {
            k[t] = variables.indexOf(z.get(t));
        }

        double v = localScore(i, k);

        return Double.isNaN(v);
    }

    @Override
    public Score defaultScore() {
        return new KimEtAlScores(covariances, 1.0);
    }

    private void setCovariances(ICovarianceMatrix covariances) {
        CorrelationMatrix correlations = new CorrelationMatrix(covariances);
        this.covariances = covariances;
//        this.covariances = covariances;

        boolean exists = false;

        for (int i = 0; i < correlations.getSize(); i++) {
            for (int j = 0; j < correlations.getSize(); j++) {
                if (i == j) continue;
                double r = correlations.getValue(i, j);
                if (abs(r) > correlationThreshold) {
                    System.out.println("Absolute correlation too high: " + r);
                    exists = true;
                }
            }
        }

        if (exists) {
            throw new IllegalArgumentException("Some correlations are too high (> " + correlationThreshold
                    + ") in absolute value.");
        }


        this.N = covariances.getSampleSize();
    }

    private static int[] append(int[] z, int x) {
        int[] _z = Arrays.copyOf(z, z.length + 1);
        _z[z.length] = x;
        return _z;
    }

    public void setRuleType(RuleType ruleType) {
        this.ruleType = ruleType;
    }

    public RuleType getRuleType() {
        return ruleType;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
    }

//    public void setPenaltyDiscount(double penaltyDiscount) {
//        this.penaltyDiscount = penaltyDiscount;
//    }

    public void setCorrelationThreshold(double correlationThreshold) {
        this.correlationThreshold = correlationThreshold;
    }

    public void setTakeLog(boolean takeLog) {
        this.takeLog = takeLog;
    }

    public void setCalculateSquareEuclideanNorms(boolean calculateSquareEuclideanNorms) {
        this.calculateSquareEuclideanNorms = calculateSquareEuclideanNorms;
    }

    public double getPenaltyDiscount() {
        return penaltyDiscount;
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }

    public enum RuleType {MANUAL, BIC, NANDY, GIC2, RIC, RICc, GIC5, GIC6}
}


