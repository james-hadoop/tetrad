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

package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.MatrixUtils;
import edu.cmu.tetrad.util.ProbUtils;
import edu.cmu.tetrad.util.TetradSerializable;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;

import static edu.cmu.tetrad.search.LinearGaussianBicScore.bStar;
import static java.lang.Math.log;

/**
 * Estimates a SemIm given a CovarianceMatrix and a SemPm. (A DataSet may be
 * substituted for the CovarianceMatrix.) Uses regression to do the estimation,
 * so this is only for DAG models. But the DAG model may be reset on the fly
 * and the estimation redone. Variables whose parents have not changed will
 * not be reestimated. Intended to speed up estimation for algorithm that
 * require repeated estimation of DAG models over the same variables.
 * Assumes all variables are measured.
 *
 * @author Joseph Ramsey
 */
public final class FmlBicScorer implements TetradSerializable, Scorer {
    static final long serialVersionUID = 23L;

    private final ICovarianceMatrix covMatrix;
    private final Matrix edgeCoef;
    private final Matrix errorCovar;
    private final List<Node> variables;
    private final Matrix sampleCovar;
    private final Map<Node, Integer> nodesHash;
    private final Map<Integer, Set<Node>> parentsMap = new HashMap<>();
    private DataSet dataSet = null;
    private Graph dag = null;
    private Matrix implCovar;
    private double logDetSample;
    private double penaltyDiscount = 1;
    private Edge changedEdge = null;

    //=============================CONSTRUCTORS============================//

    /**
     * Constructs a new SemEstimator that uses the specified optimizer.
     *
     * @param dataSet a DataSet, all of whose variables are contained in
     *                the given SemPm. (They are identified by name.)
     */
    public FmlBicScorer(DataSet dataSet, double penaltyDiscount) {
        this(new CovarianceMatrix(dataSet));
        this.dataSet = dataSet;
        this.penaltyDiscount = penaltyDiscount;
    }

    /**
     * Constructs a new SemEstimator that uses the specified optimizer.
     *
     * @param covMatrix a covariance matrix, all of whose variables are
     *                  contained in the given SemPm. (They are identified by
     *                  name.)
     */
    public FmlBicScorer(ICovarianceMatrix covMatrix) {
        if (covMatrix == null) {
            throw new NullPointerException(
                    "CovarianceMatrix must not be null.");
        }

        this.variables = covMatrix.getVariables();
        this.covMatrix = covMatrix;

        int m = getVariables().size();
        this.edgeCoef = new Matrix(m, m);
        this.errorCovar = new Matrix(m, m);
        this.sampleCovar = covMatrix.getMatrix();

        this.nodesHash = new HashMap<>();
        for (int i = 0; i < variables.size(); i++) nodesHash.put(variables.get(i), i);
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static Scorer serializableInstance() {
        return new FmlBicScorer(CovarianceMatrix.serializableInstance());
    }

    private static int[] concat(int i, int[] parents) {
        int[] all = new int[parents.length + 1];
        all[0] = i;
        System.arraycopy(parents, 0, all, 1, parents.length);
        return all;
    }

    public DataSet getDataSet() {
        return dataSet;
    }

    public int getNumFreeParams() {
        return dag.getEdges().size() + dag.getNodes().size();
    }

    public int getDof() {
        return (dag.getNodes().size() * (dag.getNodes().size() + 1)) / 2 - getNumFreeParams();
    }

    public int getSampleSize() {
        return covMatrix.getSampleSize();
    }

    public List<Node> getMeasuredNodes() {
        return this.getVariables();
    }

    public Matrix getEdgeCoef() {
        return edgeCoef;
    }

    public Matrix getErrorCovar() {
        return errorCovar;
    }

    public List<Node> getVariables() {
        return variables;
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public void resetParameters(Edge edge) {
        List<Node> variables;

        if (changedEdge != null) {
            variables = new ArrayList<>();
//            variables.add(changedEdge.getNode1());
            variables.add(changedEdge.getNode2());
        } else {
            variables = this.variables;
        }

        resetParameters(variables);
    }

    //==============================PUBLIC METHODS=========================//

    private void resetParameters(List<Node> variables) {
        for (Node node : variables) {
            int idx = nodesHash.get(node);
            List<Node> parents = dag.getParents(node);

            if (new HashSet<>(parents).equals(parentsMap.get(idx))) continue;
            parentsMap.put(idx, new HashSet<>(parents));

            for (int _j = 0; _j < getVariables().size(); _j++) {
                edgeCoef.set(_j, idx, 0);
            }

            int[] pp = new int[parents.size()];
            for (int j = 0; j < parents.size(); j++) pp[j] = nodesHash.get(parents.get(j));
            Matrix covxx = sampleCovar.getSelection(pp, pp);
            Matrix covxy = sampleCovar.getSelection(pp, new int[]{idx});
            Matrix b = (covxx.inverse().times(covxy));

            for (int i = 0; i < parents.size(); i++) {
                int idx2 = nodesHash.get(parents.get(i));
                edgeCoef.set(idx2, idx, b.get(i, 0));
            }

            int[] all = concat(idx, pp);
            Matrix cov = sampleCovar.getSelection(all, all);

            Matrix bStar = bStar(b);
            double variance = (bStar.transpose().times(cov).times(bStar).get(0, 0));

//            double variance = sampleCovar.get(idx, idx);
//            variance -= covxy.getColumn(0).dotProduct(b.getColumn(0));
            errorCovar.set(idx, idx, variance);
        }
    }

    /**
     * Runs the estimator on the data and SemPm passed in through the
     * constructor. Returns the fml score of the resulting model.
     */
    public double score(Graph dag) {
        this.dag = dag;
        this.changedEdge = null;
        resetParameters(variables);
        return getBicScore();
    }

    @Override
    public double score(Edge edge) {
        this.changedEdge = edge;
        resetParameters(edge);
        return getBicScore();
    }

    public ICovarianceMatrix getCovMatrix() {
        return covMatrix;
    }

    /**
     * @return a string representation of the Sem.
     */
    public String toString() {
        return "\nSemEstimator";
    }

    /**
     * The value of the maximum likelihood function for the getModel the model
     * (Bollen 107). To optimize, this should be minimized.
     */
    private double getFml() {
        if (changedEdge != null) {
            resetParameters(changedEdge);
        } else {
            resetParameters(variables);
        }

        Matrix implCovarMeas = implCovarMeas();

        double logDetSigma = logDet(implCovarMeas);
        double traceSSigmaInv = traceABInv(sampleCovar, implCovarMeas);
        double logDetSample = logDetSample();
        int pPlusQ = getMeasuredNodes().size();

        double fml = logDetSigma + traceSSigmaInv - logDetSample - pPlusQ;

        if (Math.abs(fml) < 0) {
            fml = 0.0;
        }

        return fml;
    }

    private Matrix implCovarMeas() {
        computeImpliedCovar();
        return this.implCovar;
    }

    /**
     * @return BIC score, calculated as chisq - dof. This is equal to getFullBicScore() up to a constant.
     */
    private double getBicScore() {
        return getChiSquare() + penaltyDiscount * getNumFreeParams() * log(getSampleSize());
    }

    /**
     * @return the chi square value for the model.
     */
    private double getChiSquare() {
        return (getSampleSize() - 1) * getFml();
    }

    /**
     * @return the p-value for the model.
     */
    public double getPValue() {
        return 1.0 - ProbUtils.chisqCdf(getChiSquare(), getDof());
    }

    //============================PRIVATE METHODS==========================//

    /**
     * Adds semantic checks to the default deserialization method. This
     * method must have the standard signature for a readObject method, and
     * the body of the method must begin with "s.defaultReadObject();".
     * Other than that, any semantic checks can be specified and do not need
     * to stay the same from version to version. A readObject method of this
     * form may be added to any class, even if Tetrad sessions were
     * previously saved out using a version of the class that didn't include
     * it. (That's what the "s.defaultReadObject();" is for. See J. Bloch,
     * Effective Java, for help.
     */
    private void readObject(ObjectInputStream s)
            throws IOException, ClassNotFoundException {
        s.defaultReadObject();

        if (getCovMatrix() == null) {
            throw new NullPointerException();
        }

    }

    /**
     * Computes the implied covariance matrices of the Sem. There are two:
     * <code>implCovar </code> contains the covariances of all the variables and
     * <code>implCovarMeas</code> contains covariance for the measured variables
     * only.
     */
    private void computeImpliedCovar() {
        this.implCovar = MatrixUtils.impliedCovar(edgeCoef().transpose(), errCovar());
    }

    private Matrix errCovar() {
        return getErrorCovar();
    }

    private Matrix edgeCoef() {
        return getEdgeCoef();
    }

    private double logDet(Matrix matrix2D) {
        return log(matrix2D.det());
    }

    private double traceABInv(Matrix A, Matrix B) {

        // Note that at this point the sem and the sample covar MUST have the
        // same variables in the same order.
        try {

            Matrix product = A.times(B.inverse());

            double trace = product.trace();

            if (trace < -1e-8) {
                throw new IllegalArgumentException("Trace was negative: " + trace);
            }

            return trace;
        } catch (Exception e) {
            System.out.println(B);
            throw new RuntimeException(e);
        }
    }

    private double logDetSample() {
        if (logDetSample == 0.0 && this.sampleCovar != null) {
            double det = this.sampleCovar.det();
            logDetSample = log(det);
        }

        return logDetSample;
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }
}


