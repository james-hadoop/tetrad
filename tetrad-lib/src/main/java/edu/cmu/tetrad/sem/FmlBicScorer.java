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
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.util.*;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private DataSet dataSet = null;
    private Graph dag = null;
    private Matrix implCovarMeasC;
    private double logDetSample;
    private double fml = Double.NaN;
    private double penaltyDiscount = 2.;
    private Map<Node, Integer> nodesHash;


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

        int m = this.getVariables().size();
        this.edgeCoef = new Matrix(m, m);
        this.errorCovar = new Matrix(m, m);
        this.sampleCovar = covMatrix.getMatrix();

        this.nodesHash = new HashMap<>();
        for (int i = 0; i < variables.size(); i++) nodesHash.put(variables.get(i), i);
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

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static Scorer serializableInstance() {
        return new FmlBicScorer(CovarianceMatrix.serializableInstance());
    }

    //==============================PUBLIC METHODS=========================//

    /**
     * Runs the estimator on the data and SemPm passed in through the
     * constructor. Returns the fml score of the resulting model.
     */
    public double score(Graph dag) {
        for (Node node : dag.getNodes()) {
            int i1 = indexOf(node);
            getErrorCovar().set(i1, i1, 0);
            for (int _j = 0; _j < getVariables().size(); _j++) {
                getEdgeCoef().set(_j, i1, 0);
            }

            if (node.getNodeType() != NodeType.MEASURED) {
                continue;
            }

            int idx = indexOf(node);
            List<Node> parents = dag.getParents(node);

            for (int i = 0; i < parents.size(); i++) {
                Node nextParent = parents.get(i);
                if (nextParent.getNodeType() == NodeType.ERROR) {
                    parents.remove(nextParent);
                    break;
                }
            }

            double variance = sampleCovar.get(idx, idx);

            if (parents.size() > 0) {
                Vector nodeParentsCov = new Vector(parents.size());
                Matrix parentsCov = new Matrix(parents.size(), parents.size());

                for (int i = 0; i < parents.size(); i++) {
                    int idx2 = indexOf(parents.get(i));
                    nodeParentsCov.set(i, sampleCovar.get(idx, idx2));

                    for (int j = i; j < parents.size(); j++) {
                        int idx3 = indexOf(parents.get(j));
                        parentsCov.set(i, j, sampleCovar.get(idx2, idx3));
                        parentsCov.set(j, i, sampleCovar.get(idx3, idx2));
                    }
                }

                Vector edges = parentsCov.inverse().times(nodeParentsCov);

                for (int i = 0; i < edges.size(); i++) {
                    int idx2 = indexOf(parents.get(i));
                    edgeCoef.set(idx2, indexOf(node), edges.get(i));
                }

                variance -= nodeParentsCov.dotProduct(edges);
            }

            errorCovar.set(i1, i1, variance);
        }


        this.dag = dag;
        this.fml = Double.NaN;

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
    public double getFml() {
        if (!Double.isNaN(this.fml)) {
            return this.fml;
        }

        Matrix implCovarMeas; // Do this once.

        try {
            implCovarMeas = implCovarMeas();
        } catch (Exception e) {
            e.printStackTrace();
            return Double.NaN;
        }

        Matrix sampleCovar = this.sampleCovar;

        double logDetSigma = logDet(implCovarMeas);
        double traceSSigmaInv = traceABInv(sampleCovar, implCovarMeas);
        double logDetSample = logDetSample();
        int pPlusQ = getMeasuredNodes().size();

        double fml = logDetSigma + traceSSigmaInv - logDetSample - pPlusQ;

        if (Math.abs(fml) < 0) {
            fml = 0.0;
        }

        this.fml = fml;
        return fml;
    }

    private Matrix implCovarMeas() {
        computeImpliedCovar();
        return this.implCovarMeasC;
    }

    /**
     * @return BIC score, calculated as chisq - dof. This is equal to getFullBicScore() up to a constant.
     */
    public double getBicScore() {
        int dof = getDof();
        return getChiSquare() - penaltyDiscount * dof * Math.log(getSampleSize());
    }

    /**
     * @return the chi square value for the model.
     */
    public double getChiSquare() {
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

        // Note. Since the sizes of the temp matrices in this calculation
        // never change, we ought to be able to reuse them.
        Matrix implCovarC = MatrixUtils.impliedCovar(edgeCoef().transpose(), errCovar());

        // Submatrix of implied covar for measured vars only.
        int size = getMeasuredNodes().size();
        this.implCovarMeasC = new Matrix(size, size);

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                this.implCovarMeasC.set(i, j, implCovarC.get(i, j));
            }
        }
    }

    private Matrix errCovar() {
        return getErrorCovar();
    }

    private Matrix edgeCoef() {
        return getEdgeCoef();
    }

    private double logDet(Matrix matrix2D) {
        return Math.log(matrix2D.det());
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
            logDetSample = Math.log(det);
        }

        return logDetSample;
    }

    private int indexOf(Node node) {
        return nodesHash.get(node);
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }
}


