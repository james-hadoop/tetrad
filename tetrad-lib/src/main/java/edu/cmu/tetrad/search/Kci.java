package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.Vector;
import edu.pitt.csb.mgm.EigenDecomposition;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well44497b;

import java.text.DecimalFormat;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.median;
import static java.lang.Math.*;

/***
 * Kernal Independence Test (KCI).
 *
 * Zhang, K., Peters, J., Janzing, D., and Schölkopf, B. (2012). Kernel-based conditional independence
 * test and application in causal discovery. arXiv preprint arXiv:1202.3775.
 *
 * Please see that paper, especially Theorem 4 and Proposition 5.
 *
 * Using optimal kernel bandwidths suggested by Bowman and Azzalini (1997):
 *
 * Bowman, A. W., and Azzalini, A. (1997). Applied smoothing techniques for data analysis: the kernel
 * approach with S-Plus illustrations (Vol. 18). OUP Oxford.
 *
 * @author Vineet Raghu on 7/3/2016
 * @author jdramsey refactoring 7/4/2018
 */
public class Kci implements IndependenceTest {

    // The supplied data set, standardized
    private final DataSet data;

    // Variables in data
    private final List<Node> variables;
    private final double[] h;
    private final double[][] _data;

    // The alpha level of the test.
    private double alpha;

    // A normal distribution with 1 degree of freedom.
    private final NormalDistribution normal = new NormalDistribution(new SynchronizedRandomGenerator(
            new Well44497b(193924L)), 0, 1);

    // True if the approximation algorithms should be used instead of Theorems 3 or 4.
    private boolean approximate;

    // Convenience map from nodes to their indices in the list of variables.
    private final Map<Node, Integer> hash;

    // Eigenvalues greater than this time the maximum will be kept.
    private double threshold = 0.01;

    // Number of bostraps for Theorem 4 and Proposition 5.
    private int numBootstraps = 5000;

    // Azzalini optimal kernel widths will be multiplied by this.
    private double widthMultiplier = 1.0;

    // Record of independence facts
    private final Map<IndependenceFact, IndependenceResult> facts = new HashMap<>();

    // Epsilon for Propositio 5.
    private double epsilon = 0.001;

    private boolean verbose;
    private IndependenceFact latestFact = null;

    /**
     * Constructor.
     *
     * @param data  The dataset to analyse. Must be continuous.
     * @param alpha The alpha value of the test.
     */
    public Kci(DataSet data, double alpha) {
        this.data = DataUtils.standardizeData(data);
        _data = data.getDoubleData().transpose().toArray();

        this.variables = data.getVariables();
        int n = this.data.getNumRows();

        this.hash = new HashMap<>();

        for (int i = 0; i < variables.size(); i++) {
            this.hash.put(variables.get(i), i);
        }

        double[][] dataCols = this.data.getDoubleData().transpose().toArray();
        this.h = new double[variables.size()];

        for (int i = 0; i < this.data.getNumColumns(); i++) {
            this.h[i] = h(variables.get(i), dataCols, hash);
        }

        Matrix Ones = new Matrix(n, 1);
        for (int j = 0; j < n; j++) Ones.set(j, 0, 1);

        this.alpha = alpha;

    }

    //====================================PUBLIC METHODS==================================//

    /**
     * Returns an Independence test for a subset of the variables.
     */
    public IndependenceTest indTestSubset(List<Node> vars) {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns true if the given independence question is judged true, false if not. The independence question is of the
     * form x _||_ y | z, z = [z1,...,zn], where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    public IndependenceResult checkIndependence(Node x, Node y, List<Node> z) {
        if (Thread.currentThread().isInterrupted()) {
            return new IndependenceResult(new IndependenceFact(x, y, z),
                    true, Double.NaN);
        }

//        Node _x = data.getVariable(x.getName());
//
//        if (_x == null) {
//            throw new NullPointerException("The x variable, " + x + " was not in the data.");
//        }
//
//        Node _y = data.getVariable(y.getName());
//
//        if (_y == null) {
//            throw new NullPointerException("The y variable, " + y + " was not in the data.");
//        }
//
//        List<Node> _z = new ArrayList<>();
//        for (Node v : z) {
//            Node variable = data.getVariable(v.getName());
//
//            if (variable == null) {
//                throw new NullPointerException("A conditioning variable, " + v + " was not in the data.");
//            }
//
//            _z.add(variable);
//        }

        List<Node> allVars = new ArrayList<>();
        allVars.add(x);
        allVars.add(y);
        allVars.addAll(z);

        IndependenceFact fact = new IndependenceFact(x, y, z);
        this.latestFact = fact;

        if (facts.containsKey(fact)) {
            IndependenceResult result = facts.get(fact);

            if (verbose) {
                double p = result.getPValue();

                if (result.independent()) {
                    TetradLogger.getInstance().forceLogMessage(fact + " INDEPENDENT p = " + p);
                } else {
                    TetradLogger.getInstance().forceLogMessage(fact + " dependent p = " + p);
                }
            }

            return new IndependenceResult(fact, result.independent(), result.getPValue());
        } else {
            List<Integer> rows = getRows(this.data);

            int[] _cols = new int[allVars.size()];

            for (int i = 0; i < allVars.size(); i++) {
                Node key = allVars.get(i);
                _cols[i] = this.hash.get(key);
            }

            int[] _rows = new int[rows.size()];
            for (int i = 0; i < rows.size(); i++) _rows[i] = rows.get(i);

            DataSet data = this.data.subsetRowsColumns(_rows, _cols);
//            data = DataUtils.standardizeData(data);
            double[][] _data = data.getDoubleData().transpose().toArray();

            Map<Node, Integer> hash = new HashMap<>();
            for (int i = 0; i < allVars.size(); i++) hash.put(allVars.get(i), i);

            int N = data.getNumRows();

            Matrix ones = new Matrix(N, 1);
            for (int j = 0; j < N; j++) ones.set(j, 0, 1);

            Matrix I = Matrix.identity(N);
            Matrix H = I.minus(ones.times(ones.transpose()).scalarMult(1.0 / N));

            double[] h = new double[allVars.size()];
            int count = 0;

            double sum = 0.0;
            for (int i = 0; i < allVars.size(); i++) {
                h[i] = this.h[this.hash.get(allVars.get(i))];

                if (h[i] != 0) {
                    sum += h[i];
                    count++;
                }
            }

            double avg = sum / count;

            for (int i = 0; i < h.length; i++) {
                if (h[i] == 0) h[i] = avg;
            }

            IndependenceResult result = facts.get(fact);

            if (this.facts.get(fact) != null) {
                return new IndependenceResult(fact, result.independent(), result.getPValue());
            } else {
                if (z.isEmpty()) {
                    result = isIndependentUnconditional(x, y, fact, _data, h, N, hash);
                } else {
                    result = isIndependentConditional(x, y, z, fact, _data, N, H, I, h, hash);
                }
            }

            if (verbose) {
                double p = result.getPValue();

                if (result.independent()) {
                    TetradLogger.getInstance().forceLogMessage(fact + " INDEPENDENT p = " + p);

                } else {
                    TetradLogger.getInstance().forceLogMessage(fact + " dependent p = " + p);
                }
            }

            return new IndependenceResult(fact, result.independent(), result.getPValue());
        }
    }

    /**
     * Returns the list of variables over which this independence checker is capable of determinining independence
     * relations.
     */
    public List<Node> getVariables() {
        return this.variables;
    }

    /**
     * Returns the variable by the given name.
     */
    public Node getVariable(String name) {
        return this.data.getVariable(name);
    }

    /**
     * Returns the list of names for the variables in getNodesInEvidence.
     */
    public List<String> getVariableNames() {
        return this.data.getVariableNames();
    }

    /**
     * Returns true if y is determined the variable in z.
     */
    public boolean determines(List<Node> z, Node y) {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the significance level of the independence test.
     *
     * @throws UnsupportedOperationException if there is no significance level.
     */
    public double getAlpha() {
        return this.alpha;
    }

    /**
     * Sets the significance level.
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    /**
     * @return a string representation of this test.
     */
    public String toString() {
        return "KCI, alpha = " + new DecimalFormat("0.0###").format(getAlpha());
    }


    /**
     * @return The data model for the independence test.
     */
    public DataModel getData() {
        return this.data;
    }


    public ICovarianceMatrix getCov() {
        throw new UnsupportedOperationException();
    }

    public List<DataSet> getDataSets() {
        LinkedList<DataSet> L = new LinkedList<>();
        L.add(this.data);
        return L;
    }

    public int getSampleSize() {
        return this.data.getNumRows();
    }

    public List<Matrix> getCovMatrices() {
        throw new UnsupportedOperationException();
    }

    public double getScore() {
        return getAlpha() - facts.get(latestFact).getPValue();
    }

    public boolean isApproximate() {
        return this.approximate;
    }

    public void setApproximate(boolean approximate) {
        this.approximate = approximate;
    }

    private double getWidthMultiplier() {
        return this.widthMultiplier;
    }

    public void setWidthMultiplier(double widthMultiplier) {
        if (widthMultiplier <= 0) throw new IllegalStateException("Width must be > 0");
        this.widthMultiplier = widthMultiplier;
    }

    private int getNumBootstraps() {
        return this.numBootstraps;
    }

    public void setNumBootstraps(int numBootstraps) {
        if (numBootstraps < 1) throw new IllegalArgumentException("Num bootstraps should be >= 1: " + numBootstraps);
        this.numBootstraps = numBootstraps;
    }


    public double getThreshold() {
        return this.threshold;
    }

    public void setThreshold(double threshold) {
        if (threshold < 0.0) throw new IllegalArgumentException("Threshold must be >= 0.0: " + threshold);
        this.threshold = threshold;
    }

    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    //====================================PRIVATE METHODS==================================//

    /**
     * KCI independence for the unconditional case. Uses Theorem 4 from the paper.
     *
     * @return true just in case independence holds.
     */
    private IndependenceResult isIndependentUnconditional(Node x, Node y, IndependenceFact fact, double[][] _data,
                                               double[] _h, int N,
                                               Map<Node, Integer> hash) {
        Matrix Ones = new Matrix(N, 1);
        for (int j = 0; j < N; j++) Ones.set(j, 0, 1);

        Matrix H = Matrix.identity(N).minus(Ones.times(Ones.transpose()).scalarMult(1.0 / N));

        Matrix kx = center(kernelMatrix(_data, x, null, getWidthMultiplier(), hash, N, _h), H);
        Matrix ky = center(kernelMatrix(_data, y, null, getWidthMultiplier(), hash, N, _h), H);

        try {
            if (isApproximate()) {
                double sta = kx.times(ky).trace();
                double mean_appr = kx.trace() * ky.trace() / N;
                double var_appr = 2 * kx.times(kx).trace() * ky.times(ky).trace() / (N * N);
                double k_appr = mean_appr * mean_appr / var_appr;
                double theta_appr = var_appr / mean_appr;
                double p = 1.0 - new GammaDistribution(k_appr, theta_appr).cumulativeProbability(sta);
                boolean indep = p > getAlpha();
                IndependenceResult result = new IndependenceResult(fact, indep, p);
                this.facts.put(fact, result);
                return result;
            } else {
                return theorem4(kx, ky, fact, N);
            }
        } catch (Exception e) {
            e.printStackTrace();
            IndependenceResult result = new IndependenceResult(fact, false, 0.0);
            this.facts.put(fact, result);
            return result;
        }
    }

    /**
     * KCI independence for the conditional case. Uses Theorem 3 from the paper.
     *
     * @return true just in case independence holds.
     */
    private IndependenceResult isIndependentConditional(Node x, Node y, List<Node> z, IndependenceFact fact, double[][] _data,
                                             int N, Matrix H, Matrix I, double[] _h, Map<Node, Integer> hash) {
        Matrix kx;
        Matrix ky;

        try {
            Matrix KXZ = center(kernelMatrix(_data, x, z, getWidthMultiplier(), hash, N, _h), H);
            Matrix Ky = center(kernelMatrix(_data, y, null, getWidthMultiplier(), hash, N, _h), H);
            Matrix KZ = center(kernelMatrix(_data, null, z, getWidthMultiplier(), hash, N, _h), H);

            Matrix Rz = (KZ.plus(I.scalarMult(this.epsilon)).inverse().scalarMult(this.epsilon));

            kx = symmetrized(Rz.times(KXZ).times(Rz.transpose()));
            ky = symmetrized(Rz.times(Ky).times(Rz.transpose()));

            return proposition5(kx, ky, fact, N);
        } catch (Exception e) {
            e.printStackTrace();
            boolean indep = false;
            IndependenceResult result = new IndependenceResult(fact, indep, 0.0);
            this.facts.put(fact, result);
            return result;
        }
    }

    private IndependenceResult theorem4(Matrix kx, Matrix ky, IndependenceFact fact, int N) {

        double T = (1.0 / N) * (kx.times(ky).trace());

        // Eigen decomposition of kx and ky.
        Eigendecomposition eigendecompositionx = new Eigendecomposition(kx).invoke();
        List<Double> evx = eigendecompositionx.getTopEigenvalues();

        Eigendecomposition eigendecompositiony = new Eigendecomposition(ky).invoke();
        List<Double> evy = eigendecompositiony.getTopEigenvalues();

        // Calculate formula (9).
        int sum = 0;

        for (int j = 0; j < getNumBootstraps(); j++) {
            double tui = 0.0;

            for (double lambdax : evx) {
                for (double lambday : evy) {
                    tui += lambdax * lambday * getChisqSample();
                }
            }

            tui /= N * N;

            if (tui > T) sum++;
        }

        // Calculate p.
        double p = sum / (double) getNumBootstraps();
        boolean indep = p > getAlpha();
        IndependenceResult result = new IndependenceResult(fact, indep, p);
        this.facts.put(fact, result);
        return result;
    }

    private IndependenceResult proposition5(Matrix kx, Matrix ky, IndependenceFact fact, int N) {
        double T = (1.0 / N) * kx.times(ky).trace();

        Eigendecomposition eigendecompositionx = new Eigendecomposition(kx).invoke();
        Matrix vx = eigendecompositionx.getV();
        Matrix dx = eigendecompositionx.getD();

        Eigendecomposition eigendecompositiony = new Eigendecomposition(ky).invoke();
        Matrix vy = eigendecompositiony.getV();
        Matrix dy = eigendecompositiony.getD();

        // VD
        Matrix vdx = vx.times(dx);
        Matrix vdy = vy.times(dy);

        int prod = vx.columns() * vy.columns();
        Matrix UU = new Matrix(N, prod);

        // stack
        for (int i = 0; i < vx.columns(); i++) {
            for (int j = 0; j < vy.columns(); j++) {
                for (int k = 0; k < N; k++) {
                    UU.set(k, i * dy.columns() + j, vdx.get(k, i) * vdy.get(k, j));
                }
            }
        }

        Matrix uuprod = prod > N ? UU.times(UU.transpose()) : UU.transpose().times(UU);

        if (isApproximate()) {
            double sta = kx.times(ky).trace();
            double mean_appr = uuprod.trace();
            double var_appr = 2.0 * uuprod.times(uuprod).trace();
            double k_appr = mean_appr * mean_appr / var_appr;
            double theta_appr = var_appr / mean_appr;
            double p = 1.0 - new GammaDistribution(k_appr, theta_appr).cumulativeProbability(sta);
            boolean indep = p > getAlpha();
            IndependenceResult result = new IndependenceResult(fact, indep, p);
            this.facts.put(fact, result);
            return result;
        } else {

            // Get top eigenvalues of that.
            Eigendecomposition eigendecompositionu = new Eigendecomposition(uuprod).invoke();
            List<Double> eigenu = eigendecompositionu.getTopEigenvalues();

            // Calculate formulas (13) and (14).
            int sum = 0;

            for (int j = 0; j < getNumBootstraps(); j++) {
                double s = 0.0;

                for (double lambdaStar : eigenu) {
                    s += lambdaStar * getChisqSample();
                }

                s *= 1.0 / N;

                if (s > T) sum++;
            }

            double p = sum / (double) getNumBootstraps();
            boolean indep = p > getAlpha();
            IndependenceResult result = new IndependenceResult(fact, indep, p);
            this.facts.put(fact, result);
            return result;
        }
    }

    private List<Integer> series(int size) {
        List<Integer> series = new ArrayList<>();
        for (int i = 0; i < size; i++) series.add(i);
        return series;
    }

    private Matrix center(Matrix K, Matrix H) {
        return H.times(K).times(H);
    }

    private double getChisqSample() {
        double z = this.normal.sample();
        return z * z;
    }

    // Optimal bandwidth qsuggested by Bowman and Azzalini (1997) q.31,
    // using MAD.
    private double h(Node x, double[][] _data, Map<Node, Integer> hash) {
        double[] xCol = _data[hash.get(x)];
        double[] g = new double[xCol.length];
        double median = median(xCol);
        for (int j = 0; j < xCol.length; j++) g[j] = abs(xCol[j] - median);
        double mad = median(g);
        return (1.4826 * mad) * pow((4.0 / 3.0) / xCol.length, 0.2);
    }

    private List<Integer> getTopIndices(List<Double> prod, List<Integer> allIndices, double threshold) {
        double maxEig = prod.get(allIndices.get(0));

        List<Integer> indices = new ArrayList<>();

        for (int i : allIndices) {
            if (prod.get(i) > maxEig * threshold) {
                indices.add(i);
            }
        }

        return indices;
    }

    private Matrix symmetrized(Matrix kx) {
        return (kx.plus(kx.transpose())).scalarMult(0.5);
    }

    private Matrix kernelMatrix(double[][] _data, Node x, List<Node> z, double widthMultiplier,
                                Map<Node, Integer> hash,
                                int N, double[] _h) {

        List<Integer> _z = new ArrayList<>();

        if (x != null) {
            _z.add(hash.get(x));
        }

        if (z != null) {
            for (Node z2 : z) {
                _z.add(hash.get(z2));
            }
        }

        double h = getH(_z, _h);

        Matrix result = new Matrix(N, N);

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double d = distance(_data, _z, i, j);
                double k = kernelGaussian(d, widthMultiplier * h);
                result.set(i, j, k);
                result.set(j, i, k);
            }
        }

        double k = kernelGaussian(0, widthMultiplier * h);

        for (int i = 0; i < N; i++) {
            result.set(i, i, k);
        }

        return result;
    }

    private double getH(List<Integer> _z, double[] _h) {
        double h = 0;

        for (int c : _z) {
            if (_h[c] > h) {
                h = _h[c];
            }
        }

        h *= sqrt(_z.size());
        return h;
    }

    private double kernelGaussian(double z, double width) {
        z /= width;
        return exp(-z * z);
    }

    // Euclidean distance.
    private double distance(double[][] data, List<Integer> cols, int i, int j) {
        double sum = 0.0;

        for (int col : cols) {
            double d = (data[col][i] - data[col][j]) / 2;

            if (!Double.isNaN(d)) {
                sum += d * d;
            }
        }

        return sqrt(sum);
    }

    @Override
    public boolean isVerbose() {
        return this.verbose;
    }

    @Override
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    private class Eigendecomposition {
        private final Matrix k;
        private Matrix D;
        private Matrix V;
        private List<Double> topEigenvalues;

        public Eigendecomposition(Matrix k) {
            if (k.rows() == 0 || k.columns() == 0) {
                throw new IllegalArgumentException("Empty matrix to decompose. Please don't do that to me.");
            }

            this.k = k;
        }

        public Matrix getD() {
            return this.D;
        }

        public Matrix getV() {
            return this.V;
        }

        public List<Double> getTopEigenvalues() {
            return this.topEigenvalues;
        }

        public Eigendecomposition invoke() {
            List<Integer> topIndices;

            EigenDecomposition ed = new EigenDecomposition(new BlockRealMatrix(this.k.toArray()));

            double[] arr = ed.getRealEigenvalues();

            List<Double> evxAll = new ArrayList<>();
            for (double v : arr) evxAll.add(v);

            List<Integer> indx = series(evxAll.size()); // 1 2 3...
            topIndices = getTopIndices(evxAll, indx, getThreshold());

            this.D = new Matrix(topIndices.size(), topIndices.size());

            for (int i = 0; i < topIndices.size(); i++) {
                this.D.set(i, i, sqrt(evxAll.get(topIndices.get(i))));
            }

            this.topEigenvalues = new ArrayList<>();

            for (int t : topIndices) {
                getTopEigenvalues().add(evxAll.get(t));
            }

            this.V = new Matrix(ed.getEigenvector(0).getDimension(), topIndices.size());

            for (int i = 0; i < topIndices.size(); i++) {
                RealVector t = ed.getEigenvector(topIndices.get(i));
                this.V.assignColumn(i, new Vector(t.toArray()));
            }
//            } else {
//                SingularValueDecomposition svd = new SingularValueDecomposition(new BlockRealMatrix(k.toArray()));
//
//                List<Double> evxAll = asList(svd.getSingularValues());
//
//                List<Integer> indx = series(evxAll.size()); // 1 2 3...
//                topIndices = getTopIndices(evxAll, indx, getThreshold());
//
//                D = new Matrix(topIndices.size(), topIndices.size());
//
//                for (int i = 0; i < topIndices.size(); i++) {
//                    D.set(i, i, Math.sqrt(evxAll.get(topIndices.get(i))));
//                }
//
//                RealMatrix V0 = svd.getV();
//
//                V = new Matrix(V0.getRowDimension(), topIndices.size());
//
//                for (int i = 0; i < V.columns(); i++) {
//                    double[] t = V0.getColumn(topIndices.get(i));
//                    V.assignColumn(i, new Vector(t));
//                }
//
//                topEigenvalues = new ArrayList<>();
//
//                for (int t : topIndices) {
//                    getTopEigenvalues().add(evxAll.get(t));
//                }
//
//            }

            return this;
        }
    }


    private List<Integer> getRows(DataSet dataSet) {
        List<Integer> rows = new ArrayList<>();

        for (int k = 0; k < dataSet.getNumRows(); k++) {
            rows.add(k);
        }

        return rows;
    }
}