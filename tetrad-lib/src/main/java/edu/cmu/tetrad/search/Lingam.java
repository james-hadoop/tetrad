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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.PermutationGenerator;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.linear.QRDecomposition;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static edu.cmu.tetrad.util.StatUtils.*;
import static edu.cmu.tetrad.util.TetradMatrix.eye;
import static java.lang.StrictMath.abs;
import static java.lang.StrictMath.scalb;

/**
 * Implements the LiNGAM algorithm in Shimizu, Hoyer, Hyvarinen, and Kerminen, A linear nongaussian acyclic model for
 * causal discovery, JMLR 7 (2006). Largely follows the Matlab code.
 * <p>
 * We use FGES with knowledge of causal order for the pruning step.
 *
 * @author Joseph Ramsey
 */
public class Lingam {
    private double fastIcaA = 2;
    private int fastIcaMaxIter = 200;
    private double fastIcaTolerance = 1e-6;
    private double pruneFactor = 1;

    //================================CONSTRUCTORS==========================//

    /**
     * Constructs a new LiNGAM algorithm with the given alpha level (used for pruning).
     */
    public Lingam() {
    }

    //================================PUBLIC METHODS========================//

    public Graph search(DataSet data) {
        int parallel = FastIca.PARALLEL;

        TetradMatrix X = DataUtils.center(data).getDoubleData().transpose();
        FastIca fastIca = new FastIca(X, X.rows());
        fastIca.setVerbose(false);
        fastIca.setMaxIterations(fastIcaMaxIter);
        fastIca.setAlgorithmType(parallel);
        fastIca.setTolerance(fastIcaTolerance);
        fastIca.setFunction(FastIca.LOGCOSH);
        fastIca.setRowNorm(false);
        fastIca.setAlpha(fastIcaA);
        FastIca.IcaResult result = fastIca.findComponents();

        TetradMatrix W = result.getW();
        TetradMatrix A = W.inverse();
        TetradMatrix S = result.getS();

        // X = AS, S consists of independent rows.

        TetradMatrix sCov = S.times(S.transpose()).scalarMult(1.0 / S.columns());
//        System.out.println(sCov);

        System.out.println("cov(S) = I? " + sCov.equals(eye(sCov.rows()), 1e-3));

//        System.out.println("X = " + X);
//        System.out.println("AS = " + A.times(S));

        System.out.println("X = AS? " + X.equals(A.times(S), 1e-3));

        PermutationGenerator gen1 = new PermutationGenerator(W.rows());
        int[] r = null;
        double sum1 = Double.POSITIVE_INFINITY;
        int[] p;

        while ((p = gen1.next()) != null) {
            TetradMatrix Wt2 = W.getSelection(p, p);

            double sum = 0.0;

            for (int i = 0; i < Wt2.rows(); i++) {
                sum += 1.0 / abs(Wt2.get(i, i));
            }

            if (sum < sum1) {
                sum1 = sum;
                r = copy(p);
            }
        }

        if (r == null) throw new IllegalStateException();

        TetradMatrix Wt = W.getSelection(r, r);
        TetradMatrix Wtp = Wt.copy();

        for (int i = 0; i < Wtp.rows(); i++) {
            for (int j = 0; j < Wtp.columns(); j++) {
                Wtp.set(i, j, Wtp.get(i, j) / Wtp.get(i, i));
            }
        }

        final int m = X.rows();
        TetradMatrix BHat = eye(m).minus(Wtp);

        PermutationGenerator gen2 = new PermutationGenerator(BHat.rows());
        int[] k = null;
        double sum2 = Double.POSITIVE_INFINITY;
        int[] q;

        while ((q = gen2.next()) != null) {
            TetradMatrix Bt = BHat.getSelection(q, q);

            double sum = 0.0;

            for (int i = 0; i < Bt.rows(); i++) {
                for (int j = i; j < Bt.columns(); j++) {
                    final double c = Bt.get(i, j);
                    sum += abs(c);
                }
            }

            if (sum < sum2) {
                sum2 = sum;
                k = copy(q);
            }
        }

        if (k == null) throw new IllegalStateException();

        final List<Node> variables = data.getVariables();

        System.out.println("Variables = " + variables);

        FasStable fas = new FasStable(new IndTestFisherZ(data, 0.01));
        final Graph graph = fas.search();

        List<Node> order = new ArrayList<>();
        for (int i = 0; i < variables.size(); i++) order.add(variables.get(k[i]));

        System.out.println("order = " + order);

        for (Edge edge : new ArrayList<>(graph.getEdges())) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            if (order.indexOf(y) > order.indexOf(x)) {
                graph.setEndpoint(x, y, Endpoint.ARROW);
            } else {
                graph.setEndpoint(y, x, Endpoint.ARROW);
            }
        }

////        int[] k = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//
//        Graph graph = new EdgeListGraph(variables);
//
//        for (int i = 0; i < k.length; i++) {
//            System.out.println("k[" + i + "] = " + k[i] + " " + variables.get(k[i]));
//        }
//
//        TetradMatrix bFinal = pruneEdgesByResampling(X, k);
//
//        for (int i = 0; i < bFinal.rows(); i++) {
//            for (int j = 0; j < bFinal.columns(); j++) {
//                if (i == j) continue;
//                if (bFinal.get(i, j) != 0.0) {
//                    graph.addDirectedEdge(variables.get(j), variables.get(i));
//                }
//            }
//        }

        return graph;
    }

    // Inverse permutation, restores original order.
    private int[] inverse(int[] k) {
        int[] ki = new int[k.length];

        for (int i = 0; i < ki.length; i++) {
            ki[k[i]] = i;
        }

//        int[] kj = new int[k.length];
//
//        for (int i = 0; i < k.length; i++) {
//            kj[i] = ki[k[i]];
//        }

        return ki;
    }

    public void setFastIcaA(double fastIcaA) {
        this.fastIcaA = fastIcaA;
    }

    public void setFastMaxIter(int maxIter) {
        this.fastIcaMaxIter = maxIter;
    }

    public void setFastIcaTolerance(double tolerance) {
        this.fastIcaTolerance = tolerance;
    }

    public void setPruneFactor(double pruneFactor) {
        this.pruneFactor = pruneFactor;
    }

    public void setFastIcaMaxIter(int fastIcaMaxIter) {
        this.fastIcaMaxIter = fastIcaMaxIter;
    }

    //============================PRIVATE=========================//

    private int[] range(int n) {
        int[] cols = new int[n];
        for (int i = 0; i < n; i++) cols[i] = i;
        return cols;
    }

    private int[] copy(int[] choice1) {
        return Arrays.copyOf(choice1, choice1.length);
    }

    /**
     * This is the method used in Patrik's code.
     */
    private TetradMatrix pruneEdgesByResampling(TetradMatrix X, int[] k) {
        int npieces = 20;
        int piecesize = (int) Math.floor(X.columns() / (double) npieces);
        int[] ki = inverse(k);

        List<TetradMatrix> bpieces = new ArrayList<>();

        for (int p = 0; p < npieces; p++) {
            TetradMatrix Xp = X.getSelection(k, range((p) * piecesize, (p + 1) * piecesize - 1));

//            System.out.println("Xp = " + Xp);

            Xp = DataUtils.centerData(Xp);
            TetradMatrix cov = Xp.times(Xp.transpose()).scalarMult(1.0 / Xp.columns());

            TetradMatrix invSqrt;

            try {
                invSqrt = cov.sqrt().inverse();
            } catch (Exception e) {
                continue;
            }

            QRDecomposition qr = new QRDecomposition(invSqrt.getRealMatrix());
            TetradMatrix R = new TetradMatrix(qr.getR());//.transpose();

            for (int s = 0; s < Xp.rows(); s++) {
                for (int t = 0; t < Xp.rows(); t++) {
                    R.set(s, t, R.get(s, t) / R.get(s, s));
                }
            }

            if (checkNaN(R)) continue;

            TetradMatrix bnewest = eye(Xp.rows()).minus(R);
            bpieces.add(bnewest.getSelection(ki, ki));

            System.out.println("piece = " + bnewest);
        }

        TetradMatrix means = new TetradMatrix(X.rows(), X.rows());
        TetradMatrix stds = new TetradMatrix(X.rows(), X.rows());

        TetradMatrix BFinal = new TetradMatrix(X.rows(), X.rows());

        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.rows(); j++) {
                double[] b = new double[bpieces.size()];

                for (int y = 0; y < bpieces.size(); y++) {
                    b[y] = abs(bpieces.get(y).get(i, j));
                }

                means.set(i, j, mean(b));

//                if (means.get(i, j) != 0) {
                stds.set(i, j, sd(b));

                if (abs(means.get(i, j)) < pruneFactor * stds.get(i, j)) {
                    BFinal.set(i, j, means.get(i, j));
                }
//                }
            }
        }

        System.out.println("means = " + means);
        System.out.println("stds = " + stds);

        System.out.println("BFinal = " + BFinal);

        return BFinal;
    }

    private boolean checkNaN(TetradMatrix r) {
        for (int i = 0; i < r.rows(); i++) {
            for (int j = 0; j < r.rows(); j++) {
                if (Double.isNaN(r.get(i, j))) return true;
            }
        }
        return false;
    }

    private int[] range(int i1, int i2) {
        if (i2 < i1) throw new IllegalArgumentException("i2 must be >=  i2 " + i1 + ", " + i2);
        int[] series = new int[i2 - i1 + 1];
        for (int i = 0; i <= i2 - i1; i++) {
            series[i] = i + i1;
        }
        return series;
    }

}

