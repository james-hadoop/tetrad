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
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.FmlBicScorer;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.PermutationGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.StrictMath.abs;

/**
 * Implements the LiNGAM algorithm in Shimizu, Hoyer, Hyvarinen, and Kerminen, A linear nongaussian acyclic model for
 * causal discovery, JMLR 7 (2006). Largely follows the Matlab code.
 * <p>
 * We use FGES with knowledge of causal order for the pruning step.
 *
 * @author Joseph Ramsey
 */
public class Lingam {
    private double penaltyDiscount = 2;
    private double fastIcaA = 1.1;
    private int fastIcaMaxIter = 2000;
    private double fastIcaTolerance = 1e-6;
//    private final double pruneFactor = 2;

    private IcaAlgorithm algorithm = IcaAlgorithm.DEFLATION;
    private IcaFunction function = IcaFunction.LOGCOSH;

    //================================CONSTRUCTORS==========================//

    /**
     * Constructs a new LiNGAM algorithm with the given alpha level (used for pruning).
     */
    public Lingam() {
    }

    public Graph search(DataSet data) {
        for (int j = 0; j < data.getNumColumns(); j++) {
            for (int i = 0; i < data.getNumRows(); i++) {
                if (Double.isNaN(data.getDouble(i, j))) {
                    throw new IllegalArgumentException("Please remove or impute missing values.");
                }
            }
        }


        Matrix X = data.getDoubleData();
        X = DataUtils.centerData(X).transpose();
        FastIca fastIca = new FastIca(X, X.rows());
        fastIca.setVerbose(false);
        fastIca.setMaxIterations(fastIcaMaxIter);

        if (algorithm == IcaAlgorithm.DEFLATION) {
            fastIca.setAlgorithmType(FastIca.DEFLATION);
        } else {
            fastIca.setAlgorithmType(FastIca.PARALLEL);
        }

        fastIca.setTolerance(fastIcaTolerance);

        if (function == IcaFunction.LOGCOSH) {
            fastIca.setFunction(FastIca.LOGCOSH);
        } else {
            fastIca.setFunction(FastIca.EXP);
        }

        fastIca.setRowNorm(false);
        fastIca.setAlpha(fastIcaA);
        FastIca.IcaResult result11 = fastIca.findComponents();
        Matrix W = result11.getW();

        PermutationGenerator gen1 = new PermutationGenerator(W.columns());
        int[] perm1 = new int[0];
        double sum1 = Double.NEGATIVE_INFINITY;
        int[] choice1;

        while ((choice1 = gen1.next()) != null) {
            double sum = 0.0;

            for (int i = 0; i < W.columns(); i++) {
                final double wii = W.get(choice1[i], i);
                sum += abs(wii);
            }

            if (sum > sum1) {
                sum1 = sum;
                perm1 = Arrays.copyOf(choice1, choice1.length);
            }
        }

        int[] cols = new int[W.columns()];
        for (int i = 0; i < cols.length; i++) cols[i] = i;

        Matrix WTilde = W.getSelection(perm1, cols);

        Matrix WPrime = WTilde.copy();

        for (int i = 0; i < WPrime.rows(); i++) {
            for (int j = 0; j < WPrime.columns(); j++) {
                WPrime.assignRow(i, WTilde.getRow(i).scalarMult(1.0 / WTilde.get(i, i)));
            }
        }

        final int m = data.getNumColumns();
        Matrix BHat = Matrix.identity(m).minus(WPrime);

        PermutationGenerator gen2 = new PermutationGenerator(BHat.rows());
        int[] perm2 = new int[0];
        double sum2 = Double.NEGATIVE_INFINITY;
        int[] choice2;

        while ((choice2 = gen2.next()) != null) {
            double sum = 0.0;

            for (int i = 0; i < W.rows(); i++) {
                for (int j = 0; j < i; j++) {
                    final double c = BHat.get(choice2[i], choice2[j]);
                    sum += abs(c);
                }
            }

            if (sum > sum2) {
                sum2 = sum;
                perm2 = Arrays.copyOf(choice2, choice2.length);
            }
        }

        final SemBicScore score = new SemBicScore(new CovarianceMatrix(data));
        score.setPenaltyDiscount(penaltyDiscount);

        if (true) {
            Scorer scorer = new FmlBicScorer(data, 2);

            ForwardScoreSearch gss = new ForwardScoreSearch(scorer);

            List<Node> order = new ArrayList<>();
            final List<Node> variables = data.getVariables();
            for (int i = 0; i < variables.size(); i++) {
                order.add(variables.get(perm2[i]));
            }

            return gss.search(order);
        } else {
            Fges fges = new Fges(score);

            IKnowledge knowledge = new Knowledge2();
            final List<Node> variables = data.getVariables();

            for (int i = 0; i < variables.size(); i++) {
                knowledge.addToTier(i, variables.get(perm2[i]).getName());
            }

            fges.setKnowledge(knowledge);

            final Graph graph = fges.search();
            System.out.println("graph Returning this graph: " + graph);

//        Graph graph2 = new EdgeListGraph(variables);
//
//        Matrix bFinal = pruneEdgesByResampling(X, perm2, 20, 2);
//
//        for (int i = 0; i < bFinal.rows(); i++) {
//            for (int j = 0; j < bFinal.columns(); j++) {
//                if (i == j) continue;
//                if (abs(bFinal.get(i, j)) > 0) {
//                    graph2.addDirectedEdge(variables.get(j), variables.get(i));
//                }
//            }
//        }

            return graph;
        }
    }

//    // Inverse permutation, restores original order.
//    private int[] inverse(int[] k) {
//        int[] ki = new int[k.length];
//
//        for (int i = 0; i < ki.length; i++) {
//            ki[k[i]] = i;
//        }
//
//        return ki;
//    }

    //================================PUBLIC METHODS========================//

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
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

    public void setFunction(IcaFunction function) {
        this.function = function;
    }

    public void setAlgorithm(IcaAlgorithm algorithm) {
        this.algorithm = algorithm;
    }

//    /**
//     * This is the method used in Patrik's code.
//     */
//    private Matrix pruneEdgesByResampling(Matrix X, int[] k, int npieces, double pruneFactor) {
//        int piecesize = (int) Math.floor(X.columns() / (double) npieces);
//        int[] ki = inverse(k);
//
//        List<Matrix> bpieces = new ArrayList<>();
//
//        for (int p = 0; p < npieces; p++) {
//            Matrix Xp = X.getSelection(k, range((p) * piecesize, (p + 1) * piecesize - 1));
//
//            Xp = DataUtils.centerData(Xp.transpose()).transpose();
//            Matrix cov = Xp.times(Xp.transpose()).scalarMult(1.0 / Xp.columns());
//
//            Matrix invSqrt;
//
//            try {
//                invSqrt = cov.sqrt().inverse();
//            } catch (Exception e) {
//                continue;
//            }
//
//            QRDecomposition qr = new QRDecomposition(new BlockRealMatrix(invSqrt.toArray()));
//            Matrix L = new Matrix(qr.getR().getData()).transpose();
//
//            for (int s = 0; s < Xp.rows(); s++) {
//                for (int t = 0; t < Xp.rows(); t++) {
//                    L.set(s, t, L.get(s, t) / L.get(s, s));
//                }
//            }
//
//            if (checkNaN(L)) continue;
//
//            Matrix bnewest = Matrix.identity(Xp.rows()).minus(L);
//            bpieces.add(bnewest.getSelection(ki, ki));
//
//            System.out.println("piece = " + bnewest);
//        }
//
//        Matrix means = new Matrix(X.rows(), X.rows());
//        Matrix stds = new Matrix(X.rows(), X.rows());
//
//        Matrix BFinal = new Matrix(X.rows(), X.rows());
//
//        for (int i = 0; i < X.rows(); i++) {
//            for (int j = 0; j < X.rows(); j++) {
//                double[] b = new double[bpieces.size()];
//
//                for (int y = 0; y < bpieces.size(); y++) {
//                    b[y] = (bpieces.get(y).get(i, j));
//                }
//
//                means.set(i, j, StatUtils.mean(b));
//
//                stds.set(i, j, StatUtils.sd(b));
//
//                if (abs(means.get(i, j)) < pruneFactor * stds.get(i, j)) {
//                    BFinal.set(i, j, means.get(i, j));
//                }
//            }
//        }
//
//        System.out.println("means = " + means);
//        System.out.println("stds = " + stds);
//
//        System.out.println("BFinal = " + BFinal);
//
//        return BFinal;
//    }

//    private boolean checkNaN(Matrix r) {
//        for (int i = 0; i < r.rows(); i++) {
//            for (int j = 0; j < r.rows(); j++) {
//                if (Double.isNaN(r.get(i, j))) return true;
//            }
//        }
//        return false;
//    }

//    private int[] range(int i1, int i2) {
//        if (i2 < i1) throw new IllegalArgumentException("i2 must be >=  i2 " + i1 + ", " + i2);
//        int[] series = new int[i2 - i1 + 1];
//        for (int i = 0; i <= i2 - i1; i++) {
//            series[i] = i + i1;
//        }
//        return series;
//    }

    public enum IcaAlgorithm {DEFLATION, PARALLEL}
    public enum IcaFunction {LOGCOSH, EXP}
}

