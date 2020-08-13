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

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Edges;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.PermutationGenerator;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;
import org.apache.commons.math3.linear.QRDecomposition;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
    private double fastIcaA = .9;
    private int fastIcaMaxIter = 100;
    private double fastIcaTolerance = 1e-6;

    private double pruneFactor = 3;

    //================================CONSTRUCTORS==========================//

    /**
     * Constructs a new LiNGAM algorithm with the given alpha level (used for pruning).
     */
    public Lingam() {
    }

    public Graph search(DataSet data) {
        fastIcaMaxIter = 2000;
        int parallel = FastIca.DEFLATION;
        double fastIcaTolerance = 0.0001;

        data = DataUtils.center(data);

        TetradMatrix X = data.getDoubleData().transpose();
        FastIca fastIca = new FastIca(X, X.rows());
        fastIca.setVerbose(false);
        fastIca.setMaxIterations(fastIcaMaxIter);
        fastIca.setAlgorithmType(parallel);
        fastIca.setTolerance(fastIcaTolerance);
        fastIca.setFunction(FastIca.EXP);
        fastIca.setRowNorm(false);
        fastIca.setAlpha(fastIcaA);
        FastIca.IcaResult result11 = fastIca.findComponents();

        TetradMatrix A = result11.getA();
//        TetradMatrix S = result11.getS();

//        System.out.println("X = " + X);
//        System.out.println("AS = " + A.times(S));
//        System.out.println("S = " + S);

        TetradMatrix W = A.inverse();

//        System.out.println("W = " + W);

        PermutationGenerator gen1 = new PermutationGenerator(W.rows());
        int[] perm1 = new int[0];
        double sum1 = Double.POSITIVE_INFINITY;
        int[] choice1;

        while ((choice1 = gen1.next()) != null) {
            double sum = 0.0;

            for (int i = 0; i < W.rows(); i++) {
                final double wii = W.get(choice1[i], i);
                sum += abs(1.0 / wii);
            }

            if (sum < sum1) {
                sum1 = sum;
                perm1 = Arrays.copyOf(choice1, choice1.length);
            }
        }

        int[] cols = new int[W.columns()];
        for (int i = 0; i < cols.length; i++) cols[i] = i;

        // zeroless diagonal.
        TetradMatrix WTilde = W.getSelection(perm1, cols);

        System.out.println("WTilde = " + WTilde);

//        TetradMatrix WPrime = WTilde.copy();

        for (int i = 0; i < WTilde.rows(); i++) {
            WTilde.assignRow(i, WTilde.getRow(i).scalarMult(1.0 / WTilde.get(i, i)));
        }

        // scaled
        System.out.println("WPrime = " + WTilde);

        final int m = data.getNumColumns();
        TetradMatrix BHat = TetradMatrix.identity(m).minus(WTilde);

        // beta weights
        System.out.println("BHat = " + BHat);


        PermutationGenerator gen2 = new PermutationGenerator(BHat.rows());
        int[] k = new int[0];
        double sum2 = Double.POSITIVE_INFINITY;
        int[] choice2;

        while ((choice2 = gen2.next()) != null) {
            double sum = 0.0;

            for (int j = 0; j < W.rows(); j++) {
                for (int i = 0; i <= j; i++) {
                    final double c = BHat.get(choice2[i], choice2[j]);
                    sum += c * c;
                }
            }

            if (sum < sum2) {
                sum2 = sum;
                k = Arrays.copyOf(choice2, choice2.length);
            }
        }

//        TetradMatrix BTilde = BHat.getSelection(k, k);
//
//        System.out.println("BTilde = " + BTilde);

        final SemBicScore score = new SemBicScore(new CovarianceMatrix(data));
        score.setPenaltyDiscount(penaltyDiscount);

//        Fges fges = new Fges(score);

//        IKnowledge knowledge = new Knowledge2();
        final List<Node> variables = data.getVariables();
//
        for (int i = 0; i < variables.size(); i++) {
//            knowledge.addToTier(i, variables.get(k[i]).getName());
            System.out.println("i = " + (i + 1) + " " + variables.get(k[i]));
        }

//        FasStable pc = new FasStable(new IndTestFisherZ(data, 0.01));
//        final Graph graph = pc.search();
//
//        List<Node> order = new ArrayList<>();
//
//        for (int i = 0; i < variables.size(); i++) {
//            order.add(variables.get(k[i]));
//        }
//
//        for (Edge edge : new ArrayList<>(graph.getEdges())) {
//            Node x = edge.getNode1();
//            Node y = edge.getNode2();
//
//            graph.removeEdge(edge);
//
//            if (order.indexOf(y) > order.indexOf(x)) {
//                graph.addDirectedEdge(x, y);
//            } else if (order.indexOf(y) < order.indexOf(x)) {
//                graph.addDirectedEdge(y, x);
//            }
//        }

        TetradMatrix bFinal = pruneEdgesByResampling(X, k);

        System.out.println("bFinal again " + bFinal);

        Graph graph = new EdgeListGraph(variables);

        for (int j = 0; j < X.rows(); j++) {
            for (int i = 0; i < X.rows(); i++) {
                if (bFinal.get(i, j) != 0.0) {
                    graph.addDirectedEdge(variables.get(j), variables.get(i));
                }
            }
        }

//        for (Node x : graph.getNodes()) {
//            for (Node y : graph.getNodes()) {
//                if (graph.getEdges(x, y).size() > 1) {
//                    graph.removeEdges(x, y);
//                }
//            }
//        }

        System.out.println("graph Returning this graph: " + graph);
        return graph;
    }

    /**
     * This is the method used in Patrik's code.
     */
    public TetradMatrix pruneEdgesByResampling(TetradMatrix X, int[] k) {
        if (k.length != X.rows()) {
            throw new IllegalArgumentException("Expecting a permutation.");
        }

        int npieces = 20;
        int piecesize = (int) Math.floor(X.columns() / (double) npieces);

        List<TetradMatrix> bpieces = new ArrayList<>();
//        List<Vector> diststdpieces = new ArrayList<Vector>();
//        List<Vector> cpieces = new ArrayList<Vector>();

        PIECES:
        for (int p = 0; p < npieces; p++) {

//          % Select subset of data, and permute the variables to the causal order
//          Xp = X(k,((p-1)*piecesize+1):(p*piecesize));
            TetradMatrix Xp = X.getSelection(k, range((p) * piecesize, (p + 1) * piecesize - 1));
//
//            List<Integer> allCols = new ArrayList<>();
//            for (int i = 0; i < X.columns(); i++) allCols.add(i);
//            Collections.shuffle(allCols);
//            int[] myCols = new int[(int) (allCols.size() * 0.7)];
//            for (int i = 0; i < myCols.length; i++) myCols[i] = allCols.get(i);
//            int[] allRows = new int[X.rows()];
//            for (int i = 0; i < allRows.length; i++) allRows[i] = i;
//
//            TetradMatrix Xp = X.getSelection(allRows, myCols);
//

            Xp = DataUtils.centerData(Xp);
            TetradMatrix cov = Xp.times(Xp.transpose()).scalarMult(1.0 / Xp.columns());

//          % Do QL decomposition on the inverse square root of cov
//          [Q,R] = tridecomp(cov^(-0.5),'ql');

            TetradMatrix invSqrt = cov.sqrt().inverse();

            QRDecomposition qr = new QRDecomposition(invSqrt.getRealMatrix());
            TetradMatrix R = new TetradMatrix(qr.getR());

//          % The estimated disturbance-stds are one over the abs of the diag of R
//          newestdisturbancestd = 1./diag(abs(R));

            TetradVector newestdisturbancestd = new TetradVector(X.rows());

            for (int t = 0; t < X.rows(); t++) {
                if (R.get(t, t) == 0) continue PIECES;
                newestdisturbancestd.set(t, 1.0 / Math.abs(R.get(t, t)));
            }

//          % Normalize rows of R to unit diagonal
//          R = R./(diag(R)*ones(1,dims));
//
            for (int s = 0; s < Xp.rows(); s++) {
                for (int t = 0; t < Xp.rows(); t++) {
                    if (R.get(s, s) == 0) continue PIECES;
                    R.set(s, t, R.get(s, t) / R.get(s, s));
                }
            }

//          % Calculate corresponding B
//          bnewest = eye(dims)-R;
            TetradMatrix bnewest = TetradMatrix.identity(X.rows()).minus(R);

//          % Also calculate constants
//          cnewest = R*Xpm;

//            Vector cnewest = new DenseVector(rows);
//            cnewest = R.mult(new DenseVector(Xpm), cnewest);

//          % Permute back to original variable order
//          ik = iperm(k);
//          bnewest = bnewest(ik, ik);
//          newestdisturbancestd = newestdisturbancestd(ik);
//          cnewest = cnewest(ik);

//            int[] ik = iperm(k);

//            System.out.println("ik = " + Arrays.toString(ik));

//            bnewest = bnewest.getSelection(ik, ik);
//            newestdisturbancestd = Matrices.getSubVector(newestdisturbancestd, ik);
//            cnewest = Matrices.getSubVector(cnewest, ik);

//          % Save results
//          Bpieces(:,:,p) = bnewest;
//          diststdpieces(:,p) = newestdisturbancestd;
//          cpieces(:,p) = cnewest;

            bpieces.add(bnewest);
//            diststdpieces.add(newestdisturbancestd);
//            cpieces.add(cnewest);

//
//        end

        }


//
//        for i=1:dims,
//          for j=1:dims,
//
//            themean = mean(Bpieces(i,j,:));
//            thestd = std(Bpieces(i,j,:));
//            if abs(themean)<prunefactor*thestd,
//          Bfinal(i,j) = 0;
//            else
//          Bfinal(i,j) = themean;
//            end
//
//          end
//        end

        TetradMatrix means = new TetradMatrix(X.rows(), X.rows());
        TetradMatrix stds = new TetradMatrix(X.rows(), X.rows());

        TetradMatrix BFinal = new TetradMatrix(X.rows(), X.rows());

        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.rows(); j++) {
//                double sum = 0.0;
//
//                for (int y = 0; y < npieces; y++) {
//                    sum += bpieces.get(y).get(i, j);
//                }
//
//                double themean = sum / (npieces);
//
//                double sumVar = 0.0;
//
//                for (int y = 0; y < npieces; y++) {
//                    sumVar += Math.pow((bpieces.get(y).get(i, j)) - themean, 2);
//                }
//
//                double thestd = Math.sqrt(sumVar / (npieces));
//
                double[] b = new double[bpieces.size()];

                for (int y = 0; y < bpieces.size(); y++) {
                    System.out.println("piece = " + bpieces.get(y));

                    b[y] = bpieces.get(y).get(i, j);
                }

                means.set(i, j, StatUtils.mean(b));
                stds.set(i, j, StatUtils.sd(b));

                if (Math.abs(StatUtils.mean(b)) < getPruneFactor() * StatUtils.sd(b)) {
                    BFinal.set(i, j, 0);
                } else {
                    BFinal.set(i, j, StatUtils.mean(b));
                }
            }
        }

        System.out.println("BFinal = " + BFinal);

//
//        diststdfinal = mean(diststdpieces,2);
//        cfinal = mean(cpieces,2);
//
//        % Finally, rename all the variables to the way we defined them
//        % in the function definition
//
//        Bpruned = Bfinal;
//        stde = diststdfinal;
//        ci = cfinal;

        return BFinal;
    }

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

    public double getPruneFactor() {
        return pruneFactor;
    }

    public void setPruneFactor(double pruneFactor) {
        this.pruneFactor = pruneFactor;
    }

    private int[] range(int i1, int i2) {
        if (i2 < i1) throw new IllegalArgumentException("i2 must be >=  i2 " + i1 + ", " + i2);
        int[] series = new int[i2 - i1 + 1];
        for (int i = 0; i <= i2 - i1; i++) {
            series[i] = i + i1;
        }
        try {
            return series;
        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }

    public int[] iperm(int[] k) {
        int[] ik = new int[k.length];

//        for (int i = 0; i < k.length; i++) {
//            ik[k[i]] = i;
//        }
//
        for (int i = 0; i < k.length; i++) {
            for (int j = 0; j < k.length; j++) {
                if (k[i] == j) {
                    ik[j] = i;
                }
            }
        }

        return ik;
    }
}

