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

package edu.cmu.tetrad.sem;

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.SemGraph;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.TetradLogger;

/**
 * Optimizes a SEM using RICF (see that class).
 *
 * @author Joseph Ramsey
 */
public class SemOptimizerRicf implements SemOptimizer {
    static final long serialVersionUID = 23L;
    private int numRestarts = 1;

    //=============================CONSTRUCTORS=========================//

    /**
     * Blank constructor.
     */
    public SemOptimizerRicf() {
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static SemOptimizerRicf serializableInstance() {
        return new SemOptimizerRicf();
    }

    //==============================PUBLIC METHODS========================//

    /**
     * Optimizes the fitting function of the given Sem using the Powell method
     * from Numerical Recipes by adjusting the freeParameters of the Sem.
     */
    public void optimize(SemIm semIm) {
        if (this.numRestarts < 1) this.numRestarts = 1;

        if (this.numRestarts != 1) {
            throw new IllegalArgumentException("Number of restarts must be 1 for this method.");
        }

        Matrix sampleCovar = semIm.getSampleCovar();

        if (sampleCovar == null) {
            throw new NullPointerException("Sample covar has not been set.");
        }

        if (DataUtils.containsMissingValue(sampleCovar)) {
            throw new IllegalArgumentException("Please remove or impute missing values.");
        }

        if (DataUtils.containsMissingValue(sampleCovar)) {
            throw new IllegalArgumentException("Please remove or impute missing values.");
        }

        TetradLogger.getInstance().log("info", "Trying EM...");
//        new SemOptimizerEm().optimize(semIm);

        CovarianceMatrix cov = new CovarianceMatrix(semIm.getMeasuredNodes(),
                sampleCovar, semIm.getSampleSize());

        SemGraph graph = semIm.getSemPm().getGraph();
        Ricf.RicfResult result = new Ricf().ricf(graph, cov, 0.001);

        Matrix bHat = new Matrix(result.getBhat().toArray());
        Matrix lHat = new Matrix(result.getLhat().toArray());
        Matrix oHat = new Matrix(result.getOhat().toArray());

        for (Parameter param : semIm.getFreeParameters()) {
            if (param.getType() == ParamType.COEF) {
                int i = semIm.getSemPm().getVariableNodes().indexOf(param.getNodeA());
                int j = semIm.getSemPm().getVariableNodes().indexOf(param.getNodeB());
                semIm.setEdgeCoef(param.getNodeA(), param.getNodeB(), -bHat.get(j, i));
            }

            if (param.getType() == ParamType.VAR) {
                int i = semIm.getSemPm().getVariableNodes().indexOf(param.getNodeA());
                if (lHat.get(i, i) != 0) {
                    semIm.setErrVar(param.getNodeA(), lHat.get(i, i));
                } else if (oHat.get(i, i) != 0) {
                    semIm.setErrVar(param.getNodeA(), oHat.get(i, i));
                }
            }

            if (param.getType() == ParamType.COVAR) {
                int i = semIm.getSemPm().getVariableNodes().indexOf(param.getNodeA());
                int j = semIm.getSemPm().getVariableNodes().indexOf(param.getNodeB());
                if (lHat.get(i, i) != 0) {
                    semIm.setErrCovar(param.getNodeA(), param.getNodeB(), lHat.get(j, i));
                } else if (oHat.get(i, i) != 0) {
                    semIm.setErrCovar(param.getNodeA(), param.getNodeB(), oHat.get(j, i));
                }
            }
        }

        System.out.println(result);
        System.out.println(semIm);
    }

    @Override
    public void setNumRestarts(int numRestarts) {
        this.numRestarts = numRestarts;
    }

    @Override
    public int getNumRestarts() {
        return this.numRestarts;
    }

    public String toString() {
        return "Sem Optimizer RICF";
    }
}



