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
import edu.cmu.tetrad.graph.Node;

import java.util.List;

/**
 * Implements a mixed variable polynomial BIC score for fGES.
 *
 * @author Bryan Andrews
 */
public class MVPScore implements Score {

    private final DataSet dataSet;

    // The variables of the continuousData set.
    private final List<Node> variables;

    // Likelihood function
    private final MVPLikelihood likelihood;

    // Log number of instances
    private final double logn;

    public MVPScore(DataSet dataSet, double structurePrior, int fDegree, boolean discretize) {

        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.dataSet = dataSet;
        this.variables = dataSet.getVariables();
        this.likelihood = new MVPLikelihood(dataSet, structurePrior, fDegree, discretize);
        this.logn = Math.log(dataSet.getNumRows());
    }

    public double localScore(int i, int... parents) {

        double lik = this.likelihood.getLik(i, parents);
        double dof = this.likelihood.getDoF(i, parents);
        double sp = this.likelihood.getStructurePrior(parents.length);

        if (sp > 0) {
            sp = -2 * dof * sp;
        }

        double score = 2.0 * lik - dof * this.logn + sp;

        if (Double.isNaN(score) || Double.isInfinite(score)) {
            return Double.NaN;
        } else {
            return score;
        }
    }

    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScore(y, x) - localScore(y);
    }

    private int[] append(int[] parents, int extra) {
        int[] all = new int[parents.length + 1];
        System.arraycopy(parents, 0, all, 0, parents.length);
        all[parents.length] = extra;
        return all;
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

    public int getSampleSize() {
        return this.dataSet.getNumRows();
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    @Override
    public List<Node> getVariables() {
        return this.variables;
    }

    @Override
    public Node getVariable(String targetName) {

        for (Node node : this.variables) {
            if (node.getName().equals(targetName)) {
                return node;
            }
        }

        return null;
    }

    @Override
    public int getMaxDegree() {
        return (int) Math.ceil(Math.log(this.dataSet.getNumRows()));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        return false;
    }

}



