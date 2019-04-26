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
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements a Stochastic Complexity BIC score for FGS.
 *
 * @author Bryan Andrews
 */
public class StochasticComplexityScore implements Score {

    // The data
    private DataSet dataSet;

    // The number of instances
    private int N;

    // C2N table
    private double[] C2N;

    // The AD Tree used to count discrete cells.
    private AdLeafTree adTree;

    // The variables
    private List<Node> variables;

    // fNML:false vs qFML:true
    private boolean scoreEq = false;

    /**
     * Constructs the score
     */
    public StochasticComplexityScore(DataSet dataSet, boolean scoreEq) {
        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.dataSet = dataSet;
        this.N = dataSet.getNumRows();
        this.variables = dataSet.getVariables();
        this.adTree = new AdLeafTree(dataSet);
        this.scoreEq = scoreEq;

        if (!scoreEq) {
            this.C2N = new double[this.N];
            this.C2N[0] = 1;
            for (int n = 1; n < this.N; n++) {
                this.C2N[n] = (n*Math.PI/2) * Math.exp(Math.sqrt(8/(9*n*Math.PI)) + (3*Math.PI-16)/(36*n*Math.PI));
            }
        }

    }

    private double nCr(int n, int r) {
        int rfact=1, nfact=1, nrfact=1,temp1 = n-r ,temp2 = r;
        if(r>n-r)
        {
            temp1 =r;
            temp2 =n-r;
        }
        for(int i=1;i<=n;i++)
        {
            if(i<=temp2)
            {
                rfact *= i;
                nrfact *= i;
            }
            else if(i<=temp1)
            {
                nrfact *= i;
            }
            nfact *= i;
        }
        return nfact/(double)(rfact*nrfact);
    }

    /**
     * Calculates the sample likelihood and BIC score for i given its parents in a simple SEM model
     */
    public double localScore(int i, int... parents) {
        List<DiscreteVariable> A = new ArrayList<>();
        for (int j : parents)
        {
            A.add(((DiscreteVariable) this.variables.get(j)));
        }
        A.add(((DiscreteVariable) this.variables.get(i)));
        List<List<Integer>> nCells = this.adTree.getCellLeaves(A);
        A.remove(A.size()-1);
        List<List<Integer>> dCells = this.adTree.getCellLeaves(A);
        int r = ((DiscreteVariable) this.variables.get(i)).getNumCategories();
        int q = dCells.size();

        double plik = 0;
        for (List<Integer> cell : nCells) {
            if (cell.size() == 0) {
                continue;
            }
            plik += cell.size() * Math.log(((double) cell.size()) / this.N);
        }
        for (List<Integer> cell : dCells) {
            if (cell.size() == 0) {
                continue;
            }
            plik -= cell.size() * Math.log(((double) cell.size()) / this.N);
        }

        if (this.scoreEq) {
            plik -= (reg(q*r) - reg(q));
        } else {
            for (List<Integer> cell : nCells){
                int qi = cell.size();
                plik -= Math.log(Crn(r,qi));
            }
        }

        return plik;
    }

    private double reg(double n) {
        double a = n/this.N;
        double Ca = 0.5 + 0.5*Math.sqrt(1+(4/a));

        return this.N*(Math.log(a)+(a+2)*Math.log(Ca)-(1/Ca)) - 0.5*(Ca+(2/a));
    }

    private double Crn(int r, int n) {
        if (n == 0) {
            return 1;
        } else if(r == 1) {
            return 1;
        } else if(r == 2) {
            return this.C2N[n];
        } else {
            return Crn(r-1, n) + (n/(r-2))*Crn(r-2, n);
        }
    }

    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        double diff = localScore(y, x) - localScore(y);
        return diff;
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
        return dataSet.getNumRows();
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    @Override
    public List<Node> getVariables() {
        return variables;
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
        return (int) Math.ceil(Math.log(dataSet.getNumRows()));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        return false;
    }

    public boolean getScoreEq() {
        return this.scoreEq;
    }

    public void setScoreEq(boolean scoreEq) {
        this.scoreEq = scoreEq;
    }

}



