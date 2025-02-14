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
import edu.cmu.tetrad.util.Matrix;

import java.util.List;

public class PermutationMatrixPair {

    private final Matrix matrixW;
    private DataSet matrixBhat;
    private final List<Integer> permutation;

    public PermutationMatrixPair(List<Integer> permutation, Matrix matrixW) {
        this.permutation = permutation;
        this.matrixW = matrixW;
    }

    public Matrix getMatrixW() {
        return this.matrixW;
    }

    public DataSet getMatrixBhat() {
        return this.matrixBhat;
    }

    public void setMatrixBhat(DataSet matrixBhat) {
        this.matrixBhat = matrixBhat;
    }

    public List<Integer> getPermutation() {
        return this.permutation;
    }

    public String toString() {
        return "Permutation: " + this.permutation + "matrix W : " + this.matrixW;
    }
}



