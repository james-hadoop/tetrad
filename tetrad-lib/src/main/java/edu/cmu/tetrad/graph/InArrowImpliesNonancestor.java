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

package edu.cmu.tetrad.graph;


/**
 * A graph constraint that checks to make sure that in-arrows imply
 * non-ancestry-- that is, when adding an edge A *--&gt; B, where * is any
 * endpoint, B is not an ancestor of A.<p>This can be used, e.g., to guard
 * against cycles in a DAG.  If we may assume that graph G is already acyclic,
 * then we can use this constraint to check whether adding a new directed edge
 * will create a new cycle.
 *
 * @author Joseph Ramsey
 */
public final class InArrowImpliesNonancestor implements GraphConstraint {
    static final long serialVersionUID = 23L;

    //=============================CONSTRUCTORS===========================//

    public InArrowImpliesNonancestor() {

    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static InArrowImpliesNonancestor serializableInstance() {
        return new InArrowImpliesNonancestor();
    }

    //=============================PUBLIC METHODS=========================//

    /**
     * @return a string representation of the constraint.
     */
    public String toString() {
        return "<Arrow implies non-ancestor.>";
    }
}





