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

package edu.cmu.tetrad.data;

import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Matrix;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Stores time series data as a list of continuous columns.
 *
 * @author Joseph Ramsey
 */
public final class TimeSeriesData implements DataModel {
    static final long serialVersionUID = 23L;

    /**
     * @serial
     */
    private final Matrix data2;

    /**
     * @serial
     */
    private String name;

    /**
     * @serial
     */
    private final List<String> varNames;

    /**
     * @serial
     */
    private IKnowledge knowledge = new Knowledge2();

    //==============================CONSTRUCTOR===========================//

    /**
     * Constructs a new time series data contains for the given row-major data
     * array and the given list of variables. Each row of the data, data[i],
     * contains a measured for each variable (in order) for a particular time.
     * The series of times is in increasing order.
     */
    public TimeSeriesData(Matrix matrix, List<String> varNames) {
        if (matrix == null) {
            throw new NullPointerException("Data must not be null.");
        }

        if (varNames == null) {
            throw new NullPointerException("Variables must not be null.");
        }
        for (int i = 0; i < varNames.size(); i++) {
            if (varNames.get(i) == null) {
                throw new NullPointerException(
                        "Variable at index " + i + "is null.");
            }
        }
        this.data2 = matrix;
        if (varNames.size() != matrix.columns()) {
            throw new IllegalArgumentException(
                    "Number of columns in the data " +
                            "must match the number of variables.");
        }
        this.varNames = varNames;
        this.name = "Time Series Data";
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static TimeSeriesData serializableInstance() {
        List<String> varNames = new ArrayList<>();
        varNames.add("X");
        varNames.add("Y");
        return new TimeSeriesData(new Matrix(2, 2), varNames);
    }

    //=================================PUBLIC METHODS======================//

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        if (name == null) {
            throw new NullPointerException("Name must not be null.");
        }
        this.name = name;
    }

    @Override
    public boolean isContinuous() {
        return true;
    }

    @Override
    public boolean isDiscrete() {
        return false;
    }

    @Override
    public boolean isMixed() {
        return false;
    }

    @Override
    public Node getVariable(String name) {
        return null;
    }

    @Override
    public DataModel copy() {
        return null;
    }

    public List<Node> getVariables() {
        List<String> varNames = getVariableNames();
        List<Node> vars = new LinkedList<>();

        for (String varName : varNames) {
            vars.add(new ContinuousVariable(varName));
        }

        return vars;
    }

    public IKnowledge getKnowledge() {
        System.out.println();

        return this.knowledge.copy();
    }

    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) {
            throw new NullPointerException();
        }

        this.knowledge = knowledge.copy();
    }

    /**
     * @return (a copy of) the List of Variables for the data set, in the order
     * of their columns.
     */
    public List<String> getVariableNames() {
        return this.varNames;
    }

    public Matrix getData() {
        return this.data2.copy();
    }

    public int getNumTimePoints() {
        return getData().rows();
    }

    public int getNumVars() {
        return getVariableNames().size();
    }

    public double getDatum(int row, int col) {
        return this.data2.get(row, col);
    }

    /**
     * Adds semantic checks to the default deserialization method. This method
     * must have the standard signature for a readObject method, and the body of
     * the method must begin with "s.defaultReadObject();". Other than that, any
     * semantic checks can be specified and do not need to stay the same from
     * version to version. A readObject method of this form may be added to any
     * class, even if Tetrad sessions were previously saved out using a version
     * of the class that didn't include it. (That's what the
     * "s.defaultReadObject();" is for. See J. Bloch, Effective Java, for help.
     */
    private void readObject(ObjectInputStream s)
            throws IOException, ClassNotFoundException {
        s.defaultReadObject();

        if (this.name == null) {
            throw new NullPointerException();
        }

        if (this.varNames == null) {
            throw new NullPointerException();
        }

        if (this.knowledge == null) {
            throw new NullPointerException();
        }
    }
}




