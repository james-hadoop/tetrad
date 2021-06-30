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
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Implements Chickering and Meek's (2002) locally consistent score criterion.
 *
 * @author Joseph Ramsey
 */
public class GraphScore implements Score {

    private Graph dag;
    private IndependenceFacts facts;

    // The variables of the covariance matrix.
    private List<Node> variables;

    // True if verbose output should be sent to out.
    private boolean verbose = false;

    /**
     * Constructs the score using a covariance matrix.
     */
    public GraphScore(Graph dag) {
        this.dag = dag;

        this.variables = new ArrayList<>();

        for (Node node : dag.getNodes()) {
            if (node.getNodeType() == NodeType.MEASURED) {
                this.variables.add(node);
            }
        }
    }

    public GraphScore(IndependenceFacts facts) {
        this.facts = facts;

        this.variables = new ArrayList<>();

        for (Node node : facts.getVariables()) {
            if (node.getNodeType() == NodeType.MEASURED) {
                this.variables.add(node);
            }
        }
    }

    /**
     * Calculates the sample likelihood and BIC score for y given its z in a simple SEM model
     */
    public double localScore(int y, int[] z) {
        return localScoreVisit(y, z) - z.length * 0.01;
    }

    private double localScoreVisit(int y, int[] z) {
        if (z.length == 0) {
            return 0;
        }

        int zn = z[z.length - 1];
        int[] sub1 = new int[z.length - 1];
        System.arraycopy(z, 0, sub1, 0, z.length - 1);

        int[] sub2 = new int[z.length - 1];
        System.arraycopy(sub1, 0, sub2, 0, z.length - 1);
        if (sub2.length > 0) sub2[sub2.length - 1] = zn;

        return (localScoreDiff(zn, y, sub1)) + localScoreVisit(y, sub2);
    }

    private List<Node> getVariableList(int[] indices) {
        List<Node> variables = new ArrayList<>();
        for (int i : indices) {
            variables.add(this.variables.get(i));
        }
        return variables;
    }


    @Override
    public double localScoreDiff(int x, int y, int[] z) {
        return locallyConsistentScoringCriterion(x, y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScoreDiff(x, y, new int[0]);
//        return localScore(y, x) - localScore(y);
    }

    private double locallyConsistentScoringCriterion(int x, int y, int[] z) {
        Node _y = variables.get(y);
        Node _x = variables.get(x);
        List<Node> _z = getVariableList(z);

        boolean dSeparatedFrom;

        if (dag != null) {
            dSeparatedFrom = dag.isDSeparatedFrom(_x, _y, _z);
        } else if (facts != null) {
            dSeparatedFrom = facts.isIndependent(_x, _y, _z);
        } else {
            throw new IllegalStateException("Expecting either a graph or a IndependenceFacts object.");
        }

        return dSeparatedFrom ? -1.0 : 1.0;
    }

    /**
     * Specialized scoring method for a single parent. Used to speed up the effect edges search.
     */

    public double localScore(int i, int parent) {
        throw new UnsupportedOperationException();
    }

    /**
     * Specialized scoring method for no parents. Used to speed up the effect edges search.
     */
    public double localScore(int i) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    public DataSet getDataSet() {
        throw new UnsupportedOperationException();
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    @Override
    public List<Node> getVariables() {
        return variables;
    }

    public Node getVariable(String name) {
        for (Node node : variables) {
            if (node.getName().equals(name)) {
                return node;
            }
        }

        throw new IllegalArgumentException("No variable by that name: " + name);
    }

    @Override
    public int getMaxDegree() {
        return 1000;
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        return false;
    }

    public int getSampleSize() {
        return 0;
    }

    public boolean getAlternativePenalty() {
        return false;
    }

    public void setAlternativePenalty(double alpha) {
        throw new UnsupportedOperationException("No alpha can be set when searching usign d-separation.");
    }

    public Graph getDag() {
        return new EdgeListGraph(dag);
    }

    public boolean isDSeparatedFrom(Node x, Node y, List<Node> z) {
        if (dag != null) {
            return dag.isDSeparatedFrom(x, y, z);
        } else if (facts != null) {
            return facts.isIndependent(x, y, z);
        }

        throw new IllegalArgumentException("Expecting either a DAG or an IndependenceFacts object.");
    }
    public boolean isDConnectedTo(Node x, Node y, List<Node> z) {
        return !isDSeparatedFrom(x, y, z);
    }
}



