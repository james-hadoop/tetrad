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
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.util.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Checks independence facts for variables associated with the nodes in a given graph by checking d-separation facts on
 * the underlying nodes.
 *
 * @author Joseph Ramsey
 */
public class IndTestDSep implements IndependenceTest {

    private IndependenceFacts independenceFacts;

    /**
     * The graph for which this is a variable map.
     */
    private Graph graph;

    /**
     * The list of observed variables (i.e. variables for observed nodes).
     */
    private List<Node> observedVars;
    private List<Node> _observedVars;
    private List<IndependenceFact> facts;
    private boolean verbose = false;
    private double pvalue = 0;

    public IndTestDSep(Graph graph) {
        this(graph, false);
    }

    public IndTestDSep(IndependenceFacts facts, List<Node> variables) {
        this(facts, false);
        facts.setNodes(variables);
    }

    public IndTestDSep(IndependenceFacts facts) {
        this(facts, false);
    }

    /**
     * Constructs a new independence test that returns d-separation facts for the given graph as independence results.
     */
    public IndTestDSep(Graph graph, boolean keepLatents) {
        if (graph == null) {
            throw new NullPointerException();
        }

        this.graph = graph;

        this._observedVars = calcVars(graph.getNodes(), keepLatents);
        this.observedVars = new ArrayList<>(_observedVars);
    }

    public IndTestDSep(IndependenceFacts facts, boolean keepLatents) {
        if (facts == null) {
            throw new NullPointerException();
        }

        this.independenceFacts = facts;

        this._observedVars = calcVars(facts.getVariables(), keepLatents);
        this.observedVars = new ArrayList<>(_observedVars);
    }

    /**
     * Required by IndependenceTest.
     */
    public IndependenceTest indTestSubset(List<Node> vars) {
        if (vars.isEmpty()) {
            throw new IllegalArgumentException("Subset may not be empty.");
        }

        List<Node> _vars = new ArrayList<>();

        for (Node var : vars) {
            Node _var = getVariable(var.getName());

            if (_var == null) {
                throw new IllegalArgumentException(
                        "All vars must be original vars");
            }

            _vars.add(_var);
        }

        this._observedVars = _vars;
        this.observedVars = new ArrayList<>(_observedVars);

        facts = new ArrayList<>();

        return this;
    }

    /**
     * @return the list of observed nodes in the given graph.
     */
    private List<Node> calcVars(List<Node> nodes, boolean keepLatents) {
        if (keepLatents) {
            return nodes;
        } else {
            List<Node> observedVars = new ArrayList<>();

            for (Node node : nodes) {
                if (node.getNodeType() == NodeType.MEASURED) {
                    observedVars.add(node);
                }
            }

            return observedVars;
        }
    }

    /**
     * Checks the indicated d-separation fact.
     *
     * @param x one node.
     * @param y a second node.
     * @param z a List of nodes (conditioning variables)
     * @return true iff x _||_ y | z
     */
    public IndependenceResult checkIndependence(Node x, Node y, List<Node> z) {
        if (z == null) {
            throw new NullPointerException();
        }

        for (Node node : z) {
            if (node == null) {
                throw new NullPointerException();
            }
        }

        if (!observedVars.contains(x)) {
            throw new IllegalArgumentException("Not an observed variable: " + x);
        }

        if (!observedVars.contains(y)) {
            throw new IllegalArgumentException("Not an observed variable: " + y);
        }

        for (Node _z : z) {
            if (!observedVars.contains(_z)) {
                throw new IllegalArgumentException("Not an observed variable: " + _z);
            }
        }

        boolean dSeparated;

        if (graph != null) {
            dSeparated = !getGraph().isDConnectedTo(x, y, z);
        } else {
            dSeparated = independenceFacts.isIndependent(x, y, z);
        }

//        if (verbose) {
//            if (dSeparated) {
//                double pValue = 1.0;
//                TetradLogger.getInstance().log("independencies", SearchLogUtils.independenceFactMsg(x, y, z, pValue));
//                System.out.println(SearchLogUtils.independenceFactMsg(x, y, z, pValue));
//            } else {
//                double pValue = 0.0;
//                TetradLogger.getInstance().log("dependencies", SearchLogUtils.dependenceFactMsg(x, y, z, pValue));
//                System.out.println(SearchLogUtils.dependenceFactMsg(x, y, z, pValue));
//            }
//        }

        double pValue;

        if (dSeparated) {
            if (this.facts != null) {
                this.facts.add(new IndependenceFact(x, y, z));
            }

            pValue = 1.0;
        } else {
            pValue = 0.0;
        }

        this.pvalue = pValue;

        return new IndependenceResult(new IndependenceFact(x, y, z), dSeparated, pValue);
    }

    /**
     * Auxiliary method to calculate dseparation facts directly from nodes instead of from variables.
     */
    public boolean isDSeparated(Node x, Node y, List<Node> z) {
        if (z == null) {
            throw new NullPointerException();
        }

        for (Node aZ : z) {
            if (aZ == null) {
                throw new NullPointerException();
            }
        }

        return getGraph().isDSeparatedFrom(x, y, z);
    }

    /**
     * Needed for IndependenceTest interface. P value is not meaningful here.
     */
    public double getPValue() {
        return this.pvalue;
    }

    /**
     * @return the list of TetradNodes over which this independence checker is capable of determinine independence
     * relations-- that is, all the variables in the given graph or the given data set.
     */
    public List<Node> getVariables() {
        return Collections.unmodifiableList(_observedVars);
    }

    /**
     * @return the list of variable varNames.
     */
    public List<String> getVariableNames() {
        List<Node> nodes = _observedVars;
        List<String> nodeNames = new ArrayList<>();
        for (Node var : nodes) {
            nodeNames.add(var.getName());
        }
        return nodeNames;
    }

    public boolean determines(List<Node> z, Node x1) {
        return false;
    }

    public double getAlpha() {
        return 0.5;
    }

    public void setAlpha(double alpha) {
        //
    }

    public Node getVariable(String name) {
        for (Node variable : observedVars) {
            if (variable.getName().equals(name)) {
                return variable;
            }
        }

        return null;
    }

    /**
     * @return the underlying graph.
     */
    public Graph getGraph() {
        return this.graph;
    }

    public void setGraph(Graph graph) {
        this.graph = graph;
    }

    public String toString() {
        return "D-separation";
    }

    public DataSet getData() {
        return null;
    }

    @Override
    public ICovarianceMatrix getCov() {
        return null;
    }

    @Override
    public List<DataSet> getDataSets() {
        return null;
    }

    @Override
    public int getSampleSize() {
        return 0;
    }

    @Override
    public List<Matrix> getCovMatrices() {
        return null;
    }

    @Override
    public double getScore() {
        return getPValue() == 1 ? -1 : 1;
    }

    public List<IndependenceFact> getFacts() {
        return facts;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }


}





