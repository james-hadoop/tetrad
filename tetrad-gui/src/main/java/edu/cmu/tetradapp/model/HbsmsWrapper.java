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

package edu.cmu.tetradapp.model;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.util.*;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Joseph Ramsey
 */

public class HbsmsWrapper extends AbstractAlgorithmRunner implements GraphSource {
    static final long serialVersionUID = 23L;

    public enum AlgorithmType {
        BEAM, FGES
    }

    private AlgorithmType algorithmType = AlgorithmType.BEAM;

    private String name;
    private Graph initialGraph;
    private Graph graph;
    private transient List<PropertyChangeListener> listeners;

    private Parameters params2;

    private SemIm originalSemIm;
    private SemIm newSemIm;

    //============================CONSTRUCTORS============================//

    public HbsmsWrapper(DataWrapper dataWrapper,
                        Parameters params, KnowledgeBoxModel knowledgeBoxModel) {
        super(dataWrapper, params, knowledgeBoxModel);
        this.params2 = params;
        this.graph = new EdgeListGraph(dataWrapper.getSelectedDataModel().getVariables());
    }

    public HbsmsWrapper(DataWrapper dataWrapper,
                        Parameters params) {
        super(dataWrapper, params, null);
        this.params2 = params;
        setGraph(new EdgeListGraph(dataWrapper.getSelectedDataModel().getVariables()));
    }

    private void setGraph(EdgeListGraph graph) {
        this.graph = new EdgeListGraph(graph);
        this.initialGraph = new EdgeListGraph(graph);
    }

    public HbsmsWrapper(GraphSource graphWrapper,
                        DataWrapper dataWrapper,
                        Parameters params, KnowledgeBoxModel knowledgeBoxModel) {
        super(dataWrapper, params, knowledgeBoxModel);
        this.params2 = params;
        setGraph(new EdgeListGraph(graphWrapper.getGraph()));
    }

    public HbsmsWrapper(GraphSource graphWrapper,
                        DataWrapper dataWrapper,
                        Parameters params) {
        super(dataWrapper, params);
        this.params2 = params;
        setGraph(new EdgeListGraph(graphWrapper.getGraph()));
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     *
     * @see TetradSerializableUtils
     */
    public static PcRunner serializableInstance() {
        return PcRunner.serializableInstance();
    }

    //============================PUBLIC METHODS==========================//


    public AlgorithmType getAlgorithmType() {
        return algorithmType;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return this.name;
    }

    public boolean isShuffleMoves() {
        return false;
    }

    /**
     * Executes the algorithm, producing (at least) a result workbench. Must be implemented in the extending class.
     */

    public void execute() {
        DataModel dataModel = getDataModel();
        assert (dataModel != null);

        IKnowledge knowledge = (IKnowledge) params2.get("knowledge", new Knowledge2());

        if (initialGraph == null) {
            initialGraph = new EdgeListGraph(dataModel.getVariables());
        }

        Graph graph2 = new EdgeListGraph(initialGraph);
        graph2 = GraphUtils.replaceNodes(graph2, dataModel.getVariables());

        Hbsms search;

        if (dataModel instanceof DataSet) {
            DataSet dataSet = (DataSet) dataModel;

            if (getAlgorithmType() == AlgorithmType.BEAM) {
                search = new HbsmsBeam(graph2, dataSet, knowledge);
            } else if (getAlgorithmType() == AlgorithmType.FGES) {
                search = new BffGes(graph2, dataSet);
                search.setKnowledge(knowledge);
            } else {
                throw new IllegalStateException();
            }
        }
        else if (dataModel instanceof CovarianceMatrix) {
            throw new IllegalArgumentException("HBSMS requires tabular data.");
        }
        else {
            throw new IllegalStateException();
        }

        Parameters params = getParams();

        search.setAlpha(params.getDouble("alpha", 0.001));
        search.setBeamWidth(params.getInt("beamWidth", 5));
        search.setHighPValueAlpha(params.getDouble("zeroEdgeP", 0.05));
        this.graph = search.search();

//        this.graph = search.getNewSemIm().getSemPm().getGraph();

        setOriginalSemIm(search.getOriginalSemIm());
        this.newSemIm = search.getNewSemIm();
        fireGraphChange(graph);

        if (getSourceGraph() != null) {
            GraphUtils.arrangeBySourceGraph(graph, getSourceGraph());
        } else if (knowledge.isDefaultToKnowledgeLayout()) {
            SearchGraphUtils.arrangeByKnowledgeTiers(graph, knowledge);
        } else {
            GraphUtils.circleLayout(graph, 200, 200, 150);
        }

        setResultGraph(SearchGraphUtils.patternForDag(graph, knowledge));
    }

    public boolean supportsKnowledge() {
        return true;
    }

    public ImpliedOrientation getMeekRules() {
        return new MeekRules();
    }

    @Override
    public String getAlgorithmName() {
        return "BFF";
    }

    public void setAlgorithmType(AlgorithmType algorithmType) {
        this.algorithmType = algorithmType;
    }

    public void addPropertyChangeListener(PropertyChangeListener l) {
        if (!getListeners().contains(l)) getListeners().add(l);
    }

    private void fireGraphChange(Graph graph) {
        for (PropertyChangeListener l : getListeners()) {
            l.propertyChange(new PropertyChangeEvent(this, "graph", null, graph));
        }
    }

    public Graph getGraph() {
        return getResultGraph();
    }

    /**
     * @return the names of the triple classifications. Coordinates with
     */
    public List<String> getTriplesClassificationTypes() {
        return new LinkedList<>();
    }

    /**
     * @param node The node that the classifications are for. All triple from adjacencies to this node to adjacencies to
     *             this node through the given node will be considered.
     * @return the list of triples corresponding to <code>getTripleClassificationNames</code> for the given node.
     */
    public List<List<Triple>> getTriplesLists(Node node) {
        return new LinkedList<>();
    }


    private List<PropertyChangeListener> getListeners() {
        if (listeners == null) {
            listeners = new ArrayList<>();
        }
        return listeners;
    }


    private void setOriginalSemIm(SemIm originalSemIm) {
        if (this.originalSemIm == null) {
            this.originalSemIm = originalSemIm;
        }
    }


    /**
     * Adds semantic checks to the default deserialization method. This method must have the standard signature for a
     * readObject method, and the body of the method must begin with "s.defaultReadObject();". Other than that, any
     * semantic checks can be specified and do not need to stay the same from version to version. A readObject method of
     * this form may be added to any class, even if Tetrad sessions were previously saved out using a version of the
     * class that didn't include it. (That's what the "s.defaultReadObject();" is for. See J. Bloch, Effective Java, for
     * help.
     */
    private void readObject(ObjectInputStream s)
            throws IOException, ClassNotFoundException {
        s.defaultReadObject();

        if (params2 == null) {
            params2 = new Parameters();
        }
    }

    public SemIm getOriginalSemIm() {
        return originalSemIm;
    }

    public SemIm getNewSemIm() {
        return newSemIm;
    }

    public void setNewSemIm(SemIm newSemIm) {
        this.newSemIm = newSemIm;
    }
}



