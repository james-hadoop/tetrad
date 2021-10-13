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
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fask;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TetradSerializableUtils;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.abs;

/**
 * @author kaalpurush
 */
public class PairwiseForbiddenGraphModel extends KnowledgeBoxModel {

    static final long serialVersionUID = 23L;

    /**
     * Constructor from dataWrapper edge
     */
    public PairwiseForbiddenGraphModel(KnowledgeBoxInput input, Parameters params) {
        super(new KnowledgeBoxInput[]{input}, params);
//
//        if (input == null) {
//            throw new NullPointerException();
//        }
//
//
//        createKnowledge();
//
//        TetradLogger.getInstance().log("info", "Knowledge");
//
//        // This is a conundrum. At this point I dont know whether I am in a
//        // simulation or not. If in a simulation, I should print the knowledge.
//        // If not, I should wait for resetParams to be called. For now I'm
//        // printing the knowledge if it's not empty.
//        if (!((IKnowledge) params.get("knowledge", new Knowledge2())).isEmpty()) {
//            TetradLogger.getInstance().log("knowledge", params.get("knowledge", new Knowledge2()).toString());
//        }
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     *
     * @see TetradSerializableUtils
     */
    public static PairwiseForbiddenGraphModel serializableInstance() {
        return new PairwiseForbiddenGraphModel(GraphWrapper.serializableInstance(), new Parameters());
    }

    private void createKnowledge() {
        IKnowledge knwl = getKnowledge();
        if (knwl == null) {
            return;
        }

//        DataSet dataSet = (DataSet) (((DataWrapper) getKnowledgeBoxInput()).getDataModelList().get(0));
//
//        List<Node> vars = dataSet.getVariables();
//        List<Node> contVars = new ArrayList<>();
//
//        for (Node node : vars) {
//            if (node instanceof ContinuousVariable) contVars.add(node);
//        }
//
//        DataSet cont = dataSet.subsetColumns(contVars);
//
//        Fask fask = new Fask(DataUtils.getContinuousDataSet(cont),
//                new IndTestFisherZ(cont, 0.01));
//        fask.search();
//
//        knwl.clear();
//
//        List<Node> nodes = dataSet.getVariables();
//
//        Matrix _data = dataSet.getDoubleData();
//
//        int numOfNodes = nodes.size();
//
//        for (int i = 0; i < numOfNodes; i++) {
//            for (int j = i + 1; j < numOfNodes; j++) {
//                Node n1 = nodes.get(i);
//                Node n2 = nodes.get(j);
//
//                if (!(n1 instanceof ContinuousVariable && n2 instanceof ContinuousVariable)) {
//                    continue;
//                }
//
//                if (n1.getName().startsWith("E_") || n2.getName().startsWith("E_")) {
//                    continue;
//                }
//
//                double[] x1 = _data.getColumn(i).toArray();
//                double[] x2 = _data.getColumn(j).toArray();
//
//                double a1 = new AndersonDarlingTest(x1).getASquaredStar();
//                double a2 = new AndersonDarlingTest(x2).getASquaredStar();
//
//                double v = fask.leftRight(n1, n2);
//
//                if (a1 > 1 || a2 > 1) {
//                    if (v > 0) {
//                        knwl.setForbidden(n2.getName(), n1.getName());
//                    } else {
//                        knwl.setForbidden(n1.getName(), n2.getName());
//                    }
//                } else if (a1 < 1 && a2 < 1) {
//                    if (abs(v) > 0.1) {
//                        if (v > 0) {
//                            knwl.setForbidden(n2.getName(), n1.getName());
//                        } else {
//                            knwl.setForbidden(n1.getName(), n2.getName());
//                        }
//                    }
//                }
//            }
//        }
    }

    public Graph getResultGraph() {
        return null;
    }

}
