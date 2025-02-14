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

package edu.cmu.tetradapp.model;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.Triple;
import edu.cmu.tetrad.regression.Regression;
import edu.cmu.tetrad.regression.RegressionCovariance;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.search.ImpliedOrientation;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TetradSerializableUtils;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;

/**
 * Extends AbstractAlgorithmRunner to produce a wrapper for the Regression
 * algorithm.
 *
 * @author Frank Wimberly after Joe Ramsey's PcRunner
 */
public class RegressionRunner implements AlgorithmRunner, RegressionModel {
    static final long serialVersionUID = 23L;
    private List<String> regressorNames;
    private String name;
    private final Parameters params;
    private String targetName;
    private final DataModelList dataModels;
    private Graph outGraph;
    private RegressionResult result;
    private Map<String, String> allParamsSettings;
    private final List<String> variableNames;

    private int numModels = 1;
    private int modelIndex;
    private String modelSourceName;

    //=========================CONSTRUCTORS===============================//

    /**
     * Constructs a wrapper for the given DataWrapper. The DataWrapper must
     * contain a DataSet that is either a DataSet or a DataSet or a DataList
     * containing either a DataSet or a DataSet as its selected model.
     */
    public RegressionRunner(DataWrapper dataWrapper, Parameters params) {
        if (dataWrapper == null) {
            throw new NullPointerException();
        }

        if (params == null) {
            throw new NullPointerException();
        }

        if (dataWrapper instanceof Simulation) {
            Simulation simulation = (Simulation) dataWrapper;
            this.numModels = dataWrapper.getDataModelList().size();
            this.modelIndex = 0;
            this.modelSourceName = simulation.getName();
        }

        this.params = params;

        DataModel dataModel = dataWrapper.getSelectedDataModel();

        if (dataModel instanceof DataSet) {
            DataSet _dataSet = (DataSet) dataModel;
            if (!_dataSet.isContinuous()) {
                throw new IllegalArgumentException("Data set must be continuous.");
            }
        }

        this.dataModels = dataWrapper.getDataModelList();

        this.variableNames = dataModel.getVariableNames();
        this.targetName = null;
        this.regressorNames = new ArrayList<>();

        TetradLogger.getInstance().log("info", "Linear Regression");

        if (this.result == null) {
            TetradLogger.getInstance().log("info", "Please double click this regression node to run the regession.");
        } else {
            TetradLogger.getInstance().log("result", "\n" + this.result.getResultsTable().toString());
        }
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     *
     * @see TetradSerializableUtils
     */
    public static RegressionRunner serializableInstance() {
        List<Node> variables = new LinkedList<>();
        ContinuousVariable var1 = new ContinuousVariable("X");
        ContinuousVariable var2 = new ContinuousVariable("Y");

        variables.add(var1);
        variables.add(var2);
        DataSet _dataSet = new BoxDataSet(new DoubleDataBox(3, variables.size()), variables);
        double[] col1data = {0.0, 1.0, 2.0};
        double[] col2data = {2.3, 4.3, 2.5};

        for (int i = 0; i < 3; i++) {
            _dataSet.setDouble(i, 0, col1data[i]);
            _dataSet.setDouble(i, 1, col2data[i]);
        }

        DataWrapper dataWrapper = new DataWrapper(_dataSet);
        return new RegressionRunner(dataWrapper, new Parameters());
    }

    //===========================PUBLIC METHODS============================//

    public DataModel getDataModel() {
        //return (DataModel) this.dataWrapper.getDataModelList().get(0);
        return this.dataModels.get(getModelIndex());
    }

    public Parameters getParams() {
        return this.params;
    }

    public Graph getResultGraph() {
        return this.outGraph;
    }

    private void setResultGraph(Graph graph) {
        this.outGraph = graph;
    }

    public Graph getSourceGraph() {
        return null;
    }
    //=================PUBLIC METHODS OVERRIDING ABSTRACT=================//

    /**
     * Executes the algorithm, producing (at least) a result workbench. Must be
     * implemented in the extending class.
     */
    public void execute() {
        if (this.regressorNames.size() == 0 || this.targetName == null) {
            this.outGraph = new EdgeListGraph();
            return;
        }

        if (this.regressorNames.contains(this.targetName)) {
            this.outGraph = new EdgeListGraph();
            return;
        }

        Regression regression;
        Node target;
        List<Node> regressors;

        if (getDataModel() instanceof DataSet) {
            DataSet _dataSet = (DataSet) getDataModel();
            regression = new RegressionDataset(_dataSet);
            target = _dataSet.getVariable(this.targetName);
            regressors = new LinkedList<>();

            for (String regressorName : this.regressorNames) {
                regressors.add(_dataSet.getVariable(regressorName));
            }

            double alpha = this.params.getDouble("alpha", 0.001);
            regression.setAlpha(alpha);

            this.result = regression.regress(target, regressors);
            this.outGraph = regression.getGraph();
        } else if (getDataModel() instanceof ICovarianceMatrix) {
            ICovarianceMatrix covariances = (ICovarianceMatrix) getDataModel();
            regression = new RegressionCovariance(covariances);
            target = covariances.getVariable(this.targetName);
            regressors = new LinkedList<>();

            for (String regressorName : this.regressorNames) {
                regressors.add(covariances.getVariable(regressorName));
            }

            double alpha = this.params.getDouble("alpha", 0.001);
            regression.setAlpha(alpha);

            this.result = regression.regress(target, regressors);
            this.outGraph = regression.getGraph();
        }

        setResultGraph(this.outGraph);
    }

    public boolean supportsKnowledge() {
        return false;
    }

    public ImpliedOrientation getMeekRules() {
        throw new UnsupportedOperationException();
    }

    public void setExternalGraph(Graph graph) {
    }

    public Graph getExternalGraph() {
        return null;
    }

    @Override
    public String getAlgorithmName() {
        return "Regression";
    }

    public RegressionResult getResult() {
        return this.result;
    }

    public Graph getOutGraph() {
        return this.outGraph;
    }

    @Override
    public List<String> getVariableNames() {
        return this.variableNames;
    }

    @Override
    public List<String> getRegressorNames() {
        return this.regressorNames;
    }

    @Override
    public void setRegressorName(List<String> predictors) {
        this.regressorNames = predictors;
    }

    public String getTargetName() {
        return this.targetName;
    }

    @Override
    public void setTargetName(String target) {
        this.targetName = target;
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

        if (this.params == null) {
            throw new NullPointerException();
        }

    }

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Graph getGraph() {
        return this.outGraph;
    }

    /**
     * @return the names of the triple classifications. Coordinates with
     */
    public List<String> getTriplesClassificationTypes() {
        return new LinkedList<>();
    }

    /**
     * @param node The node that the classifications are for. All triple from adjacencies to this
     *             node to adjacencies to this node through the given node will be considered.
     * @return the list of triples corresponding to <code>getTripleClassificationNames</code>
     * for the given node.
     */
    public List<List<Triple>> getTriplesLists(Node node) {
        return new LinkedList<>();
    }

    @Override
    public Map<String, String> getParamSettings() {
        Map<String, String> paramSettings = new HashMap<>();
        paramSettings.put("Algorithm", "Regression");
        return paramSettings;
    }


    @Override
    public void setAllParamSettings(Map<String, String> paramSettings) {
        this.allParamsSettings = paramSettings;
    }

    @Override
    public Map<String, String> getAllParamSettings() {
        return this.allParamsSettings;
    }

    public int getNumModels() {
        return this.numModels;
    }

    public int getModelIndex() {
        return this.modelIndex;
    }

    public String getModelSourceName() {
        return this.modelSourceName;
    }

    public void setModelIndex(int modelIndex) {
        this.modelIndex = modelIndex;
    }

    @Override
    public List<Graph> getGraphs() {
        return null;
    }
}





