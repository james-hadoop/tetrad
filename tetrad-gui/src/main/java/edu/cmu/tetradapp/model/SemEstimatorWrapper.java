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
import edu.cmu.tetrad.graph.Dag;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.sem.ParamType;
import edu.cmu.tetrad.sem.Parameter;
import edu.cmu.tetrad.sem.ScoreType;
import edu.cmu.tetrad.sem.SemEstimator;
import edu.cmu.tetrad.sem.LinearSemIm;
import edu.cmu.tetrad.sem.SemOptimizer;
import edu.cmu.tetrad.sem.SemOptimizerEm;
import edu.cmu.tetrad.sem.SemOptimizerPowell;
import edu.cmu.tetrad.sem.SemOptimizerRegression;
import edu.cmu.tetrad.sem.SemOptimizerRicf;
import edu.cmu.tetrad.sem.SemOptimizerScattershot;
import edu.cmu.tetrad.sem.LinearSemPm;
import edu.cmu.tetrad.session.SessionModel;
import edu.cmu.tetrad.util.JOptionUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TetradSerializableUtils;
import edu.cmu.tetrad.util.Unmarshallable;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import javax.swing.JOptionPane;

/**
 * Wraps a SemEstimator for use in the Tetrad application.
 *
 * @author Joseph Ramsey
 */
public class SemEstimatorWrapper implements SessionModel, Unmarshallable {

    static final long serialVersionUID = 23L;

    /**
     * @serial Can be null.
     */
    private String name;

    /**
     * @serial Cannot be null.
     */
    private SemEstimator semEstimator;
    private LinearSemPm linearSemPm;

    private boolean multipleResults = false;

    private List<SemEstimator> multipleResultList = new ArrayList<>();
    private Parameters params;

    //==============================CONSTRUCTORS==========================//
    /**
     * Private constructor for serialization only. Problem is, for the real
     * constructors, I'd like to call the degrees of freedom check, which pops
     * up a dialog. This is irritating when running unit tests. jdramsey 8/29/07
     */
    public SemEstimatorWrapper(DataSet dataSet, LinearSemPm linearSemPm, Parameters params) {
        this.params = params;
        this.linearSemPm = linearSemPm;

        DataModelList dataModel = new DataModelList();
        dataModel.add(dataSet);

        this.params = params;

        for (DataModel model : dataModel) {
            if (model instanceof DataSet) {
                this.linearSemPm = linearSemPm;
                SemEstimator estimator = new SemEstimator(dataSet, linearSemPm, getOptimizer());
                estimator.setNumRestarts(getParams().getInt("numRestarts", 1));
                estimator.setScoreType((ScoreType) getParams().get("scoreType", ScoreType.Fgls));
                if (!degreesOfFreedomCheck(linearSemPm));
                estimator.estimate();

                getMultipleResultList().add(estimator);
            } else if (model instanceof ICovarianceMatrix) {
                ICovarianceMatrix covMatrix = new CovarianceMatrix((ICovarianceMatrix) model);
                this.linearSemPm = linearSemPm;
                SemEstimator estimator = new SemEstimator(covMatrix, linearSemPm, getOptimizer());
                estimator.setNumRestarts(getParams().getInt("numRestarts", 1));
                estimator.setScoreType((ScoreType) getParams().get("scoreType", ScoreType.SemBic));
                if (!degreesOfFreedomCheck(linearSemPm));
                estimator.estimate();

                getMultipleResultList().add(estimator);
            } else {
                throw new IllegalArgumentException("Data must consist of continuous data sets or covariance matrices.");
            }
        }

        if (dataModel != null) {
            multipleResults = true;

            this.semEstimator = getMultipleResultList().get(0);
        } else {
            throw new IllegalArgumentException("Data must consist of continuous data sets or covariance matrices.");
        }

        this.semEstimator = new SemEstimator(dataSet, linearSemPm, getOptimizer());
    }

    public SemEstimatorWrapper(DataWrapper dataWrapper,
                               LinearSemPmWrapper semPmWrapper, Parameters params) {
        if (dataWrapper == null) {
            throw new NullPointerException("Data wrapper must not be null.");
        }

        if (semPmWrapper == null) {
            throw new NullPointerException(
                    "OldSem PM Wrapper must not be null.");
        }

        DataModel dataModel = dataWrapper.getDataModelList();

        this.params = params;

        if (dataModel != null) {
            multipleResults = true;

            if (setParams(semPmWrapper, (DataModelList) dataModel)) {
                return;
            }

            this.semEstimator = getMultipleResultList().get(0);
        } else {
            throw new IllegalArgumentException("Data must consist of continuous data sets or covariance matrices.");
        }

        log();
    }

    private boolean setParams(LinearSemPmWrapper semPmWrapper, DataModelList dataModel) {
        for (DataModel model : dataModel) {
            if (model instanceof DataSet) {
                DataSet dataSet = (DataSet) model;
                LinearSemPm linearSemPm = semPmWrapper.getSemPm();
                this.linearSemPm = linearSemPm;
                SemEstimator estimator = new SemEstimator(dataSet, linearSemPm, getOptimizer());
                estimator.setNumRestarts(getParams().getInt("numRestarts", 1));
                estimator.setScoreType((ScoreType) getParams().get("scoreType", ScoreType.Fgls));
                if (!degreesOfFreedomCheck(linearSemPm)) {
                    return true;
                }
                estimator.estimate();

                getMultipleResultList().add(estimator);
            } else if (model instanceof ICovarianceMatrix) {
                ICovarianceMatrix covMatrix = new CovarianceMatrix((ICovarianceMatrix) model);
                LinearSemPm linearSemPm = semPmWrapper.getSemPm();
                this.linearSemPm = linearSemPm;
                SemEstimator estimator = new SemEstimator(covMatrix, linearSemPm, getOptimizer());
                estimator.setNumRestarts(getParams().getInt("numRestarts", 1));
                estimator.setScoreType((ScoreType) getParams().get("scoreType", ScoreType.Fgls));

                if (!degreesOfFreedomCheck(linearSemPm)) {
                    return true;
                }
                estimator.estimate();

                getMultipleResultList().add(estimator);
            } else {
                throw new IllegalArgumentException("Data must consist of continuous data sets or covariance matrices.");
            }
        }
        return false;
    }

    public SemEstimatorWrapper(DataWrapper dataWrapper,
                               LinearSemImWrapper semImWrapper, Parameters params) {
    	this(dataWrapper, new LinearSemPmWrapper(semImWrapper), params);
    }
    
    private boolean degreesOfFreedomCheck(LinearSemPm linearSemPm) {
        if (linearSemPm.getDof() < 1) {
            int ret = JOptionPane.showConfirmDialog(JOptionUtils.centeringComp(),
                    "This model has nonpositive degrees of freedom (DOF = "
                    + linearSemPm.getDof() + "). "
                    + "\nEstimation will be uninformative. Are you sure you want to proceed?",
                    "Please confirm", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);

            if (ret != JOptionPane.YES_OPTION) {
                return false;
            }
        }

        return true;
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     *
     * @see TetradSerializableUtils
     */
    public static SemEstimatorWrapper serializableInstance() {
        List<Node> variables = new LinkedList<>();
        ContinuousVariable x = new ContinuousVariable("X");
        variables.add(x);
        DataSet dataSet = new BoxDataSet(new VerticalDoubleDataBox(10, variables.size()), variables);

        for (int i = 0; i < dataSet.getNumRows(); i++) {
            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                dataSet.setDouble(i, j, RandomUtil.getInstance().nextDouble());
            }
        }
        Dag dag = new Dag();
        dag.addNode(x);
        LinearSemPm pm = new LinearSemPm(dag);
        Parameters params1 = new Parameters();
        return new SemEstimatorWrapper(dataSet, pm, params1);
    }

    //============================PUBLIC METHODS=========================//
    public SemEstimator getSemEstimator() {
        return this.semEstimator;
    }

    public void setSemEstimator(SemEstimator semEstimator) {
        this.semEstimator = semEstimator;
    }

    public LinearSemIm getEstimatedSemIm() {
        return semEstimator.getEstimatedSem();
    }

    public String getSemOptimizerType() {
        return getParams().getString("semOptimizerType", "Regression");
    }

    public void setSemOptimizerType(String type) {
        getParams().set("semOptimizerType", type);
    }

    public Graph getGraph() {
        return semEstimator.getEstimatedSem().getSemPm().getGraph();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    //=============================== Private methods =======================//
    private void log() {
        TetradLogger.getInstance().log("info", "SEM Estimator:");
        TetradLogger.getInstance().log("im", "" + getEstimatedSemIm());
        TetradLogger.getInstance().log("details", "ChiSq = " + getEstimatedSemIm().getChiSquare());
        TetradLogger.getInstance().log("details", "DOF = " + getEstimatedSemIm().getSemPm().getDof());
        TetradLogger.getInstance().log("details", "P = " + getEstimatedSemIm().getPValue());
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
     *
     * @throws java.io.IOException
     * @throws ClassNotFoundException
     */
    private void readObject(ObjectInputStream s)
            throws IOException, ClassNotFoundException {
        s.defaultReadObject();

//        if (semEstimator == null) {
//            throw new NullPointerException();
//        }
    }

    public boolean isMultipleResults() {
        return multipleResults;
    }

    public List<SemEstimator> getMultipleResultList() {
        return multipleResultList;
    }

    public void setMultipleResultList(List<SemEstimator> multipleResultList) {
        this.multipleResultList = multipleResultList;
    }

    public Parameters getParams() {
        return params;
    }

    private SemOptimizer getOptimizer() {
        SemOptimizer optimizer;
        String type = getParams().getString("semOptimizerType", "Regression");

        if ("Regression".equals(type)) {
            SemOptimizer defaultOptimization = getDefaultOptimization();

            if (!(defaultOptimization instanceof SemOptimizerRegression)) {
                optimizer = defaultOptimization;
                type = getType(defaultOptimization);
                getParams().set("semOptimizerType", type);
            } else {
                optimizer = new SemOptimizerRegression();
            }
        } else if ("EM".equals(type)) {
            optimizer = new SemOptimizerEm();
        } else if ("Powell".equals(type)) {
            optimizer = new SemOptimizerPowell();
        } else if ("Random Search".equals(type)) {
            optimizer = new SemOptimizerScattershot();
        } else if ("RICF".equals(type)) {
            optimizer = new SemOptimizerRicf();
        } else if ("Powell".equals(type)) {
            optimizer = new SemOptimizerPowell();
        } else {
            if (linearSemPm != null) {
                optimizer = getDefaultOptimization();

                String _type = getType(optimizer);

                if (_type != null) {
                    getParams().set("semOptimizerType", _type);
                }
            } else {
                optimizer = null;
            }
        }

        return optimizer;
    }

    private String getType(SemOptimizer optimizer) {
        String _type = null;

        if (optimizer instanceof SemOptimizerRegression) {
            _type = "Regression";
        } else if (optimizer instanceof SemOptimizerEm) {
            _type = "EM";
        } else if (optimizer instanceof SemOptimizerPowell) {
            _type = "Powell";
        } else if (optimizer instanceof SemOptimizerScattershot) {
            _type = "Random Search";
        } else if (optimizer instanceof SemOptimizerRicf) {
            _type = "RICF";
        }

        return _type;
    }

    private boolean containsFixedParam(LinearSemPm linearSemPm) {
        return new LinearSemIm(linearSemPm).getNumFixedParams() > 0;
    }

    private static boolean containsCovarParam(LinearSemPm linearSemPm) {
        boolean containsCovarParam = false;
        List<Parameter> params = linearSemPm.getParameters();

        for (Parameter param : params) {
            if (param.getType() == ParamType.COVAR) {
                containsCovarParam = true;
                break;
            }
        }
        return containsCovarParam;
    }

    public ScoreType getScoreType() {
        return (ScoreType) params.get("scoreType", ScoreType.SemBic);
    }

    public void setScoreType(ScoreType scoreType) {
        params.set("scoreType", scoreType);
    }

    public void setNumRestarts(int numRestarts) {
        getParams().set("numRestarts", numRestarts);
    }

    public int getNumRestarts() {
        return getParams().getInt("numRestarts", 1);
    }

    private SemOptimizer getDefaultOptimization() {
        if (linearSemPm == null) {
            throw new NullPointerException(
                    "Sorry, I didn't see a SEM PM as parent to the estimator; perhaps the parents are wrong.");
        }

        boolean containsLatent = false;

        for (Node node : linearSemPm.getGraph().getNodes()) {
            if (node.getNodeType() == NodeType.LATENT) {
                containsLatent = true;
                break;
            }
        }

        SemOptimizer optimizer;

        if (containsFixedParam(linearSemPm) || linearSemPm.getGraph().existsDirectedCycle()
                || containsCovarParam(linearSemPm)) {
            optimizer = new SemOptimizerPowell();
        } else if (containsLatent) {
            optimizer = new SemOptimizerEm();
        } else {
            optimizer = new SemOptimizerRegression();
        }

        optimizer.setNumRestarts(getParams().getInt("numRestarts", 1));

        return optimizer;

//        optimizer.optimize(semIm);
//        this.semOptimizer = optimizer;
    }
}
