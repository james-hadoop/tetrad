package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.utils.HasParameterValues;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataReader;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import edu.pitt.dbmi.data.reader.tabular.MixedTabularDatasetFileReader;
import edu.pitt.dbmi.data.reader.tabular.MixedTabularDatasetReader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
public class LoadContinuousDataAndSingleGraph implements Simulation, HasParameterValues {
    static final long serialVersionUID = 23L;
    private String path;
    private Graph graph = null;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<String> usedParameters = new ArrayList<>();
    private Parameters parametersValues = new Parameters();

    public LoadContinuousDataAndSingleGraph(String path) {
        this.path = path;
        String structure = new File(path).getName();
        parametersValues.set("Structure", structure);
    }

    @Override
    public void createData(Parameters parameters) {
        this.dataSets = new ArrayList<>();

        File dir = new File(path + "/data");

        if (dir.exists()) {
            File[] files = dir.listFiles();

            File file = files[0];

            for (File _file : files) {
                if (_file.getName().startsWith(".")) continue;
                file = _file;
                break;
            }

            System.out.println("Loading data from " + file.getAbsolutePath());
            try {
                DataSet dataSet = (DataSet) DataConvertUtils.toDataModel(
                        new ContinuousTabularDatasetFileReader(file.toPath(), Delimiter.TAB).readInData());
                dataSets.add(dataSet);
            } catch (Exception e) {
                System.out.println("Couldn't parse " + file.getAbsolutePath());
            }
        }

        File dir2 = new File(path + "/graph");

        if (dir2.exists()) {
            File[] files = dir2.listFiles();

            File file = files[0];

            for (File _file : files) {
                if (_file.getName().startsWith(".")) continue;
                file = _file;
                break;
            }

            System.out.println("Loading graph from " + file.getAbsolutePath());
            this.graph = GraphUtils.loadGraphTxt(file);

            GraphUtils.circleLayout(this.graph, 225, 200, 150);
        }

//        for (Node node : dataSets.get(0).getVariables()) {
//            if (graph.getNode(node.getName()) == null) {
//                graph.removeNode(graph.getNode(node.getName()));
//            }
//        }
//
//        for (Node node : graph.getNodes()) {
//            if (dataSets.get(0).getVariable(node.getName()) == null) {
//                dataSets.get(0).removeColumn(dataSets.get(0).getVariable(node.getName()));
//            }
//        }

        if (parameters.get("numRuns") != null) {
            parameters.set("numRuns", parameters.get("numRuns"));
        } else {
            parameters.set("numRuns", dataSets.size());
        }

        System.out.println();
    }

    @Override
    public Graph getTrueGraph(int index) {
        return this.graph;
    }

    @Override
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }

    public String getDescription() {
        try {
            StringBuilder b = new StringBuilder();
            b.append("Load data sets and graphs from a directory.").append("\n\n");
            return b.toString();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public List<String> getParameters() {
        return usedParameters;
    }

    @Override
    public int getNumDataModels() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public Parameters getParameterValues() {
        return parametersValues;
    }
}
