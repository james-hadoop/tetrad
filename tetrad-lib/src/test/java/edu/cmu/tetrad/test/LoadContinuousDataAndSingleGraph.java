package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.simulation.Simulation;
import edu.cmu.tetrad.algcomparison.utils.HasParameterValues;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.data.reader.Delimiter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
public class LoadContinuousDataAndSingleGraph implements Simulation, HasParameterValues {
    static final long serialVersionUID = 23L;
    private final String path;
    private final String subdir;
    private Graph graph;
    private List<DataSet> dataSets = new ArrayList<>();
    private final List<String> usedParameters = new ArrayList<>();
    private final Parameters parametersValues = new Parameters();

    public LoadContinuousDataAndSingleGraph(String path, String subdir) {
        this.path = path;
        this.subdir = subdir;
        String structure = new File(path).getName();
        this.parametersValues.set("Structure", structure);
    }

    @Override
    public void createData(Parameters parameters, boolean newModel) {
//        if (!newModel && !dataSets.isEmpty()) return;
        if (!this.dataSets.isEmpty()) return;

        this.dataSets = new ArrayList<>();

        File dir = new File(this.path + "/" + this.subdir);

        if (dir.exists()) {
            File[] files = dir.listFiles();

            for (File file : files) {
                if (!file.getName().endsWith(".txt")) continue;
                System.out.println("Loading data from " + file.getAbsolutePath());

                try {
                    DataSet dataSet = DataUtils.loadContinuousData(file, "//", '\"',
                            "*", true, Delimiter.TAB);
                    this.dataSets.add(dataSet);

                    if (!(dataSet.isContinuous())) {
                        throw new IllegalArgumentException("Not a continuous data set: " + dataSet.getName());
                    }
                } catch (Exception e) {
                    System.out.println("Couldn't parse " + file.getAbsolutePath());
                }
            }
        }

        File dir2 = new File(this.path + "/graph");

        if (dir2.exists()) {
            File[] files = dir2.listFiles();

            if (files.length != 1) {
                throw new IllegalArgumentException("Expecting exactly one graph file.");
            }

            File file = files[0];

            System.out.println("Loading graph from " + file.getAbsolutePath());
            this.graph = GraphUtils.loadGraphTxt(file);

            GraphUtils.circleLayout(this.graph, 225, 200, 150);
        }

        if (parameters.get("numRuns") != null) {
            parameters.set("numRuns", parameters.get("numRuns"));
        } else {
            parameters.set("numRuns", this.dataSets.size());
        }

        System.out.println();
    }

    @Override
    public Graph getTrueGraph(int index) {
        return this.graph;
    }

    @Override
    public DataModel getDataModel(int index) {
        return this.dataSets.get(index);
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
        return this.usedParameters;
    }

    @Override
    public int getNumDataModels() {
        return this.dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public Parameters getParameterValues() {
        return this.parametersValues;
    }
}
