package edu.cmu.tetradapp.model;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;

/**
 * @author jdramsey
 */
public interface GraphSettable {
    Parameters getParameters();
    void setGraph(Graph newValue);

    void setModelIndex(int index);

    int getNumModels();

    String getModelSourceName();

    int getModelIndex();
}
