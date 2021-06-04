package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.List;

public interface FastForward {
    Graph search(List<Node> order);
    double score(List<Node> order);

    boolean isAssociated(Node w, Node v);
}
