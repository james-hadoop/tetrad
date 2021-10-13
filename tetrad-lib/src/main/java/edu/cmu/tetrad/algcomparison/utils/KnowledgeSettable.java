package edu.cmu.tetrad.algcomparison.utils;

import edu.cmu.tetrad.data.IKnowledge;

/**
 * Stores a knowledge object.
 *
 * @author dmalinsky
 */
public interface KnowledgeSettable {

    /**
     * @return a knowledge object.
     */
    IKnowledge getKnowledge();

    /**
     * Sets a knowledge object.
     */
    void setKnowledge(IKnowledge knowledge);
}
