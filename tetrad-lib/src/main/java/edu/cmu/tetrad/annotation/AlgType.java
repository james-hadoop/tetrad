package edu.cmu.tetrad.annotation;

/**
 * Author : Jeremy Espino MD Created 6/30/17 10:36 AM
 */
public enum AlgType {
    forbid_latent_common_causes, // PC_All, PcStableMax, FGES, IMaGES_Discrete, IMaGES_Continuous, FANG, EFANG
    allow_latent_common_causes, // FCI, RFCI, GFCI, SVARFCI, SvarGFCI
    /*DAG, */
    search_for_Markov_blankets, // FGES-MB, MBFS
    produce_undirected_graphs, // FAS, MGM, GLASSO
    orient_pairwise, // R3, RSkew, Skew
    search_for_structure_over_latents // BPC, FOFC, FTFC
}
