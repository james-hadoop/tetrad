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

package edu.cmu.tetrad.study;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.algorithm.multi.Fask;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.*;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Cpc;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.score.DSeparationScore;
import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.LinearFisherModel;
import edu.cmu.tetrad.algcomparison.simulation.LoadContinuousDataAndSingleGraph;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;

/**
 * An example script to simulate data and run a comparison analysis on it.
 *
 * @author jdramsey
 */
public class ExampleCompareSimulationContinuousPag {
    enum Type {LinearGaussian, LinearNongaussian}

    public static void main1(String... args) {
        Parameters parameters = new Parameters();
        ExampleCompareSimulationContinuousPag.Type type = Type.LinearNongaussian;

        parameters.set("numRuns", 5);
        parameters.set("numMeasures", 10);
        parameters.set("numLatents", 0);
        parameters.set("avgDegree", 4);
        parameters.set("sampleSize", 1000); // This varies.
        parameters.set("differentGraphs", true);

        parameters.set("alpha", 0.01);

        parameters.set("coefLow", 0.3);
        parameters.set("coefHigh", 0.9);

        parameters.set("varLow", 0.3);
        parameters.set("varHigh", 0.9);

        parameters.set("maxDegree", 12);
//        parameters.set("maxIndegree", 3);

        parameters.set("completeRuleSetUsed", true);

        parameters.set("depth", 8);

        parameters.set("penaltyDiscount", 2, 4, 6);
        parameters.set("faithfulnessAssumed", true, false);
        parameters.set("symmetricFirstStep", true, false);

        parameters.set("twoCycleAlpha", 0);
        parameters.set("extraEdgeThreshold", 1.0);

        parameters.set("useFasAdjacencies", true);
        parameters.set("useCorrDiffAdjacencies", false);

        Statistics statistics = new Statistics();

//        statistics.add(new ParameterColumn("numMeasures"));
        statistics.add(new ParameterColumn("avgDegree"));
//        statistics.add(new ParameterColumn("maxDegree"));
//        statistics.add(new ParameterColumn("sampleSize"));
        statistics.add(new ParameterColumn("alpha"));
        statistics.add(new ParameterColumn("penaltyDiscount"));
//        statistics.add(new ParameterColumn("maxDegree"));

        statistics.add(new ParameterColumn("faithfulnessAssumed"));
        statistics.add(new ParameterColumn("symmetricFirstStep"));
        statistics.add(new ParameterColumn("maxDegree"));


        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecall());
//        statistics.add(new MathewsCorrAdj());
//        statistics.add(new MathewsCorrArrow());
//        statistics.add(new F1Adj());
//        statistics.add(new F1Arrow());
//        statistics.add(new SHD());
        statistics.add(new ElapsedTime());

        statistics.setWeight("AP", 1.0);
        statistics.setWeight("AR", 1.0);
        statistics.setWeight("AHPC", 1.0);
//        statistics.setWeight("AHR", 1.0);
//        statistics.setWeight("SHD", 1.0);

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new PcAll(new FisherZ()));
//        algorithms.add(new R3(new FAS(new FisherZ())));
//        algorithms.add(new Fask(new FAS(new FisherZ())));
////
        algorithms.add(new Fci(new SemBicTest()));
        algorithms.add(new Fci(new FisherZ()));
        algorithms.add(new Rfci(new FisherZ()));
        algorithms.add(new FciMax(new FisherZ()));
        algorithms.add(new CFCI(new FisherZ()));
        algorithms.add(new Gfci(new FisherZ(), new SemBicScore()));
        algorithms.add(new FciNg(new FisherZ()));
        algorithms.add(new GfciNg(new FisherZ(), new SemBicScore()));


        Simulations simulations = new Simulations();

        if (type == ExampleCompareSimulationContinuousPag.Type.LinearGaussian) {
            simulations.add(new SemSimulation(new RandomForward()));
        } else if (type == ExampleCompareSimulationContinuousPag.Type.LinearNongaussian) {
            simulations.add(new LinearFisherModel(new RandomForward()));
            parameters.set("errorsNormal", false);
        }

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
//        comparison.setSortByUtility(true);
//        comparison.setShowUtilities(true);
//        comparison.setParallelized(true);

        comparison.setSaveGraphs(true);
        comparison.setSavePags(true);
        comparison.setSavePatterns(true);
//        comparison.setSaveTrueDags(true);

        comparison.setComparisonGraph(Comparison.ComparisonGraph.PAG_of_the_true_DAG);

        comparison.compareFromSimulations("comparison.pag.simulations", simulations, "comparison", algorithms, statistics, parameters);
    }

    public static void main3(String... args) {
        Parameters parameters = new Parameters();
        int sampleSize = 1000;
        int maxDegree = 12;
        ExampleCompareSimulationContinuousPag.Type type = Type.LinearNongaussian;

        parameters.set("numRuns", 1);
//        parameters.set("sampleSize", sampleSize); // This varies.
        parameters.set("differentGraphs", true);

        parameters.set("penaltyDiscount", 1, 2, 4, 6); // tookout 1

        parameters.set("coefLow", 0.3);
        parameters.set("coefHigh", 0.9);

        parameters.set("varLow", 0.3);
        parameters.set("varHigh", 0.9);

        parameters.set("maxDegree", -1);
//        parameters.set("maxIndegree", 3);

        parameters.set("completeRuleSetUsed", false);

        parameters.set("alpha", 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4/*, 0.5*/);//, 0.6, 0.7, 0.8, 0.9, 1.0);

        parameters.set("depth", -1);
        parameters.set("twoCycleAlpha", 0.000000);
//        parameters.set("extraEdgeThreshold", 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, .6, .7, 0.8, 0.9, 1.0);
        parameters.set("extraEdgeThreshold", 1.0);

        parameters.set("useFasAdjacencies", true);
        parameters.set("useCorrDiffAdjacencies", true);

//        parameters.set("alpha");
        parameters.set("kernelRegressionSampleSize", 50, 100, 200);

        parameters.set("verbose", false);


        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn("avgDegree"));
        statistics.add(new ParameterColumn("alpha"));
        statistics.add(new ParameterColumn("penaltyDiscount"));
        statistics.add(new ParameterColumn("extraEdgeThreshold"));
        statistics.add(new ParameterColumn("kernelRegressionSampleSize"));


        statistics.add(new NumEdgeTrueGraph());
        statistics.add(new NumEdgeEstGraph());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new AncestorPrecision());
        statistics.add(new AncestorRecall());
        statistics.add(new DescentantsOfMTorTp());
        statistics.add(new DescentantsOfMTorFp());
        statistics.add(new DescentantsOfMTorTn());
        statistics.add(new DescentantsOfMTorFn());
        statistics.add(new DescentantsOfMTorFpr());
        statistics.add(new DescentantsOfMTorTpr());
        statistics.add(new DescentantsOfMTorRatioTprToFpr());

//        statistics.setWeight("mtorDescFpr", 1.0);
//        statistics.setWeight("mtorDescTpr", 1.0);
        statistics.setWeight("mtorDescRatioTprToFpr", 1.0);
//        statistics.setWeight("NumEdgesEst", 1.0);

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new Pc(new FisherZ()));
        algorithms.add(new Cpc(new FisherZ()));
        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new Fci(new FisherZ()));
//        algorithms.add(new Rfci(new FisherZ()));
        algorithms.add(new FciMax(new FisherZ()));
        algorithms.add(new CFCI(new FisherZ()));
//        algorithms.add(new Gfci(new FisherZ(), new SemBicScore()));
        algorithms.add(new Fask(new FisherZ()));
//        algorithms.add(new FciNg(new FisherZ()));
        algorithms.add(new GfciNg(new FisherZ(), new SemBicScore()));


        Simulations simulations = new Simulations();

//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.all.cyclic.ground.truth"));

        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.uaac812.cyclic.ground.truth"));
//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.mcf7.cyclic.ground.truth"));
//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.bt20.cyclic.ground.truth"));
//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.bt549.cyclic.ground.truth"));

//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/hill.egf.cyclic.ground.truth"));

//        simulations.add(new LoadContinuousDataAndSingleGraph("/Users/user/Box/data/4cellLineData/sessions.for.algcomparison/sachs"));

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);

        comparison.setSaveGraphs(true);
        comparison.setSavePags(false);
        comparison.setSavePatterns(false);

        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        comparison.compareFromSimulations("comparison.hill", simulations, "comparison.hill", algorithms, statistics, parameters);

    }

    public static void main(String... args) {
        Parameters parameters = new Parameters();

        parameters.set("numRuns", 1);
        parameters.set("numMeasures", 20);
        parameters.set("avgDegree", 4);
        parameters.set("sampleSize", 50000);
        parameters.set("differentGraphs", true);
        parameters.set("coefLow", 0.3);
        parameters.set("coefHigh", 0.9);
        parameters.set("varLow", 1);
        parameters.set("varHigh", 3);
        parameters.set("includePositiveCoefs", true);
        parameters.set("includeNegativeCoefs", false);
        parameters.set("depth", -1);
//        parameters.set("faithfulnessAssumed", true, false);
//        parameters.set("symmetricFirstStep", true, false);
        parameters.set("maxDegree", 100);

        // False for non-Gausian.
        parameters.set("errorsNormal", true);
        parameters.set("twoCycleAlpha", 0);
        parameters.set("extraEdgeThreshold", 1.0);

        parameters.set("useFasAdjacencies", true);
        parameters.set("useCorrDiffAdjacencies", false);

        parameters.set("stableFAS", true);
        parameters.set("concurrentFAS", false);
        parameters.set("colliderDiscoveryRule", 2);
        parameters.set("conflictRule", 3);
        parameters.set("depth", -1);
        parameters.set("useMaxPOrientationHeuristic", true);
//        parameters.set("maxPOrientationMaxPathLength");

        parameters.set("structurePrior", 0);


        Statistics statistics = new Statistics();

        parameters.set("penaltyDiscount", 1);
        parameters.set("alpha", 0.00001);

        statistics.add(new ParameterColumn("avgDegree"));
        statistics.add(new ParameterColumn("alpha"));
        statistics.add(new ParameterColumn("penaltyDiscount"));
//        statistics.add(new ParameterColumn("structurePrior"));
//
//        statistics.add(new ParameterColumn("faithfulnessAssumed"));
//        statistics.add(new ParameterColumn("symmetricFirstStep"));
//        statistics.add(new ParameterColumn("maxDegree"));

        statistics.add(new NumEdgeEstGraph());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new Fask(new FisherZ()));
//        algorithms.add(new PcAll(new FisherZ()));
        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new Fask(new FAS(new FisherZ())));
//        algorithms.add(new R3(new FAS(new FisherZ())));

        Simulations simulations = new Simulations();

        simulations.add(new LinearFisherModel(new RandomForward()));

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(false);
        comparison.setShowSimulationIndices(false);

        comparison.setSaveGraphs(false);
        comparison.setSavePags(false);
        comparison.setSavePatterns(false);

        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);

        comparison.compareFromSimulations("comparison.fges", simulations, "comparison", algorithms, statistics, parameters);
    }
}




