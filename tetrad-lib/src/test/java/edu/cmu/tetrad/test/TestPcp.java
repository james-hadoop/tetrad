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

package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.algorithm.continuous.dag.Lingam;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.PermutationGenerator;

/**
 * Tests the PC search.
 *
 * @author Joseph Ramsey
 */
public class TestPcp {

    private void test1() {
        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Algorithms algorithms = new Algorithms();
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new Pcp(new FisherZ()));
        algorithms.add(new Lingam());
//        algorithms.add(new Fask(new FisherZ()));
//        algorithms.add(new PcStableMax(new FisherZ()));

        Statistics statistics = new Statistics();
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ElapsedTime());

        Parameters parameters = new Parameters();
        parameters.set(Params.NUM_RUNS, 1);
        parameters.set(Params.NUM_MEASURES, 10);
        parameters.set(Params.FDR_Q, 0.05);
        parameters.set(Params.COEF_LOW, 0.2);
        parameters.set(Params.COEF_HIGH, 0.7);
        parameters.set(Params.ALPHA, 0.01);
        parameters.set(Params.SAMPLE_SIZE, 1000);
        parameters.set(Params.AVG_DEGREE, 2);
        parameters.set(Params.COLLIDER_DISCOVERY_RULE, 3);

        parameters.set(Params.ERRORS_NORMAL, false);
        parameters.set(Params.PENALTY_DISCOUNT, 4);
        parameters.set(Params.SKEW_EDGE_THRESHOLD, 0.3);

        Comparison comparison = new Comparison();

        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);
        comparison.compareFromSimulations("comparison3", simulations, algorithms,
                statistics, parameters);
    }

    public void test2() {
        PermutationGenerator.testPrint(5);
    }

    public static  void main(String...args) {
        new TestPcp().test1();
    }
}




