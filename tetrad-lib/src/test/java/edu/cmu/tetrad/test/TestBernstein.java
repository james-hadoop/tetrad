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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.*;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore;
import edu.cmu.tetrad.algcomparison.simulation.LinearSemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import org.junit.Test;

/**
 * Tests to make sure the DelimiterType enumeration hasn't been tampered with.
 *
 * @author Joseph Ramsey
 */
@SuppressWarnings("ALL")
public final class TestBernstein {

    @Test
    public void test1() {
        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 10);
        params.set(Params.NUM_LATENTS, 3);
        params.set(Params.AVG_DEGREE, 3);
        params.set(Params.SAMPLE_SIZE, 1000, 10000);
        params.set(Params.NUM_RUNS, 10);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.ALPHA, 0.001);//1e-10, 1e-5, 0.001, 0.01, 0.2, 0.3, 0.4, 0.5);
        params.set(Params.PENALTY_DISCOUNT, 1, 2, 3, 4);
        params.set(Params.COEF_LOW, 0.25);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VAR_LOW, 1);
        params.set(Params.VAR_HIGH, 1);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.BREAK_TIES, true);
        params.set(Params.GRASP_USE_SCORE, true);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new Fci(new FisherZ()));
        algorithms.add(new FciMax(new FisherZ()));
        algorithms.add(new Rfci(new FisherZ()));
        algorithms.add(new Gfci(new LinearGaussianBicScore(), new FisherZ()));
        algorithms.add(new GRASPFCI(new LinearGaussianBicScore(), new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new LinearSemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.ALPHA));
        statistics.add(new ParameterColumn(Params.PENALTY_DISCOUNT));
        statistics.add(new ParameterColumn(Params.GRASP_USE_SCORE));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new AdjacencyTPR());
        statistics.add(new AdjacencyFPR());
        statistics.add(new SHD());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setSaveData(true);
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.True_PAG);

        comparison.compareFromSimulations(
                "/Users/josephramsey/tetrad/bernstein",
                "Comparison1.txt",
                simulations, algorithms, statistics, params);
    }

}





