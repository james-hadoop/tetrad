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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.BOSS;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.GSP;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;

/**
 * Tests to make sure the DelimiterType enumeration hasn't been tampered with.
 *
 * @author Joseph Ramsey
 */
public final class TestBoss {


    @Test
    public void testBoss() {
        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 10);
        params.set(Params.AVG_DEGREE, 1, 2, 4, 5, 6, 7, 8);
        params.set(Params.SAMPLE_SIZE, 100000);
        params.set(Params.NUM_RUNS, 50);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0.25);
        params.set(Params.COEF_HIGH, 1.0);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
        algorithms.add(new GSP(new SemBicScore()));
//        algorithms.add(new Fges(new SemBicScore()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.DEPTH));
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveData(false);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss", simulations, algorithms, statistics, params);
    }

    @Test
    public void testBoss2() {
        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 20);
        params.set(Params.AVG_DEGREE, 4);
        params.set(Params.SAMPLE_SIZE, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000,
                50000, 100000, 200000);
        params.set(Params.NUM_RUNS, 30);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0.2);
        params.set(Params.COEF_HIGH, 0.8);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
        algorithms.add(new Fges(new SemBicScore()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();
        statistics.add(new ParameterColumn(Params.SAMPLE_SIZE));
        statistics.add(new ParameterColumn(Params.AVG_DEGREE));
        statistics.add(new CorrectSkeleton());
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new SHD());
        statistics.add(new F1All());
        statistics.add(new ElapsedTime());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.Pattern_of_the_true_DAG);
        comparison.setSaveData(false);

        comparison.compareFromSimulations("/Users/josephramsey/tetrad/boss", simulations, algorithms, statistics, params);
    }
}





