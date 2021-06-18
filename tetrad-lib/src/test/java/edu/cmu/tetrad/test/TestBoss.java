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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.EbicScore;
import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.K3;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

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
        params.set(Params.NUM_MEASURES, 5);
        params.set(Params.AVG_DEGREE, 0, 1, 2, 3, 4);
        params.set(Params.SAMPLE_SIZE, 100000);
        params.set(Params.NUM_RUNS, 50);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0.25);
        params.set(Params.COEF_HIGH, 1.0);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
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
        params.set(Params.VERBOSE, false);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
        algorithms.add(new Fges(new SemBicScore()));
        algorithms.add(new PcAll(new FisherZ()));

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

    @Test
    public void testBoss3() {
//        RandomUtil.getInstance().setSeed(386829384L);

        Parameters params = new Parameters();
        params.set(Params.NUM_MEASURES, 60);
        params.set(Params.AVG_DEGREE, 5);
        params.set(Params.SAMPLE_SIZE, 1000);
        params.set(Params.NUM_RUNS, 1);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.COEF_LOW, 0.1);
        params.set(Params.COEF_HIGH, 0.9);
        params.set(Params.VERBOSE, false);
        params.set(Params.CACHE_SCORES, true);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new EbicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
//        algorithms.add(new Fges(new SemBicScore()));
//        algorithms.add(new PcAll(new FisherZ()));

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

    @Test
    public void test3() {
        int numNodes = 7;
        int numEdges = 10;
        int sampleSize = 1000;

        Graph graph = GraphUtils.randomGraph(numNodes, 0,
                numEdges, 100, 100, 100, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet data = im.simulateData(sampleSize, false);

        search(4, 6, data);

    }

    public void search(int i, int j, DataSet data) {
        if (!(j > i)) throw new IllegalArgumentException("Need to fetch a node " +
                "to the right of the prefix.");

        Score score = new edu.cmu.tetrad.search.SemBicScore(data);
        K3 k3 = new K3(score);
        LinkedList<Node> br = new LinkedList<>(data.getVariables());

        Node fetched = br.get(j);

        System.out.println("Fetched node = " + fetched);

        br.remove(fetched);
        br.add(0, fetched);

        System.out.println("br = " + br);

        double[] scores = new double[i];

        for (int p = 0; p < i; p++) {
            scores[p] = k3.getScore(new K3.Subproblem(br.get(p), new HashSet<>(getPrefix(br, p))));
        }

        double _score3 = 0;
        for (double v : scores) _score3 += v;
        double min = _score3;
        int bestw = -0;

        for (int w = 1; w < i; w++) {
            swap(w - 1, w, scores, br, i, k3);

            double _score2 = 0;
            for (double v : scores) _score2 += v;

            if (_score2 < min) {
                bestw = w;
                min = _score2;
            }
        }

        br.remove(fetched);
        br.add(bestw, fetched);

        System.out.println("br = " + br);
        System.out.println();
    }

    private void swap(int vv, int ww, double[] scores, LinkedList<Node> br, int prefixSize, K3 k3) {
        if (vv < 0) throw new IllegalArgumentException("vv < 0");
        if (vv >= prefixSize) throw new IllegalArgumentException("vv > prefixSize");
        if (!(ww >= vv)) throw new IllegalArgumentException("w < prefixSize");
        if (ww != vv + 1) throw new IllegalArgumentException("ww is not directly subsequent to vv");

        Node v = br.get(vv);

        br.remove(v);
        br.add(ww, v);

        for (int r = 0; r < prefixSize; r++) {
            if (r == vv || r == ww) {
                double score1 = k3.getScore(new K3.Subproblem(br.get(r), new HashSet<>(getPrefix(br, r))));
                scores[r] = score1;
            }
        }
    }

    public List<Node> getPrefix(List<Node> order, int i) {
        List<Node> prefix = new ArrayList<>();
        for (int j = 0; j < i; j++) {
            prefix.add(order.get(j));
        }
        return prefix;
    }
}





