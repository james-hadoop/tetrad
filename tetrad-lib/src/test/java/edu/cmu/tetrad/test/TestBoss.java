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
import edu.cmu.tetrad.algcomparison.simulation.SemSimulationTrueModel;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.GraphNode;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;
import org.junit.Test;

import java.util.*;

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
        params.set(Params.AVG_DEGREE, 0, 1, 2, 4, 5, 6, 7, 8);
        params.set(Params.SAMPLE_SIZE, 10000);
        params.set(Params.NUM_RUNS, 20);
        params.set(Params.RANDOMIZE_COLUMNS, false);
        params.set(Params.PENALTY_DISCOUNT, 2);
        params.set(Params.COEF_LOW, 0.25);
        params.set(Params.COEF_HIGH, 1.0);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);


        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new SemBicScore()));
//        algorithms.add(new GSP(new SemBicScore()));
        algorithms.add(new Fges(new SemBicScore()));

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
        params.set(Params.NUM_MEASURES, 100);
        params.set(Params.AVG_DEGREE, 10);
        params.set(Params.SAMPLE_SIZE, 1000);
        params.set(Params.NUM_RUNS, 1);
        params.set(Params.RANDOMIZE_COLUMNS, true);
        params.set(Params.PENALTY_DISCOUNT, 1);
        params.set(Params.COEF_LOW, 0);
        params.set(Params.COEF_HIGH, 1);
        params.set(Params.VERBOSE, false);
        params.set(Params.COLLIDER_DISCOVERY_RULE, 2, 3);
        params.set(Params.CACHE_SCORES, true);
        params.set(Params.NUM_STARTS, 1);
        params.set(Params.CACHE_SCORES, false);

        Algorithms algorithms = new Algorithms();
        algorithms.add(new BOSS(new EbicScore()));
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

    @Test
    public void testFromIndependence() {
        Node x = new GraphNode("X");
        Node y = new GraphNode("Y");
        Node z = new GraphNode("Z");
        Node w = new GraphNode("W");

        IndependenceFacts facts = new IndependenceFacts();

        facts.add(new IndependenceFact(x, z, list(y)));
        facts.add(new IndependenceFact(y, w, list(x, z)));
        facts.add(new IndependenceFact(z, x, list(y)));
//        facts.add(new IndependenceFact(x, w, list(y))); // unfaithful.

//        IndTestIndependenceFacts test = new IndTestIndependenceFacts(facts);
//
//        Pc pc = new Pc(test);
//        pc.setVerbose(true);
//        System.out.println(pc.search());

        ScoreIndependenceFacts score = new ScoreIndependenceFacts(facts);
        TeyssierScorer scorer = new TeyssierScorer(score);
        scorer.score(score.getVariables());

//        System.out.println("Initial score = " + -scorer.score());

        BestOrderScoreSearch boss = new BestOrderScoreSearch(score);
        boss.setCachingScores(false);
        boss.bossSearchPromotion(scorer);
    }

    private List<Node> list(Node...nodes) {
        List<Node> list = new ArrayList<>();
        Collections.addAll(list, nodes);
        return list;
    }
}





