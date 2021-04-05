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

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.sem.DagScorer;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.sem.SemIm;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Best Fit Finder using the GES algorithm.
 * </p>
 * Improves the P value of a SEM IM by adding, removing, or reversing single edges.
 *
 * @author Joseph Ramsey
 */

public final class HbsmsGes implements Hbsms {
    private final DataSet data;
    private IKnowledge knowledge = new Knowledge2();
    private final Graph graph;
    private double alpha = 0.05;
    private final Set<GraphWithPValue> significantModels = new HashSet<>();
    private SemIm originalSemIm;
    private SemIm newSemIm;
    private final Scorer scorer;
    private final Scorer scorer2;

    public HbsmsGes(Graph graph, DataSet data) {
        if (graph == null) throw new NullPointerException("Graph not specified.");
        this.data = DataUtils.getBootstrapSample(data, (int) (data.getNumRows() * 0.8));

        DagInPatternIterator iterator = new DagInPatternIterator(graph, getKnowledge(), true,
                true);
        graph = iterator.next();
        graph = SearchGraphUtils.patternForDag(graph);

        if (GraphUtils.containsBidirectedEdge(graph)) {
            throw new IllegalArgumentException("Contains bidirected edge.");
        }

        this.graph = graph;

        List<DataSet> split = DataUtils.split(data, 0.8);

        this.scorer = new DagScorer(split.get(0));
        this.scorer2 = new DagScorer(split.get(1));
    }

    private void saveModelIfSignificant(Graph graph) {
        double pValue = scoreGraph(graph, scorer).getPValue();

        if (pValue > alpha) {
            getSignificantModels().add(new GraphWithPValue(graph, pValue));
        }
    }

    public static class GraphWithPValue {
        private final Graph graph;
        private final double pValue;

        public GraphWithPValue(Graph graph, double pValue) {
            this.graph = graph;
            this.pValue = pValue;
        }

        public Graph getGraph() {
            return graph;
        }

        public double getPValue() {
            return pValue;
        }

        public int hashCode() {
            return 17 * graph.hashCode();
        }

        public boolean equals(Object o) {
            if (o == null) return false;
            if (!(o instanceof GraphWithPValue)) return false;
            GraphWithPValue p = (GraphWithPValue) o;
            return (p.graph.equals(graph));
        }
    }

    public Score scoreGraph(Graph graph, Scorer scorer) {
        Graph dag = SearchGraphUtils.dagFromPattern(graph, getKnowledge());

        if (dag == null) {
            return Score.negativeInfinity();
        }

        scorer.score(dag);
        return new Score(this.scorer);
    }

    public Graph getGraph() {
        return graph;
    }

    public SemIm getOriginalSemIm() {
        return originalSemIm;
    }

    public SemIm getNewSemIm() {
        return newSemIm;
    }

    public void setHighPValueAlpha(double highPValueAlpha) {
    }

    public Score scoreDag(Graph dag) {
        scorer.score(dag);
        return new Score(scorer);
    }

    public Graph search() {
        Score score1 = scoreGraph(getGraph(), scorer);
        double score = score1.getScore();
        System.out.println(getGraph());
        System.out.println(score);

        originalSemIm = score1.getEstimatedSem();

        Fges fges = new Fges(new SemBicScore(data));
        fges.setInitialGraph(graph);
        fges.setKnowledge(knowledge);

        Graph model = SearchGraphUtils.dagFromPattern(fges.search());

        saveModelIfSignificant(model);

        Score _score = scoreGraph(getGraph(), scorer2);
        newSemIm = _score.getEstimatedSem();

        return new EdgeListGraph(getGraph());
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setBeamWidth(int beamWidth) {
//        if (beamWidth < 1) throw new IllegalArgumentException();
        // Do nothing. We don't care about beam width.
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public Set<GraphWithPValue> getSignificantModels() {
        return significantModels;
    }

    public static class Score {
        private final Scorer scorer;
        private final double pValue;
        private final double fml;
        private final double chisq;
        private double bic;
        //        private double aic;
        private int dof;

        public Score(Scorer scorer) {
            this.scorer = scorer;
            this.pValue = scorer.getPValue();
            this.fml = scorer.getFml();
            this.chisq = scorer.getChiSquare();
            this.bic = scorer.getBicScore();
            this.dof = scorer.getDof();
        }

        private Score() {
            this.scorer = null;
            this.pValue = 0.0;
            this.fml = Double.POSITIVE_INFINITY;
            this.chisq = 0.0;
        }

        public SemIm getEstimatedSem() {
            return scorer.getEstSem();
        }

        public double getPValue() {
            return pValue;
        }

        public double getScore() {
            return -bic;
        }

        public double getFml() {
            return fml;
        }

        public static Score negativeInfinity() {
            return new Score();
        }

        public int getDof() {
            return dof;
        }

        public double getChiSquare() {
            return chisq;
        }

        public double getBic() {
            return bic;
        }
    }
}


