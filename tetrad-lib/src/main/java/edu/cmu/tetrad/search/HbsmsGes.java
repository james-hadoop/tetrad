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
import edu.cmu.tetrad.sem.*;

import java.util.HashSet;
import java.util.Set;

/**
 * Heuristic Best Significant Model Search using the GES algorithm.
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

    public HbsmsGes(Graph graph, DataSet data) {
        if (graph == null) throw new NullPointerException("Graph not specified.");
        this.data = data;
        this.graph = SearchGraphUtils.patternForDag(graph);
        this.scorer = new FmlBicScorer(data, 2);
    }

    public Graph search() {
        Score score1 = scoreGraph(getGraph(), scorer);

        Fges fges = new Fges(new SemBicScore(data));
        fges.setInitialGraph(graph);
        fges.setKnowledge(knowledge);

        Graph best = SearchGraphUtils.dagFromCpdag(fges.search(), knowledge);

        SemPm pm = new SemPm(best);
        SemEstimator est = new SemEstimator(data, pm);
        this.newSemIm = est.estimate();

        return new EdgeListGraph(getGraph());
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
        Graph dag = SearchGraphUtils.dagFromCpdag(graph, knowledge);

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

    public Score scoreDag(Graph dag) {
        scorer.score(dag);
        return new Score(scorer);
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
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public Set<GraphWithPValue> getSignificantModels() {
        return significantModels;
    }

    public static class Score {
        private final Scorer scorer;

        public Score(Scorer scorer) {
            this.scorer = scorer;
        }
    }
}


