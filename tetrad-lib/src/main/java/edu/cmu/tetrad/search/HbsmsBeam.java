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

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.sem.DagScorer;
import edu.cmu.tetrad.sem.Scorer;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.util.ChoiceGenerator;

import static java.lang.Math.log;

/**
 * Heuristic Best Significant Model Search using a beam search.
 *
 * @author Joseph Ramsey
 */
public final class HbsmsBeam implements Hbsms {
    private final Graph initialGraph;
    private final Scorer scorer;
    private GlobalImprovementSearch gis;
    private final IKnowledge knowledge = new Knowledge2();
    private SemIm newSemIm;

    public HbsmsBeam(Graph graph, DataSet data, IKnowledge knowledge) {
        gis = new GlobalImprovementSearch(data);
        gis.setKnowledge(knowledge);
        this.initialGraph = new EdgeListGraph(graph);
        CovarianceMatrix cov = new CovarianceMatrix(data);
        this.scorer = new DagScorer(cov);
    }

    public Graph search() {
        Graph best = new EdgeListGraph(initialGraph);
        best = gis.improve(best);
        Score score = scoreGraph(best, scorer);
        this.newSemIm = score.getEstimatedSem();
        return best;
    }

    @Override
    public SemIm getOriginalSemIm() {
        throw new UnsupportedOperationException();
    }

    @Override
    public SemIm getNewSemIm() {
        return newSemIm;
    }

    public Score scoreGraph(Graph graph, Scorer scorer) {
        Graph dag = new EdgeListGraph(graph);// SearchGraphUtils.dagFromPattern(graph, getKnowledge());
        scorer.score(dag);
        return new Score(scorer);
    }

    public void setAlpha(double alpha) {
        throw new UnsupportedOperationException();
    }

    public void setBeamWidth(int beamWidth) {
        throw new UnsupportedOperationException();
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.gis.setKnowledge(knowledge);
    }

    public static class Score {
        private final Scorer scorer;

        public Score(Scorer scorer) {
            this.scorer = scorer;
        }

        public SemIm getEstimatedSem() {
            return scorer.getEstSem();
        }

        public double getPValue() {
            return scorer.getPValue();
        }

        public double getScore() {
            double gamma = 1;

            return getChiSquare() - (scorer.getDof() * log(scorer.getSampleSize())
                    + 2 * gamma * ChoiceGenerator.logCombinations(scorer.getDof() - 1, scorer.getNumFreeParams()));
        }

        public double getFml() {
            return scorer.getFml();
        }

        public int getDof() {
            return scorer.getDof();
        }

        public double getChiSquare() {
            return (scorer.getSampleSize() - 1) * getFml();
        }

        public double getBic() {
            return getChiSquare() - scorer.getDof() * Math.log(scorer.getSampleSize());
        }
    }
}


