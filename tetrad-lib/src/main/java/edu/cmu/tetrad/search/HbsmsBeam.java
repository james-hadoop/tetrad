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
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.sem.*;

/**
 * Heuristic Best Significant Model Search using a beam search.
 *
 * @author Joseph Ramsey
 */
public final class HbsmsBeam implements Hbsms {
    private final Graph initialGraph;
    private final Scorer scorer;
    private final CovarianceMatrix cov;
    private ForwardScoreSearch gis;
    private final IKnowledge knowledge = new Knowledge2();
    private SemIm newSemIm;

    public HbsmsBeam(Graph graph, DataSet data, IKnowledge knowledge) {
        graph = GraphUtils.replaceNodes(graph, data.getVariables());
        CovarianceMatrix cov = new CovarianceMatrix(data);
        this.cov = cov;
        this.scorer = new FmlBicScorer(cov);
//        this.scorer = new FgesScorer(data);
        gis = new ForwardScoreSearch(this.scorer);
//        gis.setKnowledge(knowledge);
        this.initialGraph = new EdgeListGraph(graph);
    }

    public Graph search() {


        Graph best = gis.search(scorer.getVariables());
        SemPm pm = new SemPm(best);
        SemEstimator est = new SemEstimator(cov, pm);
        this.newSemIm = est.estimate();
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

    public double scoreGraph(Graph graph, Scorer scorer) {
        Graph dag = new EdgeListGraph(graph);
        return scorer.score(dag);
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
//        this.gis.setKnowledge(knowledge);
    }
}


