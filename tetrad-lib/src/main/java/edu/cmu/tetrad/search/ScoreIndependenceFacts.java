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

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IndependenceFacts;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DepthChoiceGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Checks conditional independence against a list of conditional independence facts, manually entered.
 *
 * @author Joseph Ramsey
 * @see ChiSquareTest
 */
public final class ScoreIndependenceFacts implements Score {

    private final IndependenceFacts facts;
    private boolean verbose = false;

    public ScoreIndependenceFacts(IndependenceFacts facts) {
        this.facts = facts;
    }

    public IndependenceTest indTestSubset(List<Node> vars) {
        throw new UnsupportedOperationException();
    }

    @Override
    public double localScore(int i, int... parents) {
        List<Node> nodes = facts.getVariables();

        Node y = nodes.get(i);

        List<Node> z = new ArrayList<>();
        for (int pi : parents) z.add(nodes.get(pi));

//        System.out.println("Scoring " + y + " | " + z);

        double score = 0;

        z.remove(y);

        for (Node p : new ArrayList<>(z)) {
            List<Node> _z = new ArrayList<>(z);
            _z.remove(p);

            boolean indep = false;

            DepthChoiceGenerator gen = new DepthChoiceGenerator(_z.size(), _z.size());
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> cond = GraphUtils.asList(choice, _z);

                if (facts.isIndependent(y, p, cond)) {
                    indep = true;
                    break;
                }
            }

            if (!indep) {
                score++;
            }
        }

//        score -= parents.length * 0.1;

        System.out.println("SCORE " + y + " | " + z + " = " + score);

        return (score + 1) / (parents.length + 1);
    }

    @Override
    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    private static int[] append(int[] z, int x) {
        int[] _z = Arrays.copyOf(z, z.length + 1);
        _z[z.length] = x;
        return _z;
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScoreDiff(x, y, new int[0]);
    }

    @Override
    public double localScore(int i, int parent) {
        return localScore(i, new int[]{parent});
    }

    @Override
    public double localScore(int i) {
        return localScore(i, new int[0]);
    }

    public List<Node> getVariables() {
        return facts.getVariables();
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return false;
    }

    public Node getVariable(String name) {
        if (name == null) throw new NullPointerException();

        List<Node> variables = facts.getVariables();

        for (Node node : variables) {
            if (name.equals(node.getName())) {
                return node;
            }
        }

        return null;
    }

    @Override
    public int getMaxDegree() {
        return 0;
    }

    public List<String> getVariableNames() {
        return facts.getVariableNames();
    }

    public boolean determines(List<Node> z, Node y) {
        return false;
    }

    public double getAlpha() {
        return Double.NaN;
    }

    public void setAlpha(double alpha) {
        throw new UnsupportedOperationException();
    }

    public DataModel getData() {
        return facts;
    }

    @Override
    public int getSampleSize() {
        return 0;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
}





