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

import edu.cmu.tetrad.data.Clusters;
import edu.cmu.tetrad.graph.*;

import java.rmi.MarshalledObject;
import java.util.*;

/**
 * Bron-Kerbosch Degeneracy.
 *
 * @author Bryan Andrews
 */
public class BKD {

    private Graph graph;

    private List maximalCliques;

    public BKD(Graph graph){
        this.graph = graph;
        this.maximalCliques = new ArrayList();
    }

    public List getCliques(boolean usesDegeneracyOrder) {
        Set<Node> P;
        Set<Node> R;
        Set<Node> X;
        if (usesDegeneracyOrder){
            List<Node> order = getDegeneracyOrder();
            Set<Node> before = new HashSet<>();
            Set<Node> after = new HashSet<>(order);
            for (int i = 0; i < order.size(); i++) {
                Node v = order.get(i);
                after.remove(v);
                P = new HashSet<>(this.graph.getAdjacentNodes(v));
                R = new HashSet<>();
                X = new HashSet<>(this.graph.getAdjacentNodes(v));
                P.retainAll(after);
                R.add(v);
                X.retainAll(before);
                BronKerboschPivot(P, R, X);
                before.add(v);
            }
        } else {
            P = new HashSet<>(this.graph.getNodes());
            R = new HashSet<>();
            X = new HashSet<>();
            BronKerboschPivot(P, R, X);
        }

        return this.maximalCliques;
    }

    private void BronKerboschPivot(Set<Node> P, Set<Node> R, Set<Node> X) {
        if (P.isEmpty() && X.isEmpty()) {
            this.maximalCliques.add(R);
        } else {
            List<Node> V;
            Set<Node> P_;
            Set<Node> R_;
            Set<Node> X_;
            Node u = getPivot(P, X);
            List<Node> PU = new ArrayList<>(P);
            PU.removeAll(this.graph.getAdjacentNodes(u));
            for (Node v : PU) {
                V = this.graph.getAdjacentNodes(v);
                P_ = new HashSet<>(P);
                R_ = new HashSet<>(R);
                X_ = new HashSet<>(X);
                P_.retainAll(V);
                R_.add(v);
                X_.retainAll(V);
                BronKerboschPivot(P_,R_,X_);
                P.remove(v);
                X.add(v);
            }
        }
    }

    private Node getPivot(Set<Node> P, Set<Node> X){
        Set<Node> V = new HashSet();
        V.addAll(P);
        V.addAll(X);
        int max = -1;
        Node u = null;
        List<Node> U;
        for (Node v : V) {
            U = this.graph.getAdjacentNodes(v);
            U.retainAll(P);
            if (U.size() > max) {
                max = U.size();
                u = v;
            }
        }

        return u;
    }

    private List<Node> getDegeneracyOrder() {
        Graph subgraph = new EdgeListGraph(this.graph);
        List<Node> nodes = subgraph.getNodes();
        List<Node> order = new ArrayList<>();
        int i = 0;
        int min = nodes.size();
        Node remove = null;
        while (subgraph.getNumNodes() > 0) {
            if (i >= nodes.size()) {
                if (min > 0) {
                    order.add(remove);
                    nodes.remove(remove);
                    subgraph.removeNode(remove);
                }
                i = 0;
                min = nodes.size();
                remove = null;
            }
            Node node = nodes.get(i);
            if (subgraph.getDegree(node) <= 1) {
                    order.add(node);
                    nodes.remove(node);
                    subgraph.removeNode(node);
                    min = 0;
            } else {
                if (subgraph.getDegree(node) < min) {
                    min = subgraph.getDegree(node);
                    remove = node;
                }
                i++;
            }
        }

        return order;
    }

}





