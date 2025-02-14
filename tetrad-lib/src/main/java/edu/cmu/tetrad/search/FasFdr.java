///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015, 2022 by Peter Spirtes, Richard        //
// Scheines, Joseph Ramsey, and Clark Glymour.                               //
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

import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.*;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static java.lang.Math.sqrt;

/**
 * Implements the "fast adjacency search" used in several causal algorithm in this package. In the fast adjacency
 * search, at a given stage of the search, an edge X*-*Y is removed from the graph if X _||_ Y | S, where S is a subset
 * of size d either of adj(X) or of adj(Y), where d is the depth of the search. The fast adjacency search performs this
 * procedure for each pair of adjacent edges in the graph and for each depth d = 0, 1, 2, ..., d1, where d1 is either
 * the maximum depth or else the first such depth at which no edges can be removed. The interpretation of this adjacency
 * search is different for different algorithm, depending on the assumptions of the algorithm. A mapping from {x, y} to
 * S({x, y}) is returned for edges x *-* y that have been removed.
 * <p>
 * This variant does each depth twice, gathering up the p values in the first round, using FDR to estimate a cutoff
 * for acceptance, and rerunning using the specified cutoff.
 *
 * @author Joseph Ramsey.
 */
public class FasFdr implements IFas {

    private final Matrix cov;

    private final double alpha;
    /**
     * The search graph. It is assumed going in that all of the true adjacencies of x are in this graph for every node
     * x. It is hoped (i.e. true in the large sample limit) that true adjacencies are never removed.
     */
    private final Graph graph;

    /**
     * The independence test. This should be appropriate to the types
     */
    private final IndependenceTest test;

    /**
     * Specification of which edges are forbidden or required.
     */
    private IKnowledge knowledge = new Knowledge2();

    /**
     * The maximum number of variables conditioned on in any conditional independence test. If the depth is -1, it will
     * be taken to be the maximum value, which is 1000. Otherwise, it should be set to a non-negative integer.
     */
    private int depth = 1000;

    /**
     * The number of independence tests.
     */
    private final int numIndependenceTests;


    /**
     * The logger, by default the empty logger.
     */
    private final TetradLogger logger = TetradLogger.getInstance();

    /**
     * The number of dependence judgements. Temporary.
     */
    private final int numDependenceJudgement;

    /**
     * The sepsets found during the search.
     */
    private SepsetMap sepset = new SepsetMap();

    //    private List<Double> pValues = new ArrayList<Double>();

    private final NumberFormat nf = new DecimalFormat("0.00E0");

    /**
     * True iff verbose output should be printed.
     */
    private boolean verbose;
    private final List<Double> pValueList = new ArrayList<>();

    private PrintStream out = System.out;

    //==========================CONSTRUCTORS=============================//

    /**
     * Constructs a new FastAdjacencySearch.
     */
    public FasFdr(IndependenceTest test, int numIndependenceTests, int numDependenceJudgement) {
        this.graph = new EdgeListGraph(test.getVariables());
        this.test = test;
        this.alpha = test.getAlpha();
        this.cov = test.getCov().getMatrix();
        this.numIndependenceTests = numIndependenceTests;
        this.numDependenceJudgement = numDependenceJudgement;
    }

    //==========================PUBLIC METHODS===========================//

    /**
     * Discovers all adjacencies in data.  The procedure is to remove edges in the graph which connect pairs of
     * variables which are independent conditional on some other set of variables in the graph (the "sepset"). These are
     * removed in tiers.  First, edges which are independent conditional on zero other variables are removed, then edges
     * which are independent conditional on one other variable are removed, then two, then three, and so on, until no
     * more edges can be removed from the graph.  The edges which remain in the graph after this procedure are the
     * adjacencies in the data.
     *
     * @return a SepSet, which indicates which variables are independent conditional on which other variables
     */
    public Graph search() {
        this.logger.log("info", "Starting Fast Adjacency Search.");
        this.graph.removeEdges(this.graph.getEdges());

        this.sepset = new SepsetMap();

        int _depth = this.depth;

        if (_depth == -1) {
            _depth = 1000;
        }


        List<Node> nodes = this.graph.getNodes();
        Map<Node, Set<Node>> adjacencies = emptyGraph(nodes);

        searchICov(nodes, this.test, adjacencies, true);
        searchiCovAll(nodes, this.test, adjacencies);

        for (int d = 0; d <= _depth; d++) {
            searchAtDepth(nodes, this.test, adjacencies, d);

            if (!(freeDegree(nodes, adjacencies) > this.depth)) {
                break;
            }
        }

        this.pValueList.clear();

        for (int d = 0; d <= _depth; d++) {
            this.test.setAlpha(this.alpha);
            Map<Node, Set<Node>> _adjacencies = copy(adjacencies);
            searchAtDepth(nodes, this.test, adjacencies, d);
            double cutoff = StatUtils.fdrCutoff(this.test.getAlpha(), this.pValueList, false);
            adjacencies = _adjacencies;
            this.test.setAlpha(cutoff);
            boolean more = searchAtDepth(nodes, this.test, adjacencies, d);

            if (!more) {
                break;
            }
        }

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node x = nodes.get(i);
                Node y = nodes.get(j);

                if (adjacencies.get(x).contains(y)) {
                    this.graph.addUndirectedEdge(x, y);
                }
            }
        }

        this.logger.log("info", "Finishing Fast Adjacency Search.");

        return this.graph;
    }

    private Map<Node, Set<Node>> emptyGraph(List<Node> nodes) {
        Map<Node, Set<Node>> adjacencies = new HashMap<>();

        for (Node node : nodes) {
            adjacencies.put(node, new TreeSet<>());
        }
        return adjacencies;
    }

    private void searchiCovAll(List<Node> nodes, IndependenceTest test, Map<Node, Set<Node>> adjacencies) {
        boolean removed;

        do {
            removed = false;

            for (Node x : nodes) {
                List<Node> adjx = new ArrayList<>(adjacencies.get(x));

                for (Node y : adjx) {
                    if (!adjacencies.get(x).contains(y)) continue;
                    List<Node> adjy = new ArrayList<>(adjacencies.get(y));
                    List<Node> adj = new ArrayList<>(adjx);
                    for (Node node : adjy) if (!adj.contains(node)) adj.add(node);
                    removed = removed || searchICov(adj, test, adjacencies, false);
                }
            }
        } while (removed);
    }


    private Map<Node, Set<Node>> copy(Map<Node, Set<Node>> adjacencies) {
        Map<Node, Set<Node>> copy = new HashMap<>();

        for (Node node : adjacencies.keySet()) {
            copy.put(node, new HashSet<>(adjacencies.get(node)));
        }

        return copy;
    }

    public int getDepth() {
        return this.depth;
    }

    public void setDepth(int depth) {
        if (depth < -1) {
            throw new IllegalArgumentException(
                    "Depth must be -1 (unlimited) or >= 0.");
        }

        this.depth = depth;
    }

    public IKnowledge getKnowledge() {
        return this.knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) {
            throw new NullPointerException("Cannot set knowledge to null");
        }
        this.knowledge = knowledge;
    }

    //==============================PRIVATE METHODS======================/

    private boolean searchICov(List<Node> nodes, IndependenceTest test, Map<Node, Set<Node>> adjacencies,
                               boolean addDependencies) {
        if (nodes.size() < 2) return false;

        boolean removed = false;

        int[] n = new int[nodes.size()];
        List<Node> variables = test.getVariables();

        for (int i = 0; i < nodes.size(); i++) {
            n[i] = variables.indexOf(nodes.get(i));
        }

        Matrix inv = this.cov.getSelection(n, n).inverse();
        int sampleSize = test.getCov().getSampleSize();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node x = nodes.get(i);
                Node y = nodes.get(j);

                double r = -inv.get(i, j) / sqrt(inv.get(i, i) * inv.get(j, j));

                double fisherZ = sqrt(sampleSize - (nodes.size() - 2) - 3.0) *
                        0.5 * (Math.log(1.0 + r) - Math.log(1.0 - r));
                double pvalue = 2.0 * (1.0 - RandomUtil.getInstance().normalCdf(0, 1, Math.abs(fisherZ)));

                boolean independent = pvalue > test.getAlpha();

                if (addDependencies) {
                    if (independent) {
                        List<Node> theRest = new ArrayList<>();

                        for (Node node : nodes) {
                            if (node != x && node != y) theRest.add(node);
                        }

                        getSepsets().set(x, y, theRest);

                        removed = true;
                    } else if (!forbiddenEdge(x, y)) {
                        adjacencies.get(x).add(y);
                        adjacencies.get(y).add(x);

                    }
                } else {
                    if (independent) {
                        if (!adjacencies.get(x).contains(y)) continue;

                        List<Node> theRest = new ArrayList<>();

                        for (Node node : nodes) {
                            if (node != x && node != y) theRest.add(node);
                        }

                        adjacencies.get(x).remove(y);
                        adjacencies.get(y).remove(x);

                        getSepsets().set(x, y, theRest);

                        if (this.verbose) {
                            IndependenceResult result = this.test.checkIndependence(x, y, theRest);
                            this.out.println(x + " _||_ " + y + " | the rest" + " p = " +
                                    this.nf.format(result.getPValue()));
                        }

                        removed = true;
                    }
                }
            }
        }

        return removed;
    }

    private int freeDegree(List<Node> nodes, Map<Node, Set<Node>> adjacencies) {
        int max = 0;

        for (Node x : nodes) {
            Set<Node> opposites = adjacencies.get(x);

            for (Node y : opposites) {
                Set<Node> adjx = new HashSet<>(opposites);
                adjx.remove(y);

                if (adjx.size() > max) {
                    max = adjx.size();
                }
            }
        }

        return max;
    }

    private boolean forbiddenEdge(Node x, Node y) {
        String name1 = x.getName();
        String name2 = y.getName();

        if (this.knowledge.isForbidden(name1, name2) &&
                this.knowledge.isForbidden(name2, name1)) {
            this.logger.log("edgeRemoved", "Removed " + Edges.undirectedEdge(x, y) + " because it was " +
                    "forbidden by background knowledge.");

            return true;
        }

        return false;
    }

    private boolean searchAtDepth(List<Node> nodes, IndependenceTest test, Map<Node, Set<Node>> adjacencies, int depth) {
        int count = 0;

        for (Node x : nodes) {
            if (++count % 100 == 0) this.out.println("count " + count + " of " + nodes.size());

            List<Node> adjx = new ArrayList<>(adjacencies.get(x));

            EDGE:
            for (Node y : adjx) {
                List<Node> _adjx = new ArrayList<>(adjacencies.get(x));
                _adjx.remove(y);
                List<Node> ppx = possibleParents(x, _adjx, this.knowledge);

                if (ppx.size() >= depth) {
                    ChoiceGenerator cg = new ChoiceGenerator(ppx.size(), depth);
                    int[] choice;

                    while ((choice = cg.next()) != null) {
                        List<Node> condSet = GraphUtils.asList(choice, ppx);
                        IndependenceResult result;

                        try {
                            result = test.checkIndependence(x, y, condSet);
                            this.pValueList.add(result.getPValue());
                        } catch (Exception e) {
                            result = new IndependenceResult(new IndependenceFact(x, y, condSet),
                                    false, Double.NaN);
                        }

                        boolean noEdgeRequired =
                                this.knowledge.noEdgeRequired(x.getName(), y.getName());

                        if (result.independent() && noEdgeRequired) {
                            adjacencies.get(x).remove(y);
                            adjacencies.get(y).remove(x);
                            getSepsets().set(x, y, condSet);

                            if (this.verbose) {
                                this.out.println(SearchLogUtils.independenceFact(x, y, condSet) + " p = " +
                                        this.nf.format(result.getPValue()));
                            }
                            continue EDGE;
                        }

                    }
                }
            }
        }

        return freeDegree(nodes, adjacencies) > depth;
    }

    private List<Node> possibleParents(Node x, List<Node> adjx,
                                       IKnowledge knowledge) {
        List<Node> possibleParents = new LinkedList<>();
        String _x = x.getName();

        for (Node z : adjx) {
            String _z = z.getName();

            if (possibleParentOf(_z, _x, knowledge)) {
                possibleParents.add(z);
            }
        }

        return possibleParents;
    }

    private boolean possibleParentOf(String z, String x, IKnowledge knowledge) {
        return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
    }

    public int getNumIndependenceTests() {
        return this.numIndependenceTests;
    }

    public int getNumDependenceJudgments() {
        return this.numDependenceJudgement;
    }

    public SepsetMap getSepsets() {
        return this.sepset;
    }

    public boolean isVerbose() {
        return this.verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    @Override
    public boolean isAggressivelyPreventCycles() {
        return false;
    }

    @Override
    public IndependenceTest getIndependenceTest() {
        return null;
    }

    @Override
    public Graph search(List<Node> nodes) {
        return null;
    }

    @Override
    public long getElapsedTime() {
        return 0;
    }

    @Override
    public List<Node> getNodes() {
        return null;
    }

    @Override
    public List<Triple> getAmbiguousTriples(Node node) {
        return null;
    }

    @Override
    public void setOut(PrintStream out) {
        this.out = out;
    }
}


