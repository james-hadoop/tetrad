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
package edu.cmu.tetrad.graph;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Edge.Property;
import edu.cmu.tetrad.graph.EdgeTypeProbability.EdgeType;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import nu.xom.*;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.RecursiveTask;
import java.util.regex.Matcher;

import static java.lang.Math.min;
import static java.util.Collections.shuffle;

/**
 * Basic graph utilities.
 *
 * @author Joseph Ramsey
 */
public final class GraphUtils {

    /**
     * Arranges the nodes in the graph in a circle.
     *
     * @param radius  The radius of the circle in pixels; a good default is 150.
     * @param centerx The x coordinate for the center of the layout.
     * @param centery The y coordinate for the center of the layout.
     * @param radius  The radius of the layout.
     */
    public static void circleLayout(Graph graph, int centerx, int centery, int radius) {
        if (graph == null) {
            return;
        }
        List<Node> nodes = graph.getNodes();
        Collections.sort(nodes);

        double rad = 6.28 / nodes.size();
        double phi = .75 * 6.28;    // start from 12 o'clock.

        for (Node node : nodes) {
            int centerX = centerx + (int) (radius * Math.cos(phi));
            int centerY = centery + (int) (radius * Math.sin(phi));

            node.setCenterX(centerX);
            node.setCenterY(centerY);

            phi += rad;
        }
    }

    public static void kamadaKawaiLayout(Graph graph, boolean randomlyInitialized, double naturalEdgeLength, double springConstant, double stopEnergy) {
        KamadaKawaiLayout layout = new KamadaKawaiLayout(graph);
        layout.setRandomlyInitialized(randomlyInitialized);
        layout.setNaturalEdgeLength(naturalEdgeLength);
        layout.setSpringConstant(springConstant);
        layout.setStopEnergy(stopEnergy);
        layout.doLayout();
    }

    public static void fruchtermanReingoldLayout(Graph graph) {
        FruchtermanReingoldLayout layout = new FruchtermanReingoldLayout(graph);
        layout.doLayout();
    }

    public static Graph randomDag(int numNodes, int numLatentConfounders, int maxNumEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        return GraphUtils.randomDag(nodes, numLatentConfounders, maxNumEdges, maxDegree, maxIndegree, maxOutdegree, connected);
    }

    public static Dag randomDag(List<Node> nodes, int numLatentConfounders, int maxNumEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {
        return new Dag(GraphUtils.randomGraph(nodes, numLatentConfounders, maxNumEdges, maxDegree, maxIndegree, maxOutdegree, connected));
    }

    public static Graph randomGraph(int numNodes, int numLatentConfounders, int numEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        return GraphUtils.randomGraph(nodes, numLatentConfounders, numEdges, maxDegree, maxIndegree, maxOutdegree, connected);
    }

    public static Graph randomGraph(List<Node> nodes, int numLatentConfounders, int maxNumEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {

        // It is still unclear whether we should use the random forward edges method or the
        // random uniform method to create random DAGs, hence this method.
        // jdramsey 12/8/2015
        return GraphUtils.randomGraphRandomForwardEdges(nodes, numLatentConfounders, maxNumEdges, maxDegree, maxIndegree, maxOutdegree, connected, true);
//        return randomGraphUniform(nodes, numLatentConfounders, maxNumEdges, maxDegree, maxIndegree, maxOutdegree, connected);
    }

    public static Graph randomGraphUniform(List<Node> nodes, int numLatentConfounders, int maxNumEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {
        int numNodes = nodes.size();

        if (numNodes <= 0) {
            throw new IllegalArgumentException("NumNodes most be > 0: " + numNodes);
        }

        if (maxNumEdges < 0 || maxNumEdges > numNodes * (numNodes - 1)) {
            throw new IllegalArgumentException("NumEdges must be " + "at least 0 and at most (#nodes)(#nodes - 1) / 2: " + maxNumEdges);
        }

        if (numLatentConfounders < 0 || numLatentConfounders > numNodes) {
            throw new IllegalArgumentException("Max # latent confounders must be " + "at least 0 and at most the number of nodes: " + numLatentConfounders);
        }

        for (Node node : nodes) {
            node.setNodeType(NodeType.MEASURED);
        }

        UniformGraphGenerator generator;

        if (connected) {
            generator = new UniformGraphGenerator(UniformGraphGenerator.CONNECTED_DAG);
        } else {
            generator = new UniformGraphGenerator(UniformGraphGenerator.ANY_DAG);
        }

        generator.setNumNodes(numNodes);
        generator.setMaxEdges(maxNumEdges);
        generator.setMaxDegree(maxDegree);
        generator.setMaxInDegree(maxIndegree);
        generator.setMaxOutDegree(maxOutdegree);
        generator.generate();
        Graph dag = generator.getDag(nodes);

        // Create a list of nodes. Add the nodes in the list to the
        // dag. Arrange the nodes in a circle.
        GraphUtils.fixLatents1(numLatentConfounders, dag);

        GraphUtils.circleLayout(dag, 200, 200, 150);

        return dag;
    }

    private static List<Node> getCommonCauses(Graph dag) {
        List<Node> commonCauses = new ArrayList<>();
        List<Node> nodes = dag.getNodes();

        for (Node node : nodes) {
            List<Node> children = dag.getChildren(node);

            if (children.size() >= 2) {
                commonCauses.add(node);
            }
        }

        return commonCauses;
    }

    public static Graph randomGraphRandomForwardEdges(int numNodes, int numLatentConfounders, int numEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {

        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        return GraphUtils.randomGraph(nodes, numLatentConfounders, numEdges, maxDegree, maxIndegree, maxOutdegree, connected);
    }

    public static Graph randomGraphRandomForwardEdges(List<Node> nodes, int numLatentConfounders, int numEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected) {
        return GraphUtils.randomGraphRandomForwardEdges(nodes, numLatentConfounders, numEdges, maxDegree, maxIndegree, maxOutdegree, connected, true);
    }

    public static Graph randomGraphRandomForwardEdges(List<Node> nodes, int numLatentConfounders, int numEdges, int maxDegree, int maxIndegree, int maxOutdegree, boolean connected, boolean layoutAsCircle) {
        if (nodes.size() <= 0) {
            throw new IllegalArgumentException("NumNodes most be > 0");
        }

        // Believe it or not ths is needed.
        long size = nodes.size();

        if (numEdges < 0 || numEdges > size * (size - 1)) {
            throw new IllegalArgumentException("numEdges must be " + "greater than 0 and <= (#nodes)(#nodes - 1) / 2: " + numEdges);
        }

        if (numLatentConfounders < 0 || numLatentConfounders > nodes.size()) {
            throw new IllegalArgumentException("MaxNumLatents must be " + "greater than 0 and less than the number of nodes: " + numLatentConfounders);
        }

        LinkedList<List<Integer>> allEdges = new LinkedList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                List<Integer> pair = new ArrayList<>(2);
                pair.add(i);
                pair.add(j);
                allEdges.add(pair);
            }
        }

        Graph dag;

        int numTriesForGraph = 0;

        do {
            dag = new EdgeListGraph(nodes);

            shuffle(allEdges);

            while (!allEdges.isEmpty() && dag.getNumEdges() < numEdges) {
                List<Integer> e = allEdges.removeFirst();

                Node n1 = nodes.get(e.get(0));
                Node n2 = nodes.get(e.get(1));

                if (dag.getIndegree(n2) >= maxIndegree) {
                    continue;
                }

                if (dag.getOutdegree(n1) >= maxOutdegree) {
                    continue;
                }

                if (dag.getIndegree(n1) + dag.getOutdegree(n1) >= maxDegree) {
                    continue;
                }

                if (dag.getIndegree(n2) + dag.getOutdegree(n2) >= maxDegree) {
                    continue;
                }

                dag.addDirectedEdge(n1, n2);
            }
        } while (++numTriesForGraph < 1000 && connected && (connectedComponents(dag).size() != 1));

        GraphUtils.fixLatents4(numLatentConfounders, dag);

        if (layoutAsCircle) {
            GraphUtils.circleLayout(dag, 200, 200, 150);
        }

        return dag;
    }

    private static boolean existsDisconnectedNode(Graph dag) {
        for (Node node : dag.getNodes()) {
            if (dag.getDegree(node) == 0) {
                return true;
            }
        }

        return false;
    }

    public static Graph scaleFreeGraph(int numNodes, int numLatentConfounders, double alpha, double beta, double delta_in, double delta_out) {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        return GraphUtils.scaleFreeGraph(nodes, numLatentConfounders, alpha, beta, delta_in, delta_out);
    }

    private static Graph scaleFreeGraph(List<Node> _nodes, int numLatentConfounders, double alpha, double beta, double delta_in, double delta_out) {

        if (alpha + beta >= 1) {
            throw new IllegalArgumentException("For the Bollobas et al. algorithm," + "\napha + beta + gamma = 1, so alpha + beta must be < 1.");
        }

//        System.out.println("# nodes = " + _nodes.size() + " latents = " + numLatentConfounders
//                + "  alpha = " + alpha + " beta = " + beta + " delta_in = " + delta_in + " delta_out = " + delta_out);
//        >>> print inspect.getsource(nx.scale_free_graph)
//        def scale_free_graph(n,
//                alpha=0.41,
//                beta=0.54,
//                gamma=0.05,
//                delta_in=0.2,
//                delta_out=0,
//                create_using=None,
//                seed=None):
//        """Return a scale free directed graph.
//
//        Parameters
//                ----------
//        n : integer
//        Number of nodes in graph
//        alpha : float
//        Probability for adding a new node connected to an existing node
//        chosen randomly according to the in-degree distribution.
//                beta : float
//        Probability for adding an edge between two existing nodes.
//                One existing node is chosen randomly according the in-degree
//        distribution and the other chosen randomly according to the out-degree
//        distribution.
//                gamma : float
//        Probability for adding a new node conecgted to an existing node
//        chosen randomly according to the out-degree distribution.
//                delta_in : float
//        Bias for choosing ndoes from in-degree distribution.
//                delta_out : float
//        Bias for choosing ndoes from out-degree distribution.
//                create_using : graph, optional (default MultiDiGraph)
//        Use this graph instance to start the process (default=3-cycle).
//        seed : integer, optional
//        Seed for random number generator
//
//                Examples
//        --------
//                >>> G=nx.scale_free_graph(100)
//
//        Notes
//                -----
//                The sum of alpha, beta, and gamma must be 1.
//
//        References
//                ----------
        //        .. [1] B. Bollob{\'a}s, C. Borgs, J. Chayes, and O. Riordan,
        //            Directed scale-free graphs,
        //                    Proceedings of the fourteenth annual ACM-SIAM symposium on
        //            Discrete algorithm, 132--139, 2003.
        //            """
//
//            def _choose_node(G,distribution,delta):
//            cumsum=0.0
//            # normalization
//            psum=float(sum(distribution.values()))+float(delta)*len(distribution)
//            r=random.random()
//            for i in range(0,len(distribution)):
//            cumsum+=(distribution[i]+delta)/psum
//            if r < cumsum:
//            break
//            return i
//
//            if create_using is None:
//            # start with 3-cycle
//            G = nx.MultiDiGraph()
//            G.add_edges_from([(0,1),(1,2),(2,0)])
//            else:
//            # keep existing graph structure?
//            G = create_using
//            if not (G.is_directed() and G.is_multigraph()):
//            raise nx.NetworkXError(\
//                    "MultiDiGraph required in create_using")
//
//            if alpha <= 0:
//            raise ValueError('alpha must be >= 0.')
//            if beta <= 0:
//            raise ValueError('beta must be >= 0.')
//            if gamma <= 0:
//            raise ValueError('beta must be >= 0.')
//
//            if alpha+beta+gamma !=1.0:
//            raise ValueError('alpha+beta+gamma must equal 1.')
//
//            G.name="directed_scale_free_graph(%s,alpha=%s,beta=%s,gamma=%s,delta_in=%s,delta_out=%s)"%(n,alpha,beta,gamma,delta_in,delta_out)
//
//            # seed random number generated (uses None as default)
//                random.seed(seed)
//
//                while len(G)<n:
//                r = random.random()
//                # random choice in alpha,beta,gamma ranges
//                if r<alpha:
//                # alpha
//                # add new node v
//                        v = len(G)
//                # choose w according to in-degree and delta_in
//                w = _choose_node(G, G.in_degree(),delta_in)
//                elif r < alpha+beta:
//                # beta
//                # choose v according to out-degree and delta_out
//                v = _choose_node(G, G.out_degree(),delta_out)
//                # choose w according to in-degree and delta_in
//                w = _choose_node(G, G.in_degree(),delta_in)
//                else:
//                # gamma
//                # choose v according to out-degree and delta_out
//                v = _choose_node(G, G.out_degree(),delta_out)
//                # add new node w
//                        w = len(G)
//                G.add_edge(v,w)
//
//                return G
        shuffle(_nodes);

        LinkedList<Node> nodes = new LinkedList<>();
        nodes.add(_nodes.get(0));

        Graph G = new EdgeListGraph(_nodes);

        if (alpha <= 0) {
            throw new IllegalArgumentException("alpha must be > 0.");
        }
        if (beta <= 0) {
            throw new IllegalArgumentException("beta must be > 0.");
        }

        double gamma = 1.0 - alpha - beta;

        if (gamma <= 0) {
            throw new IllegalArgumentException("alpha + beta must be < 1.");
        }

        if (delta_in <= 0) {
            throw new IllegalArgumentException("delta_in must be > 0.");
        }
        if (delta_out <= 0) {
            throw new IllegalArgumentException("delta_out must be > 0.");
        }

        Map<Node, Set<Node>> parents = new HashMap<>();
        Map<Node, Set<Node>> children = new HashMap<>();
        parents.put(_nodes.get(0), new HashSet<>());
        children.put(_nodes.get(0), new HashSet<>());

        while (nodes.size() < _nodes.size()) {
            double r = RandomUtil.getInstance().nextDouble();
            int v, w;

            if (r < alpha) {
                v = nodes.size();
                w = GraphUtils.chooseNode(GraphUtils.indegrees(nodes, parents), delta_in);
                Node m = _nodes.get(v);
                nodes.addFirst(m);
                parents.put(m, new HashSet<>());
                children.put(m, new HashSet<>());
                v = 0;
                w++;
            } else if (r < alpha + beta) {
                v = GraphUtils.chooseNode(GraphUtils.outdegrees(nodes, children), delta_out);
                w = GraphUtils.chooseNode(GraphUtils.indegrees(nodes, parents), delta_in);
                if (!(w > v)) {
                    continue;
                }
            } else {
                v = GraphUtils.chooseNode(GraphUtils.outdegrees(nodes, children), delta_out);
                w = nodes.size();
                Node m = _nodes.get(w);
                nodes.addLast(m);
                parents.put(m, new HashSet<>());
                children.put(m, new HashSet<>());
            }

            if (G.isAdjacentTo(nodes.get(v), nodes.get(w))) {
                continue;
            }

            G.addDirectedEdge(nodes.get(v), nodes.get(w));

            parents.get(nodes.get(w)).add(nodes.get(v));
            children.get(nodes.get(v)).add(nodes.get(w));
        }

        GraphUtils.fixLatents1(numLatentConfounders, G);

        GraphUtils.circleLayout(G, 200, 200, 150);

        return G;
    }

    private static int chooseNode(int[] distribution, double delta) {
        double cumsum = 0.0;
        double psum = GraphUtils.sum(distribution) + delta * distribution.length;
        double r = RandomUtil.getInstance().nextDouble();

        for (int i = 0; i < distribution.length; i++) {
            cumsum += (distribution[i] + delta) / psum;

            if (r < cumsum) {
                return i;
            }
        }

        throw new IllegalArgumentException("Didn't pick a node.");
    }

    private static int sum(int[] distribution) {
        int sum = 0;
        for (int w : distribution) {
            sum += w;
        }
        return sum;
    }

    private static int[] indegrees(List<Node> nodes, Map<Node, Set<Node>> parents) {
        int[] indegrees = new int[nodes.size()];

        for (int i = 0; i < nodes.size(); i++) {
            indegrees[i] = parents.get(nodes.get(i)).size();
        }

        return indegrees;
    }

    private static int[] outdegrees(List<Node> nodes, Map<Node, Set<Node>> children) {
        int[] outdegrees = new int[nodes.size()];

        for (int i = 0; i < nodes.size(); i++) {
            outdegrees[i] = children.get(nodes.get(i)).size();
        }

        return outdegrees;
    }

    public static void fixLatents1(int numLatentConfounders, Graph graph) {
        List<Node> commonCauses = GraphUtils.getCommonCauses(graph);
        int index = 0;

        while (index++ < numLatentConfounders) {
            if (commonCauses.size() == 0) {
                break;
            }
            int i = RandomUtil.getInstance().nextInt(commonCauses.size());
            Node node = commonCauses.get(i);
            node.setNodeType(NodeType.LATENT);
            commonCauses.remove(i);
        }
    }

    // JMO's method for fixing latents
    public static void fixLatents4(int numLatentConfounders, Graph graph) {
        if (numLatentConfounders == 0) {
            return;
        }

        List<Node> commonCausesAndEffects = GraphUtils.getCommonCausesAndEffects(graph);
        int index = 0;

        while (index++ < numLatentConfounders) {
            if (commonCausesAndEffects.size() == 0) {
                index--;
                break;
            }
            int i = RandomUtil.getInstance().nextInt(commonCausesAndEffects.size());
            Node node = commonCausesAndEffects.get(i);
            node.setNodeType(NodeType.LATENT);
            commonCausesAndEffects.remove(i);
        }

        List<Node> nodes = graph.getNodes();
        while (index++ < numLatentConfounders) {
            int r = RandomUtil.getInstance().nextInt(nodes.size());
            if (nodes.get(r).getNodeType() == NodeType.LATENT) {
                index--;
            } else {
                nodes.get(r).setNodeType(NodeType.LATENT);
            }
        }
    }

    //Helper method for fixLatents4
    //Common effects refers to common effects with at least one child
    private static List<Node> getCommonCausesAndEffects(Graph dag) {
        List<Node> commonCausesAndEffects = new ArrayList<>();
        List<Node> nodes = dag.getNodes();

        for (Node node : nodes) {
            List<Node> children = dag.getChildren(node);

            if (children.size() >= 2) {
                commonCausesAndEffects.add(node);
            } else {
                List<Node> parents = dag.getParents(node);
                if (parents.size() >= 2 && children.size() >= 1) {
                    commonCausesAndEffects.add(node);
                }
            }
        }

        return commonCausesAndEffects;
    }

    /**
     * Makes a cyclic graph by repeatedly adding cycles of length of 3, 4, or 5
     * to the graph, then finally adding two cycles.
     */
    public static Graph cyclicGraph2(int numNodes, int numEdges, int maxDegree) {

        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        Graph graph = new EdgeListGraph(nodes);

        LOOP:
        while (graph.getEdges().size() < numEdges /*&& ++count1 < 100*/) {
//            int cycleSize = RandomUtil.getInstance().nextInt(2) + 4;
            int cycleSize = RandomUtil.getInstance().nextInt(3) + 3;

            // Pick that many nodes randomly
            List<Node> cycleNodes = new ArrayList<>();
            int count2 = -1;

            for (int i = 0; i < cycleSize; i++) {
                Node node = nodes.get(RandomUtil.getInstance().nextInt(nodes.size()));

                if (cycleNodes.contains(node)) {
                    i--;
                    ++count2;
                    if (count2 < 10) {
                        continue;
                    }
                }

                cycleNodes.add(node);
            }

            for (Node cycleNode : cycleNodes) {
                if (graph.getDegree(cycleNode) >= maxDegree) {
                    continue LOOP;
                }
            }

            Edge edge;

            // Make sure you won't created any two cycles (this will be done later, explicitly)
            for (int i = 0; i < cycleNodes.size() - 1; i++) {
                edge = Edges.directedEdge(cycleNodes.get(i + 1), cycleNodes.get(i));

                if (graph.containsEdge(edge)) {
                    continue LOOP;
                }
            }

            edge = Edges.directedEdge(cycleNodes.get(0), cycleNodes.get(cycleNodes.size() - 1));

            if (graph.containsEdge(edge)) {
                continue;
            }

            for (int i = 0; i < cycleNodes.size() - 1; i++) {
                edge = Edges.directedEdge(cycleNodes.get(i), cycleNodes.get(i + 1));

                if (!graph.containsEdge(edge)) {
                    graph.addEdge(edge);

                    if (graph.getNumEdges() == numEdges) {
                        break LOOP;
                    }
                }
            }

            edge = Edges.directedEdge(cycleNodes.get(cycleNodes.size() - 1), cycleNodes.get(0));

            if (!graph.containsEdge(edge)) {
                graph.addEdge(edge);

                if (graph.getNumEdges() == numEdges) {
                    break;
                }
            }
        }

        GraphUtils.circleLayout(graph, 200, 200, 150);

        return graph;
    }

    /**
     * Makes a cyclic graph by repeatedly adding cycles of length of 3, 4, or 5
     * to the graph, then finally adding two cycles.
     */
    public static Graph cyclicGraph3(int numNodes, int numEdges, int maxDegree, double probCycle, double probTwoCycle) {

        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + (i + 1)));
        }

        Graph graph = new EdgeListGraph(nodes);

        LOOP:
        while (graph.getEdges().size() < numEdges /*&& ++count1 < 100*/) {
//            int cycleSize = RandomUtil.getInstance().nextInt(2) + 4;
            int cycleSize = RandomUtil.getInstance().nextInt(3) + 3;

            // Pick that many nodes randomly
            List<Node> cycleNodes = new ArrayList<>();
            int count2 = -1;

            for (int i = 0; i < cycleSize; i++) {
                Node node = nodes.get(RandomUtil.getInstance().nextInt(nodes.size()));

                if (cycleNodes.contains(node)) {
                    i--;
                    ++count2;
                    if (count2 < 10) {
                        continue;
                    }
                }

                cycleNodes.add(node);
            }

            for (Node cycleNode : cycleNodes) {
                if (graph.getDegree(cycleNode) >= maxDegree) {
                    continue LOOP;
                }
            }

            Edge edge;

            // Make sure you won't created any two cycles (this will be done later, explicitly)
            for (int i = 0; i < cycleNodes.size() - 1; i++) {
                edge = Edges.directedEdge(cycleNodes.get(i + 1), cycleNodes.get(i));

                if (graph.containsEdge(edge)) {
                    continue LOOP;
                }
            }

            edge = Edges.directedEdge(cycleNodes.get(0), cycleNodes.get(cycleNodes.size() - 1));

            if (graph.containsEdge(edge)) {
                continue;
            }

            for (int i = 0; i < cycleNodes.size() - 1; i++) {
                edge = Edges.directedEdge(cycleNodes.get(i), cycleNodes.get(i + 1));

                if (!graph.containsEdge(edge)) {
                    graph.addEdge(edge);

                    if (graph.getNumEdges() == numEdges) {
                        break LOOP;
                    }
                }
            }

            if (RandomUtil.getInstance().nextDouble() < probCycle) {
                edge = Edges.directedEdge(cycleNodes.get(cycleNodes.size() - 1), cycleNodes.get(0));
            } else {
                edge = Edges.directedEdge(cycleNodes.get(0), cycleNodes.get(cycleNodes.size() - 1));
            }

            if (!graph.containsEdge(edge)) {
                graph.addEdge(edge);

                if (graph.getNumEdges() == numEdges) {
                    break;
                }
            }
        }

        Set<Edge> edges = graph.getEdges();

        for (Edge edge : edges) {
            Node a = edge.getNode1();
            Node b = edge.getNode2();
            if (RandomUtil.getInstance().nextDouble() < probTwoCycle) {
                graph.removeEdges(a, b);
                graph.addEdge(Edges.directedEdge(a, b));
                graph.addEdge(Edges.directedEdge(b, a));
            }
        }

        GraphUtils.circleLayout(graph, 200, 200, 150);

        return graph;
    }

    public static void addTwoCycles(Graph graph, int numTwoCycles) {
        List<Edge> edges = new ArrayList<>(graph.getEdges());
        shuffle(edges);

        for (int i = 0; i < min(numTwoCycles, edges.size()); i++) {
            Edge edge = edges.get(i);
            Edge reversed = Edges.directedEdge(edge.getNode2(), edge.getNode1());

            if (graph.containsEdge(reversed)) {
                i--;
                continue;
            }

            graph.addEdge(reversed);
        }
    }

    /**
     * Arranges the nodes in the result graph according to their positions in
     * the source graph.
     *
     * @return true if all of the nodes were arranged, false if not.
     */
    public static boolean arrangeBySourceGraph(Graph resultGraph, Graph sourceGraph) {
        if (resultGraph == null) {
            throw new IllegalArgumentException("Graph must not be null.");
        }

        if (sourceGraph == null) {
            GraphUtils.circleLayout(resultGraph, 200, 200, 150);
            return true;
        }

        boolean arrangedAll = true;

        // There is a source graph. Position the nodes in the
        // result graph correspondingly.
        for (Node o : resultGraph.getNodes()) {
            String name = o.getName();
            Node sourceNode = sourceGraph.getNode(name);

            if (sourceNode == null) {
                arrangedAll = false;
                continue;
            }

            o.setCenterX(sourceNode.getCenterX());
            o.setCenterY(sourceNode.getCenterY());
        }

        return arrangedAll;
    }

    public static void arrangeByLayout(Graph graph, HashMap<String, PointXy> layout) {
        for (Node node : graph.getNodes()) {
            PointXy point = layout.get(node.getName());
            node.setCenter(point.getX(), point.getY());
        }
    }

    /**
     * @return the node associated with a given error node. This should be the
     * only child of the error node, E --&gt; N.
     */
    public static Node getAssociatedNode(Node errorNode, Graph graph) {
        if (errorNode.getNodeType() != NodeType.ERROR) {
            throw new IllegalArgumentException("Can only get an associated node " + "for an error node: " + errorNode);
        }

        List<Node> children = graph.getChildren(errorNode);

        if (children.size() != 1) {
            System.out.println("children of " + errorNode + " = " + children);
            System.out.println(graph);

            throw new IllegalArgumentException("An error node should have only " + "one child, which is its associated node: " + errorNode);
        }

        return children.get(0);
    }

    /**
     * @return true if <code>set</code> is a clique in <code>graph</code>.
     * R. Silva, June 2004
     */
    public static boolean isClique(Collection<Node> set, Graph graph) {
        List<Node> setv = new LinkedList<>(set);
        for (int i = 0; i < setv.size() - 1; i++) {
            for (int j = i + 1; j < setv.size(); j++) {
                if (!graph.isAdjacentTo(setv.get(i), setv.get(j))) {
                    return false;
                }
            }
        }
        return true;
    }

    public static Graph loadRSpecial(File file) {
        DataSet eg = null;

        try {
            ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(file.toPath(), Delimiter.COMMA);
            reader.setHasHeader(false);
            Data data = reader.readInData();
            eg = (DataSet) DataConvertUtils.toDataModel(data);
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }

        if (eg == null) throw new NullPointerException();

        List<Node> vars = eg.getVariables();

        Graph graph = new EdgeListGraph(vars);

        for (int i = 0; i < vars.size(); i++) {
            for (int j = 0; j < vars.size(); j++) {
                if (i == j) continue;
                if (eg.getDouble(i, j) == 1 && eg.getDouble(j, i) == 1) {
                    if (!graph.isAdjacentTo(vars.get(i), vars.get(j))) {
                        graph.addUndirectedEdge(vars.get(i), vars.get(j));
                    }
                } else if (eg.getDouble(i, j) == 1 && eg.getDouble(j, i) == 0) {
                    graph.addDirectedEdge(vars.get(i), vars.get(j));
                }
            }
        }

        return graph;
    }

    /**
     * Calculates the Markov blanket of a target in a DAG. This includes the
     * target, the parents of the target, the children of the target, the
     * parents of the children of the target, edges from parents to target,
     * target to children, parents of children to children, and parent to
     * parents of children. (Edges among children are implied by the inclusion
     * of edges from parents of children to children.) Edges among parents and
     * among parents of children not explicitly included above are not included.
     * (Joseph Ramsey 8/6/04)
     *
     * @param target a node in the given DAG.
     * @param dag    the DAG with respect to which a Markov blanket DAG is to to be
     *               calculated. All of the nodes and edges of the Markov Blanket DAG are in
     *               this DAG.
     */
    public static Dag markovBlanketDag(Node target, Graph dag) {
        if (dag.getNode(target.getName()) == null) {
            throw new NullPointerException("Target node not in graph: " + target);
        }

        Graph blanket = new EdgeListGraph();
        blanket.addNode(target);

        // Add parents of target.
        List<Node> parents = dag.getParents(target);
        for (Node parent1 : parents) {
            blanket.addNode(parent1);

            blanket.addDirectedEdge(parent1, target);
        }

        // Add children of target and parents of children of target.
        List<Node> children = dag.getChildren(target);
        List<Node> parentsOfChildren = new LinkedList<>();
        for (Node child : children) {
            if (!blanket.containsNode(child)) {
                blanket.addNode(child);
            }

            blanket.addDirectedEdge(target, child);

            List<Node> parentsOfChild = dag.getParents(child);
            parentsOfChild.remove(target);
            for (Node aParentsOfChild : parentsOfChild) {
                if (!parentsOfChildren.contains(aParentsOfChild)) {
                    parentsOfChildren.add(aParentsOfChild);
                }

                if (!blanket.containsNode(aParentsOfChild)) {
                    blanket.addNode(aParentsOfChild);
                }

                blanket.addDirectedEdge(aParentsOfChild, child);
            }
        }

        // Add in edges connecting parents and parents of children.
        parentsOfChildren.removeAll(parents);

        for (Node parent2 : parents) {
            for (Node aParentsOfChildren : parentsOfChildren) {
                Edge edge1 = dag.getEdge(parent2, aParentsOfChildren);
                Edge edge2 = blanket.getEdge(parent2, aParentsOfChildren);

                if (edge1 != null && edge2 == null) {
                    Edge newEdge = new Edge(parent2, aParentsOfChildren, edge1.getProximalEndpoint(parent2), edge1.getProximalEndpoint(aParentsOfChildren));

                    blanket.addEdge(newEdge);
                }
            }
        }

        // Add in edges connecting children and parents of children.
        for (Node aChildren1 : children) {

            for (Node aParentsOfChildren : parentsOfChildren) {
                Edge edge1 = dag.getEdge(aChildren1, aParentsOfChildren);
                Edge edge2 = blanket.getEdge(aChildren1, aParentsOfChildren);

                if (edge1 != null && edge2 == null) {
                    Edge newEdge = new Edge(aChildren1, aParentsOfChildren, edge1.getProximalEndpoint(aChildren1), edge1.getProximalEndpoint(aParentsOfChildren));

                    blanket.addEdge(newEdge);
                }
            }
        }

        return new Dag(blanket);
    }

    /**
     * @return the connected components of the given graph, as a list of lists
     * of nodes.
     */
    public static List<List<Node>> connectedComponents(Graph graph) {
        graph = new EdgeListGraph(graph);
        List<List<Node>> components = new LinkedList<>();
        LinkedList<Node> unsortedNodes = new LinkedList<>(graph.getNodes());

        while (!unsortedNodes.isEmpty()) {
            Node seed = unsortedNodes.removeFirst();
            Set<Node> component = new ConcurrentSkipListSet<>();
            GraphUtils.collectComponentVisit(seed, component, graph, unsortedNodes);
            components.add(new ArrayList<>(component));
        }

        return components;
    }

    /**
     * Assumes node should be in component.
     */
    private static void collectComponentVisit(Node node, Set<Node> component, Graph graph, List<Node> unsortedNodes) {
        if (TaskManager.getInstance().isCanceled()) {
            return;
        }

        component.add(node);
        unsortedNodes.remove(node);
        List<Node> adj = graph.getAdjacentNodes(node);

        for (Node anAdj : adj) {
            if (!component.contains(anAdj)) {
                GraphUtils.collectComponentVisit(anAdj, component, graph, unsortedNodes);
            }
        }
    }

    /**
     * @param graph The graph in which a directed cycle is sought.
     * @return the first directed cycle encountered in <code>graph</code>.
     */
    public static List<Node> directedCycle(Graph graph) {
        for (Node node : graph.getNodes()) {
            List<Node> path = GraphUtils.directedPathFromTo(graph, node, node);

            if (path != null) {
                return path;
            }
        }

        return null;
    }

    /**
     * @param graph The graph in which a directed path is sought.
     * @param node1 The 'from' node.
     * @param node2 The 'to'node.
     * @return A path from <code>node1</code> to <code>node2</code>, or null if
     * there is no path.
     */
    private static List<Node> directedPathFromTo(Graph graph, Node node1, Node node2) {
        return GraphUtils.directedPathVisit(graph, node1, node2, new LinkedList<>());
    }

    /**
     * @return the path of the first directed path found from node1 to node2, if
     * any.
     */
    private static List<Node> directedPathVisit(Graph graph, Node node1, Node node2, LinkedList<Node> path) {
        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverseDirected(node1, edge);

            if (child == null) {
                continue;
            }

            if (child == node2) {
                return path;
            }

            if (path.contains(child)) {
                continue;
            }

            if (GraphUtils.directedPathVisit(graph, child, node2, path) != null) {
                return path;
            }
        }

        path.removeLast();
        return null;
    }

    //all adjancencies are directed <=> there is no uncertainty about who the parents of 'node' are.
    public static boolean allAdjacenciesAreDirected(Node node, Graph graph) {
        List<Edge> nodeEdges = graph.getEdges(node);
        for (Edge edge : nodeEdges) {
            if (!edge.isDirected()) {
                return false;
            }
        }
        return true;
    }

    public static Graph removeBidirectedOrientations(Graph estCpdag) {
        estCpdag = new EdgeListGraph(estCpdag);

        // Make bidirected edges undirected.
        for (Edge edge : estCpdag.getEdges()) {
            if (Edges.isBidirectedEdge(edge)) {
                estCpdag.removeEdge(edge);
                estCpdag.addUndirectedEdge(edge.getNode1(), edge.getNode2());
            }
        }

        return estCpdag;
    }

    public static Graph removeBidirectedEdges(Graph estCpdag) {
        estCpdag = new EdgeListGraph(estCpdag);

        // Remove bidirected edges altogether.
        for (Edge edge : new ArrayList<>(estCpdag.getEdges())) {
            if (Edges.isBidirectedEdge(edge)) {
                estCpdag.removeEdge(edge);
            }
        }

        return estCpdag;
    }

    public static Graph undirectedGraph(Graph graph) {
        Graph graph2 = new EdgeListGraph(graph.getNodes());

        for (Edge edge : graph.getEdges()) {
            if (!graph2.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                graph2.addUndirectedEdge(edge.getNode1(), edge.getNode2());
            }
        }

        return graph2;
    }

    public static Graph completeGraph(Graph graph) {
        Graph graph2 = new EdgeListGraph(graph.getNodes());

        graph2.removeEdges(new ArrayList<>(graph2.getEdges()));

        List<Node> nodes = graph2.getNodes();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                Node node1 = nodes.get(i);
                Node node2 = nodes.get(j);
                graph2.addUndirectedEdge(node1, node2);
            }
        }

        return graph2;
    }

    public static List<List<Node>> directedPathsFromTo(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.directedPathsFromToVisit(graph, node1, node2, new LinkedList<>(), paths, maxLength);
        return paths;
    }

    private static void directedPathsFromToVisit(Graph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        if (maxLength != -1 && path.size() > maxLength - 2) {
            return;
        }

        int witnessed = 0;

        for (Node node : path) {
            if (node == node1) {
                witnessed++;
            }
        }

        if (witnessed > 1) {
            return;
        }

        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverseDirected(node1, edge);

            if (child == null) {
                continue;
            }

            if (path.contains(child)) {
                continue;
            }

            if (child == node2) {
                LinkedList<Node> _path = new LinkedList<>(path);
                _path.add(child);
                paths.add(_path);
                continue;
            }

            GraphUtils.directedPathsFromToVisit(graph, child, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    public static List<List<Node>> semidirectedPathsFromTo(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.semidirectedPathsFromToVisit(graph, node1, node2, new LinkedList<>(), paths, maxLength);
        return paths;
    }

    private static void semidirectedPathsFromToVisit(Graph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        if (maxLength != -1 && path.size() > maxLength - 2) {
            return;
        }

        int witnessed = 0;

        for (Node node : path) {
            if (node == node1) {
                witnessed++;
            }
        }

        if (witnessed > 1) {
            return;
        }

        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverseSemiDirected(node1, edge);

            if (child == null) {
                continue;
            }

            if (child == node2) {
                LinkedList<Node> _path = new LinkedList<>(path);
                _path.add(child);
                paths.add(_path);
                continue;
            }

            if (path.contains(child)) {
                continue;
            }

            GraphUtils.semidirectedPathsFromToVisit(graph, child, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    public static List<List<Node>> allPathsFromTo(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.allPathsFromToVisit(graph, node1, node2, new LinkedList<>(), paths, maxLength);
        return paths;
    }

    private static void allPathsFromToVisit(Graph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        path.addLast(node1);

        if (path.size() > (maxLength == -1 ? 1000 : maxLength)) {
            return;
        }

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverse(node1, edge);

            if (child == null) {
                continue;
            }

            if (path.contains(child)) {
                continue;
            }

            if (child == node2) {
                List<Node> _path = new ArrayList<>(path);
                _path.add(child);
                paths.add(_path);
                continue;
            }

            GraphUtils.allPathsFromToVisit(graph, child, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    public static List<List<Node>> allDirectedPathsFromTo(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.allDirectedPathsFromToVisit(graph, node1, node2, new LinkedList<>(), paths, maxLength);
        return paths;
    }

    private static void allDirectedPathsFromToVisit(Graph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        path.addLast(node1);

        if (path.size() > (maxLength == -1 ? 1000 : maxLength)) {
            return;
        }

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverseDirected(node1, edge);

            if (child == null) {
                continue;
            }

            if (path.contains(child)) {
                continue;
            }

            if (child == node2) {
                List<Node> _path = new ArrayList<>(path);
                _path.add(child);
                paths.add(_path);
                continue;
            }

            GraphUtils.allDirectedPathsFromToVisit(graph, child, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    public static List<List<Node>> treks(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.treks(graph, node1, node2, new LinkedList<>(), paths, maxLength);
        return paths;
    }

    private static void treks(Graph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        if (path.size() > (maxLength == -1 ? 1000 : maxLength - 2)) {
            return;
        }

        if (path.contains(node1)) {
            return;
        }

        if (node1 == node2) {
            return;
        }

        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node next = Edges.traverse(node1, edge);

            // Must be a directed edge.
            if (!edge.isDirected()) {
                continue;
            }

            // Can't have any colliders on the path.
            if (path.size() > 1) {
                Node node0 = path.get(path.size() - 2);

                if (next == node0) {
                    continue;
                }

                if (graph.isDefCollider(node0, node1, next)) {
                    continue;
                }
            }

            // Found a path.
            if (next == node2) {
                LinkedList<Node> _path = new LinkedList<>(path);
                _path.add(next);
                paths.add(_path);
                continue;
            }

            // Nodes may only appear on the path once.
            if (path.contains(next)) {
                continue;
            }

            GraphUtils.treks(graph, next, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    public static List<List<Node>> treksIncludingBidirected(SemGraph graph, Node node1, Node node2) {
        List<List<Node>> paths = new LinkedList<>();
        GraphUtils.treksIncludingBidirected(graph, node1, node2, new LinkedList<>(), paths);
        return paths;
    }

    private static void treksIncludingBidirected(SemGraph graph, Node node1, Node node2, LinkedList<Node> path, List<List<Node>> paths) {
        if (!graph.isShowErrorTerms()) {
            throw new IllegalArgumentException("The SEM Graph must be showing its error terms; this method " + "doesn't traverse two edges between the same nodes well.");
        }

        if (path.contains(node1)) {
            return;
        }

        if (node1 == node2) {
            return;
        }

        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node next = Edges.traverse(node1, edge);

            // Must be a directed edge or a bidirected edge.
            if (!(edge.isDirected() || Edges.isBidirectedEdge(edge))) {
                continue;
            }

            // Can't have any colliders on the path.
            if (path.size() > 1) {
                Node node0 = path.get(path.size() - 2);

                if (next == node0) {
                    continue;
                }

                if (graph.isDefCollider(node0, node1, next)) {
                    continue;
                }
            }

            // Found a path.
            if (next == node2) {
                LinkedList<Node> _path = new LinkedList<>(path);
                _path.add(next);
                paths.add(_path);
                continue;
            }

            // Nodes may only appear on the path once.
            if (path.contains(next)) {
                continue;
            }

            GraphUtils.treksIncludingBidirected(graph, next, node2, path, paths);
        }

        path.removeLast();
    }

    /**
     * @return the edges that are in <code>graph1</code> but not in
     * <code>graph2</code>, as a list of undirected edges..
     */
    public static List<Edge> adjacenciesComplement(Graph graph1, Graph graph2) {
        List<Edge> edges = new ArrayList<>();

        for (Edge edge1 : graph1.getEdges()) {
            String name1 = edge1.getNode1().getName();
            String name2 = edge1.getNode2().getName();

            Node node21 = graph2.getNode(name1);
            Node node22 = graph2.getNode(name2);

            if (node21 == null || node22 == null || !graph2.isAdjacentTo(node21, node22)) {
                edges.add(Edges.nondirectedEdge(edge1.getNode1(), edge1.getNode2()));
            }
        }

        return edges;
    }

    /**
     * @return a new graph in which the bidirectred edges of the given graph
     * have been changed to undirected edges.
     */
    public static Graph bidirectedToUndirected(Graph graph) {
        Graph newGraph = new EdgeListGraph(graph);

        for (Edge edge : newGraph.getEdges()) {
            if (Edges.isBidirectedEdge(edge)) {
                newGraph.removeEdge(edge);
                newGraph.addUndirectedEdge(edge.getNode1(), edge.getNode2());
            }
        }

        return newGraph;
    }

    /**
     * @return a new graph in which the undirectred edges of the given graph
     * have been changed to bidirected edges.
     */
    public static Graph undirectedToBidirected(Graph graph) {
        Graph newGraph = new EdgeListGraph(graph);

        for (Edge edge : newGraph.getEdges()) {
            if (Edges.isUndirectedEdge(edge)) {
                newGraph.removeEdge(edge);
                newGraph.addBidirectedEdge(edge.getNode1(), edge.getNode2());
            }
        }

        return newGraph;
    }

    public static String pathString(List<Node> path, Graph graph) {
        return GraphUtils.pathString(graph, path, new LinkedList<>());
    }

    public static String pathString(Graph graph, Node... x) {
        List<Node> path = new ArrayList<>();
        Collections.addAll(path, x);
        return GraphUtils.pathString(graph, path, new LinkedList<>());
    }

    private static String pathString(Graph graph, List<Node> path, List<Node> conditioningVars) {
        StringBuilder buf = new StringBuilder();

        if (path.size() < 2) {
            return "NO PATH";
        }

        buf.append(path.get(0).toString());

        if (conditioningVars.contains(path.get(0))) {
            buf.append("(C)");
        }

        for (int m = 1; m < path.size(); m++) {
            Node n0 = path.get(m - 1);
            Node n1 = path.get(m);

            Edge edge = graph.getEdge(n0, n1);


            if (edge == null) {
                buf.append("(-)");
            } else {
                Endpoint endpoint0 = edge.getProximalEndpoint(n0);
                Endpoint endpoint1 = edge.getProximalEndpoint(n1);

                if (endpoint0 == Endpoint.ARROW) {
                    buf.append("<");
                } else if (endpoint0 == Endpoint.TAIL) {
                    buf.append("-");
                } else if (endpoint0 == Endpoint.CIRCLE) {
                    buf.append("o");
                }

                buf.append("-");

                if (endpoint1 == Endpoint.ARROW) {
                    buf.append(">");
                } else if (endpoint1 == Endpoint.TAIL) {
                    buf.append("-");
                } else if (endpoint1 == Endpoint.CIRCLE) {
                    buf.append("o");
                }
            }

            buf.append(n1.toString());

            if (conditioningVars.contains(n1)) {
                buf.append("(C)");
            }
        }
        return buf.toString();
    }

    /**
     * Converts the given graph, <code>originalGraph</code>, to use the new
     * variables (with the same names as the old).
     *
     * @param originalGraph The graph to be converted.
     * @param newVariables  The new variables to use, with the same names as the
     *                      old ones.
     * @return A new, converted, graph.
     */
    public static Graph replaceNodes(Graph originalGraph, List<Node> newVariables) {
        Map<String, Node> newNodes = new HashMap<>();

        for (Node node : newVariables) {
            newNodes.put(node.getName(), node);
        }

        Graph convertedGraph = new EdgeListGraph(newVariables);

        for (Edge edge : originalGraph.getEdges()) {
            Node node1 = newNodes.get(edge.getNode1().getName());
            Node node2 = newNodes.get(edge.getNode2().getName());

            if (node1 == null) {
                node1 = edge.getNode1();

            }
            if (!convertedGraph.containsNode(node1)) {
                convertedGraph.addNode(node1);
            }

            if (node2 == null) {
                node2 = edge.getNode2();
            }

            if (!convertedGraph.containsNode(node2)) {
                convertedGraph.addNode(node2);
            }

            if (!convertedGraph.containsNode(node1)) {
                convertedGraph.addNode(node1);
            }

            if (!convertedGraph.containsNode(node2)) {
                convertedGraph.addNode(node2);
            }

            Endpoint endpoint1 = edge.getEndpoint1();
            Endpoint endpoint2 = edge.getEndpoint2();
            Edge newEdge = new Edge(node1, node2, endpoint1, endpoint2);
            convertedGraph.addEdge(newEdge);
        }

        for (Triple triple : originalGraph.getUnderLines()) {
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            convertedGraph.addUnderlineTriple(convertedGraph.getNode(triple.getX().getName()), convertedGraph.getNode(triple.getY().getName()), convertedGraph.getNode(triple.getZ().getName()));
        }

        for (Triple triple : originalGraph.getDottedUnderlines()) {
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            convertedGraph.addDottedUnderlineTriple(convertedGraph.getNode(triple.getX().getName()), convertedGraph.getNode(triple.getY().getName()), convertedGraph.getNode(triple.getZ().getName()));
        }

        for (Triple triple : originalGraph.getAmbiguousTriples()) {
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            convertedGraph.addAmbiguousTriple(convertedGraph.getNode(triple.getX().getName()), convertedGraph.getNode(triple.getY().getName()), convertedGraph.getNode(triple.getZ().getName()));
        }

        return convertedGraph;
    }

    public static Graph restrictToMeasured(Graph graph) {
        graph = new EdgeListGraph(graph);

        for (Node node : graph.getNodes()) {
            if (node.getNodeType() == NodeType.LATENT) {
                graph.removeNode(node);
            }
        }

        return graph;
    }

    /**
     * Converts the given list of nodes, <code>originalNodes</code>, to use the
     * new variables (with the same names as the old).
     *
     * @param originalNodes The list of nodes to be converted.
     * @param newNodes      A list of new nodes, containing as a subset nodes with
     *                      the same names as those in <code>originalNodes</code>. the old ones.
     * @return The converted list of nodes.
     */
    public static List<Node> replaceNodes(List<Node> originalNodes, List<Node> newNodes) {
        List<Node> convertedNodes = new LinkedList<>();

        for (Node node : originalNodes) {
            if (node == null) {
                throw new NullPointerException("Null node among original nodes.");
            }

            for (Node _node : newNodes) {
                if (_node == null) {
                    throw new NullPointerException("Null node among new nodes.");
                }

                if (node.getName().equals(_node.getName())) {
                    convertedNodes.add(_node);
                    break;
                }
            }
        }

        return convertedNodes;
    }

    /**
     * Counts the adjacencies that are in graph1 but not in graph2.
     *
     * @throws IllegalArgumentException if graph1 and graph2 are not namewise
     *                                  isomorphic.
     */
    public static int countAdjErrors(Graph graph1, Graph graph2) {
        if (graph1 == null) {
            throw new NullPointerException("The reference graph is missing.");
        }

        if (graph2 == null) {
            throw new NullPointerException("The target graph is missing.");
        }

        graph2 = GraphUtils.replaceNodes(graph2, graph1.getNodes());

        graph1 = GraphUtils.undirectedGraph(graph1);
        graph2 = GraphUtils.undirectedGraph(graph2);

        int count = 0;

        Set<Edge> edges1 = graph1.getEdges();

        for (Edge edge : edges1) {
            if (!graph2.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                ++count;
            }
        }

        return count;
    }

    /**
     * Counts the arrowpoints that are in graph1 but not in graph2.
     */
    public static int countArrowptErrors(Graph graph1, Graph graph2) {
        if (graph1 == null) {
            throw new NullPointerException("The reference graph is missing.");
        }

        if (graph2 == null) {
            throw new NullPointerException("The target graph is missing.");
        }

        graph2 = GraphUtils.replaceNodes(graph2, graph1.getNodes());

        assert graph2 != null;

        int count = 0;

        for (Edge edge1 : graph1.getEdges()) {
            Node node1 = edge1.getNode1();
            Node node2 = edge1.getNode2();

            Edge edge2 = graph2.getEdge(node1, node2);

            if (edge1.getEndpoint1() == Endpoint.ARROW) {
                if (edge2 == null) {
                    count++;
                } else if (edge2.getProximalEndpoint(edge1.getNode1()) != Endpoint.ARROW) {
                    count++;
                }
            }

            if (edge1.getEndpoint2() == Endpoint.ARROW) {
                if (edge2 == null) {
                    count++;
                } else if (edge2.getProximalEndpoint(edge1.getNode2()) != Endpoint.ARROW) {
                    count++;
                }
            }
        }

        for (Edge edge1 : graph2.getEdges()) {
            Node node1 = edge1.getNode1();
            Node node2 = edge1.getNode2();

            Edge edge2 = graph1.getEdge(node1, node2);

            if (edge1.getEndpoint1() == Endpoint.ARROW) {
                if (edge2 == null) {
                    count++;
                } else if (edge2.getProximalEndpoint(edge1.getNode1()) != Endpoint.ARROW) {
                    count++;
                }
            }

            if (edge1.getEndpoint2() == Endpoint.ARROW) {
                if (edge2 == null) {
                    count++;
                } else if (edge2.getProximalEndpoint(edge1.getNode2()) != Endpoint.ARROW) {
                    count++;
                }
            }
        }

        return count;
    }

    public static int getNumCorrectArrowpts(Graph correct, Graph estimated) {
        correct = GraphUtils.replaceNodes(correct, estimated.getNodes());

        Set<Edge> edges = estimated.getEdges();
        int numCorrect = 0;

        for (Edge estEdge : edges) {
            assert correct != null;
            Edge correctEdge = correct.getEdge(estEdge.getNode1(), estEdge.getNode2());
            if (correctEdge == null) {
                continue;
            }

            if (estEdge.getProximalEndpoint(estEdge.getNode1()) == Endpoint.ARROW && correctEdge.getProximalEndpoint(estEdge.getNode1()) == Endpoint.ARROW) {
                numCorrect++;
            }

            if (estEdge.getProximalEndpoint(estEdge.getNode2()) == Endpoint.ARROW && correctEdge.getProximalEndpoint(estEdge.getNode2()) == Endpoint.ARROW) {
                numCorrect++;
            }
        }

        return numCorrect;
    }

    /**
     * Converts the given list of nodes, <code>originalNodes</code>, to use the
     * replacement nodes for them by the same name in the given
     * <code>graph</code>.
     *
     * @param originalNodes The list of nodes to be converted.
     * @param graph         A graph to be used as a source of new nodes.
     * @return A new, converted, graph.
     */
    public static List<Node> replaceNodes(List<Node> originalNodes, Graph graph) {
        List<Node> convertedNodes = new LinkedList<>();

        for (Node node : originalNodes) {
            convertedNodes.add(graph.getNode(node.getName()));
        }

        return convertedNodes;
    }

    /**
     * @return an empty graph with the given number of nodes.
     */
    public static Graph emptyGraph(int numNodes) {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < numNodes; i++) {
            nodes.add(new GraphNode("X" + i));
        }

        return new EdgeListGraph(nodes);
    }

    /**
     * Converts a graph to a Graphviz .dot file
     */
    public static String graphToDot(Graph graph) {
        StringBuilder builder = new StringBuilder();
        builder.append("digraph g {\n");
        for (Edge edge : graph.getEdges()) {
            String n1 = edge.getNode1().getName();
            String n2 = edge.getNode2().getName();

            Endpoint end1 = edge.getEndpoint1();
            Endpoint end2 = edge.getEndpoint2();

            if (n1.compareTo(n2) > 0) {
                String temp = n1;
                n1 = n2;
                n2 = temp;

                Endpoint tmp = end1;
                end1 = end2;
                end2 = tmp;
            }
            builder.append(" \"").append(n1).append("\" -> \"").append(n2).append("\" [");

            if (end1 != Endpoint.TAIL) {
                builder.append("dir=both, ");
            }

            builder.append("arrowtail=");
            if (end1 == Endpoint.ARROW) {
                builder.append("normal");
            } else if (end1 == Endpoint.TAIL) {
                builder.append("none");
            } else if (end1 == Endpoint.CIRCLE) {
                builder.append("odot");
            }
            builder.append(", arrowhead=");
            if (end2 == Endpoint.ARROW) {
                builder.append("normal");
            } else if (end2 == Endpoint.TAIL) {
                builder.append("none");
            } else if (end2 == Endpoint.CIRCLE) {
                builder.append("odot");
            }

            // Bootstrapping
            List<EdgeTypeProbability> edgeTypeProbabilities = edge.getEdgeTypeProbabilities();
            if (edgeTypeProbabilities != null && !edgeTypeProbabilities.isEmpty()) {
                StringBuilder label = new StringBuilder(n1 + " - " + n2);
                for (EdgeTypeProbability edgeTypeProbability : edgeTypeProbabilities) {
                    EdgeType edgeType = edgeTypeProbability.getEdgeType();
                    double probability = edgeTypeProbability.getProbability();
                    if (probability > 0) {
                        StringBuilder edgeTypeString = new StringBuilder();
                        switch (edgeType) {
                            case nil:
                                edgeTypeString = new StringBuilder("no edge");
                                break;
                            case ta:
                                edgeTypeString = new StringBuilder("-->");
                                break;
                            case at:
                                edgeTypeString = new StringBuilder("&lt;--");
                                break;
                            case ca:
                                edgeTypeString = new StringBuilder("o->");
                                break;
                            case ac:
                                edgeTypeString = new StringBuilder("&lt;-o");
                                break;
                            case cc:
                                edgeTypeString = new StringBuilder("o-o");
                                break;
                            case aa:
                                edgeTypeString = new StringBuilder("&lt;-&gt;");
                                break;
                            case tt:
                                edgeTypeString = new StringBuilder("---");
                                break;
                        }

                        List<Property> properties = edgeTypeProbability.getProperties();
                        if (properties != null && properties.size() > 0) {
                            for (Property property : properties) {
                                edgeTypeString.append(" ").append(property.toString());
                            }
                        }

                        label.append("\\n[").append(edgeTypeString).append("]:").append(edgeTypeProbability.getProbability());
                    }
                }
                builder.append(", label=\"").append(label).append("\", fontname=courier");
            }

            builder.append("]; \n");
        }
        builder.append("}");

        return builder.toString();
    }

    public static void graphToDot(Graph graph, File file) {
        try {
            Writer writer = new FileWriter(file);
            writer.write(graphToDot(graph));
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * @return an XML element representing the given graph. (Well, only a basic
     * graph for now...)
     */
    public static Element convertToXml(Graph graph) {
        Element element = new Element("graph");

        Element variables = new Element("variables");
        element.appendChild(variables);

        for (Node node : graph.getNodes()) {
            Element variable = new Element("variable");
            Text text = new Text(node.getName());
            variable.appendChild(text);
            variables.appendChild(variable);
        }

        Element edges = new Element("edges");
        element.appendChild(edges);

        for (Edge edge : graph.getEdges()) {
            Element _edge = new Element("edge");
            Text text = new Text(edge.toString());
            _edge.appendChild(text);
            edges.appendChild(_edge);
        }

        Set<Triple> ambiguousTriples = graph.getAmbiguousTriples();

        if (!ambiguousTriples.isEmpty()) {
            Element underlinings = new Element("ambiguities");
            element.appendChild(underlinings);

            for (Triple triple : ambiguousTriples) {
                Element underlining = new Element("ambiguities");
                Text text = new Text(niceTripleString(triple));
                underlining.appendChild(text);
                underlinings.appendChild(underlining);
            }
        }

        Set<Triple> underlineTriples = graph.getUnderLines();

        if (!underlineTriples.isEmpty()) {
            Element underlinings = new Element("underlines");
            element.appendChild(underlinings);

            for (Triple triple : underlineTriples) {
                Element underlining = new Element("underline");
                Text text = new Text(niceTripleString(triple));
                underlining.appendChild(text);
                underlinings.appendChild(underlining);
            }
        }

        Set<Triple> dottedTriples = graph.getDottedUnderlines();

        if (!dottedTriples.isEmpty()) {
            Element dottedUnderlinings = new Element("dottedUnderlines");
            element.appendChild(dottedUnderlinings);

            for (Triple triple : dottedTriples) {
                Element dottedUnderlining = new Element("dottedUnderline");
                Text text = new Text(niceTripleString(triple));
                dottedUnderlining.appendChild(text);
                dottedUnderlinings.appendChild(dottedUnderlining);
            }
        }

        return element;
    }

    private static String niceTripleString(Triple triple) {
        return triple.getX() + ", " + triple.getY() + ", " + triple.getZ();
    }

    public static String graphToXml(Graph graph) {
        Document document = new Document(convertToXml(graph));
        OutputStream out = new ByteArrayOutputStream();
        Serializer serializer = new Serializer(out);
        serializer.setLineSeparator("\n");
        serializer.setIndent(2);

        try {
            serializer.write(document);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return out.toString();
    }

    public static Graph parseGraphXml(Element graphElement, Map<String, Node> nodes) throws ParsingException {
        if (!"graph".equals(graphElement.getLocalName())) {
            throw new IllegalArgumentException("Expecting graph element: " + graphElement.getLocalName());
        }

        if (!("variables".equals(graphElement.getChildElements().get(0).getLocalName()))) {
            throw new ParsingException("Expecting variables element: " + graphElement.getChildElements().get(0).getLocalName());
        }

        Element variablesElement = graphElement.getChildElements().get(0);
        Elements variableElements = variablesElement.getChildElements();
        List<Node> variables = new ArrayList<>();

        for (int i = 0; i < variableElements.size(); i++) {
            Element variableElement = variableElements.get(i);

            if (!("variable".equals(variablesElement.getChildElements().get(i).getLocalName()))) {
                throw new ParsingException("Expecting variable element.");
            }

            String value = variableElement.getValue();

            if (nodes == null) {
                variables.add(new GraphNode(value));
            } else {
                variables.add(nodes.get(value));
            }
        }

        Graph graph = new EdgeListGraph(variables);

//        graphNotes.add(noteAttribute.getValue());
        if (!("edges".equals(graphElement.getChildElements().get(1).getLocalName()))) {
            throw new ParsingException("Expecting edges element.");
        }

        Element edgesElement = graphElement.getChildElements().get(1);
        Elements edgesElements = edgesElement.getChildElements();

        for (int i = 0; i < edgesElements.size(); i++) {
            Element edgeElement = edgesElements.get(i);

            if (!("edge".equals(edgeElement.getLocalName()))) {
                throw new ParsingException("Expecting edge element: " + edgeElement.getLocalName());
            }

            String value = edgeElement.getValue();

            final String regex = "([A-Za-z0-9_-]*:?[A-Za-z0-9_-]*) ?(.)-(.) ?([A-Za-z0-9_-]*:?[A-Za-z0-9_-]*)";
//            String regex = "([A-Za-z0-9_-]*) ?([<o])-([o>]) ?([A-Za-z0-9_-]*)";

            java.util.regex.Pattern pattern = java.util.regex.Pattern.compile(regex);
            Matcher matcher = pattern.matcher(value);

            if (!matcher.matches()) {
                throw new ParsingException("Edge doesn't match pattern.");
            }

            String var1 = matcher.group(1);
            String leftEndpoint = matcher.group(2);
            String rightEndpoint = matcher.group(3);
            String var2 = matcher.group(4);

            Node node1 = graph.getNode(var1);
            Node node2 = graph.getNode(var2);
            Endpoint endpoint1;

            switch (leftEndpoint) {
                case "<":
                    endpoint1 = Endpoint.ARROW;
                    break;
                case "o":
                    endpoint1 = Endpoint.CIRCLE;
                    break;
                case "-":
                    endpoint1 = Endpoint.TAIL;
                    break;
                default:
                    throw new IllegalStateException("Expecting an endpoint: " + leftEndpoint);
            }

            Endpoint endpoint2;

            switch (rightEndpoint) {
                case ">":
                    endpoint2 = Endpoint.ARROW;
                    break;
                case "o":
                    endpoint2 = Endpoint.CIRCLE;
                    break;
                case "-":
                    endpoint2 = Endpoint.TAIL;
                    break;
                default:
                    throw new IllegalStateException("Expecting an endpoint: " + rightEndpoint);
            }

            Edge edge = new Edge(node1, node2, endpoint1, endpoint2);
            graph.addEdge(edge);
        }

        int size = graphElement.getChildElements().size();
        if (2 >= size) {
            return graph;
        }

        int p = 2;

        if ("ambiguities".equals(graphElement.getChildElements().get(p).getLocalName())) {
            Element ambiguitiesElement = graphElement.getChildElements().get(p);
            Set<Triple> triples = parseTriples(variables, ambiguitiesElement, "ambiguity");
            graph.setAmbiguousTriples(triples);
            p++;
        }

        if (p >= size) {
            return graph;
        }

        if ("underlines".equals(graphElement.getChildElements().get(p).getLocalName())) {
            Element ambiguitiesElement = graphElement.getChildElements().get(p);
            Set<Triple> triples = parseTriples(variables, ambiguitiesElement, "underline");
            graph.setUnderLineTriples(triples);
            p++;
        }

        if (p >= size) {
            return graph;
        }

        if ("dottedunderlines".equals(graphElement.getChildElements().get(p).getLocalName())) {
            Element ambiguitiesElement = graphElement.getChildElements().get(p);
            Set<Triple> triples = parseTriples(variables, ambiguitiesElement, "dottedunderline");
            graph.setDottedUnderLineTriples(triples);
        }

        return graph;
    }

    /**
     * A triples element has a list of three (comman separated) nodes as text.
     */
    private static Set<Triple> parseTriples(List<Node> variables, Element triplesElement, String s) {
        Elements elements = triplesElement.getChildElements(s);

        Set<Triple> triples = new HashSet<>();

        for (int q = 0; q < elements.size(); q++) {
            Element tripleElement = elements.get(q);
            String value = tripleElement.getValue();

            String[] tokens = value.split(",");

            if (tokens.length != 3) {
                throw new IllegalArgumentException("Expecting a triple: " + value);
            }

            String x = tokens[0].trim();
            String y = tokens[1].trim();
            String z = tokens[2].trim();

            Node _x = getNode(variables, x);
            Node _y = getNode(variables, y);
            Node _z = getNode(variables, z);

            Triple triple = new Triple(_x, _y, _z);
            triples.add(triple);
        }
        return triples;
    }

    private static Node getNode(List<Node> nodes, String x) {
        for (Node node : nodes) {
            if (x.equals(node.getName())) {
                return node;
            }
        }

        return null;
    }

    public static Element getRootElement(File file) throws ParsingException, IOException {
        Builder builder = new Builder();
        Document document = builder.build(file);
        return document.getRootElement();
    }

    /**
     * @param graph The graph to be saved.
     * @param file  The file to save it in.
     * @param xml   True if to be saved in XML, false if in text.
     * @return I have no idea whey I'm returning this; it's already closed...
     */
    public static PrintWriter saveGraph(Graph graph, File file, boolean xml) {
        PrintWriter out;

        try {
            out = new PrintWriter(new FileOutputStream(file));
//            out.print(graph);

            if (xml) {
                out.print(graphToXml(graph));
            } else {
                out.print(graph);
            }
            out.flush();
            out.close();
        } catch (IOException e1) {
            throw new IllegalArgumentException("Output file could not " + "be opened: " + file);
        }
        return out;
    }

    public static Graph loadGraph(File file) {

        Element root;
        Graph graph;

        try {
            root = getRootElement(file);
            graph = parseGraphXml(root, null);
        } catch (ParsingException e1) {
            throw new IllegalArgumentException("Could not parse " + file, e1);
        } catch (IOException e1) {
            throw new IllegalArgumentException("Could not read " + file, e1);
        }

        return graph;
    }

    public static Graph loadGraphTxt(File file) {
        try {
            Reader in1 = new FileReader(file);
            return readerToGraphTxt(in1);

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException();
        }
    }

    public static Graph loadGraphRuben(File file) {
        try {
            final String commentMarker = "//";
            final char quoteCharacter = '"';
            final String missingValueMarker = "*";
            final boolean hasHeader = false;

            DataSet dataSet = DataUtils.loadContinuousData(file, commentMarker, quoteCharacter, missingValueMarker, hasHeader, Delimiter.TAB);

            List<Node> nodes = dataSet.getVariables();
            Graph graph = new EdgeListGraph(nodes);

            for (int i = 0; i < nodes.size(); i++) {
                for (int j = i + 1; j < nodes.size(); j++) {
                    if (dataSet.getDouble(i, j) != 0D) {
                        graph.addDirectedEdge(nodes.get(i), nodes.get(j));
                    }
                }
            }

            return graph;

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException();
        }
    }

    public static Graph loadGraphJson(File file) {
        try {
            Reader in1 = new FileReader(file);
            return readerToGraphJson(in1);
        } catch (Exception e) {
            e.printStackTrace();
        }

        throw new IllegalStateException();
    }

    public static Graph readerToGraphTxt(String graphString) throws IOException {
        return readerToGraphTxt(new CharArrayReader(graphString.toCharArray()));
    }

    public static Graph readerToGraphTxt(Reader reader) throws IOException {
        Graph graph = new EdgeListGraph();
        try (BufferedReader in = new BufferedReader(reader)) {
            for (String line = in.readLine(); line != null; line = in.readLine()) {
                line = line.trim();
                switch (line) {
                    case "Graph Nodes:":
                        extractGraphNodes(graph, in);
                        break;
                    case "Graph Edges:":
                        extractGraphEdges(graph, in);
                        break;
                }
            }
        }

        return graph;
    }

    public static Graph readerToGraphRuben(Reader reader) throws IOException {
        Graph graph = new EdgeListGraph();
        try (BufferedReader in = new BufferedReader(reader)) {
            for (String line = in.readLine(); line != null; line = in.readLine()) {
                line = line.trim();
                switch (line) {
                    case "Graph Nodes:":
                        extractGraphNodes(graph, in);
                        break;
                    case "Graph Edges:":
                        extractGraphEdges(graph, in);
                        break;
                }
            }
        }

        return graph;
    }

    private static void extractGraphEdges(Graph graph, BufferedReader in) throws IOException {
        for (String line = in.readLine(); line != null; line = in.readLine()) {
            line = line.trim();

            if (line.isEmpty()) {
                continue;
            }

            String[] tokens = line.split("\\s+");

            String number = tokens[0];

            String[] tokensa = number.split("\\.");

            int fromIndex;

            try {
                Integer.parseInt(tokensa[0]);
                fromIndex = 1;
            } catch (NumberFormatException e) {
                fromIndex = 0;
            }

            String from = tokens[fromIndex];

            line = line.substring(line.indexOf(from) + from.length()).trim();
            tokens = line.split("\\s+");

            String edge = tokens[0];

            if ("Attributes:".equals(edge)) break;

            line = line.substring(line.indexOf(edge) + edge.length()).trim();
            tokens = line.split("\\s+");

            String to = tokens[0];
            line = line.substring(line.indexOf(to) + to.length()).trim();

            Node _from = graph.getNode(from);
            Node _to = graph.getNode(to);

            if (_from == null) {
                graph.addNode(new GraphNode(from));
                _from = graph.getNode(from);
            }

            if (_to == null) {
                graph.addNode(new GraphNode(to));
                _to = graph.getNode(to);
            }

            char end1 = edge.charAt(0);
            char end2 = edge.charAt(2);

            Endpoint _end1;
            Endpoint _end2;

            if (end1 == '<') {
                _end1 = Endpoint.ARROW;
            } else if (end1 == 'o') {
                _end1 = Endpoint.CIRCLE;
            } else if (end1 == '-') {
                _end1 = Endpoint.TAIL;
            } else {
                throw new IllegalArgumentException("Unrecognized endpoint: " + end1 + ", for edge " + edge);
            }

            if (end2 == '>') {
                _end2 = Endpoint.ARROW;
            } else if (end2 == 'o') {
                _end2 = Endpoint.CIRCLE;
            } else if (end2 == '-') {
                _end2 = Endpoint.TAIL;
            } else {
                throw new IllegalArgumentException("Unrecognized endpoint: " + end2 + ", for edge " + edge);
            }

            Edge _edge = new Edge(_from, _to, _end1, _end2);

            //Bootstrapping
            if (line.contains("[no edge]") || line.contains(" --> ") || line.contains(" &lt;-- ") || line.contains(" o-> ") || line.contains(" &lt;-o ") || line.contains(" o-o ") || line.contains(" &lt;-&gt; ") || line.contains(" --- ")) {

                // String bootstrap_format = "[no edge]:0.0000;[n1 --> n2]:0.0000;[n1 &lt;-- n2]:0.0000;[n1 o-> n2]:0.0000;[n1 &lt;-o n2]:0.0000;[n1 o-o n2]:0.0000;[n1 &lt;-> n2]:0.0000;[n1 --- n2]:0.0000;";
                int last_semicolon = line.lastIndexOf(";");
                String bootstraps;
                if (last_semicolon != -1) {
                    bootstraps = line.substring(0, last_semicolon + 1);
                } else {
                    bootstraps = line;
                }

                line = line.substring(bootstraps.length()).trim();

                String[] bootstrap = bootstraps.split(";");
                for (String s : bootstrap) {
                    String[] token = s.split(":");
                    if (token.length < 2) {
                        continue;
                    }

                    String orient = token[0];
                    double prob = Double.parseDouble(token[1]);

                    if (orient.equalsIgnoreCase("[no edge]")) {
                        _edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeType.nil, prob));
                    } else {
                        orient = orient.replace("[", "").replace("]", "");
                        EdgeTypeProbability etp;
                        if (orient.contains(" --> ")) {
                            etp = new EdgeTypeProbability(EdgeType.ta, prob);
                        } else if (orient.contains(" &lt;-- ")) {
                            etp = new EdgeTypeProbability(EdgeType.at, prob);
                        } else if (orient.contains(" o-> ")) {
                            etp = new EdgeTypeProbability(EdgeType.ca, prob);
                        } else if (orient.contains(" &lt;-o ")) {
                            etp = new EdgeTypeProbability(EdgeType.ac, prob);
                        } else if (orient.contains(" o-o ")) {
                            etp = new EdgeTypeProbability(EdgeType.cc, prob);
                        } else if (orient.contains(" &lt;-> ")) {
                            etp = new EdgeTypeProbability(EdgeType.aa, prob);
                        } else {// [n1 --- n2]
                            etp = new EdgeTypeProbability(EdgeType.tt, prob);
                        }
                        String[] _edge_property = orient.trim().split("\\s+");
                        if (_edge_property.length > 3) {
                            for (int j = 3; j < _edge_property.length; j++) {
                                etp.addProperty(Property.valueOf(_edge_property[j]));
                            }
                        }
                        _edge.addEdgeTypeProbability(etp);
                    }

                }
            }

            if (line.length() > 0) {
                tokens = line.split("\\s+");

                for (String token : tokens) {
                    _edge.addProperty(Property.valueOf(token));
                }
            }

            graph.addEdge(_edge);
        }
    }

    private static void extractGraphNodes(Graph graph, BufferedReader in) throws IOException {
        for (String line = in.readLine(); line != null; line = in.readLine()) {
            line = line.trim();
            if (line.isEmpty()) {
                break;
            }

            Arrays.stream(line.split("[,;]")).map(GraphNode::new).forEach(graph::addNode);
        }
    }

    public static Graph readerToGraphJson(Reader reader) throws IOException {
        BufferedReader in = new BufferedReader(reader);

        StringBuilder json = new StringBuilder();
        String line;

        while ((line = in.readLine()) != null) {
            json.append(line.trim());
        }

        return JsonUtils.parseJSONObjectToTetradGraph(json.toString());
    }

    public static Graph loadGraphGcpCausaldag(File file) {
        System.out.println("KK " + file.getAbsolutePath());
        File parentFile = file.getParentFile().getParentFile();
        parentFile = new File(parentFile, "data");
        File dataFile = new File(parentFile, file.getName().replace("causaldag.gsp", "data"));

        System.out.println(dataFile.getAbsolutePath());

        List<Node> variables = null;

        try {
            ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(dataFile.toPath(), Delimiter.TAB);
            Data data = reader.readInData();

            DataSet dataSet = (DataSet) DataConvertUtils.toDataModel(data);

            variables = dataSet.getVariables();
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            Reader in1 = new FileReader(file);
            return GraphUtils.readerToGraphCausaldag(in1, variables);

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException();
        }
    }

    public static Graph readerToGraphCausaldag(Reader reader, List<Node> variables) throws IOException {
        Graph graph = new EdgeListGraph(variables);
        try (BufferedReader in = new BufferedReader(reader)) {
            for (String line = in.readLine(); line != null; line = in.readLine()) {
                line = line.trim();

                String[] tokens = line.split("[\\[\\]]");

                for (String t : tokens) {
//                    System.out.println(t);

                    String[] tokens2 = t.split("[,|]");

                    if (tokens2[0].isEmpty()) continue;

                    Node x = variables.get(Integer.parseInt(tokens2[0]));

                    for (int j = 1; j < tokens2.length; j++) {
                        if (tokens2[j].isEmpty()) continue;

                        Node y = variables.get(Integer.parseInt(tokens2[j]));

                        graph.addDirectedEdge(y, x);
                    }
                }
            }
        }

        return graph;
    }

    public static HashMap<String, PointXy> grabLayout(List<Node> nodes) {
        HashMap<String, PointXy> layout = new HashMap<>();

        for (Node node : nodes) {
            layout.put(node.getName(), new PointXy(node.getCenterX(), node.getCenterY()));
        }

        return layout;
    }

    /**
     * @return A list of triples of the form X*-&gt;Y&lt;-*Z.
     */
    public static List<Triple> getCollidersFromGraph(Node node, Graph graph) {
        List<Triple> colliders = new ArrayList<>();

        List<Node> adj = graph.getAdjacentNodes(node);
        if (adj.size() < 2) {
            return new LinkedList<>();
        }

        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
        int[] choice;

        while ((choice = gen.next()) != null) {
            Node x = adj.get(choice[0]);
            Node z = adj.get(choice[1]);

            Endpoint endpt1 = graph.getEdge(x, node).getProximalEndpoint(node);
            Endpoint endpt2 = graph.getEdge(z, node).getProximalEndpoint(node);

            if (endpt1 == Endpoint.ARROW && endpt2 == Endpoint.ARROW) {
                colliders.add(new Triple(x, node, z));
            }
        }

        return colliders;
    }

    /**
     * @return A list of triples of the form X, Y, Z, where X, Y, Z is a
     * definite noncollider in the given graph.
     */
    public static List<Triple> getNoncollidersFromGraph(Node node, Graph graph) {
        List<Triple> noncolliders = new ArrayList<>();

        List<Node> adj = graph.getAdjacentNodes(node);
        if (adj.size() < 2) {
            return new LinkedList<>();
        }

        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
        int[] choice;

        while ((choice = gen.next()) != null) {
            Node x = adj.get(choice[0]);
            Node z = adj.get(choice[1]);

            Endpoint endpt1 = graph.getEdge(x, node).getProximalEndpoint(node);
            Endpoint endpt2 = graph.getEdge(z, node).getProximalEndpoint(node);

            if (endpt1 == Endpoint.ARROW && endpt2 == Endpoint.TAIL || endpt1 == Endpoint.TAIL && endpt2 == Endpoint.ARROW || endpt1 == Endpoint.TAIL && endpt2 == Endpoint.TAIL) {
                noncolliders.add(new Triple(x, node, z));
            }
        }

        return noncolliders;
    }

    /**
     * @return A list of triples of the form &lt;X, Y, Z&gt;, where &lt;X, Y, Z&gt; is a
     * definite noncollider in the given graph.
     */
    public static List<Triple> getAmbiguousTriplesFromGraph(Node node, Graph graph) {
        List<Triple> ambiguousTriples = new ArrayList<>();

        List<Node> adj = graph.getAdjacentNodes(node);
        if (adj.size() < 2) {
            return new LinkedList<>();
        }

        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
        int[] choice;

        while ((choice = gen.next()) != null) {
            Node x = adj.get(choice[0]);
            Node z = adj.get(choice[1]);

            if (graph.isAmbiguousTriple(x, node, z)) {
                ambiguousTriples.add(new Triple(x, node, z));
            }
        }

        return ambiguousTriples;
    }

    /**
     * @return A list of triples of the form &lt;X, Y, Z&gt;, where &lt;X, Y, Z&gt; is a
     * definite noncollider in the given graph.
     */
    public static List<Triple> getUnderlinedTriplesFromGraph(Node node, Graph graph) {
        List<Triple> underlinedTriples = new ArrayList<>();
        Set<Triple> allUnderlinedTriples = graph.getUnderLines();

        List<Node> adj = graph.getAdjacentNodes(node);
        if (adj.size() < 2) {
            return new LinkedList<>();
        }

        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
        int[] choice;

        while ((choice = gen.next()) != null) {
            Node x = adj.get(choice[0]);
            Node z = adj.get(choice[1]);

            if (allUnderlinedTriples.contains(new Triple(x, node, z))) {
                underlinedTriples.add(new Triple(x, node, z));
            }
        }

        return underlinedTriples;
    }

    /**
     * @return A list of triples of the form &lt;X, Y, Z&gt;, where &lt;X, Y, Z&gt; is a
     * definite noncollider in the given graph.
     */
    public static List<Triple> getDottedUnderlinedTriplesFromGraph(Node node, Graph graph) {
        List<Triple> dottedUnderlinedTriples = new ArrayList<>();
        Set<Triple> allDottedUnderlinedTriples = graph.getDottedUnderlines();

        List<Node> adj = graph.getAdjacentNodes(node);
        if (adj.size() < 2) {
            return new LinkedList<>();
        }

        ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
        int[] choice;

        while ((choice = gen.next()) != null) {
            Node x = adj.get(choice[0]);
            Node z = adj.get(choice[1]);

            if (allDottedUnderlinedTriples.contains(new Triple(x, node, z))) {
                dottedUnderlinedTriples.add(new Triple(x, node, z));
            }
        }

        return dottedUnderlinedTriples;
    }

    /**
     * A standard matrix graph representation for directed graphs. a[i][j] = 1
     * is j-->i and -1 if i-->j.
     *
     * @throws IllegalArgumentException if <code>graph</code> is not a directed
     *                                  graph.
     */
    private static int[][] incidenceMatrix(Graph graph) throws IllegalArgumentException {
        List<Node> nodes = graph.getNodes();
        int[][] m = new int[nodes.size()][nodes.size()];

        for (Edge edge : graph.getEdges()) {
            if (!Edges.isDirectedEdge(edge)) {
                throw new IllegalArgumentException("Not a directed graph.");
            }
        }

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                Node x1 = nodes.get(i);
                Node x2 = nodes.get(j);
                Edge edge = graph.getEdge(x1, x2);

                if (edge == null) {
                    m[i][j] = 0;
                } else if (edge.getProximalEndpoint(x1) == Endpoint.ARROW) {
                    m[i][j] = 1;
                } else if (edge.getProximalEndpoint(x1) == Endpoint.TAIL) {
                    m[i][j] = -1;
                }
            }
        }

        return m;
    }

    public static String loadGraphRMatrix(Graph graph) throws IllegalArgumentException {
        int[][] m = GraphUtils.incidenceMatrix(graph);

        TextTable table = new TextTable(m[0].length + 1, m.length + 1);

        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                table.setToken(i + 1, j + 1, m[i][j] + "");
            }
        }

        for (int i = 0; i < m.length; i++) {
            table.setToken(i + 1, 0, (i + 1) + "");
        }

        List<Node> nodes = graph.getNodes();

        for (int j = 0; j < m[0].length; j++) {
            table.setToken(0, j + 1, nodes.get(j).getName());
        }

        return table.toString();
    }

    public static Graph loadGraphPcAlgMatrix(DataSet dataSet) {
        List<Node> vars = dataSet.getVariables();

        Graph graph = new EdgeListGraph(vars);

        for (int i = 0; i < dataSet.getNumRows(); i++) {
            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                if (i == j) {
                    continue;
                }
                int g = dataSet.getInt(i, j);
                int h = dataSet.getInt(j, i);

                if (g == 1 && h == 1 && !graph.isAdjacentTo(vars.get(i), vars.get(j))) {
                    graph.addUndirectedEdge(vars.get(i), vars.get(j)); //
                } else if (g == 1 && h == 0) {
                    graph.addDirectedEdge(vars.get(j), vars.get(i));
                }
            }
        }

        return graph;
    }

    public static Graph loadGraphBNTPcMatrix(List<Node> vars, DataSet dataSet) {
        Graph graph = new EdgeListGraph(vars);

        for (int i = 0; i < dataSet.getNumRows(); i++) {
            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                int g = dataSet.getInt(i, j);
                int h = dataSet.getInt(j, i);

                if (g == 1 && h == 1 && !graph.isAdjacentTo(vars.get(i), vars.get(j))) {
                    graph.addUndirectedEdge(vars.get(i), vars.get(j));
                } else if (g == -1 && h == 0) {
                    graph.addDirectedEdge(vars.get(i), vars.get(j));
                }
            }
        }

        return graph;
    }

    public static String graphRMatrixTxt(Graph graph) throws IllegalArgumentException {
        int[][] m = GraphUtils.incidenceMatrix(graph);

        TextTable table = new TextTable(m[0].length + 1, m.length + 1);

        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                table.setToken(i + 1, j + 1, m[i][j] + "");
            }
        }

        for (int i = 0; i < m.length; i++) {
            table.setToken(i + 1, 0, (i + 1) + "");
        }

        List<Node> nodes = graph.getNodes();

        for (int j = 0; j < m[0].length; j++) {
            table.setToken(0, j + 1, nodes.get(j).getName());
        }

        return table.toString();

    }

    public static boolean containsBidirectedEdge(Graph graph) {
        boolean containsBidirected = false;

        for (Edge edge : graph.getEdges()) {
            if (Edges.isBidirectedEdge(edge)) {
                containsBidirected = true;
                break;
            }
        }
        return containsBidirected;
    }

    public static boolean existsDirectedPathFromTo(Node node1, Node node2, Graph graph) {
        return GraphUtils.existsDirectedPathVisit(node1, node2, new LinkedList<>(), -1, graph);
    }

    public static boolean existsDirectedPathFromTo(Node node1, Node node2, int depth, Graph graph) {
        return node1 == node2 || GraphUtils.existsDirectedPathVisit(node1, node2, new LinkedList<>(), depth, graph);
    }

    private static boolean existsDirectedPathVisit(Node node1, Node node2, LinkedList<Node> path, int depth, Graph graph) {
        path.addLast(node1);

        if (depth == -1) depth = Integer.MAX_VALUE;

        if (path.size() >= depth) {
            return false;
        }

        for (Edge edge : graph.getEdges(node1)) {
            Node child = Edges.traverseDirected(node1, edge);

            if (graph.getEdges(node1, child).size() == 2) {
                return true;
            }

            if (child == null) {
                continue;
            }

            if (child == node2) {
                return true;
            }

            if (path.contains(child)) {
                continue;
            }

            if (GraphUtils.existsDirectedPathVisit(child, node2, path, depth, graph)) {
                return true;
            }
        }

        path.removeLast();
        return false;
    }

    public static LinkedList<Triple> listColliderTriples(Graph graph) {
        LinkedList<Triple> colliders = new LinkedList<>();

        for (Node node : graph.getNodes()) {
            List<Node> adj = graph.getAdjacentNodes(node);

            if (adj.size() < 2) {
                continue;
            }

            ChoiceGenerator gen = new ChoiceGenerator(adj.size(), 2);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> others = GraphUtils.asList(choice, adj);

                if (graph.isDefCollider(others.get(0), node, others.get(1))) {
                    colliders.add(new Triple(others.get(0), node, others.get(1)));
                }
            }
        }
        return colliders;
    }

    /**
     * Constructs a list of nodes from the given <code>nodes</code> list at the
     * given indices in that list.
     *
     * @param indices The indices of the desired nodes in <code>nodes</code>.
     * @param nodes   The list of nodes from which we select a sublist.
     * @return the The sublist selected.
     */
    public static List<Node> asList(int[] indices, List<Node> nodes) {
        List<Node> list = new LinkedList<>();

        for (int index : indices) {
            list.add(nodes.get(index));
        }

        return list;
    }

    public static Set<Node> asSet(int[] indices, List<Node> nodes) {
        Set<Node> set = new HashSet<>();

        for (int i : indices) {
            set.add(nodes.get(i));
        }

        return set;
    }

    public static int numDirectionalErrors(Graph result, Graph cpdag) {
        int count = 0;

        for (Edge edge : result.getEdges()) {
            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            Node _node1 = cpdag.getNode(node1.getName());
            Node _node2 = cpdag.getNode(node2.getName());

            Edge _edge = cpdag.getEdge(_node1, _node2);

            if (_edge == null) {
                continue;
            }

            if (Edges.isDirectedEdge(edge)) {
                if (_edge.pointsTowards(_node1)) {
                    count++;
                } else if (Edges.isUndirectedEdge(_edge)) {
                    count++;
                }
            }
        }

        return count;
    }

    public static int numBidirected(Graph result) {
        int numBidirected = 0;

        for (Edge edge : result.getEdges()) {
            if (Edges.isBidirectedEdge(edge)) {
                numBidirected++;
            }
        }

        return numBidirected;
    }

    public static int degree(Graph graph) {
        int maxDegree = 0;

        for (Node node : graph.getNodes()) {
            int n = graph.getEdges(node).size();
            if (n > maxDegree) {
                maxDegree = n;
            }
        }

        return maxDegree;
    }

    /**
     * Finds a causal order for the given graph that is follows the order
     * of the given initialorder as possible.
     * @param graph The graph to find a causal order for. Must be acyclic, though
     *              it need not be a DAG.
     * @param initialorder The order to try to get as close to as possible.
     * @return Such a causal order.
     */
    public static List<Node> getCausalOrdering(Graph graph, List<Node> initialorder) {
        if (graph.existsDirectedCycle()) {
            throw new IllegalArgumentException("Graph must be acyclic.");
        }

        List<Node> found = new ArrayList<>();
        boolean _found = true;

        while (_found) {
            _found = false;

            for (Node node : initialorder) {
                HashSet<Node> nodes = new HashSet<>(found);
                if (!nodes.contains(node) && nodes.containsAll(graph.getParents(node))) {
                    found.add(node);
                    _found = true;
                }
            }
        }

        return found;
    }

    public static String getIntersectionComparisonString(List<Graph> graphs) {
        if (graphs == null || graphs.isEmpty()) {
            return "";
        }

        StringBuilder b = GraphUtils.undirectedEdges(graphs);

        b.append(GraphUtils.directedEdges(graphs));

        return b.toString();
    }

    private static StringBuilder undirectedEdges(List<Graph> graphs) {
        List<Graph> undirectedGraphs = new ArrayList<>();

        for (Graph graph : graphs) {
            Graph graph2 = new EdgeListGraph(graph);
            graph2.reorientAllWith(Endpoint.TAIL);
            undirectedGraphs.add(graph2);
        }

        Map<String, Node> exemplars = new HashMap<>();

        for (Graph graph : undirectedGraphs) {
            for (Node node : graph.getNodes()) {
                exemplars.put(node.getName(), node);
            }
        }

        Set<Node> nodeSet = new HashSet<>();

        for (String s : exemplars.keySet()) {
            nodeSet.add(exemplars.get(s));
        }
        List<Node> nodes = new ArrayList<>(nodeSet);
        List<Graph> undirectedGraphs2 = new ArrayList<>();

        for (int i = 0; i < graphs.size(); i++) {
            Graph graph = GraphUtils.replaceNodes(undirectedGraphs.get(i), nodes);
            undirectedGraphs2.add(graph);
        }

        Set<Edge> undirectedEdgesSet = new HashSet<>();

        for (Graph graph : undirectedGraphs2) {
            undirectedEdgesSet.addAll(graph.getEdges());
        }

        List<Edge> undirectedEdges = new ArrayList<>(undirectedEdgesSet);

        undirectedEdges.sort((o1, o2) -> {
            String name11 = o1.getNode1().getName();
            String name12 = o1.getNode2().getName();
            String name21 = o2.getNode1().getName();
            String name22 = o2.getNode2().getName();

            int major = name11.compareTo(name21);
            int minor = name12.compareTo(name22);

            if (major == 0) {
                return minor;
            } else {
                return major;
            }
        });

        List<List<Edge>> groups = new ArrayList<>();
        for (int i = 0; i < graphs.size(); i++) {
            groups.add(new ArrayList<>());
        }

        for (Edge edge : undirectedEdges) {
            int count = 0;

            for (Graph graph : undirectedGraphs2) {
                if (graph.containsEdge(edge)) {
                    count++;
                }
            }

            if (count == 0) {
                throw new IllegalArgumentException();
            }

            groups.get(count - 1).add(edge);
        }

        StringBuilder b = new StringBuilder();

        for (int i = groups.size() - 1; i >= 0; i--) {
            b.append("\n\nIn ").append(i + 1).append(" graph").append((i > 0) ? "s" : "").append("...\n");

            for (int j = 0; j < groups.get(i).size(); j++) {
                b.append("\n").append(j + 1).append(". ").append(groups.get(i).get(j));
            }
        }

        return b;
    }

    private static StringBuilder directedEdges(List<Graph> directedGraphs) {
        Set<Edge> directedEdgesSet = new HashSet<>();

        Map<String, Node> exemplars = new HashMap<>();

        for (Graph graph : directedGraphs) {
            for (Node node : graph.getNodes()) {
                exemplars.put(node.getName(), node);
            }
        }

        Set<Node> nodeSet = new HashSet<>();

        for (String s : exemplars.keySet()) {
            nodeSet.add(exemplars.get(s));
        }

        List<Node> nodes = new ArrayList<>(nodeSet);

        List<Graph> directedGraphs2 = new ArrayList<>();

        for (Graph directedGraph : directedGraphs) {
            Graph graph = GraphUtils.replaceNodes(directedGraph, nodes);
            directedGraphs2.add(graph);
        }

        for (Graph graph : directedGraphs2) {
            directedEdgesSet.addAll(graph.getEdges());
        }

        List<Edge> directedEdges = new ArrayList<>(directedEdgesSet);

        directedEdges.sort((o1, o2) -> {
            String name11 = o1.getNode1().getName();
            String name12 = o1.getNode2().getName();
            String name21 = o2.getNode1().getName();
            String name22 = o2.getNode2().getName();

            int major = name11.compareTo(name21);
            int minor = name12.compareTo(name22);

            if (major == 0) {
                return minor;
            } else {
                return major;
            }
        });

        List<List<Edge>> groups = new ArrayList<>();
        for (int i = 0; i < directedGraphs2.size(); i++) {
            groups.add(new ArrayList<>());
        }
        Set<Edge> contradicted = new HashSet<>();
        Map<Edge, Integer> directionCounts = new HashMap<>();

        for (Edge edge : directedEdges) {
            if (!edge.isDirected()) {
                continue;
            }

            int count1 = 0;
            int count2 = 0;

            for (Graph graph : directedGraphs2) {
                if (graph.containsEdge(edge)) {
                    count1++;
                } else if (graph.containsEdge(edge.reverse())) {
                    count2++;
                }
            }

            if (count1 != 0 && count2 != 0 && !contradicted.contains(edge.reverse())) {
                contradicted.add(edge);
            }

            directionCounts.put(edge, count1);
            directionCounts.put(edge.reverse(), count2);

            if (count1 == 0) {
                groups.get(count2 - 1).add(edge);
            }

            if (count2 == 0) {
                groups.get(count1 - 1).add(edge);
            }
        }

        StringBuilder b = new StringBuilder();

        for (int i = groups.size() - 1; i >= 0; i--) {
            b.append("\n\nUncontradicted in ").append(i + 1).append(" graph").append((i > 0) ? "s" : "").append("...\n");

            for (int j = 0; j < groups.get(i).size(); j++) {
                b.append("\n").append(j + 1).append(". ").append(groups.get(i).get(j));
            }
        }

        b.append("\n\nContradicted:\n");
        int index = 1;

        for (Edge edge : contradicted) {
            b.append("\n").append(index++).append(". ").append(Edges.undirectedEdge(edge.getNode1(), edge.getNode2())).append(" (--> ").append(directionCounts.get(edge)).append(" &lt;-- ").append(directionCounts.get(edge.reverse())).append(")");
        }

        return b;
    }

    private static boolean uncontradicted(Edge edge1, Edge edge2) {
        if (edge1 == null || edge2 == null) {
            return true;
        }

        Node x = edge1.getNode1();
        Node y = edge1.getNode2();

        if (edge1.pointsTowards(x) && edge2.pointsTowards(y)) {
            return false;
        } else return !edge1.pointsTowards(y) || !edge2.pointsTowards(x);
    }

    public static String edgeMisclassifications(double[][] counts, NumberFormat nf) {
        StringBuilder builder = new StringBuilder();

        TextTable table2 = new TextTable(9, 7);

        table2.setToken(1, 0, "---");
        table2.setToken(2, 0, "o-o");
        table2.setToken(3, 0, "o->");
        table2.setToken(4, 0, "&lt;-o");
        table2.setToken(5, 0, "-->");
        table2.setToken(6, 0, "&lt;--");
        table2.setToken(7, 0, "&lt;->");
        table2.setToken(8, 0, "No Edge");
        table2.setToken(0, 1, "---");
        table2.setToken(0, 2, "o-o");
        table2.setToken(0, 3, "o->");
        table2.setToken(0, 4, "-->");
        table2.setToken(0, 5, "&lt;->");
        table2.setToken(0, 6, "No Edge");

        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 6; j++) {
                if (i == 7 && j == 5) {
                    table2.setToken(7 + 1, 5 + 1, "*");
                } else {
                    table2.setToken(i + 1, j + 1, "" + nf.format(counts[i][j]));
                }
            }
        }

        builder.append(table2);

        double correctEdges = 0;
        double estimatedEdges = 0;

        for (int i = 0; i < counts.length; i++) {
            for (int j = 0; j < counts[0].length - 1; j++) {
                if ((i == 0 && j == 0) || (i == 1 && j == 1) || (i == 2 && j == 2) || (i == 4 && j == 3) || (i == 6 && j == 4)) {
                    correctEdges += counts[i][j];
                }

                estimatedEdges += counts[i][j];
            }
        }

        NumberFormat nf2 = new DecimalFormat("0.00");

        builder.append("\nRatio correct edges to estimated edges = ").append(nf2.format((correctEdges / estimatedEdges)));

        return builder.toString();
    }

    public static String edgeMisclassifications(int[][] counts) {
        StringBuilder builder = new StringBuilder();

        TextTable table2 = new TextTable(9, 7);

        table2.setToken(1, 0, "---");
        table2.setToken(2, 0, "o-o");
        table2.setToken(3, 0, "o->");
        table2.setToken(4, 0, "&lt;-o");
        table2.setToken(5, 0, "-->");
        table2.setToken(6, 0, "<--");
        table2.setToken(7, 0, "<->");
        table2.setToken(8, 0, "No Edge");
        table2.setToken(0, 1, "---");
        table2.setToken(0, 2, "o-o");
        table2.setToken(0, 3, "o->");
        table2.setToken(0, 4, "-->");
        table2.setToken(0, 5, "<->");
        table2.setToken(0, 6, "No Edge");

        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 6; j++) {
                if (i == 7 && j == 5) {
                    table2.setToken(7 + 1, 5 + 1, "*");
                } else {
                    table2.setToken(i + 1, j + 1, "" + counts[i][j]);
                }
            }
        }

        builder.append(table2);

        int correctEdges = 0;
        int estimatedEdges = 0;

        for (int i = 0; i < counts.length; i++) {
            for (int j = 0; j < counts[0].length - 1; j++) {
                if ((i == 0 && j == 0) || (i == 1 && j == 1) || (i == 2 && j == 2) || (i == 4 && j == 3) || (i == 6 && j == 4)) {
                    correctEdges += counts[i][j];
                }

                estimatedEdges += counts[i][j];
            }
        }

        NumberFormat nf2 = new DecimalFormat("0.00");

        builder.append("\nRatio correct edges to estimated edges = ").append(nf2.format((correctEdges / (double) estimatedEdges)));

        return builder.toString();
    }

    public static void addPagColoring(Graph graph) {
        for (Edge edge : graph.getEdges()) {
            if (!Edges.isDirectedEdge(edge)) {
                continue;
            }

            Node x = Edges.getDirectedEdgeTail(edge);
            Node y = Edges.getDirectedEdgeHead(edge);

            graph.removeEdge(edge);
            graph.addEdge(edge);

            Edge xyEdge = graph.getEdge(x, y);
            graph.removeEdge(xyEdge);

            if (!GraphUtils.existsSemiDirectedPath(x, y, graph)) {
                edge.addProperty(Property.dd); // green.
            } else {
                edge.addProperty(Property.pd);
            }

            graph.addEdge(xyEdge);

            if (graph.defVisible(edge)) {
                edge.addProperty(Property.nl); // bold.
            } else {
                edge.addProperty(Property.pl);
            }
        }
    }

    // Returns true if a path consisting of undirected and directed edges toward 'to' exists of
    // length at most 'bound'. Cycle checker in other words.
    public static boolean existsSemiDirectedPath(Node from, Node to, Graph graph) {
        Queue<Node> Q = new LinkedList<>();
        Set<Node> V = new HashSet<>();

        for (Node u : graph.getAdjacentNodes(from)) {
            Edge edge = graph.getEdge(from, u);
            Node c = GraphUtils.traverseSemiDirected(from, edge);

            if (c == null) {
                continue;
            }

            if (!V.contains(c)) {
                V.add(c);
                Q.offer(c);
            }
        }

        while (!Q.isEmpty()) {
            Node t = Q.remove();

            if (t == to) {
                return true;
            }

            for (Node u : graph.getAdjacentNodes(t)) {
                Edge edge = graph.getEdge(t, u);
                Node c = GraphUtils.traverseSemiDirected(t, edge);

                if (c == null) {
                    continue;
                }

                if (!V.contains(c)) {
                    V.add(c);
                    Q.offer(c);
                }
            }
        }

        return false;
    }

    public static int[][] edgeMisclassificationCounts(Graph leftGraph, Graph topGraph, boolean print) {
//        topGraph = GraphUtils.replaceNodes(topGraph, leftGraph.getNodes());

        class CountTask extends RecursiveTask<Counts> {

            private final List<Edge> edges;
            private final Graph leftGraph;
            private final Graph topGraph;
            private final Counts counts;
            private final int[] count;
            private final int chunk;
            private final int from;
            private final int to;

            public CountTask(int chunk, int from, int to, List<Edge> edges, Graph leftGraph, Graph topGraph, int[] count) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
                this.edges = edges;
                this.leftGraph = leftGraph;
                this.topGraph = topGraph;
                this.counts = new Counts();
                this.count = count;
            }

            @Override
            protected Counts compute() {
                int range = this.to - this.from;

                if (range <= this.chunk) {
                    for (int i = this.from; i < this.to; i++) {
                        int j = ++this.count[0];

                        Edge edge = this.edges.get(i);

                        Node x = edge.getNode1();
                        Node y = edge.getNode2();

                        Edge left = this.leftGraph.getEdge(x, y);
                        Edge top = this.topGraph.getEdge(x, y);

                        int m = GraphUtils.getTypeLeft(left, top);
                        int n = GraphUtils.getTypeTop(top);

                        this.counts.increment(m, n);
                    }

                    return this.counts;
                } else {
                    int mid = (this.to + this.from) / 2;
                    CountTask left = new CountTask(this.chunk, this.from, mid, this.edges, this.leftGraph, this.topGraph, this.count);
                    CountTask right = new CountTask(this.chunk, mid, this.to, this.edges, this.leftGraph, this.topGraph, this.count);

                    left.fork();
                    Counts rightAnswer = right.compute();
                    Counts leftAnswer = left.join();

                    leftAnswer.addAll(rightAnswer);
                    return leftAnswer;
                }
            }

            public Counts getCounts() {
                return this.counts;
            }
        }

        Set<Edge> edgeSet = new HashSet<>();
        edgeSet.addAll(topGraph.getEdges());
        edgeSet.addAll(leftGraph.getEdges());

//        System.out.println("Union formed");
        if (print) {
            System.out.println("Top graph " + topGraph.getEdges().size());
            System.out.println("Left graph " + leftGraph.getEdges().size());
            System.out.println("All edges " + edgeSet.size());
        }

        List<Edge> edges = new ArrayList<>(edgeSet);

//        System.out.println("Finding pool");
        ForkJoinPoolInstance pool = ForkJoinPoolInstance.getInstance();

//        System.out.println("Starting count task");
        CountTask task = new CountTask(500, 0, edges.size(), edges, leftGraph, topGraph, new int[1]);
        Counts counts = pool.getPool().invoke(task);

//        System.out.println("Finishing count task");
        return counts.countArray();
    }

    private static Set<Edge> complement(Set<Edge> edgeSet, Graph topGraph) {
        Set<Edge> complement = new HashSet<>(edgeSet);
        complement.removeAll(topGraph.getEdges());
        return complement;
    }

    private static int getTypeTop(Edge edgeTop) {
        if (edgeTop == null) {
            return 5;
        }

        if (Edges.isUndirectedEdge(edgeTop)) {
            return 0;
        }

        if (Edges.isNondirectedEdge(edgeTop)) {
            return 1;
        }

        if (Edges.isPartiallyOrientedEdge(edgeTop)) {
            return 2;
        }

        if (Edges.isDirectedEdge(edgeTop)) {
            return 3;
        }

        if (Edges.isBidirectedEdge(edgeTop)) {
            return 4;
        }

        return 5;

//        throw new IllegalArgumentException("Unsupported edge type : " + edgeTop);
    }

    private static int getTypeLeft(Edge edgeLeft, Edge edgeTop) {
        if (edgeLeft == null) {
            return 7;
        }

        if (edgeTop == null) {
            edgeTop = edgeLeft;
        }

        if (Edges.isUndirectedEdge(edgeLeft)) {
            return 0;
        }

        if (Edges.isNondirectedEdge(edgeLeft)) {
            return 1;
        }

        Node x = edgeLeft.getNode1();
        Node y = edgeLeft.getNode2();

        if (Edges.isPartiallyOrientedEdge(edgeLeft)) {
            if ((edgeLeft.pointsTowards(x) && edgeTop.pointsTowards(y)) || (edgeLeft.pointsTowards(y) && edgeTop.pointsTowards(x))) {
                return 3;
            } else {
                return 2;
            }
        }

        if (Edges.isDirectedEdge(edgeLeft)) {
            if ((edgeLeft.pointsTowards(x) && edgeTop.pointsTowards(y)) || (edgeLeft.pointsTowards(y) && edgeTop.pointsTowards(x))) {
                return 5;
            } else {
                return 4;
            }
        }

        if (Edges.isBidirectedEdge(edgeLeft)) {
            return 6;
        }

        throw new IllegalArgumentException("Unsupported edge type : " + edgeLeft);
    }

    private static int getTypeLeft2(Edge edgeLeft) {
        if (edgeLeft == null) {
            return 7;
        }

        if (Edges.isUndirectedEdge(edgeLeft)) {
            return 0;
        }

        if (Edges.isNondirectedEdge(edgeLeft)) {
            return 1;
        }

        if (Edges.isPartiallyOrientedEdge(edgeLeft)) {
            return 2;
        }

        if (Edges.isDirectedEdge(edgeLeft)) {
            return 4;
        }

        if (Edges.isBidirectedEdge(edgeLeft)) {
            return 6;
        }

        throw new IllegalArgumentException("Unsupported edge type : " + edgeLeft);
    }

    public static Set<Set<Node>> maximalCliques(Graph graph, List<Node> nodes) {
        Set<Set<Node>> report = new HashSet<>();
        GraphUtils.brokKerbosh1(new HashSet<>(), new HashSet<>(nodes), new HashSet<>(), report, graph);
        return report;
    }

    //        BronKerbosch1(R, P, X):
    //            if P and X are both empty:
    //                   report R as a maximal clique
    //            for each vertex v in P:
    //                   BronKerbosch1(R â {v}, P â N(v), X â N(v))
    //                   P := P \ {v}
    //                   X := X â {v}
    private static void brokKerbosh1(Set<Node> R, Set<Node> P, Set<Node> X, Set<Set<Node>> report, Graph graph) {
        if (P.isEmpty() && X.isEmpty()) {
            report.add(new HashSet<>(R));
        }

        for (Node v : new HashSet<>(P)) {
            Set<Node> _R = new HashSet<>(R);
            Set<Node> _P = new HashSet<>(P);
            Set<Node> _X = new HashSet<>(X);
            _R.add(v);
            _P.retainAll(graph.getAdjacentNodes(v));
            _X.retainAll(graph.getAdjacentNodes(v));
            GraphUtils.brokKerbosh1(_R, _P, _X, report, graph);
            P.remove(v);
            X.add(v);
        }
    }

    public static String graphToText(Graph graph) {
        // add edge properties relating to edge coloring of PAGs
        if (graph.isPag()) {
            GraphUtils.addPagColoring(graph);
        }

        Formatter fmt = new Formatter();
        fmt.format("%s%n%n", GraphUtils.graphNodesToText(graph, "Graph Nodes:", ';'));
        fmt.format("%s%n", GraphUtils.graphEdgesToText(graph, "Graph Edges:"));

        // Graph Attributes
        String graphAttributes = GraphUtils.graphAttributesToText(graph, "Graph Attributes:");
        if (graphAttributes != null) {
            fmt.format("%s%n", graphAttributes);
        }

        // Nodes Attributes
        String graphNodeAttributes = GraphUtils.graphNodeAttributesToText(graph, "Graph Node Attributes:", ';');
        if (graphNodeAttributes != null) {
            fmt.format("%s%n", graphNodeAttributes);
        }

        Set<Triple> ambiguousTriples = graph.getAmbiguousTriples();
        if (!ambiguousTriples.isEmpty()) {
            fmt.format("%n%n%s", GraphUtils.triplesToText(ambiguousTriples, "Ambiguous triples (i.e. list of triples for which there is ambiguous data about whether they are colliders or not):"));
        }

        Set<Triple> underLineTriples = graph.getUnderLines();
        if (!underLineTriples.isEmpty()) {
            fmt.format("%n%n%s", GraphUtils.triplesToText(underLineTriples, "Underline triples:"));
        }

        Set<Triple> dottedUnderLineTriples = graph.getDottedUnderlines();
        if (!dottedUnderLineTriples.isEmpty()) {
            fmt.format("%n%n%s", GraphUtils.triplesToText(dottedUnderLineTriples, "Dotted underline triples:"));
        }

        return fmt.toString();
    }

    public static String graphNodeAttributesToText(Graph graph, String title, char delimiter) {
        List<Node> nodes = graph.getNodes();

        Map<String, Map<String, Object>> graphNodeAttributes = new LinkedHashMap<>();
        for (Node node : nodes) {
            Map<String, Object> attributes = node.getAllAttributes();

            if (!attributes.isEmpty()) {
                for (String key : attributes.keySet()) {
                    Object value = attributes.get(key);

                    Map<String, Object> nodeAttributes = graphNodeAttributes.get(key);
                    if (nodeAttributes == null) {
                        nodeAttributes = new LinkedHashMap<>();
                    }
                    nodeAttributes.put(node.getName(), value);

                    graphNodeAttributes.put(key, nodeAttributes);
                }
            }
        }

        if (!graphNodeAttributes.isEmpty()) {
            StringBuilder sb = (title == null || title.length() == 0) ? new StringBuilder() : new StringBuilder(String.format("%s", title));

            for (String key : graphNodeAttributes.keySet()) {
                Map<String, Object> nodeAttributes = graphNodeAttributes.get(key);
                int size = nodeAttributes.size();
                int count = 0;

                sb.append(String.format("%n%s: [", key));

                for (String nodeName : nodeAttributes.keySet()) {
                    Object value = nodeAttributes.get(nodeName);

                    sb.append(String.format("%s: %s", nodeName, value));

                    count++;

                    if (count < size) {
                        sb.append(delimiter);
                    }

                }

                sb.append("]");
            }

            return sb.toString();
        }

        return null;
    }

    public static String graphAttributesToText(Graph graph, String title) {
        Map<String, Object> attributes = graph.getAllAttributes();
        if (!attributes.isEmpty()) {
            StringBuilder sb = (title == null || title.length() == 0) ? new StringBuilder() : new StringBuilder(String.format("%s%n", title));

            for (String key : attributes.keySet()) {
                Object value = attributes.get(key);

                sb.append(key);
                sb.append(": ");
                if (value instanceof String) {
                    sb.append(value);
                } else if (value instanceof Number) {
                    sb.append(String.format("%f%n", ((Number) value).doubleValue()));
                }
            }

            return sb.toString();
        }

        return null;
    }

    public static String graphNodesToText(Graph graph, String title, char delimiter) {
        StringBuilder sb = (title == null || title.length() == 0) ? new StringBuilder() : new StringBuilder(String.format("%s%n", title));

        List<Node> nodes = graph.getNodes();
        int size = nodes.size();
        int count = 0;
        for (Node node : nodes) {
            count++;
            sb.append(node.getName());
            if (count < size) {
                sb.append(delimiter);
            }
        }

        return sb.toString();
    }

    public static String graphEdgesToText(Graph graph, String title) {
        Formatter fmt = new Formatter();

        if (title != null && title.length() > 0) {
            fmt.format("%s%n", title);
        }

        List<Edge> edges = new ArrayList<>(graph.getEdges());

        Edges.sortEdges(edges);

        int count = 0;

        for (Edge edge : edges) {
            count++;

            // We will print edge's properties in the edge (via toString() function) level.
            //List<Edge.Property> properties = edge.getProperties();
            final String f = "%d. %s";
            Object[] o = new Object[2 /*+ properties.size()*/];/*+ properties.size()*/// <- here we include its properties (nl dd pl pd)
            o[0] = count;
            o[1] = edge; // <- here we include its properties (nl dd pl pd)
            fmt.format(f, o);
            fmt.format("\n");
        }

        return fmt.toString();
    }

    public static String triplesToText(Set<Triple> triples, String title) {
        Formatter fmt = new Formatter();

        if (title != null && title.length() > 0) {
            fmt.format("%s%n", title);
        }

        int size = (triples == null) ? 0 : triples.size();
        if (size > 0) {
            int count = 0;
            for (Triple triple : triples) {
                count++;
                if (count < size) {
                    fmt.format("%s%n", triple);
                } else {
                    fmt.format("%s", triple);
                }
            }
        }

        return fmt.toString();
    }

    public static TwoCycleErrors getTwoCycleErrors(Graph trueGraph, Graph estGraph) {
        Set<Edge> trueEdges = trueGraph.getEdges();
        Set<Edge> trueTwoCycle = new HashSet<>();

        for (Edge edge : trueEdges) {
            if (!edge.isDirected()) {
                continue;
            }

            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            if (trueEdges.contains(Edges.directedEdge(node2, node1))) {
                Edge undirEdge = Edges.undirectedEdge(node1, node2);
                trueTwoCycle.add(undirEdge);
            }
        }

        Set<Edge> estEdges = estGraph.getEdges();
        Set<Edge> estTwoCycle = new HashSet<>();

        for (Edge edge : estEdges) {
            if (!edge.isDirected()) {
                continue;
            }

            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            if (estEdges.contains(Edges.directedEdge(node2, node1))) {
                Edge undirEdge = Edges.undirectedEdge(node1, node2);
                estTwoCycle.add(undirEdge);
            }
        }

        Graph trueTwoCycleGraph = new EdgeListGraph(trueGraph.getNodes());

        for (Edge edge : trueTwoCycle) {
            trueTwoCycleGraph.addEdge(edge);
        }

        Graph estTwoCycleGraph = new EdgeListGraph(estGraph.getNodes());

        for (Edge edge : estTwoCycle) {
            estTwoCycleGraph.addEdge(edge);
        }

        estTwoCycleGraph = GraphUtils.replaceNodes(estTwoCycleGraph, trueTwoCycleGraph.getNodes());

        int adjFn = GraphUtils.countAdjErrors(trueTwoCycleGraph, estTwoCycleGraph);
        int adjFp = GraphUtils.countAdjErrors(estTwoCycleGraph, trueTwoCycleGraph);

        Graph undirectedGraph = GraphUtils.undirectedGraph(estTwoCycleGraph);
        int adjCorrect = undirectedGraph.getNumEdges() - adjFp;

        return new TwoCycleErrors(adjCorrect, adjFn, adjFp);
    }

    public static boolean isDConnectedTo(Node x, Node y, List<Node> z, Graph graph) {
        return GraphUtils.isDConnectedTo1(x, y, z, graph);
    }

    // Breadth first.
    private static boolean isDConnectedTo1(Node x, Node y, List<Node> z, Graph graph) {
        class EdgeNode {

            private final Edge edge;
            private final Node node;

            public EdgeNode(Edge edge, Node node) {
                this.edge = edge;
                this.node = node;
            }

            public int hashCode() {
                return this.edge.hashCode() + this.node.hashCode();
            }

            public boolean equals(Object o) {
                if (!(o instanceof EdgeNode)) {
                    throw new IllegalArgumentException();
                }
                EdgeNode _o = (EdgeNode) o;
                return _o.edge == this.edge && _o.node == this.node;
            }
        }

        Queue<EdgeNode> Q = new ArrayDeque<>();
        Set<EdgeNode> V = new HashSet<>();

        if (x == y) {
            return true;
        }

        for (Edge edge : graph.getEdges(x)) {
            if (edge.getDistalNode(x) == y) {
                return true;
            }
            EdgeNode edgeNode = new EdgeNode(edge, x);
            Q.offer(edgeNode);
            V.add(edgeNode);
        }

        while (!Q.isEmpty()) {
            EdgeNode t = Q.poll();

            Edge edge1 = t.edge;
            Node a = t.node;
            Node b = edge1.getDistalNode(a);

            for (Edge edge2 : graph.getEdges(b)) {
                Node c = edge2.getDistalNode(b);
                if (c == a) {
                    continue;
                }

                if (GraphUtils.reachable(edge1, edge2, a, z, graph)) {
                    if (c == y) {
                        return true;
                    }

                    EdgeNode u = new EdgeNode(edge2, b);

                    if (!V.contains(u)) {
                        V.add(u);
                        Q.offer(u);
                    }
                }
            }
        }

        return false;
    }

    public static boolean isDConnectedTo(List<Node> x, List<Node> y, List<Node> z, Graph graph) {
        Set<Node> zAncestors = GraphUtils.zAncestors(z, graph);

        Queue<OrderedPair<Node>> Q = new ArrayDeque<>();
        Set<OrderedPair<Node>> V = new HashSet<>();

        for (Node _x : x) {
            for (Node node : graph.getAdjacentNodes(_x)) {
                if (y.contains(node)) {
                    return true;
                }
                OrderedPair<Node> edge = new OrderedPair<>(_x, node);
                Q.offer(edge);
                V.add(edge);
            }
        }

        while (!Q.isEmpty()) {
            OrderedPair<Node> t = Q.poll();

            Node b = t.getFirst();
            Node a = t.getSecond();

            for (Node c : graph.getAdjacentNodes(b)) {
                if (c == a) {
                    continue;
                }

                boolean collider = graph.isDefCollider(a, b, c);
                if (!((collider && zAncestors.contains(b)) || (!collider && !z.contains(b)))) {
                    continue;
                }

                if (y.contains(c)) {
                    return true;
                }

                OrderedPair<Node> u = new OrderedPair<>(b, c);
                if (V.contains(u)) {
                    continue;
                }

                V.add(u);
                Q.offer(u);
            }
        }

        return false;
    }

    public static Set<Node> getDconnectedVars(Node y, List<Node> z, Graph graph) {
        Set<Node> Y = new HashSet<>();

        class EdgeNode {

            private final Edge edge;
            private final Node node;

            public EdgeNode(Edge edge, Node node) {
                this.edge = edge;
                this.node = node;
            }

            public int hashCode() {
                return this.edge.hashCode() + this.node.hashCode();
            }

            public boolean equals(Object o) {
                if (!(o instanceof EdgeNode)) {
                    throw new IllegalArgumentException();
                }
                EdgeNode _o = (EdgeNode) o;
                return _o.edge == this.edge && _o.node == this.node;
            }
        }

        Queue<EdgeNode> Q = new ArrayDeque<>();
        Set<EdgeNode> V = new HashSet<>();

        for (Edge edge : graph.getEdges(y)) {
            EdgeNode edgeNode = new EdgeNode(edge, y);
            Q.offer(edgeNode);
            V.add(edgeNode);
            Y.add(edge.getDistalNode(y));
        }

        while (!Q.isEmpty()) {
            EdgeNode t = Q.poll();

            Edge edge1 = t.edge;
            Node a = t.node;
            Node b = edge1.getDistalNode(a);

            for (Edge edge2 : graph.getEdges(b)) {
                Node c = edge2.getDistalNode(b);
                if (c == a) {
                    continue;
                }

                if (GraphUtils.reachable(edge1, edge2, a, z, graph)) {
                    EdgeNode u = new EdgeNode(edge2, b);

                    if (!V.contains(u)) {
                        V.add(u);
                        Q.offer(u);
                        Y.add(c);
                    }
                }
            }
        }

        return Y;
    }

    // Depth first.
    public static boolean isDConnectedTo2(Node x, Node y, List<Node> z, Graph graph) {
        LinkedList<Node> path = new LinkedList<>();

        path.add(x);

        for (Node c : graph.getAdjacentNodes(x)) {
            if (GraphUtils.isDConnectedToVisit2(x, c, y, path, z, graph)) {
                return true;
            }
        }

        return false;
    }

    public static boolean isDConnectedTo2(Node x, Node y, List<Node> z, Graph graph, LinkedList<Node> path) {
        path.add(x);

        for (Node c : graph.getAdjacentNodes(x)) {
            if (GraphUtils.isDConnectedToVisit2(x, c, y, path, z, graph)) {
                return true;
            }
        }

        return false;
    }

    private static boolean isDConnectedToVisit2(Node a, Node b, Node y, LinkedList<Node> path, List<Node> z, Graph graph) {
        if (b == y) {
            path.addLast(b);
            return true;
        }

        if (path.contains(b)) {
            return false;
        }

        path.addLast(b);

        for (Node c : graph.getAdjacentNodes(b)) {
            if (a == c) {
                continue;
            }

            if (GraphUtils.reachable(a, b, c, z, graph)) {
                if (GraphUtils.isDConnectedToVisit2(b, c, y, path, z, graph)) {
//                    path.removeLast();
                    return true;
                }
            }
        }

        path.removeLast();
        return false;
    }

    public static boolean isDConnectedTo3(Node x, Node y, List<Node> z, Graph graph) {
        return GraphUtils.reachableDConnectedNodes(x, z, graph).contains(y);
    }

    private static Set<Node> reachableDConnectedNodes(Node x, List<Node> z, Graph graph) {
        Set<Node> R = new HashSet<>();
        R.add(x);

        Queue<OrderedPair<Node>> Q = new ArrayDeque<>();
        Set<OrderedPair<Node>> V = new HashSet<>();

        for (Node node : graph.getAdjacentNodes(x)) {
            OrderedPair<Node> edge = new OrderedPair<>(x, node);
            Q.offer(edge);
            V.add(edge);
            R.add(node);
        }

        while (!Q.isEmpty()) {
            OrderedPair<Node> t = Q.poll();

            Node a = t.getFirst();
            Node b = t.getSecond();

            for (Node c : graph.getAdjacentNodes(b)) {
                if (c == a) {
                    continue;
                }
                if (!GraphUtils.reachable(a, b, c, z, graph, null)) {
                    continue;
                }
                R.add(c);

                OrderedPair<Node> u = new OrderedPair<>(b, c);
                if (V.contains(u)) {
                    continue;
                }

                V.add(u);
                Q.offer(u);
            }
        }

        return R;
    }

    // Finds a sepset for x and y, if there is one; otherwise, returns null.
    public static List<Node> getSepset(Node x, Node y, Graph graph) {
        List<Node> sepset = GraphUtils.getSepsetVisit(x, y, graph);
        if (sepset == null) {
            sepset = GraphUtils.getSepsetVisit(y, x, graph);
        }
        return sepset;
    }

    private static List<Node> getSepsetVisit(Node x, Node y, Graph graph) {
        if (x == y) {
            return null;
        }

        List<Node> z = new ArrayList<>();

        List<Node> _z;

        do {
            _z = new ArrayList<>(z);

            Set<Node> path = new HashSet<>();
            path.add(x);
            Set<Triple> colliders = new HashSet<>();

            for (Node b : graph.getAdjacentNodes(x)) {
                if (GraphUtils.sepsetPathFound(x, b, y, path, z, graph, colliders, -1)) {
                    return null;
                }
            }
        } while (!new HashSet<>(z).equals(new HashSet<>(_z)));

        return z;
    }

    private static boolean sepsetPathFound(Node a, Node b, Node y, Set<Node> path, List<Node> z, Graph graph, Set<Triple> colliders, int bound) {
        if (b == y) {
            return true;
        }

        if (path.contains(b)) {
            return false;
        }

        if (path.size() > (bound == -1 ? 1000 : bound)) {
            return false;
        }

        path.add(b);

        if (b.getNodeType() == NodeType.LATENT || z.contains(b)) {
            List<Node> passNodes = GraphUtils.getPassNodes(a, b, z, graph);

            for (Node c : passNodes) {
                if (GraphUtils.sepsetPathFound(b, c, y, path, z, graph, colliders, bound)) {
                    path.remove(b);
                    return true;
                }
            }

            path.remove(b);
            return false;
        } else {
            boolean found1 = false;
            Set<Triple> _colliders1 = new HashSet<>();

            for (Node c : GraphUtils.getPassNodes(a, b, z, graph)) {
                if (GraphUtils.sepsetPathFound(b, c, y, path, z, graph, _colliders1, bound)) {
                    found1 = true;
                    break;
                }
            }

            if (!found1) {
                path.remove(b);
                colliders.addAll(_colliders1);
                return false;
            }

            z.add(b);
            boolean found2 = false;
            Set<Triple> _colliders2 = new HashSet<>();

            for (Node c : GraphUtils.getPassNodes(a, b, z, graph)) {
                if (GraphUtils.sepsetPathFound(b, c, y, path, z, graph, _colliders2, bound)) {
                    found2 = true;
                    break;
                }
            }

            if (!found2) {
                path.remove(b);
                colliders.addAll(_colliders2);
                return false;
            }

            z.remove(b);
            path.remove(b);
            return true;
        }
    }

    private static Set<Triple> colliders(Node b, Graph graph, Set<Triple> colliders) {
        Set<Triple> _colliders = new HashSet<>();

        for (Triple collider : colliders) {
            if (graph.isAncestorOf(collider.getY(), b)) {
                _colliders.add(collider);
            }
        }

        return _colliders;
    }

    private static boolean reachable(Node a, Node b, Node c, List<Node> z, Graph graph) {
        boolean collider = graph.isDefCollider(a, b, c);

        if ((!collider || graph.isUnderlineTriple(a, b, c)) && !z.contains(b)) {
            return true;
        }

        boolean ancestor = GraphUtils.isAncestor(b, z, graph);
        return collider && ancestor;
    }

    private static boolean reachable(Edge e1, Edge e2, Node a, List<Node> z, Graph graph) {
        Node b = e1.getDistalNode(a);
        Node c = e2.getDistalNode(b);

        boolean collider = e1.getProximalEndpoint(b) == Endpoint.ARROW && e2.getProximalEndpoint(b) == Endpoint.ARROW;

        if ((!collider || graph.isUnderlineTriple(a, b, c)) && !z.contains(b)) {
            return true;
        }

        boolean ancestor = GraphUtils.isAncestor(b, z, graph);
        return collider && ancestor;
    }

    private static boolean reachable(Node a, Node b, Node c, List<Node> z, Graph graph, Set<Triple> colliders) {
        boolean collider = graph.isDefCollider(a, b, c);

        if (!collider && !z.contains(b)) {
            return true;
        }

        boolean ancestor = GraphUtils.isAncestor(b, z, graph);

        boolean colliderReachable = collider && ancestor;

        if (colliders != null && collider && !ancestor) {
            colliders.add(new Triple(a, b, c));
        }

        return colliderReachable;
    }

    private static boolean isAncestor(Node b, List<Node> z, Graph graph) {
        if (z.contains(b)) {
            return true;
        }

        Queue<Node> Q = new ArrayDeque<>();
        Set<Node> V = new HashSet<>();

        for (Node node : z) {
            Q.offer(node);
            V.add(node);
        }

        while (!Q.isEmpty()) {
            Node t = Q.poll();
            if (t == b) {
                return true;
            }

            for (Node c : graph.getParents(t)) {
                if (!V.contains(c)) {
                    Q.offer(c);
                    V.add(c);
                }
            }
        }

        return false;

    }

    private static List<Node> getPassNodes(Node a, Node b, List<Node> z, Graph graph) {
        List<Node> passNodes = new ArrayList<>();

        for (Node c : graph.getAdjacentNodes(b)) {
            if (c == a) {
                continue;
            }

            if (GraphUtils.reachable(a, b, c, z, graph)) {
                passNodes.add(c);
            }
        }

        return passNodes;
    }

    private static Set<Node> zAncestors(List<Node> z, Graph graph) {
        Queue<Node> Q = new ArrayDeque<>();
        Set<Node> V = new HashSet<>();

        for (Node node : z) {
            Q.offer(node);
            V.add(node);
        }

        while (!Q.isEmpty()) {
            Node t = Q.poll();

            for (Node c : graph.getParents(t)) {
                if (!V.contains(c)) {
                    Q.offer(c);
                    V.add(c);
                }
            }
        }

        return V;
    }

    public static Set<Node> zAncestors2(List<Node> z, Graph graph) {
        Set<Node> ancestors = new HashSet<>(z);

        boolean changed = true;

        while (changed) {
            changed = false;

            for (Node n : new HashSet<>(ancestors)) {
                List<Node> parents = graph.getParents(n);

                if (!ancestors.containsAll(parents)) {
                    ancestors.addAll(parents);
                    changed = true;
                }
            }
        }

        return ancestors;
    }

    public static boolean existsInducingPath(Node x, Node y, Graph graph) {
        if (x.getNodeType() != NodeType.MEASURED) {
            throw new IllegalArgumentException();
        }
        if (y.getNodeType() != NodeType.MEASURED) {
            throw new IllegalArgumentException();
        }

        LinkedList<Node> path = new LinkedList<>();
        path.add(x);

        for (Node b : graph.getAdjacentNodes(x)) {
            if (GraphUtils.existsInducingPathVisit(graph, x, b, x, y, path)) {
                return true;
            }
        }

        return false;
    }

    // Needs to be public.
    public static boolean existsInducingPathVisit(Graph graph, Node a, Node b, Node x, Node y, LinkedList<Node> path) {
        if (path.contains(b)) {
            return false;
        }

        path.addLast(b);

        if (b == y) {
            return true;
        }

        for (Node c : graph.getAdjacentNodes(b)) {
            if (c == a) {
                continue;
            }

            if (b.getNodeType() == NodeType.MEASURED) {
                if (!graph.isDefCollider(a, b, c)) {
                    continue;
                }

            }

            if (graph.isDefCollider(a, b, c)) {
                if (!(graph.isAncestorOf(b, x) || graph.isAncestorOf(b, y))) {
                    continue;
                }
            }

            if (GraphUtils.existsInducingPathVisit(graph, b, c, x, y, path)) {
                return true;
            }
        }

        path.removeLast();
        return false;
    }

    public static Set<Node> getInducedNodes(Node x, Graph graph) {
        if (x.getNodeType() != NodeType.MEASURED) {
            throw new IllegalArgumentException();
        }

        LinkedList<Node> path = new LinkedList<>();
        path.add(x);

        Set<Node> induced = new HashSet<>();

        for (Node b : graph.getAdjacentNodes(x)) {
            GraphUtils.collectInducedNodesVisit(graph, x, b, path, induced);
        }

        return induced;
    }

    private static void collectInducedNodesVisit(Graph graph, Node x, Node b, LinkedList<Node> path, Set<Node> induced) {
        if (path.contains(b)) {
            return;
        }

        if (induced.contains(b)) {
            return;
        }

        path.addLast(b);

        if (GraphUtils.isInducingPath(graph, path)) {
            induced.add(b);
        }

        for (Node c : graph.getAdjacentNodes(b)) {
            GraphUtils.collectInducedNodesVisit(graph, x, c, path, induced);
        }

        path.removeLast();
    }

    public static boolean isInducingPath(Graph graph, LinkedList<Node> path) {
        if (path.size() < 2) {
            return false;
        }
        if (path.get(0).getNodeType() != NodeType.MEASURED) {
            return false;
        }
        if (path.get(path.size() - 1).getNodeType() != NodeType.MEASURED) {
            return false;
        }

        System.out.println("Path = " + path);

        Node x = path.get(0);
        Node y = path.get(path.size() - 1);

        for (int i = 0; i < path.size() - 2; i++) {
            Node a = path.get(i);
            Node b = path.get(i + 1);
            Node c = path.get(i + 2);

            if (b.getNodeType() == NodeType.MEASURED) {
                if (!graph.isDefCollider(a, b, c)) {
                    return false;
                }
            }

            if (graph.isDefCollider(a, b, c)) {
                if (!(graph.isAncestorOf(b, x) || graph.isAncestorOf(b, y))) {
                    return false;
                }
            }
        }

        return true;
    }

    public static List<Node> getInducingPath(Node x, Node y, Graph graph) {
        if (x.getNodeType() != NodeType.MEASURED) {
            throw new IllegalArgumentException();
        }
        if (y.getNodeType() != NodeType.MEASURED) {
            throw new IllegalArgumentException();
        }

        LinkedList<Node> path = new LinkedList<>();
        path.add(x);

        for (Node b : graph.getAdjacentNodes(x)) {
            if (GraphUtils.existsInducingPathVisit(graph, x, b, x, y, path)) {
                return path;
            }
        }

        return null;
    }

    public static List<Node> possibleDsep(Node x, Node y, Graph graph, int maxPathLength, IndependenceTest test) {
        Set<Node> dsep = new HashSet<>();

        Queue<OrderedPair<Node>> Q = new ArrayDeque<>();
        Set<OrderedPair<Node>> V = new HashSet<>();

        Map<Node, Set<Node>> previous = new HashMap<>();
        previous.put(x, null);

        OrderedPair<Node> e = null;
        int distance = 0;

        assert graph != null;
        Set<Node> adjacentNodes = new HashSet<>(graph.getAdjacentNodes(x));

        for (Node b : adjacentNodes) {
            if (b == y) {
                continue;
            }
            OrderedPair<Node> edge = new OrderedPair<>(x, b);
            if (e == null) {
                e = edge;
            }
            Q.offer(edge);
            V.add(edge);
            GraphUtils.addToSet(previous, b, x);
            dsep.add(b);
        }

        while (!Q.isEmpty()) {
            OrderedPair<Node> t = Q.poll();

            if (e == t) {
                e = null;
                distance++;
                if (distance > 0 && distance > (maxPathLength == -1 ? 1000 : maxPathLength)) {
                    break;
                }
            }

            Node a = t.getFirst();
            Node b = t.getSecond();

            if (GraphUtils.existOnePathWithPossibleParents(previous, b, x, b, graph)) {
                dsep.add(b);
            }

            for (Node c : graph.getAdjacentNodes(b)) {
                if (c == a) {
                    continue;
                }
                if (c == x) {
                    continue;
                }
                if (c == y) {
                    continue;
                }

                GraphUtils.addToSet(previous, b, c);

                if (graph.isDefCollider(a, b, c) || graph.isAdjacentTo(a, c)) {
                    OrderedPair<Node> u = new OrderedPair<>(a, c);
                    if (V.contains(u)) {
                        continue;
                    }

                    V.add(u);
                    Q.offer(u);

                    if (e == null) {
                        e = u;
                    }
                }
            }
        }

        Map<Node, Double> scores = new HashMap<>();

        for (Node node : dsep) {
            test.checkIndependence(x, y, Collections.singletonList(node));
            scores.put(node, test.getScore());
        }

        dsep.remove(x);
        dsep.remove(y);

        List<Node> _dsep = new ArrayList<>(dsep);

        Collections.sort(_dsep);
        Collections.reverse(_dsep);

        return _dsep;
    }

    private static boolean existOnePathWithPossibleParents(Map<Node, Set<Node>> previous, Node w, Node x, Node b, Graph graph) {
        if (w == x) {
            return true;
        }

        Set<Node> p = previous.get(w);
        if (p == null) {
            return false;
        }

        for (Node r : p) {
            if (r == b || r == x) {
                continue;
            }

            if ((GraphUtils.existsSemidirectedPath(r, x, graph)) || GraphUtils.existsSemidirectedPath(r, b, graph)) {
                return true;
            }
        }

        return false;
    }

    private static void addToSet(Map<Node, Set<Node>> previous, Node b, Node c) {
        previous.computeIfAbsent(c, k -> new HashSet<>());
        Set<Node> list = previous.get(c);
        list.add(b);
    }

    public static boolean existsSemidirectedPath(Node from, Node to, Graph G) {
        Queue<Node> Q = new LinkedList<>();
        Set<Node> V = new HashSet<>();

        for (Node u : G.getAdjacentNodes(from)) {
            Edge edge = G.getEdge(from, u);
            Node c = GraphUtils.traverseSemiDirected(from, edge);

            if (c == null) {
                continue;
            }

            if (!V.contains(c)) {
                V.add(c);
                Q.offer(c);
            }
        }

        while (!Q.isEmpty()) {
            Node t = Q.remove();

            if (t == to) {
                return true;
            }

            for (Node u : G.getAdjacentNodes(t)) {
                Edge edge = G.getEdge(t, u);
                Node c = GraphUtils.traverseSemiDirected(t, edge);

                if (c == null) {
                    continue;
                }

                if (!V.contains(c)) {
                    V.add(c);
                    Q.offer(c);
                }
            }
        }

        return false;
    }

    public static int getDegree(Graph graph) {
        int max = 0;

        for (Node node : graph.getNodes()) {
            if (graph.getAdjacentNodes(node).size() > max) {
                max = graph.getAdjacentNodes(node).size();
            }
        }

        return max;
    }

    public static int getIndegree(Graph graph) {
        int max = 0;

        for (Node node : graph.getNodes()) {
            if (graph.getAdjacentNodes(node).size() > max) {
                max = graph.getIndegree(node);
            }
        }

        return max;
    }

    // Returns true if a path consisting of undirected and directed edges toward 'to' exists of
    // length at most 'bound'. Cycle checker in other words.
    public static List<Node> existsUnblockedSemiDirectedPath(Node from, Node to, Set<Node> cond, int bound, Graph graph) {
        Queue<Node> Q = new LinkedList<>();
        Set<Node> V = new HashSet<>();
        Q.offer(from);
        V.add(from);
        Node e = null;
        int distance = 0;
        Map<Node, Node> back = new HashMap<>();

        while (!Q.isEmpty()) {
            Node t = Q.remove();
            if (t == to) {
                LinkedList<Node> _back = new LinkedList<>();
                _back.add(to);
                return _back;
            }

            if (e == t) {
                e = null;
                distance++;
                if (distance > (bound == -1 ? 1000 : bound)) {
                    return null;
                }
            }

            for (Node u : graph.getAdjacentNodes(t)) {
                Edge edge = graph.getEdge(t, u);
                Node c = GraphUtils.traverseSemiDirected(t, edge);
                if (c == null) {
                    continue;
                }
                if (cond.contains(c)) {
                    continue;
                }

                if (c == to) {
                    back.put(c, t);
                    LinkedList<Node> _back = new LinkedList<>();
                    _back.addLast(to);
                    Node f = to;

                    for (int i = 0; i < 10; i++) {
                        f = back.get(f);
                        if (f == null) {
                            break;
                        }
                        _back.addFirst(f);
                    }

                    return _back;
                }

                if (!V.contains(c)) {
                    back.put(c, t);
                    V.add(c);
                    Q.offer(c);

                    if (e == null) {
                        e = u;
                    }
                }
            }
        }

        return null;
    }

    // Used to find semidirected paths for cycle checking.
    public static Node traverseSemiDirected(Node node, Edge edge) {
        if (node == edge.getNode1()) {
            if (edge.getEndpoint1() == Endpoint.TAIL || edge.getEndpoint1() == Endpoint.CIRCLE) {
                return edge.getNode2();
            }
        } else if (node == edge.getNode2()) {
            if (edge.getEndpoint2() == Endpoint.TAIL || edge.getEndpoint2() == Endpoint.CIRCLE) {
                return edge.getNode1();
            }
        }
        return null;
    }

    public static Graph getComparisonGraph(Graph graph, Parameters params) {
        String type = params.getString("graphComparisonType");

        if ("DAG".equals(type)) {
            params.set("graphComparisonType", "DAG");
            return new EdgeListGraph(graph);
        } else if ("CPDAG".equals(type)) {
            params.set("graphComparisonType", "CPDAG");
            return SearchGraphUtils.cpdagForDag(graph);
        } else if ("PAG".equals(type)) {
            params.set("graphComparisonType", "PAG");
            return new DagToPag(graph).convert();
        } else {
            params.set("graphComparisonType", "DAG");
            return new EdgeListGraph(graph);
        }
    }

    /**
     * Check to see if a set of variables Z satisfies the back-door criterion
     * relative to node x and node y.
     *
     * @author Kevin V. Bui (March 2020)
     */
    public boolean isSatisfyBackDoorCriterion(Graph graph, Node x, Node y, List<Node> z) {
        Dag dag = new Dag(graph);

        // make sure no nodes in z is a descendant of x
        if (z.stream().anyMatch(zNode -> dag.isDescendentOf(zNode, x))) {
            return false;
        }

        // make sure zNodes bock every path between node x and node y that contains an arrow into node x
        List<List<Node>> directedPaths = GraphUtils.allDirectedPathsFromTo(graph, x, y, -1);
        directedPaths.forEach(nodes -> {
            // remove all variables that are not on the back-door path
            nodes.forEach(node -> {
                if (!(node == x || node == y)) {
                    dag.removeNode(node);
                }
            });
        });

        return dag.isDSeparatedFrom(x, y, z);
    }

    private static class EdgeNode {

        private final Edge edge;
        private final Node node;

        public EdgeNode(Edge edge, Node node) {
            if (edge.getNode1() == node) {
                this.edge = edge;
            } else if (edge.getNode2() == node) {
                this.edge = edge.reverse();
            } else {
                throw new IllegalArgumentException("Edge does not contain node.");
            }

            this.node = node;
        }

        public int hashCode() {
            return this.edge.hashCode() + 7 * this.node.hashCode();
        }

        public boolean equals(Object o) {
            if (o == null) {
                return false;
            }

            if (!(o instanceof EdgeNode)) {
                return false;
            }

            EdgeNode _o = (EdgeNode) o;
            return _o.edge.equals(this.edge) && _o.node.equals(this.node);
        }

        public Edge getEdge() {
            return this.edge;
        }

        public Node getNode() {
            return this.node;
        }
    }

    private static class Counts {

        private final int[][] counts;

        public Counts() {
            this.counts = new int[8][6];
        }

        public void increment(int m, int n) {
            this.counts[m][n]++;
        }

        public int getCount(int m, int n) {
            return this.counts[m][n];
        }

        public void addAll(Counts counts2) {
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 6; j++) {
                    this.counts[i][j] += counts2.getCount(i, j);
                }
            }
        }

        public int[][] countArray() {
            return this.counts;
        }
    }

    public static class GraphComparison {

        private final int[][] counts;
        private final int adjFn;
        private final int adjFp;
        private final int adjCorrect;
        private final int arrowptFn;
        private final int arrowptFp;
        private final int arrowptCorrect;

        private final double adjPrec;
        private final double adjRec;
        private final double arrowptPrec;
        private final double arrowptRec;

        private final int shd;
        private final int twoCycleFn;
        private final int twoCycleFp;
        private final int twoCycleCorrect;

        private final List<Edge> edgesAdded;
        private final List<Edge> edgesRemoved;
        private final List<Edge> edgesReorientedFrom;
        private final List<Edge> edgesReorientedTo;
        private final List<Edge> edgesAdjacencies;

        public GraphComparison(int adjFn, int adjFp, int adjCorrect, int arrowptFn, int arrowptFp, int arrowptCorrect, double adjPrec, double adjRec, double arrowptPrec, double arrowptRec, int shd, int twoCycleCorrect, int twoCycleFn, int twoCycleFp, List<Edge> edgesAdded, List<Edge> edgesRemoved, List<Edge> edgesReorientedFrom, List<Edge> edgesReorientedTo, List<Edge> edgesAdjacencies, int[][] counts) {
            this.adjFn = adjFn;
            this.adjFp = adjFp;
            this.adjCorrect = adjCorrect;
            this.arrowptFn = arrowptFn;
            this.arrowptFp = arrowptFp;
            this.arrowptCorrect = arrowptCorrect;

            this.adjPrec = adjPrec;
            this.adjRec = adjRec;
            this.arrowptPrec = arrowptPrec;
            this.arrowptRec = arrowptRec;

            this.shd = shd;
            this.twoCycleCorrect = twoCycleCorrect;
            this.twoCycleFn = twoCycleFn;
            this.twoCycleFp = twoCycleFp;
            this.edgesAdded = edgesAdded;
            this.edgesRemoved = edgesRemoved;
            this.edgesReorientedFrom = edgesReorientedFrom;
            this.edgesReorientedTo = edgesReorientedTo;
            this.edgesAdjacencies = edgesAdjacencies;

            this.counts = counts;
        }

        public int getAdjFn() {
            return this.adjFn;
        }

        public int getAdjFp() {
            return this.adjFp;
        }

        public int getAdjCor() {
            return this.adjCorrect;
        }

        public int getAhdFn() {
            return this.arrowptFn;
        }

        public int getAhdFp() {
            return this.arrowptFp;
        }

        public int getAhdCor() {
            return this.arrowptCorrect;
        }

        public int getShd() {
            return this.shd;
        }

        public int getTwoCycleFn() {
            return this.twoCycleFn;
        }

        public int getTwoCycleFp() {
            return this.twoCycleFp;
        }

        public int getTwoCycleCorrect() {
            return this.twoCycleCorrect;
        }

        public List<Edge> getEdgesAdded() {
            return this.edgesAdded;
        }

        public List<Edge> getEdgesRemoved() {
            return this.edgesRemoved;
        }

        public List<Edge> getEdgesReorientedFrom() {
            return this.edgesReorientedFrom;
        }

        public List<Edge> getEdgesReorientedTo() {
            return this.edgesReorientedTo;
        }

        public List<Edge> getCorrectAdjacencies() {
            return this.edgesAdjacencies;
        }

        public double getAdjPrec() {
            return this.adjPrec;
        }

        public double getAdjRec() {
            return this.adjRec;
        }

        public double getAhdPrec() {
            return this.arrowptPrec;
        }

        public double getAhdRec() {
            return this.arrowptRec;
        }

        public int[][] getCounts() {
            return this.counts;
        }
    }

    public static class TwoCycleErrors {

        public int twoCycCor;
        public int twoCycFn;
        public int twoCycFp;

        public TwoCycleErrors(int twoCycCor, int twoCycFn, int twoCycFp) {
            this.twoCycCor = twoCycCor;
            this.twoCycFn = twoCycFn;
            this.twoCycFp = twoCycFp;
        }

        public String toString() {
            return "2c cor = " + this.twoCycCor + "\t" + "2c fn = " + this.twoCycFn + "\t" + "2c fp = " + this.twoCycFp;
        }

    }

}
