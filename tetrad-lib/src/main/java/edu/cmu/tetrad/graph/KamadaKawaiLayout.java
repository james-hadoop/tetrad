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

import edu.cmu.tetrad.util.Matrix;
import edu.cmu.tetrad.util.NumberFormatUtil;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradAlgebra;

import javax.swing.*;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.List;

/**
 * Lays out a graph by placing springs between the nodes and letting the system
 * settle (one node at a time).
 *
 * @author Joseph Ramsey
 */
public final class KamadaKawaiLayout {

    /**
     * The graph being laid out.
     */
    private final Graph graph;

    /**
     * The list of nodes used to construct d, p, k, and l.
     */
    private List<Node> componentNodes;

    /**
     * Natural length of an edge.
     */
    private double naturalEdgeLength = 80.0;

    /**
     * Spring constant; higher for more elasticity.
     */
    private double springConstant = 0.5;

    /**
     * Node i is at (p[i][0], p[i][1]).
     */
    private double[][] p;

    /**
     * l[i][j] is the natural length of the spring between node i and node j
     * defined by L * d[i][j].
     */
    private double[][] l;

    /**
     * k[i][j] is the strength of the spring between node i and node j, defined
     * by K / (d[i][j] * d[i][j]).
     */
    private double[][] k;

    /**
     * Leftmost x coord minus 100.0 to lay out the next component.
     */
    private double leftmostX = -50.;

    /**
     * Monitors progress for the user.
     */
    private ProgressMonitor monitor;

    /**
     * True if nodes should be initialized in random locations, false if they
     * should be initialized in their getModel locations.
     */
    private boolean randomlyInitialized;

    /**
     * The max delta at which the algorithm will stop settling.
     */
    private double stopEnergy = 1.0;

    //==============================CONSTRUCTORS===========================//

    public KamadaKawaiLayout(Graph graph) {
        if (graph == null) {
            throw new NullPointerException();
        }

        this.graph = GraphUtils.undirectedGraph(graph);
    }

    //============================PUBLIC METHODS==========================//

    public void doLayout() {
        GraphUtils.circleLayout(this.graph, 300, 300, 200);

        this.monitor = new ProgressMonitor(null, "Energy settling...",
                "Energy = ?", 0, 100);
        getMonitor().setMillisToDecideToPopup(10);
        getMonitor().setMillisToPopup(0);
        getMonitor().setProgress(0);

        List<List<Node>> components =
                GraphUtils.connectedComponents(this.graph);

        components.sort((o1, o2) -> {
            int i1 = o1.size();
            int i2 = o2.size();
            return Integer.compare(i2, i1);
        });

        for (List<Node> component1 : components) {
            initialize(component1, isRandomlyInitialized());
            layoutComponent(component1);
        }

        getMonitor().setProgress(100);
    }


    private boolean isRandomlyInitialized() {
        return this.randomlyInitialized;
    }

    public void setRandomlyInitialized(boolean randomlyInitialized) {
        this.randomlyInitialized = randomlyInitialized;
    }

    private double getStopEnergy() {
        return this.stopEnergy;
    }

    public void setStopEnergy(double stopEnergy) {
        if (stopEnergy <= 0.0) {
            throw new IllegalArgumentException(
                    "Stop energy must be greater than" + " zero.");
        }

        this.stopEnergy = stopEnergy;
    }


    private double getNaturalEdgeLength() {
        return this.naturalEdgeLength;
    }

    public void setNaturalEdgeLength(double naturalEdgeLength) {
        if (naturalEdgeLength < 0.0) {
            throw new IllegalArgumentException(
                    "Natural edge length should be " + "greater than zero.");
        }

        this.naturalEdgeLength = naturalEdgeLength;
    }

    private double getSpringConstant() {
        return this.springConstant;
    }

    public void setSpringConstant(double springConstant) {
        if (springConstant < 0.0) {
            throw new IllegalArgumentException(
                    "Spring constant should be " + "greater than zero.");
        }

        this.springConstant = springConstant;
    }

    //============================PRIVATE METHODS=========================//

    /**
     * Initializes the layout for the given nodes.
     *
     * @param nodes               the nodes to be laid out.
     * @param randomlyInitialized true if the spring layout should start from a
     *                            randomlyInitialized position, false if the
     *                            spring layout should start from the getModel
     *                            positions of the nodes.
     */
    private void initialize(List<Node> nodes, boolean randomlyInitialized) {
        setComponentNodes(Collections.unmodifiableList(nodes));

        this.p = new double[nodes.size()][2];
        int[][] d;
        this.l = new double[nodes.size()][nodes.size()];
        this.k = new double[nodes.size()][nodes.size()];

        if (randomlyInitialized) {
            for (int i = 0; i < nodes.size(); i++) {
                this.p[i][0] = RandomUtil.getInstance().nextInt(600);
                this.p[i][1] = RandomUtil.getInstance().nextInt(600);
            }
        } else {
            for (int i = 0; i < nodes.size(); i++) {
                Node node = nodes.get(i);
                this.p[i][0] = node.getCenterX();
                this.p[i][1] = node.getCenterY();
            }

        }

        d = allPairsShortestPath();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    continue;
                }
                this.l[i][j] = getNaturalEdgeLength() * d[i][j];
            }
        }

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    continue;
                }
                this.k[i][j] = getSpringConstant() / (d[i][j] * d[i][j]);
            }
        }
    }

    private void layoutComponent(List<Node> componentNodes) {
        setComponentNodes(componentNodes);
        optimize(getStopEnergy());
        shiftComponentToRight(componentNodes);
    }

    private void shiftComponentToRight(List<Node> componentNodes) {
        double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE;

        for (int i = 0; i < componentNodes.size(); i++) {
            if (this.p[i][0] < minX) {
                minX = this.p[i][0];
            }
            if (this.p[i][1] < minY) {
                minY = this.p[i][1];
            }
        }

        this.leftmostX += 100.;

        for (int i = 0; i < componentNodes.size(); i++) {
            this.p[i][0] += this.leftmostX - minX;
            this.p[i][1] += 40.0 - minY;
        }

        for (int i = 0; i < componentNodes.size(); i++) {
            if (this.p[i][0] > this.leftmostX) {
                this.leftmostX = this.p[i][0];
            }
        }

        for (int i = 0; i < componentNodes.size(); i++) {
            Node node = componentNodes.get(i);
            node.setCenterX((int) this.p[i][0]);
            node.setCenterY((int) this.p[i][1]);
        }
    }

    private void optimize(double deltaCutoff) {
        NumberFormat nf = NumberFormatUtil.getInstance().getNumberFormat();

        double initialMaxDelta = -1.;
        double maxDelta;
        final int jump = 100;
        Matrix a = new Matrix(2, 2);
        Matrix b = new Matrix(2, 1);
        int oldM = -1;

        do {
            if (this.monitor.isCanceled()) {
                return;
            }

            int[] m = new int[1];
            maxDelta = maxDelta(m);

            if (initialMaxDelta == -1) {
                initialMaxDelta = maxDelta;
            }

            if (m[0] == oldM) {
                this.p[m[0]][0] += RandomUtil.getInstance().nextInt(
                        2 * jump) - jump;
                this.p[m[0]][1] += RandomUtil.getInstance().nextInt(
                        2 * jump) - jump;
                continue;
            }

            oldM = m[0];

            int progress =
                    (int) (99.0 - 98.0 * maxDelta / (0.5 * initialMaxDelta));
            if (progress < 1) {
                progress = 1;
            }
            if (progress > 99) {
                progress = 99;
            }
            getMonitor().setProgress(progress);
            getMonitor().setNote("Energy = " + nf.format(maxDelta));

            if (m[0] == -1) {
                throw new IllegalStateException();
            }

            double oldDelta = Double.NaN;
            double delta;

            while ((delta = delta(m[0])) > deltaCutoff) {
                Thread.yield();
                if (this.monitor.isCanceled()) {
                    return;
                }

                if (Math.abs(delta - oldDelta) < 0.001) {
                    this.p[m[0]][0] += RandomUtil.getInstance().nextInt(
                            2 * jump) - jump;
                    this.p[m[0]][1] += RandomUtil.getInstance().nextInt(
                            2 * jump) - jump;
                    continue;
                }

                final double h = 1.e-2;

                double partialXX = secondPartial(m[0], 0, 0);
                double partialXY = secondPartial(m[0], 0, 1);
                double partialX = firstPartial(m[0], 0, h);
                double partialYY = secondPartial(m[0], 1, 1);
                double partialY = firstPartial(m[0], 1, h);

                a.set(0, 0, partialXX);
                a.set(0, 1, partialXY);
                a.set(1, 0, partialXY);
                a.set(1, 1, partialYY);

                b.set(0, 0, -partialX);
                b.set(1, 0, -partialY);

                Matrix c;

                try {
                    c = TetradAlgebra.solve(a, b);
                } catch (Exception e) {
                    this.p[m[0]][0] += RandomUtil.getInstance().nextInt(
                            2 * jump) - jump;
                    this.p[m[0]][1] += RandomUtil.getInstance().nextInt(
                            2 * jump) - jump;
                    continue;
                }

                double dx = c.get(0, 0);
                double dy = c.get(1, 0);

                this.p[m[0]][0] += dx;
                this.p[m[0]][1] += dy;

                oldDelta = delta;
            }
        } while (maxDelta > deltaCutoff);
    }

    private double energy() {
        int n = this.p.length;
        double sum = 0.0;

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                sum += 0.5 * this.k[i][j] * Math.pow(distance(i, j) - this.l[i][j], 2.0);
            }
        }

        return sum;
    }

    private double maxDelta(int[] index) {
        double maxDelta = Double.NEGATIVE_INFINITY;
        int m = -1;

        for (int i = 0; i < getComponentNodes().size(); i++) {
            double delta = delta(i);

            if (delta == Double.NEGATIVE_INFINITY) {
                throw new IllegalStateException();
            }

            if (delta > maxDelta) {
                maxDelta = delta;
                m = i;
            }
        }

        //        System.out.println("maxDelta = " + maxDelta);
        index[0] = m;
        return maxDelta;
    }

    private double delta(int i) {
        double partialX = firstPartial(i, 0, 1.e-4);
        double partialY = firstPartial(i, 1, 1.e-4);
        return Math.sqrt(partialX * partialX + partialY * partialY);
    }

    private double firstPartial(int i, int var, double h) {
        double storedCoord = this.p[i][var];

        this.p[i][var] -= h;
        double energy1 = energy();

        this.p[i][var] += h;
        double energy2 = energy();

        this.p[i][var] = storedCoord;
        return (energy2 - energy1) / (2. * h);
    }

    private double secondPartial(int m, int i, int j) {
        double storedX = this.p[m][0];
        double storedY = this.p[m][1];

        this.p[m][i] += 0.01;
        this.p[m][j] += 0.01;
        double ff1 = energy();

        this.p[m][j] -= 2 * 0.01;
        double ff2 = energy();

        this.p[m][i] -= 2 * 0.01;
        this.p[m][j] += 2 * 0.01;
        double ff3 = energy();

        this.p[m][j] -= 2 * 0.01;
        double ff4 = energy();

        this.p[m][0] = storedX;
        this.p[m][1] = storedY;

        return (ff1 - ff2 - ff3 + ff4) / (4.0 * 0.01 * 0.01);
    }

    private double distance(int i, int j) {
        double x1 = this.p[i][0];
        double y1 = this.p[i][1];
        double x2 = this.p[j][0];
        double y2 = this.p[j][1];

        return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    }

    /**
     * Floyd's all-pairs shortest-path algorithm, restricted to integral
     * lengths. Returns an int[][] matrix I, where I[i][j] is the length of the
     * shortest path from i to j.
     */
    private int[][] allPairsShortestPath() {
        int[][] I1 = new int[getComponentNodes().size()][getComponentNodes()
                .size()];
        int[][] I2 = new int[getComponentNodes().size()][getComponentNodes()
                .size()];
        int infinity = getComponentNodes().size() * getComponentNodes().size();

        for (int i = 0; i < getComponentNodes().size(); i++) {
            for (int j = 0; j < getComponentNodes().size(); j++) {
                Node node1 = getComponentNodes().get(i);
                Node node2 = getComponentNodes().get(j);
                if (this.graph.getEdge(node1, node2) != null) {
                    I2[i][j] = 1;
                } else {
                    I2[i][j] = infinity;
                }
            }
        }

        for (int k = 0; k < getComponentNodes().size(); k++) {
            int[][] temp = I1;
            I1 = I2;
            I2 = temp;

            for (int i = 0; i < getComponentNodes().size(); i++) {
                for (int j = 0; j < getComponentNodes().size(); j++) {
                    I2[i][j] = Math.min(I1[i][j], I1[i][k] + I1[k][j]);
                }
            }

        }

        return I2;
    }

    private ProgressMonitor getMonitor() {
        return this.monitor;
    }

    private List<Node> getComponentNodes() {
        return this.componentNodes;
    }

    private void setComponentNodes(List<Node> componentNodes) {
        this.componentNodes = componentNodes;
    }
}





