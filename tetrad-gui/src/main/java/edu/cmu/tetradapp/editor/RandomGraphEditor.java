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

package edu.cmu.tetradapp.editor;

import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetradapp.util.DoubleTextField;
import edu.cmu.tetradapp.util.IntTextField;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.text.DecimalFormat;

/**
 * Edits the parameters for generating random graphs.
 *
 * @author Joseph Ramsey
 */
class RandomGraphEditor extends JPanel {
    private final Parameters parameters;
    private final IntTextField numNodesField;
    private final IntTextField numLatentsField;
    private final DoubleTextField avgDegreeField;
    private final IntTextField maxIndegreeField;
    private final IntTextField maxOutdegreeField;
    private final IntTextField maxDegreeField;
    private final JComboBox<String> connectedBox;
    private final IntTextField numTwoCyclesField;
    private final IntTextField minCycleLengthField;

    /**
     * Constructs a dialog to edit the given workbench randomization
     * parameters.
     */
    public RandomGraphEditor(boolean cyclicAllowed, Parameters parameters) {
        this(new EdgeListGraph(), cyclicAllowed, parameters);
    }

    /**
     * Constructs a dialog to edit the given workbench randomization
     * parameters.
     * //     * @param preferredNumNodes an integer which, if greater than 1, will revise the number of nodes,
     * //     * number of edges,a nd number of latent nodes. Useful if the interface suggests a number of nodes
     * //     * that overrides the number of nodes set in the preferences.
     */
    public RandomGraphEditor(Graph oldGraph, boolean cyclicAllowed, Parameters parameters) {
        if (parameters == null) {
            throw new NullPointerException();
        }

        this.parameters = parameters;

        numNodesField = new IntTextField(getNumMeasuredNodes(), 4);
        numLatentsField = new IntTextField(getNumLatents(), 4);
        avgDegreeField = new DoubleTextField(getAvgDegree(), 4, new DecimalFormat("0.0"));
        maxIndegreeField = new IntTextField(getMaxIndegree(), 4);
        maxOutdegreeField = new IntTextField(getMaxOutdegree(), 4);
        maxDegreeField = new IntTextField(getMaxDegree(), 4);
        connectedBox = new JComboBox<>(new String[]{"No", "Yes"});
        JComboBox<String> addCyclesBox = new JComboBox<String>(new String[]{"No", "Yes"});
        numTwoCyclesField = new IntTextField(getMinNumCycles(), 4);
        minCycleLengthField = new IntTextField(getMinCycleLength(), 4);

        // set up text and ties them to the parameters object being edited.
        numNodesField.setFilter((value, oldValue) -> {
            if (value == numNodesField.getValue()) {
                return oldValue;
            }

            try {
                setNumMeasuredNodes(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        numLatentsField.setFilter((value, oldValue) -> {
            if (value == numLatentsField.getValue()) {
                return oldValue;
            }

            try {
                setNumLatents(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        avgDegreeField.setFilter((value, oldValue) -> {
            if (value == avgDegreeField.getValue()) {
                return oldValue;
            }

            try {
                setAvgDegree(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        maxIndegreeField.setFilter((value, oldValue) -> {
            if (value == maxIndegreeField.getValue()) {
                return oldValue;
            }

            try {
                setMaxIndegree(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        maxOutdegreeField.setFilter((value, oldValue) -> {
            if (value == maxOutdegreeField.getValue()) {
                return oldValue;
            }

            try {
                setMaxOutdegree(value);
            } catch (Exception e) {
                // Ignore.
            }

            return value;
        });

        maxDegreeField.setFilter((value, oldValue) -> {
            if (value == maxDegreeField.getValue()) {
                return oldValue;
            }

            try {
                setMaxDegree(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        if (isConnected()) {
            connectedBox.setSelectedItem("Yes");
        } else {
            connectedBox.setSelectedItem("No");
        }

        minCycleLengthField.setEnabled(isAddCycles());

        connectedBox.setMaximumSize(connectedBox.getPreferredSize());
        connectedBox.addActionListener(e -> {
            if (!(e.getSource() instanceof JComboBox)) {
                return;
            }

            JComboBox box = (JComboBox) e.getSource();
            if ("Yes".equals(box.getSelectedItem())) {
                setConnected(true);
            } else if ("No".equals(box.getSelectedItem())) {
                setConnected(false);
            } else {
                throw new IllegalArgumentException();
            }
        });

        if (isAddCycles()) {
            addCyclesBox.setSelectedItem("Yes");
        } else {
            addCyclesBox.setSelectedItem("No");
        }

        addCyclesBox.setMaximumSize(addCyclesBox.getPreferredSize());
        addCyclesBox.addActionListener(e -> {
            JComboBox box = (JComboBox) e.getSource();
            if ("Yes".equals(box.getSelectedItem())) {
                setAddCycles(true);
                minCycleLengthField.setEnabled(true);
            } else if ("No".equals(box.getSelectedItem())) {
                setAddCycles(false);
                minCycleLengthField.setEnabled(false);
            } else {
                throw new IllegalArgumentException();
            }
        });

        numTwoCyclesField.setFilter((value, oldValue) -> {
            if (value == numTwoCyclesField.getValue()) {
                return oldValue;
            }

            try {
                setMinNumCycles(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        minCycleLengthField.setFilter((value, oldValue) -> {
            if (value == minCycleLengthField.getValue()) {
                return oldValue;
            }

            try {
                setMinCycleLength(value);
            } catch (Exception ignored) {
            }

            return value;
        });

        // construct the workbench.
        setLayout(new BorderLayout());

        Box b1 = Box.createVerticalBox();

        Box b2 = Box.createHorizontalBox();
        b2.add(new JLabel("Parameters for Random DAG:"));
        b2.add(Box.createHorizontalGlue());
        b1.add(b2);
        b1.add(Box.createVerticalStrut(5));

        Box b10 = Box.createHorizontalBox();
        b10.add(new JLabel("Number of measured nodes:"));
        b10.add(Box.createRigidArea(new Dimension(10, 0)));
        b10.add(Box.createHorizontalGlue());
        b10.add(numNodesField);
        b1.add(b10);

        Box b11 = Box.createHorizontalBox();
        b11.add(new JLabel("Max # latent confounders:"));
        b11.add(Box.createHorizontalStrut(25));
        b11.add(Box.createHorizontalGlue());
        b11.add(numLatentsField);
        b1.add(b11);
        b1.add(Box.createVerticalStrut(5));

        Box b12 = Box.createHorizontalBox();
        b12.add(new JLabel("Average Degree:"));
        b12.add(Box.createHorizontalGlue());
        b12.add(avgDegreeField);
        b1.add(b12);
        b1.add(Box.createVerticalStrut(5));

        Box b14 = Box.createHorizontalBox();
        b14.add(new JLabel("Maximum indegree:"));
        b14.add(Box.createHorizontalGlue());
        b14.add(maxIndegreeField);
        b1.add(b14);

        Box b15 = Box.createHorizontalBox();
        b15.add(new JLabel("Maximum outdegree:"));
        b15.add(Box.createHorizontalGlue());
        b15.add(maxOutdegreeField);
        b1.add(b15);

        Box b13 = Box.createHorizontalBox();
        b13.add(new JLabel("Maximum degree:"));
        b13.add(Box.createHorizontalGlue());
        b13.add(maxDegreeField);
        b1.add(b13);
        b1.add(Box.createVerticalStrut(5));

        Box b16 = Box.createHorizontalBox();
        b16.add(new JLabel("Connected:"));
        b16.add(Box.createHorizontalGlue());
        b16.add(connectedBox);
        b1.add(b16);
        b1.add(Box.createVerticalStrut(5));

        Box d = Box.createVerticalBox();
        b1.setBorder(new TitledBorder(""));
        d.add(b1);

        if (cyclicAllowed) {
            Box c1 = Box.createVerticalBox();

            Box c2 = Box.createHorizontalBox();
            c2.add(new JLabel("Create a cyclic graph?"));
            c2.add(Box.createHorizontalGlue());
            c2.add(addCyclesBox);
            c1.add(c2);
            c1.add(Box.createVerticalStrut(5));

            Box c3 = Box.createHorizontalBox();
            c3.add(new JLabel("Number of two cycles to add:"));
            c3.add(Box.createHorizontalGlue());
            c3.add(numTwoCyclesField);
            c1.add(c3);
            c1.add(Box.createVerticalStrut(5));

            c1.setBorder(new TitledBorder(""));

            d.add(Box.createVerticalStrut(5));
            d.add(c1);
        }

        add(d, BorderLayout.CENTER);
    }

    public void setEnabled(boolean enabled) {
        super.setEnabled(enabled);

        if (isChooseFixed() && enabled) {
            numNodesField.setEnabled(enabled);
            numLatentsField.setEnabled(enabled);
            avgDegreeField.setEnabled(enabled);
            maxIndegreeField.setEnabled(false);
            maxOutdegreeField.setEnabled(false);
            maxDegreeField.setEnabled(false);
            connectedBox.setEnabled(false);
        } else {
            numNodesField.setEnabled(enabled);
            numLatentsField.setEnabled(enabled);
            avgDegreeField.setEnabled(enabled);
            maxIndegreeField.setEnabled(enabled);
            maxOutdegreeField.setEnabled(enabled);
            maxDegreeField.setEnabled(enabled);
            connectedBox.setEnabled(enabled);
        }
    }

    public boolean isChooseFixed() {
        return parameters.getBoolean("graphChooseFixed", true);
    }

    public int getNumNodes() {
        return getNumMeasuredNodes() + getNumLatents();
    }

    private int getNumMeasuredNodes() {
        return parameters.getInt(Params.NUM_MEASURES, 10);
    }

    private void setNumMeasuredNodes(int numMeasuredNodes) {
        if (numMeasuredNodes + getNumLatents() < 2) {
            throw new IllegalArgumentException("Number of nodes Must be greater than or equal to 2.");
        }

        parameters.set(Params.NUM_MEASURES, numMeasuredNodes);

        if (isConnected()) {
            setAvgDegree(Math.max(getAvgDegree(), numMeasuredNodes + getNumLatents()));
        }
    }

    public int getNumLatents() {
        return parameters.getInt(Params.NUM_LATENTS, 0);
    }

    private void setNumLatents(int numLatentNodes) {
        if (numLatentNodes < 0) {
            throw new IllegalArgumentException(
                    "Max # latent confounders must be" + " >= 0: " +
                            numLatentNodes);
        }

        parameters.set(Params.NUM_LATENTS, numLatentNodes);
    }

    public double getAvgDegree() {
        return parameters.getDouble(Params.AVG_DEGREE);
    }


    private void setAvgDegree(double avgDegree) {
        parameters.set(Params.AVG_DEGREE, avgDegree);
    }

    public int getMaxDegree() {
        return parameters.getInt(Params.MAX_DEGREE, 100);
    }

    private void setMaxDegree(int maxDegree) {
        if (!isConnected() && maxDegree < 1) {
            parameters.set("maxDegree");
            return;
        }

        if (isConnected() && maxDegree < 3) {
            parameters.set(Params.MAX_DEGREE, 3);
            return;
        }

        parameters.set("maxDegree", maxDegree);
    }

    public int getMaxIndegree() {
        return parameters.getInt(Params.MAX_INDEGREE, 100);
    }

    private void setMaxIndegree(int maxIndegree) {
        if (isConnected() && maxIndegree < 2) {
            parameters.set(Params.MAX_INDEGREE);
            return;
        }

        parameters.set(Params.MAX_INDEGREE, maxIndegree);
    }

    public int getMaxOutdegree() {
        return parameters.getInt(Params.MAX_OUTDEGREE, 100);
    }

    private void setMaxOutdegree(int maxOutDegree) {
        if (!isConnected() && maxOutDegree < 1) {
            parameters.set(Params.MAX_OUTDEGREE, 1);
            return;
        }

        if (isConnected() && maxOutDegree < 2) {
            parameters.set(Params.MAX_OUTDEGREE, 2);
            return;
        }

        parameters.set(Params.MAX_OUTDEGREE, maxOutDegree);
    }

    public boolean isConnected() {
        return parameters.getBoolean(Params.CONNECTED, false);
    }

    private void setConnected(boolean connected) {
        parameters.set(Params.CONNECTED, connected);

        if (connected) {
            if (getMaxIndegree() < 2) {
                setMaxIndegree(2);
            }

            if (getMaxOutdegree() < 2) {
                setMaxOutdegree(2);
            }

            if (getMaxDegree() < 3) {
                setMaxDegree(3);
            }

            if (getAvgDegree() < getNumNodes()) {
                setAvgDegree(getNumNodes());
            }
        }
    }

    public boolean isAddCycles() {
        return parameters.getBoolean("randomGraphAddCycles", false);
    }

    private void setAddCycles(boolean addCycles) {
        parameters.set("randomGraphAddCycles", addCycles);
    }

    public int getMinNumCycles() {
        int minNumCycles = parameters.getInt("randomGraphMinNumCycles", 0);
        System.out.println("get min num cycles = " + minNumCycles);

        return minNumCycles;
    }

    private void setMinNumCycles(int minNumCycles) {

        System.out.println("set min num cycles = " + minNumCycles);

        if (minNumCycles < 0) {
            parameters.set("randomGraphMinNumCycles", 0);
            return;
        }

        parameters.set("randomGraphMinNumCycles", minNumCycles);
    }

    public int getMinCycleLength() {
        return parameters.getInt("randomGraphMinCycleLength", 2);
    }

    private void setMinCycleLength(int minCycleLength) {
        if (minCycleLength < 2) {
            parameters.set("randomGraphMinCycleLength", 2);
            return;
        }

        parameters.set("randomGraphMinCycleLength", minCycleLength);
    }
}





