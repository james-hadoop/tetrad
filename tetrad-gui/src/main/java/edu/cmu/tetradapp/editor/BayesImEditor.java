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
package edu.cmu.tetradapp.editor;

import edu.cmu.tetrad.bayes.BayesIm;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetradapp.model.BayesImWrapper;
import edu.cmu.tetradapp.workbench.GraphWorkbench;

import javax.swing.*;
import java.awt.*;

/**
 * An editor for Bayes net instantiated models. Assumes that the workbench and
 * parameterized model have been established (that is, that the nodes have been
 * identified and named and that the number and names of the values for the
 * nodes have been specified) and allows the user to set conditional
 * probabilities of node values given combinations of parent values.
 *
 * @author Aaron Powers
 * @author Joseph Ramsey jdramsey@andrew.cmu.edu
 */
public class BayesImEditor extends JPanel {

    private static final long serialVersionUID = 1L;

    private final JPanel targetPanel;
    /**
     * The wizard that allows the user to modify parameter values for this IM.
     */
    private BayesImEditorWizard wizard;
    private final BayesImWrapper wrapper;

    /**
     * Constructs a new instanted model editor from a Bayes IM.
     */
    public BayesImEditor(BayesImWrapper wrapper) {
        this.wrapper = wrapper;
        setLayout(new BorderLayout());

        this.targetPanel = new JPanel();
        this.targetPanel.setLayout(new BorderLayout());

        setEditorPanel();

        add(this.targetPanel, BorderLayout.CENTER);
        validate();

        if (wrapper.getNumModels() > 1) {
            JComboBox<Integer> comp = new JComboBox<>();

            for (int i = 0; i < wrapper.getNumModels(); i++) {
                comp.addItem(i + 1);
            }

            comp.setSelectedIndex(wrapper.getModelIndex());

            comp.addActionListener(e -> {
                wrapper.setModelIndex(comp.getSelectedIndex() - 1);
                setEditorPanel();
                validate();
            });

            comp.setMaximumSize(comp.getPreferredSize());

            Box b = Box.createHorizontalBox();
            b.add(new JLabel("Using model"));
            b.add(comp);
            b.add(new JLabel("from "));
            b.add(new JLabel(wrapper.getModelSourceName()));
            b.add(Box.createHorizontalGlue());

            add(b, BorderLayout.NORTH);
        }
    }

    private void setEditorPanel() {
        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());

        BayesIm bayesIm = this.wrapper.getBayesIm();
        Graph graph = bayesIm.getDag();
        GraphWorkbench workbench = new GraphWorkbench(graph);

        JMenuBar menuBar = new JMenuBar();
        JMenu file = new JMenu("File");
        menuBar.add(file);
        file.add(new SaveBayesImXmlAction(this));
        file.add(new LoadBayesImXmlAction(this.wrapper, this));
        file.add(new LoadBayesImXsdlXmlAction(this.wrapper, this));
        file.add(new SaveScreenshot(this, true, "Save Screenshot..."));
        file.add(new SaveComponentImage(workbench, "Save Graph Image..."));
        setLayout(new BorderLayout());
        panel.add(menuBar, BorderLayout.NORTH);

        this.wizard = new BayesImEditorWizard(bayesIm, workbench);
        this.wizard.enableEditing(false);

        this.wizard.addPropertyChangeListener(evt -> {
            if ("editorValueChanged".equals(evt.getPropertyName())) {
                firePropertyChange("modelChanged", null, null);
            }
        });

        JScrollPane workbenchScroll = new JScrollPane(workbench);
        JScrollPane wizardScroll = new JScrollPane(getWizard());

        workbenchScroll.setPreferredSize(new Dimension(450, 450));

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                workbenchScroll, wizardScroll);
        splitPane.setOneTouchExpandable(true);
        splitPane.setDividerLocation(workbenchScroll.getPreferredSize().width);
        panel.add(splitPane, BorderLayout.CENTER);

        setName("Bayes IM Editor");
        getWizard().addPropertyChangeListener(evt -> {
            if ("editorClosing".equals(evt.getPropertyName())) {
                firePropertyChange("editorClosing", null, getName());
            }

            if ("closeFrame".equals(evt.getPropertyName())) {
                firePropertyChange("closeFrame", null, null);
                firePropertyChange("editorClosing", true, true);
            }

            if ("modelChanged".equals(evt.getPropertyName())) {
                firePropertyChange("modelChanged", evt.getOldValue(),
                        evt.getNewValue());
            }
        });

        this.targetPanel.add(panel, BorderLayout.CENTER);
        revalidate();
        repaint();
    }

    /**
     * Sets the name of this editor.
     */
    public void setName(String name) {
        String oldName = getName();
        super.setName(name);
        firePropertyChange("name", oldName, getName());
    }

    /**
     * @return a reference to this editor.
     */
    public BayesImEditorWizard getWizard() {
        return this.wizard;
    }

    public void getBayesIm() {
        removeAll();
        setEditorPanel();
        revalidate();
        repaint();
        firePropertyChange("modelChanged", null, null);
    }

}
