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

import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetradapp.model.GraphWrapper;
import edu.cmu.tetradapp.model.Misclassifications;
import org.jetbrains.annotations.NotNull;

import javax.swing.*;
import java.awt.*;

/**
 * Provides a little display/editor for notes in the session workbench. This may
 * be elaborated in the future to allow marked up text.
 *
 * @author Joseph Ramsey
 */
public class MisclassificationsEditor extends JPanel {

    private static final long serialVersionUID = -5291697901326757833L;

    /**
     * The model for the note.
     */
    private final Misclassifications comparison;
    private final Parameters params;
    private JTextArea area;

    /**
     * Constructs the editor given the model
     */
    public MisclassificationsEditor(Misclassifications comparison) {
        this.comparison = comparison;
        this.params = comparison.getParams();
        setup();
    }

    //============================ Private Methods =========================//

    private void setup() {
        JTabbedPane pane = new JTabbedPane(JTabbedPane.LEFT);

        String compareString = comparison.getComparisonString();

        JPanel panel = new JPanel();

        Font font = new Font("Monospaced", Font.PLAIN, 14);
        area = new JTextArea();
        area.setText(compareString);

        area.setFont(font);

        JScrollPane scroll = new JScrollPane(area);
        scroll.setPreferredSize(new Dimension(400, 400));

        panel.add(Box.createVerticalStrut(10));

        Box box = Box.createHorizontalBox();
        panel.add(box);
        panel.add(Box.createVerticalStrut(10));

        Box box1 = Box.createHorizontalBox();
        box1.add(new JLabel("Graph Comparison: "));
        box1.add(Box.createHorizontalGlue());

        add(box1);
        setLayout(new BorderLayout());

        pane.add("Comparison", scroll);

        GraphEditor graphEditor = new GraphEditor(new GraphWrapper(comparison.getTargetGraph()));
        graphEditor.enableEditing(false);
        pane.add("Target Graph", graphEditor.getWorkbench());

        graphEditor = new GraphEditor(new GraphWrapper(comparison.getReferenceGraph()));
        graphEditor.enableEditing(false);
        pane.add("True Graph", graphEditor.getWorkbench());

        add(pane, BorderLayout.CENTER);
        add(menubar(), BorderLayout.NORTH);
    }

    @NotNull
    private JMenuBar menubar() {
        JMenuBar menubar = new JMenuBar();
        JMenu menu = new JMenu("Compare To...");
        JMenuItem graph = new JCheckBoxMenuItem("DAG");
        JMenuItem cpdag = new JCheckBoxMenuItem("CPDAG");
        JMenuItem pag = new JCheckBoxMenuItem("PAG");

        ButtonGroup group = new ButtonGroup();
        group.add(graph);
        group.add(cpdag);
        group.add(pag);

        menu.add(graph);
        menu.add(cpdag);
        menu.add(pag);

        menubar.add(menu);

        switch (params.getString("graphComparisonType")) {
            case "CPDAG":
                menu.setText("Compare to CPDAG...");
                cpdag.setSelected(true);
                break;
            case "PAG":
                menu.setText("Compare to PAG...");
                pag.setSelected(true);
                break;
            default:
                menu.setText("Compare to DAG...");
                graph.setSelected(true);
                break;
        }

        graph.addActionListener(e -> {
            params.set("graphComparisonType", "DAG");
            menu.setText("Compare to DAG...");

            area.setText(comparison.getComparisonString());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);

            area.repaint();

        });

        cpdag.addActionListener(e -> {
            params.set("graphComparisonType", "CPDAG");
            menu.setText("Compare to CPDAG...");

            area.setText(comparison.getComparisonString());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);

            area.repaint();

        });

        pag.addActionListener(e -> {
            params.set("graphComparisonType", "PAG");
            menu.setText("Compare to PAG...");

            area.setText(comparison.getComparisonString());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);
            area.repaint();
        });

        return menubar;
    }
}
