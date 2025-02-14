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

import edu.cmu.tetrad.util.NumberFormatUtil;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetradapp.util.DoubleTextField;
import edu.cmu.tetradapp.util.IntTextField;

import javax.swing.*;

/**
 * Edits the properties of a measurement params.
 *
 * @author Joseph Ramsey
 */
class DiscDetIndepParamsEditor extends JComponent {

    /**
     * The parameters object being edited.
     */
    private final Parameters params;

    /**
     * A text field to allow the user to enter the number of dishes to
     * generate.
     */
    private final DoubleTextField alphaField;

    /**
     * A text field to allow the user to enter the number of dishes to
     * generate.
     */
    private final IntTextField depthField;

    /**
     * Constructs a dialog to edit the given gene simulation parameters object.
     */
    public DiscDetIndepParamsEditor(Parameters params) {
        this.params = params;

        // set up text and ties them to the parameters object being edited.
        this.alphaField = new DoubleTextField(params.getDouble("alpha", 0.001), 5,
                NumberFormatUtil.getInstance().getNumberFormat());

        this.alphaField.setFilter((value, oldValue) -> {
            try {
                params().set("alpha", 0.001);
                return value;
            } catch (IllegalArgumentException e) {
                return oldValue;
            }
        });

        this.depthField = new IntTextField(params.getInt("depth", -1), 5);
        this.depthField.setFilter((value, oldValue) -> {
            try {
                params().set("depth", value);
                return value;
            } catch (IllegalArgumentException e) {
                return oldValue;
            }
        });

        buildGui();
    }

    /**
     * Constructs the Gui used to edit properties; called from each constructor.
     * Constructs labels and text fields for editing each property and adds
     * appropriate listeners.
     */
    private void buildGui() {
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        if (this.alphaField != null) {
            Box b1 = Box.createHorizontalBox();
            b1.add(new JLabel("Alpha:"));
            b1.add(Box.createHorizontalGlue());
            b1.add(this.alphaField);
            add(b1);
        }

        Box b2 = Box.createHorizontalBox();
        b2.add(new JLabel("Depth:"));
        b2.add(Box.createHorizontalGlue());
        b2.add(this.depthField);
        add(b2);

        add(Box.createHorizontalGlue());
    }

    /**
     * @return the mapping object being edited. (This probably should not be
     * public, but it is needed so that the textfields can edit the model.)
     */
    private Parameters params() {
        return this.params;
    }
}





