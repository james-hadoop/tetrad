package edu.cmu.tetradapp.app.hpc.editor;

import javax.swing.table.DefaultTableModel;
import java.util.Vector;

/**
 * 
 * Feb 10, 2017 1:57:41 PM
 * 
 * @author Chirayu (Kong) Wongchokprasitti, PhD
 * 
 */
public class HpcJobInfoTableModel extends DefaultTableModel {

	private static final long serialVersionUID = 7348053696995503854L;

	private final int buttonColumn;

	public HpcJobInfoTableModel(final Vector<Vector<String>> activeRowData, final Vector<String> activeColumnNames,
			final int buttonColumn) {
		super(activeRowData, activeColumnNames);
		this.buttonColumn = buttonColumn;
	}

	public boolean isCellEditable(int row, int column) {
		if (column == buttonColumn)
			return true;
		return false;
	}

}
