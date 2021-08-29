package edu.cmu.tetradapp.editor;

import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TextTable;
import edu.cmu.tetradapp.model.TabularComparison;
import org.jetbrains.annotations.NotNull;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

public class StatsListEditor extends JPanel {

    private static final long serialVersionUID = 8455624852328328919L;

    private final TabularComparison comparison;
    private final Parameters params;

    private Graph targetGraph;
    private Graph referenceGraph;
    private JTextArea area;
    private List<Graph> referenceGraphs;

    public StatsListEditor(TabularComparison comparison) {
        this.comparison = comparison;
        this.params = comparison.getParams();
        setup();
    }

    private void setup() {

        // We'll leave the underlying model the same but just complain if there's not exactly
        // one reference and one target graph.
        referenceGraphs = comparison.getReferenceGraphs();
        List<Graph> targetGraphs = comparison.getTargetGraphs();

        if (referenceGraphs.size() != 1) throw new IllegalArgumentException("Expecting one comparison graph.");
        if (targetGraphs.size() != 1) throw new IllegalArgumentException("Expecting one target graph.");

        JMenuBar menubar = menubar();
        targetGraph = targetGraphs.get(0);
        show(referenceGraphs, menubar);
    }

    private void show(List<Graph> referenceGraphs, JMenuBar menubar) {
        referenceGraph = GraphUtils.getComparisonGraph(referenceGraphs.get(0), params);
        setLayout(new BorderLayout());
        add(menubar, BorderLayout.NORTH);
        add(getTableDisplay(), BorderLayout.CENTER);
        revalidate();
        repaint();
    }

    private JComponent getTableDisplay() {
        area = new JTextArea();
        area.setText(tableTextWithHeader());
        area.moveCaretPosition(0);
        area.setSelectionStart(0);
        area.setSelectionEnd(0);

        area.setBorder(new EmptyBorder(5, 5, 5, 5));

        area.setFont(new Font(Font.MONOSPACED, Font.BOLD, 14));
        area.setPreferredSize(new Dimension(700, 1200));

        JScrollPane pane = new JScrollPane(area);
        pane.setPreferredSize(new Dimension(700, 700));

        Box b = Box.createVerticalBox();
        b.add(pane);

        return b;
    }

    @NotNull
    private String tableTextWithHeader() {
        TextTable table = tableText();
        return "Comparing target " + comparison.getTargetName() + " to reference " + comparison.getReferenceName()
                + "\n\n" + table;
    }

    @NotNull
    private TextTable tableText() {
        List<Statistic> statistics = statistics();

        TextTable table = new TextTable(statistics.size(), 3);
        NumberFormat nf = new DecimalFormat("0.###");

        for (int i = 0; i < statistics.size(); i++) {
            table.setToken(i, 0, statistics.get(i).getAbbreviation());
            table.setToken(i, 1, statistics.get(i).getDescription());
            double value = statistics.get(i).getValue(referenceGraph, targetGraph, null);
            table.setToken(i, 2, Double.isNaN(value) ? "-" : "" + nf.format(value));
        }

        table.setJustification(TextTable.LEFT_JUSTIFIED);
        return table;
    }

    @NotNull
    private List<Statistic> statistics() {
        List<Statistic> statistics = new ArrayList<>();

        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ArrowheadPrecisionCommonEdges());
        statistics.add(new ArrowheadRecallCommonEdges());
        statistics.add(new AdjacencyTN());
        statistics.add(new AdjacencyTN());
        statistics.add(new AdjacencyTP());
        statistics.add(new AdjacencyTPR());
        statistics.add(new AdjacencyFPR());
        statistics.add(new AdjacencyFN());
        statistics.add(new AdjacencyFP());
        statistics.add(new AdjacencyFN());
        statistics.add(new ArrowheadTN());
        statistics.add(new ArrowheadTP());
        statistics.add(new F1Adj());
        statistics.add(new F1All());
        statistics.add(new F1Arrow());
        statistics.add(new MathewsCorrAdj());
        statistics.add(new MathewsCorrArrow());
        statistics.add(new SHD_CPDAG());
        statistics.add(new NumAmbiguousTriples());
        statistics.add(new PercentAmbiguous());
        statistics.add(new PercentBidirectedEdges());
        statistics.add(new NumberOfEdgesEst());
        statistics.add(new NumberOfEdgesTrue());
        statistics.add(new TailPrecision());
        statistics.add(new TailRecall());
        statistics.add(new TwoCyclePrecision());
        statistics.add(new TwoCycleRecall());
        statistics.add(new TwoCycleFalsePositive());
        statistics.add(new TwoCycleFalseNegative());
        statistics.add(new TwoCycleTruePositive());
        statistics.add(new AverageDegreeEst());
        statistics.add(new AverageDegreeTrue());
        statistics.add(new DensityEst());
        statistics.add(new DensityTrue());
        return statistics;
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
            referenceGraph = GraphUtils.getComparisonGraph(referenceGraphs.get(0), params);

            area.setText(tableTextWithHeader());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);

            area.repaint();

        });

        cpdag.addActionListener(e -> {
            params.set("graphComparisonType", "CPDAG");
            menu.setText("Compare to CPDAG...");
            referenceGraph = GraphUtils.getComparisonGraph(referenceGraphs.get(0), params);

            area.setText(tableTextWithHeader());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);

            area.repaint();

        });

        pag.addActionListener(e -> {
            params.set("graphComparisonType", "PAG");
            menu.setText("Compare to PAG...");
            referenceGraph = GraphUtils.getComparisonGraph(referenceGraphs.get(0), params);

            area.setText(tableTextWithHeader());
            area.moveCaretPosition(0);
            area.setSelectionStart(0);
            area.setSelectionEnd(0);
            area.repaint();
        });

        return menubar;
    }
}
