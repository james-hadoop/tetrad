package edu.cmu.tetrad.algcomparison.statistic;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.dbmi.data.reader.Delimiter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The arrow precision. This counts arrowheads maniacally, wherever they occur in the graphs.
 * The true positives are the number of arrowheads in both the true and estimated graphs.
 * Thus, if the true contains X*->Y and estimated graph either does not contain an edge from
 * X to Y or else does not contain an arrowhead at X for an edge from X to Y, one false
 * positive is counted. Similarly for false negatives.
 *
 * @author jdramsey
 */
public class DescentantsOfMTorFpr implements Statistic {
    static final long serialVersionUID = 23L;

    @Override
    public String getAbbreviation() {
        return "mtorDescFpr";
    }

    @Override
    public String getDescription() {
        return "Descendants of mTor FPR";
    }

    @Override
    public double getValue(Graph trueGraph, Graph estGraph) {
        HillStats hillStats = new HillStats(trueGraph, estGraph);

        List<Node> all = hillStats.getAll();
        List<Node> positives = hillStats.getPositives();
        List<Node> negatives = hillStats.getNegatives();
        List<Node> tn = hillStats.getTn();
        List<Node> fn = hillStats.getFn();
        List<Node> fp = hillStats.getFp();


        System.out.println();
        System.out.println("All:");

        for (int i = 0; i < all.size(); i++) {
            System.out.println((i + 1) + ". " + all.get(i));
        }

        System.out.println();

        System.out.println("Positives:");

        for (int i = 0; i < positives.size(); i++) {
            System.out.println((i + 1) + ". " + positives.get(i));
        }

        System.out.println();

        System.out.println("Negatives:");

        for (int i = 0; i < negatives.size(); i++) {
            System.out.println((i + 1) + ". " + negatives.get(i));
        }

        System.out.println();

        System.out.println("False Negatives:");

        for (int i = 0; i < fn.size(); i++) {
            System.out.println((i + 1) + ". " + fn.get(i));
        }

        System.out.println();

        System.out.println("False Positives:");

        for (int i = 0; i < fp.size(); i++) {
            System.out.println((i + 1) + ". " + fp.get(i));
        }

        System.out.println();


        return fp.size() / (double) (fp.size() + tn.size());
    }

    @Override
    public double getNormValue(double value) {
        return 1.0 - value;
    }
}
