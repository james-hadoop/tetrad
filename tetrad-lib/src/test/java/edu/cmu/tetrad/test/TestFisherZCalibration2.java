package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.LocalConsistencyCriterion;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestDSep;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static junit.framework.TestCase.assertNotNull;

public class TestFisherZCalibration2 {

    @Test
    public void test1() {
        RandomUtil.getInstance().setSeed(105034020L);
        toTest(0.05);
    }

    private void toTest(double alpha) {
        Parameters parameters = new Parameters();
        parameters.set(Params.ALPHA, alpha);
        parameters.set(Params.DEPTH, 2);
        parameters.set(Params.PENALTY_DISCOUNT, 1);
        parameters.set(Params.STRUCTURE_PRIOR, 0);
        parameters.set(Params.COEF_LOW, .2);
        parameters.set(Params.COEF_HIGH, .7);
        int numDraws = 100;
        int sampleSize = 1000;

        System.out.println("# draws = " + numDraws + " sample size = " + sampleSize);

        Graph graph = null;
        DataSet data = null;

        try {
            graph =GraphUtils.loadGraphTxt(new File("/Users/user/Downloads/simulation.for.wayne.long2/save/graph/graph.1.txt"));

            ContinuousTabularDatasetFileReader reader = new ContinuousTabularDatasetFileReader(
                    new File("/Users/user/Downloads/simulation.for.wayne.long2/save/data/data.1.txt").toPath(),
                    Delimiter.TAB
            );

            data = (DataSet) DataConvertUtils.toDataModel(reader.readInData());
        } catch (IOException e) {
            e.printStackTrace();
        }

        assertNotNull(graph);
        assertNotNull(data);

        IndTestFisherZ test1 = (IndTestFisherZ) new FisherZ().getTest(data, parameters);
        IndependenceTest test2 = new LocalConsistencyCriterion().getTest(data, parameters);

        List<Node> variables = data.getVariables();
        graph = GraphUtils.replaceNodes(graph, variables);

        IndependenceTest dsep = new IndTestDSep(graph);

        for (int depth : new int[]{0, 1, 2, 3, 4}) {
            testOneDepth(parameters, numDraws, test1, test2, variables, dsep, depth);
        }
    }

    private void testOneDepth(Parameters parameters, int numDraws, IndTestFisherZ test1, IndependenceTest test2, List<Node> variables, IndependenceTest dsep, int depth) {
        int countSame = 0;
        int fn1 = 0;
        int fp1 = 0;

        int fn2 = 0;
        int fp2 = 0;

        int fn3 = 0;
        int fp3 = 0;

        int fn4 = 0;
        int fp4 = 0;

        int ds = 0;

        for (int i = 0; i < numDraws; i++) {
            Collections.shuffle(variables);
            Collections.shuffle(variables);
            Collections.shuffle(variables);

            Node x = variables.get(0);
            Node y = variables.get(1);

            List<Node> z = new ArrayList<>();
            for (int j = 0; j < depth; j++) {
                z.add(variables.get(j + 2));
            }

            double alpha = parameters.getDouble("alpha");

            boolean fzInd1 = test1.getPValue(x, y, z) > alpha;
            boolean fzInd2 = test1.getPValue2(x, y, z) > alpha;
            boolean fzInd3 = test1.getPValue3(x, y, z) > alpha;
            boolean fzInd4 = test1.getPValue4(x, y, z) > alpha;
            boolean _dsep = dsep.isIndependent(x, y, z);

            if (fzInd1 && !_dsep) fn1++;
            if (!fzInd1 && _dsep) fp1++;

            if (fzInd2 && !_dsep) fn2++;
            if (!fzInd2 && _dsep) fp2++;

            if (fzInd3 && !_dsep) fn3++;
            if (!fzInd3 && _dsep) fp3++;

            if (fzInd4 && !_dsep) fn4++;
            if (!fzInd4 && _dsep) fp4++;

            System.out.println(new IndependenceFact(x, y, z) + " p without 2 = " + test1.getPValue(x, y, z));

            if (_dsep) ds++;
        }

        TextTable table = new TextTable(5, 3);
        table.setToken(0, 1, "FP");
        table.setToken(0, 2, "FN");
        table.setToken(1, 0, "Fisher Z Joe without 2");
        table.setToken(2, 0, "Fisher Z Joe WITH 2");
        table.setToken(3, 0, "Fisher Z Wayne without 2");
        table.setToken(4, 0, "Fisher Z Wayne WITH 2");

        table.setToken(1, 1, "" + fp1);
        table.setToken(1, 2, "" + fn1);

        table.setToken(2, 1, "" + fp2);
        table.setToken(2, 2, "" + fn2);

        table.setToken(3, 1, "" + fp3);
        table.setToken(3, 2, "" + fn3);

        table.setToken(4, 1, "" + fp4);
        table.setToken(4, 2, "" + fn4);

        System.out.println();
        System.out.println("Depth = " + depth);
        System.out.println();
        System.out.println(table);

        System.out.println();

        double alpha = parameters.getDouble(Params.ALPHA);
        System.out.println("alpha = " + alpha);
        System.out.println("alpha^ Joe without 2 = " + fp1 / (double) ds);
        System.out.println("alpha^ Joe WITH 2 = " + fp2 / (double) ds);
        System.out.println("alpha^ Wayne without 2 = " + fp3 / (double) ds);
        System.out.println("alpha^ Wayne WITH 2 = " + fp4 / (double) ds);
    }
}
