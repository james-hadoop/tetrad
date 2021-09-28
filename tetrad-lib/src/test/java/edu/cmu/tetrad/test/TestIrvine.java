package edu.cmu.tetrad.test;

import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.*;
import edu.cmu.tetrad.algcomparison.independence.ConditionalGaussianLRT;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ConditionalGaussianBicScore;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.data.DataReader;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import edu.pitt.dbmi.data.reader.tabular.MixedTabularDatasetFileReader;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;

public class TestIrvine {


    @Test
    public void test1() {
        File dir = new File("/Users/josephramsey/Documents/GitHub/tetrad-example-analyses");

        File[] files = dir.listFiles();

        for (File dataDir : files) {
            if (dataDir.getName().startsWith("README")) continue;
            if (dataDir.getName().startsWith("LICENSE")) continue;
            if (dataDir.getName().startsWith(".")) continue;

            File dataDir2 = new File(dataDir, "tetrad.formatted.data");

            File[] obj = dataDir2.listFiles();

            if (obj == null) {
                System.out.println("Couldn't fine the directory for " + dataDir.getName());
            }

            if (obj != null) {
                for (File f : obj) {
                    if (f.getName().startsWith(".")) continue;

                    String name = f.getName();

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.PC, Algorithm.CPC, Algorithm.PCMAX};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FGES, Algorithm.BOSS};

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FCI, Algorithm.FCIMAX, Algorithm.RFCI};
                    Algorithm[] algorithms = new Algorithm[]{Algorithm.BOSS, Algorithm.BFCI};

                    IndependenceWrapper test;
                    ScoreWrapper score;

                    Parameters params = new Parameters();
                    params.set(Params.ALPHA, 0.001);
                    params.set(Params.PENALTY_DISCOUNT, 2);
                    params.set(Params.NUM_STARTS, 10);
                    params.set(Params.COMPLETE_RULE_SET_USED, true);

                    for (Algorithm algorithm : algorithms) {
                        File outfile = new File(dataDir, "results1/"
                                + algorithm.toString().toLowerCase(Locale.ROOT) + "." + name);

                        if (outfile.exists()) continue;

                        DataSet _data = null;

                        if (name.contains("continuous")) {

                            try {
                                Data data = new ContinuousTabularDatasetFileReader(f.toPath(), Delimiter.TAB).readInData();
                                _data = (DataSet) DataConvertUtils.toDataModel(data);

                                if (_data == null) {
                                    System.out.println("Couldn't find that dataset: " + f);
                                } else {
                                    _data = DataUtils.shuffleColumns(_data);
                                }

                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            assert _data != null;
                            test = new FisherZ();
                            score = new edu.cmu.tetrad.algcomparison.score.LinearGaussianBicScore();
                        } else {

//                            continue;

                            try {
                                Data data = new MixedTabularDatasetFileReader(f.toPath(), Delimiter.TAB, 20).readInData();
                                _data = (DataSet) DataConvertUtils.toDataModel(data);

                                if (_data == null) {
                                    System.out.println("Couldn't find that dataset: " + f);
                                } else {
                                    _data = DataUtils.shuffleColumns(_data);
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            assert _data != null;
                            test = new ConditionalGaussianLRT();
                            score = new ConditionalGaussianBicScore();
                        }

                        edu.cmu.tetrad.algcomparison.algorithm.Algorithm search = null;

                        if (algorithm == Algorithm.PC) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc(test);
                        } else if (algorithm == Algorithm.CPC) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Cpc(test);
                        } else if (algorithm == Algorithm.PCMAX) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.PcMax(test);
                        } else if (algorithm == Algorithm.FGES) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Fges(score);
                        } else if (algorithm == Algorithm.BOSS) {
                            search = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.BOSS(score, test);
                        } else if (algorithm == Algorithm.FCI) {
                            search = new Fci(test);
                        } else if (algorithm == Algorithm.FCIMAX) {
                            search = new FciMax(test);
                        } else if (algorithm == Algorithm.RFCI) {
                            search = new Rfci(test);
                        } else if (algorithm == Algorithm.GFCI) {
                            search = new Gfci(score, test);
                        } else if (algorithm == Algorithm.BFCI) {
                            search = new BFCI(score, test);
                        }

                        assert search != null;
                        Graph graph = search.search(_data, params, null);

                        try {
                            if (!outfile.getParentFile().mkdirs()) {
                                System.out.println("Couldn't make directories");
                            }
                            PrintWriter writer = new PrintWriter(outfile);
                            writer.println(graph);
                            writer.close();
                        } catch (FileNotFoundException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
        }
    }

    @Test
    public void test2() {
        File dir = new File("/Users/josephramsey/Documents/GitHub/tetrad-example-analyses");

        File[] files = dir.listFiles();

        for (File dataDir : files) {
            if (dataDir.getName().startsWith("README")) continue;
            if (dataDir.getName().startsWith("LICENSE")) continue;
            if (dataDir.getName().startsWith(".")) continue;

            File dataDir2 = new File(dataDir, "tetrad.formatted.data");

            File[] obj = dataDir2.listFiles();

            if (obj == null) {
                System.out.println("Couldn't fine the directory for " + dataDir.getName());
            }

            if (obj != null) {
                for (File f : obj) {
                    if (f.getName().startsWith(".")) continue;

                    String name = f.getName();

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.PC, Algorithm.CPC, Algorithm.PCMAX};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FGES, Algorithm.BOSS};

//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.FCI, Algorithm.FCIMAX, Algorithm.RFCI};
//                    Algorithm[] algorithms = new Algorithm[]{Algorithm.GFCI, Algorithm.BFCI};

                    Parameters params = new Parameters();
                    params.set(Params.ALPHA, 0.001);
                    params.set(Params.PENALTY_DISCOUNT, 2);
                    params.set(Params.NUM_STARTS, 5);
                    params.set(Params.COMPLETE_RULE_SET_USED, true);

                    for (Algorithm algorithm : Algorithm.values()) {
                        File outfile = new File(dataDir, "results1/"
                                + algorithm.toString().toLowerCase(Locale.ROOT) + "." + name);
                        File knowledgeFile = new File(dataDir, "ground.truth.knowledge/knowledge." + name);

                        if (!knowledgeFile.exists()) {
                            continue;
                        }

                        IKnowledge knowledge = null;

                        try {
                            knowledge = new DataReader().parseKnowledge(knowledgeFile);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        if (knowledge == null) {
                            continue;
                        }

//                        System.out.println(knowledge);

                        Graph graph = GraphUtils.loadGraphTxt(outfile);
                        boolean satisfies = true;

                        for (Edge edge : graph.getEdges()) {
                            Node x = edge.getNode1();
                            Node y = edge.getNode2();

                            String xn = x.getName();
                            String yn = y.getName();

                            if (graph.isDirectedFromTo(x, y) && knowledge.isForbidden(xn, yn)) {
//                                System.out.println("Violates: " + edge + " " + outfile.getName());
                                satisfies = false;
                            } else if (graph.isDirectedFromTo(y, x) && knowledge.isForbidden(yn, xn)) {
//                                System.out.println("Violates: " + edge + " " + outfile.getName());
                                satisfies = false;
                            }
//                            else if (graph.isDirectedFromTo(y, x) && knowledge.isRequired(xn, yn)) {
//                                satisfies = false;
//                            } else if (graph.isDirectedFromTo(x, y) && knowledge.isRequired(yn, xn)) {
//                                satisfies = false;
//                            }
                        }

                        if (satisfies) {
                            System.out.println("Satisfies ground truth: " + outfile.getName() + "(" + graph.getNumEdges() + " edge)");
                        }
                    }
                }
            }
        }
    }

    private enum Algorithm {PC, CPC, PCMAX, FGES, BOSS, FCI, FCIMAX, RFCI, GFCI, BFCI}
}

