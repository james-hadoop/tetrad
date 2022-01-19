import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.GRASP;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.OTHER_PERM_ALGS;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.EbicScore;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Grasp;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.sem.LinearSemIm;
import edu.cmu.tetrad.sem.LinearSemPm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.Params;
import edu.cmu.tetrad.util.RandomUtil;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.floor;

public class TestAnneAnalysis2 {


    public static void main(String... args) {
        new TestAnneAnalysis2().run1();
    }

    private void run1() {
        System.out.println("V\tPD\tN\tAFP\tAFN\tAP\tAR\tAHP\tAHR\tSHD");
        int numRuns = 50;

        for (int v : new int[]{5, 10, 20}) {

            //            for (double alpha : new double[]{0.001}) {;//1e-8, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8}) {
                for (double pd : new double[]{1}) {
                try {
                    List<Node> vars = new ArrayList<>();
                    for (int i = 0; i < v; i++) vars.add(new ContinuousVariable("x" + (i + 1)));

                    for (int n : new int[]{50, 100, 500, 1000, 5000, 10000, 50000}) {//, 100000}) {
                        ExpectedStat efp = new ExpectedStat();
                        ExpectedStat efn = new ExpectedStat();
                        ExpectedStat eap = new ExpectedStat();
                        ExpectedStat ear = new ExpectedStat();
                        ExpectedStat eahp = new ExpectedStat();
                        ExpectedStat eahr = new ExpectedStat();
                        ExpectedStat eshd = new ExpectedStat();

                        for (int i = 0; i < numRuns; i++) {
                            Graph trueCpdag;
                            Graph dag;
                            ICovarianceMatrix cov;

                            {
                                dag = getTrueG2(vars, .2, .7);
                                trueCpdag = SearchGraphUtils.cpdagForDag(dag);
                                cov = getCov(n, dag);
                            }

                            Graph estCpdag;

                            {
                                {
                                    GRASP boss = new GRASP(new EbicScore(), new FisherZ());
                                    Parameters parameters = new Parameters();
                                    parameters.set(Params.ALPHA, 0.001);
                                    parameters.set(Params.PENALTY_DISCOUNT, pd);
                                    parameters.set(Params.ZS_RISK_BOUND, 1);
                                    parameters.set(Params.VERBOSE, false);
                                    parameters.set(Params.BOSS_SCORE_TYPE, false);
                                    parameters.set(Params.BREAK_TIES, true);
                                    parameters.set(Params.OUTPUT_CPDAG, true);
                                    parameters.set(Params.GRASP_USE_SCORE, true);
                                    estCpdag = boss.search(cov, parameters, dag);
                                }

//                                {
//                                    edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc pc
//                                            = new edu.cmu.tetrad.algcomparison.algorithm.oracle.cpdag.Pc(new FisherZ());
//                                    Parameters parameters = new Parameters();
//                                    parameters.set(Params.STABLE_FAS, true);
//                                    parameters.set(Params.ALPHA, alpha);
//                                    parameters.set(Params.VERBOSE, false);
//                                    estCpdag = pc.search(cov, parameters, trueCpdag);
//                                }

//                                {
//                                    Algorithm boss = new Fges(new ZhangShenBoundScore());
//                                    Parameters parameters = new Parameters();
//                                    parameters.set(Params.PENALTY_DISCOUNT, 2);
//                                    parameters.set(Params.ZS_RISK_BOUND, 1);
//                                    parameters.set(Params.VERBOSE, false);
//                                    estCpdag = boss.search(cov, parameters, dag);
//                                }

                                estCpdag = GraphUtils.replaceNodes(estCpdag, trueCpdag.getNodes());

                                efp.add(new AdjacencyFP().getValue(trueCpdag, estCpdag, null));
                                efn.add(new AdjacencyFN().getValue(trueCpdag, estCpdag, null));
                                eap.add(new AdjacencyPrecision().getValue(trueCpdag, estCpdag, null));
                                ear.add(new AdjacencyRecall().getValue(trueCpdag, estCpdag, null));
                                eahp.add(new ArrowheadPrecision().getValue(trueCpdag, estCpdag, null));
                                eahr.add(new ArrowheadRecall().getValue(trueCpdag, estCpdag, null));
                                eshd.add(new SHD().getValue(trueCpdag, estCpdag, null));
                            }
                        }

                        NumberFormat nf = new DecimalFormat("0.00");

                        System.out.println(v + "\t" + pd + "\t" + n + "\t" + nf.format(efp.expected()) + "\t" + nf.format(efn.expected())
                                + "\t" + nf.format(eap.expected()) + "\t" + nf.format(ear.expected())
                                + "\t" + nf.format(eahp.expected()) + "\t" + nf.format(eahr.expected())
                                + "\t" + nf.format(eshd.expected()));

                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static class ExpectedStat {
        private double sum;
        private double count;

        public ExpectedStat() {
            sum = 0.0;
            count = 0.0;
        }

        public void add(double stat) {
            if (!Double.isNaN(stat)) {
                sum += stat;
                count = count + 1;
            }
        }

        public double expected() {
            return sum / count;
        }
    }

    private ICovarianceMatrix getCov(int n, Graph dag) {
        LinearSemPm pm = new LinearSemPm(dag);

        Parameters parameters = new Parameters();
        parameters.set(Params.COEF_LOW, 0.1);
        parameters.set(Params.COEF_HIGH, 0.5);
        parameters.set(Params.COEF_SYMMETRIC, true);
        parameters.set(Params.VAR_LOW, .25);
        parameters.set(Params.VAR_HIGH, 4);

        LinearSemIm im = new LinearSemIm(pm, parameters);
        DataSet dataSet = im.simulateData(n, false);
        return new CovarianceMatrix(dataSet);
    }

    private Graph getTrueG2(List<Node> vars, double from, double to) {
        double density = RandomUtil.getInstance().nextUniform(from, to);
        int v = vars.size();
        int e = (int) floor((density * (v * (v - 1.))) / 2.);

        return GraphUtils.randomGraph(vars, 0, e, 100, 100,
                100, true);
    }
}
