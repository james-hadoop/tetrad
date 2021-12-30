import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fges;
import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class TestAnneAnalysis3 {


    public static void main(String... args) {
        new TestAnneAnalysis3().run1();
    }

    private void run1() {
        double penalty = 2.0;

        for (int p : new int[]{5, 10, 20}) {
            for (int n : new int[]{50, 50, 100, 500, 1000, 5000, 10000, 50000}) {

                String nString;

                switch (n) {
                    case 50:
                        nString = "50";
                        break;
                    case 100:
                        nString = "100";
                        break;
                    case 500:
                        nString = "500";
                        break;
                    case 1000:
                        nString = "1K";
                        break;
                    case 5000:
                        nString = "5K";
                        break;
                    case 10000:
                        nString = "10K";
                        break;
                    case 50000:
                        nString = "50K";
                        break;
                    default:
                        throw new IllegalArgumentException("Don't have that sample size.");
                }

                try {
                    File filecor = new File("/Users/josephramsey/Downloads/sldisco_cormats_b5K/" +
                            "cormats_p" + p + "_n" + nString + "_b5K.txt");
                    File adjout = new File("/Users/josephramsey/Downloads/sldisco_adjout_b5K/" +
                            "adjmatsout_p" + p + "_n" + nString + "_b5K.txt");

                    BufferedReader incor = new BufferedReader(new FileReader(filecor));
                    incor.readLine();

                    PrintStream out = new PrintStream(adjout);

                    int d = 1;

                    for (int j = 0; j < p; j++) {
                        for (int i = 0; i < p; i++) {
                            out.print("X" + d++);

                            if (!(i == p - 1 && j == p - 1)) {
                                out.print(" ");
                            }
                        }
                    }

                    out.println();

                    for (int k = 0; k < 5000; k++) {
                        System.out.println("k = " + (k + 1));
                        String linecor = incor.readLine();

                        List<Node> vars = new ArrayList<>();
                        for (int i = 0; i < p; i++) vars.add(new ContinuousVariable("x" + (i + 1)));
                        CovarianceMatrix cov = getCov1(linecor, vars, n);

                        edu.cmu.tetrad.search.LinearGaussianBicScore score = new edu.cmu.tetrad.search.LinearGaussianBicScore(cov);
                        score.setPenaltyDiscount(penalty);

                        Fges fges = new Fges(score, 1);

//                            System.out.println(cov);

                        Graph estCpdag = fges.search();

//                        System.out.println(estCpdag);


                        for (int j = 0; j < vars.size(); j++) {
                            for (int i = 0; i < vars.size(); i++) {
                                if (estCpdag.isAdjacentTo(vars.get(j), vars.get(i))) {
                                    out.print("1");
                                } else {
                                    out.print("0");
                                }

                                if (!(i == vars.size() - 1 && j == vars.size() - 1)) {
                                    out.print(" ");
                                }
                            }
                        }

                        out.println();
                    }

                    incor.close();
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @NotNull
    private CovarianceMatrix getCov1(String linecor, List<Node> vars, int n) {
        int v = vars.size();
        double[][] matrix = new double[v][v];
        String[] tokenscor = linecor.split(" ");
        int colcor = 2;

        for (int j = 0; j < v; j++) {
            for (int k = 0; k < v; k++) {
                matrix[j][k] = Double.parseDouble(tokenscor[colcor++]);
            }
        }

        return new CovarianceMatrix(vars, matrix, n);
    }

    @NotNull
    private Graph getTrueG1(String lineadj, List<Node> vars) {
        int v = vars.size();
        Graph trueG = new EdgeListGraph(vars);
        int[][] adjmat = new int[v][v];

        String[] tokensadj = lineadj.split(" ");
        int coladj = 2;

        for (int j = 0; j < v; j++) {
            for (int k = 0; k < v; k++) {
                adjmat[j][k] = Integer.parseInt(tokensadj[coladj++]);
            }
        }

        for (int j = 0; j < v; j++) {
            for (int k = 0; k < v; k++) {
                if (adjmat[j][k] == 1) {
                    trueG.addUndirectedEdge(vars.get(j), vars.get(k));
                }
            }
        }

        return trueG;
    }
}
