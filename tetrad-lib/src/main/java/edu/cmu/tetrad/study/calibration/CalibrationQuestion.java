package edu.cmu.tetrad.study.calibration;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.simulation.LinearFisherModel;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.LargeScaleSimulation;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;
import edu.pitt.dbmi.data.reader.tabular.VerticalDiscreteTabularDatasetFileReader;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.graph.GraphUtils.loadGraphTxt;
import static java.lang.Math.abs;

public class CalibrationQuestion {

    public static void main(String... args) {
//        try {
        scenario9();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

    private static void scenario1() {
        int c = 0;
        int total = 0;
        int sampleSize = 1000;
        int numRuns = 1000;
        int numVars = 10;
        int avgDegree = 4;
        int numEdges = avgDegree * numVars / 2;
        double p = avgDegree / (double) (numVars - 1);

        System.out.println("p = " + p);

        for (int i = 0; i < numRuns; i++) {
            Node x = new ContinuousVariable("X");
            Node y = new ContinuousVariable("Y");
            Node z = new ContinuousVariable("Z");

            List<Node> nodes = new ArrayList<>();
            nodes.add(x);
            nodes.add(y);
            nodes.add(z);

            for (int n = 3; n <= numVars; n++) {
                nodes.add(new ContinuousVariable("V" + n));
            }

            Graph gt = GraphUtils.randomGraph(nodes, 0, numEdges, avgDegree, 100, 100, false);

            gt.removeEdge(x, y);
            gt.removeEdge(y, z);
            gt.removeEdge(x, z);


            gt.addDirectedEdge(x, y);
            gt.addDirectedEdge(y, z);

            if (RandomUtil.getInstance().nextDouble() <= p) {
                gt.addDirectedEdge(x, z);
            }

            SemPm pm = new SemPm(gt);
            SemIm im = new SemIm(pm);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));

            Graph ge = fges.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

            System.out.println("gt = " + gt + " ge = " + ge);

            if (ge.isAdjacentTo(x, y) && ge.isAdjacentTo(y, z)) {
                if (gt.isAdjacentTo(x, z)) {
                    c++;
                }

                total++;
            }

            System.out.println("c = " + c + " total = " + total + " q = " + (c / (double) total));

        }

        System.out.println("p = " + p + " q = " + (c / (double) total) + " numEdges = " + numEdges);
    }

    private static void scenario2() {
        int c = 0;
        int total = 0;
        int sampleSize = 1000;
        int numRuns = 500;
        int numVars = 80;
        int avgDegree = 4;
        int numEdges = avgDegree * numVars / 2;
        double p = avgDegree / (double) (numVars - 1);

        System.out.println("p = " + p);

        for (int i = 0; i < numRuns; i++) {
            Node x = new ContinuousVariable("X");
            Node y = new ContinuousVariable("Y");
            Node z = new ContinuousVariable("Z");

            List<Node> nodes = new ArrayList<>();
            nodes.add(x);
            nodes.add(y);
            nodes.add(z);

            for (int n = 3; n <= numVars; n++) {
                nodes.add(new ContinuousVariable("V" + n));
            }

            Graph gt = GraphUtils.randomGraph(nodes, 0, numEdges, 100, 100, 100, false);

            SemPm pm = new SemPm(gt);
            SemIm im = new SemIm(pm);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));

            Graph ge = fges.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

            if (ge.isAdjacentTo(x, y) && ge.isAdjacentTo(y, z) && ge.isAdjacentTo(x, z)) {
                if (gt.isAdjacentTo(x, z)) {
                    c++;
                }

                total++;

                System.out.println("Run " + i + " p = " + p + " c = " + c + " total = " + total + " q = " + (c / (double) total));
            }

        }

        System.out.println("p = " + p + " q = " + (c / (double) total));
    }

    private static void scenario3() {
        int c = 0;
        int d = 0;
        int totalc = 0;
        int totald = 0;
        int sampleSize = 1000;
        int numRuns = 500;
        int numVars = 30;
        int avgDegree = 4;
        int numEdges = avgDegree * numVars / 2;
        double p = avgDegree / (double) (numVars - 1);

        System.out.println("p = " + p);

        for (int i = 0; i < numRuns; i++) {
            Node x = new ContinuousVariable("X");
            Node y = new ContinuousVariable("Y");
            Node z = new ContinuousVariable("Z");

            List<Node> nodes = new ArrayList<>();
            nodes.add(x);
            nodes.add(y);
            nodes.add(z);

            for (int n = 3; n <= numVars; n++) {
                nodes.add(new ContinuousVariable("V" + n));
            }

            Graph gt = GraphUtils.randomGraph(nodes, 0, numEdges, 100, 100, 100, false);

            SemPm pm = new SemPm(gt);

            Parameters parameters = new Parameters();
            parameters.set("coefLow", 0.1);
            parameters.set("coefHigh", 0.5);

            SemIm im = new SemIm(pm, parameters);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));

            Graph ge = s.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

            {
                if (ge.isAdjacentTo(x, y) && ge.isAdjacentTo(y, z)) {
                    c++;
                }

                totalc++;
            }

            if (gt.isAdjacentTo(x, z)) {
                if (ge.isAdjacentTo(x, y) && ge.isAdjacentTo(y, z)) {
                    d++;
                }

                totald++;
            }

            System.out.println("Run " + i + " P(XYe & YZe) = " + (c / (double) totalc) + " P(XYe & YZe | XZt) = " + (d / (double) totald));
        }

        System.out.println("p = " + p + " q = " + (c / (double) totalc));
    }

    private static void scenario4() {
        int c = 0;
        int total = 0;
        int sampleSize = 1000;
        int numRuns = 1;
        int numVars = 20;
        int avgDegree = 10;
        int numEdges = avgDegree * numVars / 2;
        double p = avgDegree / (double) (numVars - 1);

        System.out.println("p = " + p);

        for (int i = 0; i < numRuns; i++) {
            Node x = new ContinuousVariable("X");
            Node y = new ContinuousVariable("Y");
            Node z = new ContinuousVariable("Z");

            List<Node> nodes = new ArrayList<>();
            nodes.add(x);
            nodes.add(y);
            nodes.add(z);

            for (int n = 3; n <= numVars; n++) {
                nodes.add(new ContinuousVariable("V" + n));
            }

            Graph gt = GraphUtils.randomGraph(nodes, 0, numEdges, 100, 100, 100, false);

            SemPm pm = new SemPm(gt);

            Parameters parameters = new Parameters();
            parameters.set("coefLow", 0.1);
            parameters.set("coefHigh", 0.5);

            SemIm im = new SemIm(pm);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));

            Graph ge = s.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

            ChoiceGenerator gen = new ChoiceGenerator(numVars, 3);
            int[] choice;
            int t = 0;

            while ((choice = gen.next()) != null) {
                List<Node> v = GraphUtils.asList(choice, nodes);

                Node v1 = v.get(0);
                Node v2 = v.get(1);
                Node v3 = v.get(2);

                if (ge.isAdjacentTo(v1, v2) && ge.isAdjacentTo(v2, v3) && !ge.isAdjacentTo(v1, v3)) {
                    if (gt.isAdjacentTo(v1, v3)) {
                        c++;
                    }

                    total++;

                    System.out.println("Triple " + ++t + " p = " + p + " c = " + c + " total = " + total + " q = " + (c / (double) total));

                }

            }
        }

        System.out.println("p = " + p + " q = " + (c / (double) total));
    }

    // Make GT, simulate data, run FGES, yielding GE. Then find a nonredundant set of possible false negative shields ~XZ
    // touching all UTFP legs in GE and count for how many of these XZt is in GT.
    private static void scenario5() {
        System.out.println("SEED = " + RandomUtil.getInstance().getSeed());

        int sampleSize = 1000;
        int[] numVars = new int[]{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
        int[] avgDegree = new int[]{2, 4, 6, 8};

        NumberFormat nf = new DecimalFormat("0.00");

        String[] algorithms = {"FGES", "PC", "CPC", "PC-Max"};

        PrintStream out = null;
        PrintStream rOut = null;

        for (String algorithm : algorithms) {

            System.out.println("\n============================");
            System.out.println("Algorithm = " + algorithm);
            System.out.println("============================\n");


            try {
                out = new PrintStream(new File("/Users/user/Tetrad/tetrad-lib/src/main/" +
                        "java/edu/cmu/tetrad/study/calibration/data.for.calibration." + algorithm + ".txt"));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            try {
                rOut = new PrintStream(new File("/Users/user/Tetrad/tetrad-lib/src/main/" +
                        "java/edu/cmu/tetrad/study/calibration/data.for.calibration.rOut." + algorithm + ".txt"));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            if (out == null) throw new NullPointerException("out not initialized");
            if (rOut == null) throw new NullPointerException("rOut not initialized");

            out.println("AvgDeg\t#Vars\tL\tP\tDensity\tSparsity\tR2\tBound\tAHP\tAHPC");
            rOut.println("R2");


            for (int _numVars : numVars) {
                for (int _avgDegree : avgDegree) {

                    double numEdges = (double) _avgDegree * _numVars / 2.;
                    double density = _avgDegree / (double) (_numVars - 1);

                    Graph G2 = GraphUtils.randomGraph(_numVars, 0, (int) numEdges, 100, 100, 100, false);

                    SemPm pm = new SemPm(G2);

                    Parameters parameters = new Parameters();
                    parameters.set("coefLow", 0.);
                    parameters.set("coefHigh", 0.7);

                    SemIm im = new SemIm(pm);

                    DataSet data = im.simulateData(sampleSize, false);

                    GraphSearch s = null;

                    switch (algorithm) {
                        case "FGES":
                            s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
                            break;
                        case "PC":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.FAS_SEPSETS);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                        case "CPC":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.CONSERVATIVE);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                        case "PC-Max":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.MAX_P);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                    }

                    if (s == null) {
                        throw new NullPointerException("Unrecognized algorthm type: " + algorithm);
                    }

                    Graph R = s.search();
                    R = GraphUtils.replaceNodes(R, G2.getNodes());

                    List<Node> nodes = R.getNodes();

                    ChoiceGenerator gen = new ChoiceGenerator(nodes.size(), 3);
                    int[] choice;

                    Set<Edge> L = new HashSet<>();
                    Set<Edge> M = new HashSet<>();

                    while ((choice = gen.next()) != null) {
                        List<Node> v = GraphUtils.asList(choice, nodes);

                        Node v1 = v.get(0);
                        Node v2 = v.get(1);
                        Node v3 = v.get(2);

                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v1, v3, v2);
                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v1, v2, v3);
                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v2, v1, v3);
                    }

                    Set<Edge> L1 = new HashSet<>();
                    Set<Edge> S1 = new HashSet<>();

                    for (int i = 0; i < nodes.size() - 1; i++) {
                        List<Node> adj = R.getAdjacentNodes(nodes.get(i));

                        for (int j = 1; j < adj.size(); j++) {

                            for (int k = j + 1; k < adj.size(); k++) {
                                boolean b1 = L1.contains(Edges.undirectedEdge(nodes.get(i), adj.get(j)));
                                boolean b2 = L1.contains(Edges.undirectedEdge(nodes.get(i), adj.get(k)));

                                if (!R.isAdjacentTo(adj.get(j), adj.get(k)) && !(b1 && b2)) {
                                    L1.add(Edges.undirectedEdge(nodes.get(i), adj.get(j)));
                                    L1.add(Edges.undirectedEdge(nodes.get(i), adj.get(k)));
                                    S1.add(Edges.undirectedEdge(adj.get(j), adj.get(k)));
                                }
                            }
                        }
                    }

                    int A = 0;
                    int AStar = 0;

                    for (Edge e2 : R.getEdges()) {
                        Edge e1 = G2.getEdge(e2.getNode1(), e2.getNode2());

                        if (e1 == null) {
                            continue;
                        }

                        Node n1 = e1.getNode1();
                        Node n2 = e1.getNode2();

                        if (e2.getProximalEndpoint(n1) == Endpoint.ARROW) {
                            if (L1.contains(Edges.undirectedEdge(e2.getNode1(), e2.getNode2()))) {
                                AStar++;
                            }

                            A++;
                        }

                        if (e2.getProximalEndpoint(n2) == Endpoint.ARROW) {
                            if (L1.contains(Edges.undirectedEdge(e2.getNode1(), e2.getNode2()))) {
                                AStar++;
                            }

                            A++;
                        }
                    }

                    Set<Edge> S2 = new HashSet<>();

                    for (int i = 0; i < nodes.size() - 1; i++) {
                        List<Node> adj = R.getAdjacentNodes(nodes.get(i));

                        for (int j = 1; j < adj.size(); j++) {

                            for (int k = j + 1; k < adj.size(); k++) {
                                if (!S1.contains(Edges.undirectedEdge(adj.get(j), adj.get(k)))) continue;

                                if (G2.isAdjacentTo(adj.get(j), adj.get(k))) {
                                    S2.add(Edges.undirectedEdge(adj.get(j), adj.get(k)));
                                }
                            }
                        }
                    }

                    UtRStatistic utr = new UtRStatistic();
                    double r2 = utr.getValue(G2, R, data);

                    int tp = 0;

                    for (Edge e2 : R.getEdges()) {
                        Edge e1 = G2.getEdge(e2.getNode1(), e2.getNode2());

                        if (e1 == null) {
                            continue;
                        }

                        Node n1 = e1.getNode1();
                        Node n2 = e1.getNode2();

                        if (e1.getProximalEndpoint(n1) == Endpoint.ARROW
                                && e2.getProximalEndpoint(n1) == Endpoint.ARROW) {
                            tp++;
                        }

                        if (e1.getProximalEndpoint(n2) == Endpoint.ARROW
                                && e2.getProximalEndpoint(n2) == Endpoint.ARROW) {
                            tp++;
                        }
                    }

                    Statistic ahp = new ArrowheadPrecision();

                    double ahp2 = ahp.getValue(G2, R, data);
                    double ahpc2 = tp / (double) A;

                    double rho = 0.5;

                    double bound = 1. - 2 * rho * (A / (double) AStar) * density;

                    // Dom't divide by zero anywhere.

                    if (S2.isEmpty()) {
                        continue;
                    }

                    if (AStar == 0) {
                        continue;
                    }

                    double d = _avgDegree / (double) (_numVars - 1);

                    System.out.println(
                            "S2/S1 = " + (getFormat(nf, ((double) S2.size()) / (S1.size()))
                                    + " d = " + getFormat(nf, d))
                                    + " Avg degree = " + _avgDegree
                                    + " num vars = " + _numVars
                                    + " R.numedges = " + R.getNumEdges()
                                    + " L1 = " + L1.size()
                                    + " A = " + A
                                    + " S2 = " + S2.size()
                    );

                    out.println(
                            _avgDegree + "\t" + _numVars
                                    + "\t" + L.size()
                                    + "\t" + A
                                    + "\t" + getFormat(nf, density)
                                    + "\t" + getFormat(nf, 1.0 - density)
                                    + "\t" + getFormat(nf, r2)
                                    + "\t" + getFormat(nf, bound)
                                    + "\t" + getFormat(nf, ahp2)
                                    + "\t" + getFormat(nf, ahpc2)
                    );


                    if (Double.isNaN(r2)) continue;

                    rOut.println(getFormat(nf, r2));
                }
            }

            out.close();
            rOut.close();
        }

    }

    private static void scenario6() {
        System.out.println("SEED = " + RandomUtil.getInstance().getSeed());

        int sampleSize = 1000;
        int[] numVars = new int[]{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
        int[] avgDegree = new int[]{2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

        NumberFormat nf = new DecimalFormat("0.00");

        String[] algorithms = {
                "FGES", "PC", "CPC", "PC-Max",
                "FASK"};

        PrintStream out = null;

        for (String algorithm : algorithms) {
            System.out.println("\n============================");
            System.out.println("Algorithm = " + algorithm);
            System.out.println("============================\n");

            try {
                out = new PrintStream(new File("/Users/user/Tetrad/tetrad-lib/src/main/" +
                        "java/edu/cmu/tetrad/study/calibration/data.for.calibration." + algorithm + ".txt"));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            if (out == null) throw new NullPointerException("out not initialized");

            out.println("AvgDeg\tVars\tDensity\tSparsity\tR2\tAP\tAR\tAHP\tAHPC\tAHR\tAHRC\tA2\tU2\tF2\tf2\tS2\tSh1\tSh2\tE");

            for (int _numVars : numVars) {
                for (int _avgDegree : avgDegree) {
                    double numEdges = (double) _avgDegree * _numVars / 2.;

                    if (numEdges < 0 || numEdges > _numVars * (_numVars - 1) / 2.) {
                        continue;
                    }

                    double density = _avgDegree / (double) (_numVars - 1);

                    Graph G2 = GraphUtils.randomGraph(_numVars, 0, (int) numEdges, 100, 100, 100, false);

                    LargeScaleSimulation lgs = new LargeScaleSimulation(G2);

                    lgs.setCoefRange(0.0, 1.0);
                    lgs.setIncludeNegativeCoefs(true);
                    lgs.setVarRange(1, 2);
                    lgs.setIncludePositiveCoefs(true);

                    lgs.setErrorsNormal(false);

                    DataSet data = lgs.simulateDataFisher(100, sampleSize, 1e-3, false);

                    data = DataUtils.shuffleColumns(data);

                    GraphSearch s = null;

                    switch (algorithm) {
                        case "FGES":
                            s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
                            break;
                        case "PC":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.FAS_SEPSETS);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                        case "CPC":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.CONSERVATIVE);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                        case "PC-Max":
                            s = new PcAll(new IndTestFisherZ(data, 0.001), null);
                            ((PcAll) s).setColliderDiscovery(PcAll.ColliderDiscovery.MAX_P);
                            ((PcAll) s).setConflictRule(PcAll.ConflictRule.PRIORITY);
                            ((PcAll) s).setFasType(PcAll.FasType.STABLE);
                            ((PcAll) s).setConcurrent(PcAll.Concurrent.NO);
                            break;
                        case "FASK":
                            Parameters parameters = new Parameters();
                            IndependenceWrapper test = new FisherZ();
                            s = new Fask(data, test.getTest(data, parameters));
                            ((Fask) s).setSkewEdgeThreshold(0.5);
                            break;
                    }

                    if (s == null) {
                        throw new NullPointerException("Unrecognized algorthm type: " + algorithm);
                    }

                    long start = System.currentTimeMillis();
                    Graph R = s.search();
                    long stop = System.currentTimeMillis();
                    long elapsed = stop - start;

                    R = GraphUtils.replaceNodes(R, G2.getNodes());

                    Statistic ap = new AdjacencyPrecision();
                    double ap2 = ap.getValue(G2, R, data);

                    Statistic ar = new ArrowheadRecall();
                    double ar2 = ar.getValue(G2, R, data);

                    Statistic ahp = new ArrowheadPrecision();
                    double ahp2 = ahp.getValue(G2, R, data);

                    Statistic ahr = new ArrowheadRecall();
                    double ahr2 = ahr.getValue(G2, R, data);

                    R = getCommonGraph(R, G2);  // for AHPC

                    Statistic ahpc = new ArrowheadPrecision();
                    double ahpc2 = ahpc.getValue(G2, R, data);

                    Statistic ahrc = new ArrowheadRecall();
                    double ahrc2 = ahrc.getValue(G2, R, data);


                    List<Node> nodes = R.getNodes();

                    ChoiceGenerator gen = new ChoiceGenerator(nodes.size(), 3);
                    int[] choice;

                    Set<Edge> L = new HashSet<>();
                    Set<Edge> M = new HashSet<>();

                    while ((choice = gen.next()) != null) {
                        List<Node> v = GraphUtils.asList(choice, nodes);

                        Node v1 = v.get(0);
                        Node v2 = v.get(1);
                        Node v3 = v.get(2);

                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v1, v3, v2);
                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v1, v2, v3);
                        collectUnshieldedTripleLegsAndShieldsInR(R, L, M, v2, v1, v3);
                    }


                    int Ut = 0;
                    Set<Edge> record1 = new HashSet<>();
                    Set<Edge> record2 = new HashSet<>();
                    Map<Edge, Set<Triple>> sharedTriples1 = new HashMap<>();
                    Map<Edge, Set<Triple>> sharedTriples2 = new HashMap<>();

                    for (Node node : nodes) {
                        List<Node> adj = R.getAdjacentNodes(node);

                        for (int j = 0; j < adj.size(); j++) {
                            for (int k = j + 1; k < adj.size(); k++) {
                                if (!R.isAdjacentTo(adj.get(j), adj.get(k))) {
                                    Edge edge1 = Edges.undirectedEdge(node, adj.get(j));
                                    Edge edge2 = Edges.undirectedEdge(node, adj.get(k));
                                    if (record1.contains(edge1) && record1.contains(edge2)) {
                                        continue;
                                    }

                                    record1.add(edge1);
                                    record1.add(edge2);

                                    sharedTriples1.computeIfAbsent(edge1, k1 -> new HashSet<>());
                                    sharedTriples1.computeIfAbsent(edge2, k1 -> new HashSet<>());

                                    sharedTriples1.get(edge1).add(new Triple(adj.get(j), node, adj.get(k)));
                                    sharedTriples1.get(edge2).add(new Triple(adj.get(j), node, adj.get(k)));

                                    if (G2.isAdjacentTo(node, adj.get(j)) && G2.isAdjacentTo(node, adj.get(k))
                                            && !G2.isAdjacentTo(adj.get(j), adj.get(k))) {
                                        sharedTriples2.computeIfAbsent(edge1, k1 -> new HashSet<>());
                                        sharedTriples2.computeIfAbsent(edge2, k1 -> new HashSet<>());

                                        sharedTriples2.get(edge1).add(new Triple(adj.get(j), node, adj.get(k)));
                                        sharedTriples2.get(edge2).add(new Triple(adj.get(j), node, adj.get(k)));

                                    }

                                    Ut++;
                                }
                            }
                        }
                    }

                    int shared1 = 0;

                    for (Edge edge : sharedTriples1.keySet()) {
                        if (sharedTriples1.get(edge).size() > 1) shared1++;
                    }

                    int shared2 = 0;

                    for (Edge edge : sharedTriples2.keySet()) {
                        if (sharedTriples2.get(edge).size() > 1) shared2++;
                    }

                    int A = 0;

                    for (Edge e2 : R.getEdges()) {
                        if (e2.isDirected()) {
                            A++;
                        }
                    }

                    int S2 = 0;

                    for (Node node : nodes) {
                        List<Node> adj = R.getAdjacentNodes(node);

                        for (int j = 0; j < adj.size(); j++) {
                            for (int k = j + 1; k < adj.size(); k++) {
                                if (!R.isAdjacentTo(adj.get(j), adj.get(k))) {
                                    if (G2.isAdjacentTo(node, adj.get(j)) && G2.isAdjacentTo(node, adj.get(k))) {
                                        if (G2.isAdjacentTo(adj.get(j), adj.get(k))) {
                                            Edge edge1 = Edges.undirectedEdge(node, adj.get(j));
                                            Edge edge2 = Edges.undirectedEdge(node, adj.get(k));
                                            if (record2.contains(edge1) && record2.contains(edge2)) {
                                                continue;
                                            }

                                            record2.add(edge1);
                                            record2.add(edge2);

                                            S2++;

                                            sharedTriples2.computeIfAbsent(edge1, k1 -> new HashSet<>());
                                            sharedTriples2.computeIfAbsent(edge2, k1 -> new HashSet<>());

                                            sharedTriples2.get(edge1).add(new Triple(adj.get(j), node, adj.get(k)));
                                            sharedTriples2.get(edge2).add(new Triple(adj.get(j), node, adj.get(k)));

                                        }
                                    }
                                }
                            }
                        }
                    }

                    int F2 = 0;

                    for (Edge e : R.getEdges()) {
                        Node x = e.getNode1();
                        Node y = e.getNode2();

                        if (R.isDirectedFromTo(x, y) && G2.isDirectedFromTo(y, x)) {
                            F2 = F2 + 1;
                        }

                        if (R.isDirectedFromTo(y, x) && G2.isDirectedFromTo(x, y)) {
                            F2 = F2 + 1;
                        }
                    }

                    UtRStatistic utr = new UtRStatistic();
                    double r2 = utr.getValue(G2, R, data);

                    double f2 = utr.getF();

                    double d = _avgDegree / (double) (_numVars - 1);

                    if (A == 0) continue;
                    if (Double.isNaN(ahpc2)) continue;

                    System.out.println(
                            " d = " + getFormat(nf, d)
                                    + " Avg degree = " + _avgDegree
                                    + " num vars = " + _numVars
                                    + " R.numedges = " + R.getNumEdges()
                                    + " Ut = " + Ut
                                    + " A = " + A
                    );

                    out.println(
                            _avgDegree + "\t" + _numVars
                                    + "\t" + getFormat(nf, density)
                                    + "\t" + getFormat(nf, 1.0 - density)
                                    + "\t" + getFormat(nf, r2)
                                    + "\t" + getFormat(nf, ap2)
                                    + "\t" + getFormat(nf, ar2)
                                    + "\t" + getFormat(nf, ahp2)
                                    + "\t" + getFormat(nf, ahpc2)
                                    + "\t" + getFormat(nf, ahr2)
                                    + "\t" + getFormat(nf, ahrc2)
                                    + "\t" + nf.format(A)
                                    + "\t" + nf.format(Ut)
                                    + "\t" + nf.format(F2)
                                    + "\t" + nf.format(f2)
                                    + "\t" + nf.format(S2)
                                    + "\t" + nf.format(shared1)
                                    + "\t" + nf.format(shared2)
                                    + "\t" + elapsed
                    );
                }
            }

            out.close();
        }
    }

    private static String getFormat(NumberFormat nf, double x) {
        if (Double.isNaN(x)) return "-";
        if (Double.isInfinite(x)) return "-";
        return nf.format(x);
    }

    private static void scenario7() {
        String[] algorithms = {"sachs.model", "fask.5E-5", "fask.5E-2", "friedman", "aragam.discrete", "aragam.continuous",
                "henao", "desgranges", "goudet", "magliacane", "kalainathan", "fges", "pc", "cpc", "pcmax"};

        PrintStream out = null;

        Graph G2 = loadGraphTxt(new File("/Users/user/Box/data/Sachs/files.for.fask.sachs.report/txt/sachgroundtruth.txt"));

        makeBidirectedCycleUndirected(G2);


        try {
            out = new PrintStream(new File("/Users/user/Tetrad/tetrad-lib/src/main/" +
                    "java/edu/cmu/tetrad/study/calibration/data.for.calibration.sachs.txt"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        if (out == null) throw new NullPointerException("out not initialized");

        out.println("Alg & R2 & AHPC & A & S1 & r2max \\\\");

        for (String algorithm : algorithms) {

            System.out.println("\n============================");
            System.out.println("Algorithm = " + algorithm);
            System.out.println("============================\n");

            Graph R = loadGraphTxt(new File("/Users/user/Box/data/Sachs/files.for.fask.sachs.report/txt/" + algorithm + ".txt"));

            R = GraphUtils.replaceNodes(R, G2.getNodes());

            Statistic ahpc = new ArrowheadPrecisionCommonEdges();
            double ahpc2 = ahpc.getValue(G2, R, null);

            R = getCommonGraph(R, G2);

            makeBidirectedCycleUndirected(R);

            double numEdges = G2.getNumEdges();
            double numVars = G2.getNumNodes();
            double avgDegree = 2 * numEdges / numVars;

            double density = avgDegree / (numVars - 1);

            double d = 2 * numEdges / (numVars * (numVars - 1));

            List<Node> nodes = R.getNodes();

            int Ut = 0;

            for (int i = 0; i < nodes.size() - 1; i++) {
                List<Node> adj = R.getAdjacentNodes(nodes.get(i));

                for (int j = 0; j < adj.size(); j++) {
                    for (int k = j + 1; k < adj.size(); k++) {
                        if (!R.isAdjacentTo(adj.get(j), adj.get(k))) {
                            Ut++;
                        }
                    }
                }
            }

            int A = 0;

            for (Edge e2 : R.getEdges()) {
                if (e2.isDirected()) {
                    A++;
                }
            }

            UtRStatistic utr = new UtRStatistic();
            double r2 = utr.getValue(G2, R, null);

            double r2max = (1 - ahpc2) * (A / (2 * Ut * density));
            r2max = Math.min(r2max, 1.0);

            NumberFormat nf = new DecimalFormat("0.00");

            System.out.println(
                    " d = " + getFormat(nf, d)
                            + " Avg degree = " + avgDegree
                            + " num vars = " + numVars
                            + " R.numedges = " + R.getNumEdges()
                            + " A = " + A
            );

            out.println(
                    algorithm
                            + " & " + getFormat(nf, r2)
                            + " & " + getFormat(nf, ahpc2)
                            + " & " + nf.format(A)
                            + " & " + nf.format(Ut)
                            + " & " + getFormat(nf, r2max)
                            + " \\\\ "
            );
        }

        out.close();
    }

    private static void makeBidirectedCycleUndirected(Graph r) {
        for (Edge edge : r.getEdges()) {
            if (Edges.isBidirectedEdge(edge)) {
                r.removeEdge(edge);
                r.addUndirectedEdge(edge.getNode1(), edge.getNode2());
            }

            if (r.getEdges(edge.getNode1(), edge.getNode2()).size() > 1) {
                r.removeEdges(edge.getNode1(), edge.getNode2());
                r.addUndirectedEdge(edge.getNode1(), edge.getNode2());
            }
        }
    }

    // Returns g1 restricted to the adjacencies of g2.
    private static Graph getCommonGraph(Graph g1, Graph g2) {
        g1 = GraphUtils.replaceNodes(g1, g2.getNodes());
        Graph g1b = new EdgeListGraph(g2.getNodes());

        for (Edge e : g1.getEdges()) {
            if (g2.isAdjacentTo(e.getNode1(), e.getNode2())) {
                g1b.addEdge(e);
            }
        }

        return g1b;
    }

    private static void scenario8() throws IOException {
        int maxN = 100000;
        boolean useWeightsFromFile = false;
        int initialSegment = 100;
        boolean verbose = true;

        List<DataSet> dataSets = new ArrayList<>();

        NumberFormat nf = new DecimalFormat("0000");

        for (int i = 1; i <= initialSegment; i++) {
            File data = new File("/Users/user/Box/data/pairs 4/pair" + nf.format(i) + ".txt");
            System.out.println(data.getAbsolutePath());

            DataSet dataSet = loadContinuousData(data, false, Delimiter.WHITESPACE);
            writeDataSet(new File("/Users/user/Box/data/pairs/data"), i, dataSet);

            if (dataSet.getNumRows() > maxN) {
                dataSet = DataUtils.getBootstrapSample(dataSet, maxN);
            }

//            dataSet = DataUtils.standardizeData(dataSet);
//            if (dataSet.getNumRows() > maxN) dataSet = DataUtils.getBootstrapSample(dataSet, maxN);
            dataSets.add(dataSet);

//            writeResData(dataSet, i);
        }

        File gtFile = new File(new File("/Users/user/Box/data/pairs/"), "Readme3.txt");
        DataSet groundTruthData = loadDiscreteData(gtFile, false, Delimiter.TAB);


        List<Node> v = new ArrayList<>();
        v.add(new ContinuousVariable("Cutoff"));
        v.add(new ContinuousVariable("FPR"));
        v.add(new ContinuousVariable("TPR"));
        v.add(new ContinuousVariable("PREC"));
        v.add(new ContinuousVariable("NPV"));
        v.add(new ContinuousVariable("REC"));
        v.add(new ContinuousVariable("FDR"));
        v.add(new ContinuousVariable("FRAC"));
        v.add(new ContinuousVariable("ACC"));
        v.add(new ContinuousVariable("AUC"));
        v.add(new ContinuousVariable("TP"));
        v.add(new ContinuousVariable("FP"));
        v.add(new ContinuousVariable("FN"));
        v.add(new ContinuousVariable("TN"));

        int numRows = 20;

        TextTable tabulated = new TextTable(numRows + 2, v.size());

        for (int j = 0; j < v.size(); j++) {
            tabulated.setToken(0, j, v.get(j).getName());
        }

        for (int e = 0; e <= numRows; e++) {

            int[] discrete = {47, 70, 71, 85, 107};
            int[] vector = {52, 53, 54, 55, 71, 105};
            int[] interpolatedValues = {81, 82, 83};

            List<Set<Integer>> selected = new ArrayList<>();
            Set<Integer> omitted = new TreeSet<>();

            for (int i = 0; i < 2; i++) {
                selected.add(new TreeSet<>());
            }

            // Counts
            double tp = 0;
            double fp = 0;
            double tn = 0;
            double fn = 0;

            List<Integer> fps = new ArrayList<>();
            List<Integer> fns = new ArrayList<>();

            long start = System.currentTimeMillis();

            if (verbose) {
                System.out.println("i\tTrue\tEst");
            }

            double cutoff = e / 80.;

            List<Boolean> inCategory = new ArrayList<>();
            List<Double> scores = new ArrayList<>();

            for (int i = 1; i <= initialSegment; i++) {
                if (verbose) {
                    System.out.print(i);
                }

                if (Arrays.binarySearch(discrete, i) > -1) {
                    if (verbose) {
                        System.out.println(" DISCRETE");
                    }
                }

                if (Arrays.binarySearch(vector, i) > -1) {
                    if (verbose) {
                        System.out.println(" VECTOR");
                    }
                    continue;
                }

                DataSet dataSet = dataSets.get(i - 1);

                int x0 = 0;
                int y0 = 1;

                if (i == 52) {
                    x0 = 0;
                    y0 = 4;
                }

                if (i == 53) {
                    x0 = 0;
                    y0 = 3;
                }

                if (i == 54) {
                    x0 = 1;
                    y0 = 0;
                }

                if (i == 55) {
                    x0 = 0;
                    y0 = 16;
                }

                if (i == 71) {
                    x0 = 2;
                    y0 = 7;
                }

                if (i == 105) {
                    x0 = 1;
                    y0 = 9;
                }

                List<Node> gtNodes = groundTruthData.getVariables();
                DiscreteVariable c4 = (DiscreteVariable) gtNodes.get(4);
                DiscreteVariable c5 = (DiscreteVariable) gtNodes.get(5);

                String category = c4.getCategory(groundTruthData.getInt(i - 1, 4));

                double weight;

                if (useWeightsFromFile) {
                    weight = Double.parseDouble(c5.getCategory(groundTruthData.getInt(i - 1, 5)));
                } else {
                    weight = 1;
                }

                boolean groundTruthDirection = category.equals("->");

                if (!groundTruthDirection) {
                    int z0 = x0;
                    x0 = y0;
                    y0 = z0;

                    groundTruthDirection = !groundTruthDirection;
                }

                if (RandomUtil.getInstance().nextDouble() > 0.5) {
                    int z0 = x0;
                    x0 = y0;
                    y0 = z0;

                    groundTruthDirection = !groundTruthDirection;
                }

                int estLeftRight = 0;
                Graph g = new EdgeListGraph(dataSet.getVariables());
                List<Node> nodes = dataSet.getVariables();
                g.addUndirectedEdge(nodes.get(x0), nodes.get(y0));

                TetradLogger.getInstance().setLogging(false);

                Fask fask = new Fask(dataSet, g);
                fask.setTwoCycleThreshold(0);
                fask.setLinearityAssumed(true);
                fask.setUseFasAdjacencies(false);
                fask.setSkewEdgeThreshold(0.0);
                Graph out = fask.search();
                double lr = fask.getLr();

                if (abs(lr) < cutoff) {
                    omitted.add(i);
                }

                if (out.getEdges(nodes.get(x0), nodes.get(y0)).size() == 2) {
                    if (verbose) {
                        System.out.println(" 2-CYCLE");
                    }
//                    omitted.add(i);
                } else {
                    Edge edge = out.getEdge(nodes.get(x0), nodes.get(y0));
                    if (edge != null) {
                        boolean _estLeftRight = edge.pointsTowards(nodes.get(y0));
                        estLeftRight = _estLeftRight ? 1 : -1;
                    } else {
                        omitted.add(i);
                    }
                }

//                if (Arrays.binarySearch(discrete, i) > 0) {
//                    estLeftRight = -estLeftRight;
//                    lr = -lr;
//                }

                if (omitted.contains(i)) continue;

                // FASK assumptions hold.
                boolean isADog = groundTruthDirection;

                // FASK says left or right *correctly*.
                boolean thinkItsADog = (estLeftRight == 1);

                if (isADog && thinkItsADog) {
                    tp += weight;
                } else if (!isADog && thinkItsADog) {
                    fp += weight;
                    fps.add(i);
                } else if (isADog && !thinkItsADog) {
                    fn += weight;
                    fns.add(i);
                } else if (!isADog && !thinkItsADog) {
                    tn += weight;
                }

                inCategory.add(isADog);
                scores.add(lr);


                if (verbose) {
                    if (groundTruthDirection) System.out.print("\t-->");
                    else System.out.print("\t<--");
                }

                if (verbose) {
                    if (estLeftRight == 1) System.out.print("\t-->");
                    else if (estLeftRight == -1) System.out.print("\t<--");
                    else System.out.print("\t");
                }

                if ((isADog && thinkItsADog) || (!isADog && !thinkItsADog)) {
                    selected.get(0).add(i);
                } else if ((isADog && !thinkItsADog) || (!isADog && thinkItsADog)) {
                    selected.get(1).add(i);
                }

                if (verbose) {
                    System.out.println();
                }
            }

            long stop = System.currentTimeMillis();

            NumberFormat nf2 = new DecimalFormat("0.00");

            double tpr = tp / (tp + fn);
            double fpr = fp / (fp + tn);
            double precision = tp / (tp + fp);
            double npv = tn / (tn + fn);
            double recall = tp / (tp + fn);
            double fdr = fp / (tp + fp); // false positives over positives
            double acc = (tp + tn) / (tp + fp + tn + fn);
            double fracDec = (initialSegment - omitted.size()) / (double) initialSegment;

            boolean[] _inCategory = new boolean[inCategory.size()];
            for (int i = 0; i < inCategory.size(); i++) _inCategory[i] = inCategory.get(i);

            double[] _scores = new double[scores.size()];
            for (int i = 0; i < scores.size(); i++) _scores[i] = scores.get(i);

            RocCalculator roc = new RocCalculator(_scores, _inCategory, RocCalculator.ASCENDING);
            double auc = roc.getAuc();

            System.out.println("\nSummary:\n");
            System.out.println("Accuracy = " + nf2.format(acc));
            System.out.println("Total correct = " + selected.get(0));
            System.out.println("Total incorrect: " + selected.get(1));
            System.out.println("TPR: " + tpr);
            System.out.println("FPR: " + fpr);
            System.out.println("PREC: " + precision);
            System.out.println("NPV: " + npv);
            System.out.println("REC: " + recall);
            System.out.println("Didn't classify: " + omitted);
            System.out.println("Fraction of Decisions: " + fracDec);
            System.out.println("Elapsed time = " + ((stop - start) / (double) 1000) + "s");

            System.out.println("fps = " + fps);
            System.out.println("fns = " + fns);

            int l = 0;

            if (verbose) {
                System.out.println();
                System.out.println("PER RUN:");
                System.out.println();
                System.out.println("RELEVANT CATEGORY = --> IN GROUND TRUTH");
                System.out.println("SELECTED CATEGORY = RESIDUAL FASK JUDGED -->");

                if (useWeightsFromFile) {
                    System.out.println("USING WEIGHTS FOR CASES FROM TABLE 4 IN MOOIJ ET AL.");
                }

                System.out.println();
                System.out.println("TPR = TP / (TP + FN) = RECALL");
                System.out.println("FPR = FP / (FP + TN)");
                System.out.println("PREC = TP / (TP + FP)");
                System.out.println("NPV = TN / (TN + FN)");
                System.out.println("REC = TP / (TP + FN) = TPR");
                System.out.println("FDR = FP / (TP + FP)");
                System.out.println("ACC = (TP + TN) / (TP + FP + TN + FN)");
                System.out.println("AUC = AREA UNDER ROC CURVE");
            }

            NumberFormat nf3 = new DecimalFormat("0.00");
            NumberFormat nf4 = new DecimalFormat("0.0");

            tabulated.setToken(e + 1, l++, "" + nf3.format(cutoff));
            tabulated.setToken(e + 1, l++, "" + nf3.format(fpr));
            tabulated.setToken(e + 1, l++, "" + nf3.format(tpr));
            tabulated.setToken(e + 1, l++, "" + nf3.format(precision));
            tabulated.setToken(e + 1, l++, "" + nf3.format(npv));
            tabulated.setToken(e + 1, l++, "" + nf3.format(recall));
            tabulated.setToken(e + 1, l++, "" + nf3.format(fdr));
            tabulated.setToken(e + 1, l++, "" + nf3.format(fracDec));
            tabulated.setToken(e + 1, l++, "" + nf3.format(acc));
            tabulated.setToken(e + 1, l++, "" + nf3.format(auc));
            tabulated.setToken(e + 1, l++, "" + nf4.format(tp));
            tabulated.setToken(e + 1, l++, "" + nf4.format(fp));
            tabulated.setToken(e + 1, l++, "" + nf4.format(fn));
            tabulated.setToken(e + 1, l++, "" + nf4.format(tn));
        }

        System.out.println("\n" + tabulated);
    }

    private static DataSet loadContinuousData(File data, boolean hasHeader, Delimiter delimiter) throws IOException {
        try {
            ContinuousTabularDatasetFileReader dataReader
                    = new ContinuousTabularDatasetFileReader(data.toPath(), delimiter);
            dataReader.setHasHeader(hasHeader);
            dataReader.setMissingDataMarker("*");

            Data _data = dataReader.readInData();
            return (DataSet) DataConvertUtils.toDataModel(_data);
        } catch (IOException e) {
            e.printStackTrace();
            throw e;
        }

    }

    private static DataSet loadDiscreteData(File data, boolean hasHeader, Delimiter delimiter) throws IOException {
        try {
            VerticalDiscreteTabularDatasetFileReader dataReader
                    = new VerticalDiscreteTabularDatasetFileReader(data.toPath(), delimiter);
            dataReader.setHasHeader(hasHeader);
            dataReader.setMissingDataMarker("*");
            dataReader.setQuoteCharacter('"');

            Data _data = dataReader.readInData();
            return (DataSet) DataConvertUtils.toDataModel(_data);
        } catch (IOException e) {
            e.printStackTrace();
            throw e;
        }

    }

    private static void writeDataSet(File dir, int i, DataSet d2) {
        try {
            DataWriter.writeRectangularData(d2, new PrintWriter(
                    new File(dir, "pair2." + i + ".txt")), '\t');
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void collectUnshieldedTripleLegsAndShieldsInR(Graph R, Set<Edge> L, Set<Edge> M, Node
            v1, Node v2, Node v3) {
        if (R.isAdjacentTo(v1, v2) && R.isAdjacentTo(v2, v3) && !R.isAdjacentTo(v1, v3)) {
            M.add(Edges.undirectedEdge(v1, v3));
            L.add(Edges.undirectedEdge(v2, v1));
            L.add(Edges.undirectedEdge(v2, v3));
        }
    }

    public static DataSet filter(DataSet data) {
        List<Node> variables = data.getVariables();
        int numRows = 0;

        ROWS:
        for (int row = 0; row < data.getNumRows(); row++) {
            for (int col = 0; col < data.getNumColumns(); col++) {
                Node variable = data.getVariable(col);
                if (((Variable) variable).isMissingValue(data.getObject(row, col))) {
                    continue ROWS;
                }
            }

            numRows++;
        }

        DataSet newDataSet = new BoxDataSet(new DoubleDataBox(numRows, variables.size()), variables);
        int newRow = 0;

        ROWS:
        for (int row = 0; row < data.getNumRows(); row++) {
            for (int col = 0; col < data.getNumColumns(); col++) {
                Node variable = data.getVariable(col);
                if (((Variable) variable).isMissingValue(data.getObject(row, col))) {
                    continue ROWS;
                }
            }

            for (int col = 0; col < data.getNumColumns(); col++) {
                newDataSet.setObject(newRow, col, data.getObject(row, col));
            }

            newRow++;
        }

        return newDataSet;
    }

    private static void writeResData(DataSet dataSet, int i) {
        double[] x = dataSet.copy().getDoubleData().transpose().toArray()[0];
        double[] y = dataSet.copy().getDoubleData().transpose().toArray()[1];

        double[] r = Fask.residuals(y, x, Fask.RegressionType.NONLINEAR);

//        for (int j = 0; j < x.length; j++) x[j] -= r[j];

//        RegressionResult result = RegressionDataset.regress(y, new double[][]{x});
//        r = result.getResiduals().toArray();

        double[][] res = new double[][]{x, r};


        List<Node> vars = new ArrayList<>();
        vars.add(new ContinuousVariable("C1"));
        vars.add(new ContinuousVariable("C2"));

        DataSet dataSet1 = new BoxDataSet(new VerticalDoubleDataBox(res), vars);

        writeDataSet(new File("/Users/user/Box/data/pairs/resdata"), i, dataSet1);
    }

    private static double cxy(double[] x, double[] y) {
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < x.length; i++) indices.add(i);

        Collections.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(x[o1], x[o2]);
            }
        });

        double[] xSorted = new double[x.length];
        double[] ySorted = new double[y.length];

        for (int i = 0; i < x.length; i++) {
            xSorted[i] = x[indices.get(i)];
            ySorted[i] = y[indices.get(i)];
        }


        double sum = 0.0;

        int N = xSorted.length;

        for (int j = 0; j < N - 1; j++) {
            sum += abs(ySorted[j + 1] - ySorted[j]) / (xSorted[j + 1] - xSorted[j]);
        }

        return sum / (N - 1);
    }

    public static void scenario9() {

        Parameters parameters = new Parameters();

        parameters.set(Params.NUM_MEASURES, 50);
        parameters.set(Params.AVG_DEGREE, 4);
        parameters.set(Params.NUM_RUNS, 10);
//        parameters.set(Params.DEPTH, 3);
        parameters.set(Params.ERRORS_NORMAL, true, false);
        parameters.set(Params.SKEW_EDGE_THRESHOLD, 0.3);
        parameters.set(Params.TWO_CYCLE_THRESHOLD, 0);
        parameters.set(Params.COLLIDER_DISCOVERY_RULE, 1, 2, 3);
        parameters.set(Params.FASK_LINEARITY_ASSUMED, false);
        parameters.set(Params.ALPHA, 0.001);
        parameters.set(Params.SAMPLE_SIZE, 1000);

        parameters.set(Params.COEF_LOW, 0.2);
        parameters.set(Params.COEF_HIGH, 1.5);
        parameters.set(Params.VAR_LOW, 1);
        parameters.set(Params.VAR_HIGH, 2);


        parameters.set(Params.INTERVAL_BETWEEN_SHOCKS, 100);
        parameters.set(Params.INTERVAL_BETWEEN_RECORDINGS, 100);


        Algorithms algorithms = new Algorithms();
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Pcp(new FisherZ()));
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll(new FisherZ()));
//        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.multi.Fask(new FisherZ()));

        Simulations simulations = new Simulations();
        simulations.add(new SemSimulation(new RandomForward()));
        simulations.add(new SemSimulation(new RandomForward()));

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn("errorsNormal"));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
        comparison.setSaveSemModels(true);

        comparison.compareFromSimulations("comparison2", simulations, algorithms, statistics, parameters);


    }

}
