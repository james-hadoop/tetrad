package edu.cmu.tetrad.study.calibration;

import edu.cmu.tetrad.algcomparison.algorithm.pairwise.Skew;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.sem.LargeScaleSimulation;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.reader.Data;
import edu.pitt.dbmi.data.reader.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDatasetFileReader;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import static edu.cmu.tetrad.graph.GraphUtils.loadGraphTxt;

public class CalibrationQuestion {

    public static void main(String... args) {
        scenario8();
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

            if (ge.isAdjacentTo(x, y) && ge.isAdjacentTo(y, z)//) {// && gt.isAdjacentTo(x, y) && gt.isAdjacentTo(y, z)) {
//                    && !ge.isAdjacentTo(x, z)
            ) {
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

//            if (gt.isAdjacentTo(x, y) && gt.isAdjacentTo(y, z)) {

            SemPm pm = new SemPm(gt);
            SemIm im = new SemIm(pm);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges fges = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
//            edu.cmu.tetrad.search.Pc fges = new edu.cmu.tetrad.search.Pc(new edu.cmu.tetrad.search.IndTestFisherZ(data, 0.1));

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
//        }

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

//            if (!gt.isAdjacentTo(x, y)) gt.addDirectedEdge(x, y);
//            if (!gt.isAdjacentTo(y, z)) gt.addDirectedEdge(y, z);
//            if (!gt.isAdjacentTo(x, z)) gt.addDirectedEdge(x, z);

//            gt.removeEdge(x, z);


//            if (gt.isAdjacentTo(x, y) && gt.isAdjacentTo(y, z)) {

            SemPm pm = new SemPm(gt);

            Parameters parameters = new Parameters();
            parameters.set("coefLow", 0.1);
            parameters.set("coefHigh", 0.5);

            SemIm im = new SemIm(pm, parameters);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
//                edu.cmu.tetrad.search.Pc s = new edu.cmu.tetrad.search.Pc(new edu.cmu.tetrad.search.IndTestFisherZ(data, 0.1));

            Graph ge = s.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

//            if (!ge.isAdjacentTo(x, y)) ge.addDirectedEdge(x, y);
//            if (!ge.isAdjacentTo(y, z)) ge.addDirectedEdge(y, z);


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

//            System.out.println("Run " + i + " p = " + p + " c = " + c + " totalc = " + totalc + " d = " + d
//                    + " totald = " + totald + " qc = " + (c / (double) totalc) + " qd = " + (d / (double) totald));


            System.out.println("Run " + i + " P(XYe & YZe) = " + (c / (double) totalc) + " P(XYe & YZe | XZt) = " + (d / (double) totald));
        }
//        }

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

//            if (!gt.isAdjacentTo(x, y)) gt.addDirectedEdge(x, y);
//            if (!gt.isAdjacentTo(y, z)) gt.addDirectedEdge(y, z);
//            if (!gt.isAdjacentTo(x, z)) gt.addDirectedEdge(x, z);

//            gt.removeEdge(x, z);


//            if (gt.isAdjacentTo(x, y) && gt.isAdjacentTo(y, z)) {

            SemPm pm = new SemPm(gt);

            Parameters parameters = new Parameters();
            parameters.set("coefLow", 0.1);
            parameters.set("coefHigh", 0.5);

            SemIm im = new SemIm(pm);

            DataSet data = im.simulateData(sampleSize, false);

            edu.cmu.tetrad.search.Fges s = new edu.cmu.tetrad.search.Fges(new edu.cmu.tetrad.search.SemBicScore(data));
//                edu.cmu.tetrad.search.Pc s = new edu.cmu.tetrad.search.Pc(new edu.cmu.tetrad.search.IndTestFisherZ(data, 0.1));

            Graph ge = s.search();
            ge = GraphUtils.replaceNodes(ge, gt.getNodes());

//            if (!ge.isAdjacentTo(x, y)) ge.addDirectedEdge(x, y);
//            if (!ge.isAdjacentTo(y, z)) ge.addDirectedEdge(y, z);

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

                    DataSet data = lgs.simulateDataFisher(100, 100, sampleSize, 1e-3, false);

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
                            ((Fask) s).setUseSkewAdjacencies(true);
                            ((Fask) s).setExtraEdgeThreshold(0.5);
                            ((Fask) s).setAlpha(0.2);
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

//        Graph G2 = loadGraphTxt(new File("/Users/user/Box/data/Sachs/files.for.fask.sachs.report/txt/ground.truth.sachs.txt"));
        Graph G2 = loadGraphTxt(new File("/Users/user/Box/data/Sachs/files.for.fask.sachs.report/txt/sachgroundtruth.txt"));
//        Graph G2 = loadGraphTxt(new File("/Users/user/Box/data/Sachs/peter.ground.truth2.txt"));

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

            ChoiceGenerator gen = new ChoiceGenerator(nodes.size(), 3);
            int[] choice;

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

            int S2 = 0;

            for (int i = 0; i < nodes.size() - 1; i++) {
                List<Node> adj = R.getAdjacentNodes(nodes.get(i));

                for (int j = 0; j < adj.size(); j++) {
                    for (int k = j + 1; k < adj.size(); k++) {
                        if (!R.isAdjacentTo(adj.get(j), adj.get(k))) {
                            if (G2.isAdjacentTo(nodes.get(i), adj.get(j)) && G2.isAdjacentTo(nodes.get(i), adj.get(k))) {
                                if (G2.isAdjacentTo(adj.get(j), adj.get(k))) {
                                    S2++;
                                }
                            }
                        }
                    }
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

    private static void scenario8() {
        File dir = new File("/Users/user/Box/data/pairs");

        NumberFormat nf = new DecimalFormat("0000");

        int correct = 0;
        int total = 0;

        int skippedWeakSkewness = 0;
        int skippedSingular = 0;
        int skippedIndependent = 0;
        int excluded = 0;

        for (int i = 1; i <= 108; i++) {


//            File data = new File(dir, "pair" + nf.format(i) + ".txt");
            File data = new File(dir, "pair." + i + ".resave.txt");
            File des = new File(dir, "pair" + nf.format(i) + "_des.txt");

            File pairsMetaFile = new File(dir, "pairmeta.txt");

            ContinuousTabularDatasetFileReader dataReader
                    = new ContinuousTabularDatasetFileReader(data.toPath(), Delimiter.TAB);
            dataReader.setHasHeader(true);
            dataReader.setMissingDataMarker("*");

            ContinuousTabularDatasetFileReader dataReader2
                    = new ContinuousTabularDatasetFileReader(pairsMetaFile.toPath(), Delimiter.WHITESPACE);
            dataReader2.setHasHeader(false);
            dataReader2.setMissingDataMarker("*");


            DataSet dataSet = null;
            DataSet pairsMeta = null;

            try {
                Data _data = dataReader.readInData();
                dataSet = (DataSet) DataConvertUtils.toDataModel(_data);

                if (dataSet.getNumColumns() > 2 && Double.isNaN(dataSet.getDouble(0, dataSet.getNumColumns() - 1))) {
                    dataSet.removeColumn(dataSet.getNumColumns() - 1);
                }

//                DataWriter.writeRectangularData(dataSet, new PrintWriter(
//                        new File(dir,"pair." + i + ".resave.txt")), '\t');

                Data pairsMetaData = dataReader2.readInData();
                pairsMeta = (DataSet) DataConvertUtils.toDataModel(pairsMetaData);

//                System.out.println(pairsMeta);
//                System.out.println();
            } catch (IOException e) {
                e.printStackTrace();
            }

            assert dataSet != null;
            TetradMatrix doubles = dataSet.getDoubleData();

            double[] data0 = doubles.getColumn(0).toArray();
            double[] data1 = doubles.getColumn(1).toArray();

            double sk0 = StatUtils.skewness(data0);
            double sk1 = StatUtils.skewness(data1);

//            if (abs(sk0) < 0.0001 || abs(sk1) < 0.0001) {
//                System.out.println("Skipping " + i + " because at least one skewness is too weak.");
//                System.out.println("sk0 = " + sk0);
//                System.out.println("sk1 = " + sk1);
////                skippedWeakSkewness++;
////                continue;
//            }

//            for (int m = 0; m < dataSet.getNumRows(); m++) {
//                if (sk0 < 0) {
//                    dataSet.setDouble(m, 0, -dataSet.getDouble(m, 0));
//                }
//
//                if (sk1 < 0) {
//                    dataSet.setDouble(m, 1, -dataSet.getDouble(m, 1));
//                }
//            }

            boolean trueLeftRight = false;
            boolean trueRightLeft = false;
            boolean trueConfounded = false;
            boolean dependentError = false;
            boolean notTwoCycle = false;
            boolean symmetricLeftRight = false;
            boolean exclude = false;

            try {
                BufferedReader reader = new BufferedReader(new FileReader(des));
                String line;

                while ((line = reader.readLine()) != null) {
                    if (line.contains("#-->")) {
                        trueLeftRight = true;
                        break;
                    }
                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("#<--")) {
                        trueRightLeft = true;
                        break;
                    }
                }

                assert pairsMeta != null;

                int m = i - 1;

//                trueLeftRight = pairsMeta.getDouble(m, 1) == 1;
//                trueRightLeft = !trueLeftRight;

//                if (trueLeftRight != trueLeftRight2) {
//                    System.out.println("Discrepancy in " + i);
//                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("#<T>")) {
                        trueConfounded = true;
                        break;
                    }
                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("DependentError")) {
                        dependentError = true;
                        break;
                    }
                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("SymmetricLeftRight")) {
                        symmetricLeftRight = true;
                        break;
                    }
                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("NotTwoCycle")) {
                        notTwoCycle = true;
                        break;
                    }
                }

                reader = new BufferedReader(new FileReader(des));
                while ((line = reader.readLine()) != null) {
                    if (line.contains("Exclude")) {
                        exclude = true;
                        break;
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

//            if (exclude) {
//                excluded++;
//                continue;
//            }

            List<Node> nodes = new ArrayList<>();
            List<Node> variables = dataSet.getVariables();
            nodes.add(variables.get(0));
            nodes.add(variables.get(1));

            Graph g = new EdgeListGraph(dataSet.getVariables());
            g.addUndirectedEdge(nodes.get(0), nodes.get(1));


            List<DataSet> dataSets = new ArrayList<>();
            dataSets.add(DataUtils.getContinuousDataSet(dataSet));

            Lofs2 lofs = new Lofs2(g, dataSets);
            lofs.setRule(Lofs2.Rule.Skew);

            Graph out;


//            out = lofs.orient();

            Fask fask = new Fask(dataSet, g);
//            Fask fask = new Fask(dataSet, new IndTestFisherZ(dataSet, 0.001));
            fask.setAlpha(.00);

//            if (dependentError || notTwoCycle || symmetricLeftRight) {
//                fask.setAlpha(0);
//            }


//            try {
                out = fask.search();
                if (!out.isAdjacentTo(nodes.get(0), nodes.get(1))) {
                    System.out.println("Skipping " + i + " because the variables are judged to be independent.");
                    skippedIndependent++;
                    continue;
                }
//            } catch (Exception e) {
//                System.out.println("Rerunning " + i + " without 2-cycle check because of singularity.");
//
//                fask.setAlpha(0);
//
//                out = fask.search();
//                if (!out.isAdjacentTo(nodes.get(0), nodes.get(1))) {
//                    System.out.println("Skipping " + i + " because the variables are judged to be independent.");
//                    skippedIndependent++;
//                    continue;
//                }
//            }

            boolean estLeftRight = out.containsEdge(Edges.directedEdge(nodes.get(0), nodes.get(1)));
            boolean estRightLeft = out.containsEdge(Edges.directedEdge(nodes.get(1), nodes.get(0)));

            if (dependentError) {
                estLeftRight = !estLeftRight;
                estRightLeft = !estRightLeft;
//                continue;
            }

            if (symmetricLeftRight) {
                estLeftRight = !estLeftRight;
                estRightLeft = !estRightLeft;
            }

            if ((trueLeftRight && estLeftRight) || (trueRightLeft && estRightLeft)) continue;
//
//            if (!(estLeftRight && estRightLeft)) continue;

            System.out.println("\nindex = " + i);

            if (trueLeftRight) System.out.println("true -->");
            if (trueRightLeft) System.out.println("true <--");
            if (trueConfounded) System.out.println("True confounded");
            if (dependentError) System.out.println("Dependent error");
            if (symmetricLeftRight) System.out.println("Symmetric Left-Right");
            if (notTwoCycle) System.out.println("Not a Two Cycle");

            if (estLeftRight) System.out.println("est -->");
            if (estRightLeft) System.out.println("est <--");

            if (trueLeftRight && estLeftRight) {
                System.out.println("Correct left-right");
                correct++;
            }
            if (trueRightLeft && estRightLeft) {
                System.out.println("Correct right-left");
                correct++;
            }

            total++;

            System.out.println("Running total: correct = " + correct + " total = " + total);
        }

        NumberFormat nf2 = new DecimalFormat("0.00");

        System.out.println("\nSummary:\n");
        System.out.println("Score = " + nf2.format((correct / (double) total)));
        System.out.println("Total correct = " + correct);
        System.out.println("Total not skipped = " + total);
        System.out.println("Skipped because of weak skewnesses = " + skippedWeakSkewness);
        System.out.println("Skipped because of singularity = " + skippedSingular);
        System.out.println("Skipped because variables were judged independent = " + skippedIndependent);
        System.out.println("Skipped because explicitly excluded = " + excluded);
    }

    private static void collectUnshieldedTripleLegsAndShieldsInR(Graph R, Set<Edge> L, Set<Edge> M, Node
            v1, Node v2, Node v3) {
        if (R.isAdjacentTo(v1, v2) && R.isAdjacentTo(v2, v3) && !R.isAdjacentTo(v1, v3)) {
            M.add(Edges.undirectedEdge(v1, v3));
            L.add(Edges.undirectedEdge(v2, v1));
            L.add(Edges.undirectedEdge(v2, v3));
        }
    }
}
