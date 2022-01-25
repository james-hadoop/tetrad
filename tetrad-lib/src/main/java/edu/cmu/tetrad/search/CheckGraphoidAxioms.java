package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.PermutationGenerator;

import java.util.*;

/**
 * Checks the graphoid axioms for a set of Independence Model statements.
 *
 * @author josephramsey
 */
public class CheckGraphoidAxioms {
    static final long serialVersionUID = 23L;
    private final Set<GraphoidIndFact> facts;

    public CheckGraphoidAxioms(Set<GraphoidIndFact> facts) {
        this.facts = facts;
    }

    /**
     * If you have a list of IndependenceFacts, this will return the corresponding
     * set of GraphoidIndFacts.
     *
     * @param facts The list of IndependenceFacts.
     * @return The set of GraphoidIndFacts.
     */
    public Set<GraphoidIndFact> getGraphoidFacts(List<IndependenceFact> facts) {
        Set<GraphoidIndFact> graphoidIndFacts = new HashSet<>();

        for (IndependenceFact fact : facts) {
            List<Node> XX = Collections.singletonList(fact.getX());
            List<Node> YY = Collections.singletonList(fact.getY());
            List<Node> ZZ = new ArrayList<>(fact.getZ());
            graphoidIndFacts.add(new GraphoidIndFact(XX, YY, ZZ));
        }

        return graphoidIndFacts;
    }

    public boolean semigraphoid() {
        return symmetry() && decomposition() && weakUnion() && contraction();
    }

    public boolean graphoid() {
        return semigraphoid() && intersection();
    }

    public boolean compositionalGraphoid() {
        return graphoid() && composition();
    }

    // XX ⊥⊥ YY | Z ==> YY ⊥⊥ XX | ZZ
    private boolean symmetry() {
        for (GraphoidIndFact fact : facts) {
            if (!facts.contains(new GraphoidIndFact(fact.getYY(), fact.getXX(), fact.getZZ()))) {
                return false;
            }
        }

        return true;
    }

    // XX ⊥⊥ Y ∪ W |ZZ ==> (XX ⊥⊥ Y |ZZ) ∧ (XX ⊥⊥ W |ZZ)
    private boolean decomposition() {
        for (GraphoidIndFact fact : facts) {
            List<Node> XX = fact.getXX();
            List<Node> YY = fact.getYY();
            List<Node> ZZ = fact.getZZ();

            PermutationGenerator gen = new PermutationGenerator(YY.size());
            int[] perm;

            while ((perm = gen.next()) != null) {
                List<Node> Y = GraphUtils.asList(perm, YY);
                List<Node> W = new ArrayList<>(YY);
                W.removeAll(Y);

                if (Y.isEmpty() || W.isEmpty()) continue;

                if (!(facts.contains(new GraphoidIndFact(XX, Y, ZZ))
                        && facts.contains(new GraphoidIndFact(XX, W, ZZ)))) {
                    return false;
                }
            }
        }

        return true;
    }

    // XX _||_ Y U W | ZZ ==> X _||_ Y | ZZ U W
    private boolean weakUnion() {
        for (GraphoidIndFact fact : facts) {
            List<Node> XX = fact.getXX();
            List<Node> YY = fact.getYY();
            List<Node> ZZ = fact.getZZ();

            PermutationGenerator gen = new PermutationGenerator(YY.size());
            int[] perm;

            while ((perm = gen.next()) != null) {
                List<Node> Y = GraphUtils.asList(perm, YY);
                List<Node> W = new ArrayList<>(YY);
                W.removeAll(Y);

                if (Y.isEmpty() || W.isEmpty()) continue;

                List<Node> Z = new ArrayList<>(ZZ);
                Z.addAll(W);

                if (!(facts.contains(new GraphoidIndFact(XX, Y, Z)))) {
                    return false;
                }
            }
        }

        return true;
    }

    // (XX ⊥⊥ YY |ZZ) ∧(XX ⊥⊥ W |ZZ ∪ YY) ==> XX ⊥⊥ YY ∪ W |ZZ
    private boolean contraction() {
        for (GraphoidIndFact fact : facts) {
            List<Node> XX = fact.getXX();
            List<Node> YY = fact.getYY();
            List<Node> ZZ = fact.getZZ();

            List<Node> Z = new ArrayList<>(ZZ);
            Z.addAll(YY);

            for (GraphoidIndFact fact2 : facts) {
                if (fact2.getZZ().equals(Z)) {
                    List<Node> W = fact.getYY();

                    List<Node> Y = new ArrayList<>(YY);
                    Y.addAll(W);

                    if (!facts.contains(new GraphoidIndFact(XX, Y, ZZ))) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    // (XX ⊥⊥ YY |Z ∪ W) ∧(XX ⊥⊥ W |Z ∪ YY) ==> XX ⊥⊥ YY ∪ W |Z
    private boolean intersection() {
        for (GraphoidIndFact fact : facts) {
            List<Node> XX = fact.getXX();
            List<Node> YY = fact.getYY();
            List<Node> ZZ = fact.getZZ();

            PermutationGenerator gen = new PermutationGenerator(YY.size());
            int[] perm;

            while ((perm = gen.next()) != null) {
                List<Node> Z = GraphUtils.asList(perm, ZZ);
                List<Node> W = new ArrayList<>(YY);
                W.removeAll(Z);

                if (Z.isEmpty() || W.isEmpty()) continue;

                List<Node> Z2 = new ArrayList<>(Z);
                Z2.addAll(YY);

                if ((facts.contains(new GraphoidIndFact(XX, W, Z2)))) {
                    List<Node> Y2 = new ArrayList<>(YY);
                    YY.addAll(W);

                    if (!fact.getXX().contains(new GraphoidIndFact(XX, Y2, Z))) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    // (X ⊥⊥ YY | ZZ) ∧ (X ⊥⊥ W |ZZ) ==> X ⊥⊥ YY ∪ W |ZZ
    private boolean composition() {
        for (GraphoidIndFact fact : facts) {
            List<Node> XX = fact.getXX();
            List<Node> YY = fact.getYY();
            List<Node> ZZ = fact.getZZ();

            for (GraphoidIndFact fact2 : facts) {
                if (fact2.getXX().equals(XX) && fact2.getZZ().equals(ZZ)) {
                    List<Node> W = fact2.YY;

                    List<Node> Y = new ArrayList<>(YY);
                    Y.addAll(W);

                    if (!facts.contains(new GraphoidIndFact(XX, Y, ZZ))) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public static class GraphoidIndFact {
        private final List<Node> XX;
        private final List<Node> YY;
        private final List<Node> ZZ;

        public GraphoidIndFact(List<Node> XX, List<Node> YY, List<Node> ZZ) {
            this.XX = XX;
            this.YY = YY;
            this.ZZ = ZZ;
        }

        public List<Node> getXX() {
            return XX;
        }

        public List<Node> getYY() {
            return YY;
        }

        public List<Node> getZZ() {
            return ZZ;
        }
    }
}
