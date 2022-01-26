package edu.cmu.tetrad.search;

import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;

/**
 * Checks the graphoid axioms for a set of Independence Model statements.
 *
 * @author josephramsey
 */
public class GraphoidAxioms {
    private final Set<GraphoidIndFact> facts;
    private boolean trivialtyAssumed = false;

    /**
     * Constructor.
     *
     * @param facts A set of GraphoidIdFacts.
     */
    public GraphoidAxioms(Set<GraphoidIndFact> facts) {
        this.facts = new LinkedHashSet<>(facts);
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

    /**
     * X ⊥⊥ Y | Z ==> Y ⊥⊥ X | Z
     */
    public boolean symmetry() {
        for (GraphoidIndFact fact : facts) {
            GraphoidIndFact newFact = new GraphoidIndFact(fact.getY(), fact.getX(), fact.getZ());

            if (!facts.contains(newFact)) {
                TetradLogger.getInstance().forceLogMessage("Symmetry fails: Missing "
                        + newFact);
                return false;
            }
        }

        return true;
    }

    /**
     * X ⊥⊥ (Y ∪ W) |Z ==> (X ⊥⊥ Y |Z) ∧ (X ⊥⊥ W |Z)
     */
    public boolean decomposition() {
        for (GraphoidIndFact fact : facts) {
            Set<Node> X = fact.getX();
            Set<Node> YW = fact.getY();
            Set<Node> Z = fact.getZ();

            List<Node> YWList = new ArrayList<>(YW);

            DepthChoiceGenerator gen = new DepthChoiceGenerator(YWList.size(), YWList.size());
            int[] choice;

            while ((choice = gen.next()) != null) {
                Set<Node> Y = GraphUtils.asSet(choice, YWList);
                Set<Node> W = new HashSet<>(YW);
                W.removeAll(Y);

                if (trivialtyAssumed && X.isEmpty() || Y.isEmpty() || W.isEmpty()) continue;

                GraphoidIndFact fact1 = new GraphoidIndFact(X, Y, Z);
                GraphoidIndFact fact2 = new GraphoidIndFact(X, W, Z);

                if (!(facts.contains(fact1))) {
                    TetradLogger.getInstance().forceLogMessage("Decomposition fails:" +
                            " Have " + fact +
                            "; Missing " + fact1);
                    return false;
                }

                if (!(facts.contains(fact2))) {
                    TetradLogger.getInstance().forceLogMessage("Decomposition fails:" +
                            " Have " + fact +
                            "; Missing " + fact2);
                    return false;
                }
            }
        }

        return true;
    }


    /**
     * X _||_ Y U W | Z ==> X _||_ Y | Z U W
     */
    public boolean weakUnion() {
        for (GraphoidIndFact fact : facts) {
            Set<Node> X = fact.getX();
            Set<Node> YW = fact.getY();
            Set<Node> Z = fact.getZ();

            List<Node> YWList = new ArrayList<>(YW);

            DepthChoiceGenerator gen = new DepthChoiceGenerator(YW.size(), YW.size());
            int[] choice;

            while ((choice = gen.next()) != null) {
                Set<Node> Y = GraphUtils.asSet(choice, YWList);
                Set<Node> W = new HashSet<>(YW);
                W.removeAll(Y);

                Set<Node> ZW = new HashSet<>(Z);
                ZW.addAll(W);

                if (trivialtyAssumed && X.isEmpty() || Y.isEmpty() || W.isEmpty()) continue;

                GraphoidIndFact newFact = new GraphoidIndFact(X, Y, ZW);

                if (!(facts.contains(newFact))) {
                    TetradLogger.getInstance().forceLogMessage("Weak Union fails:" +
                            " Have " + fact +
                            "; Missing " + newFact);
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * (X ⊥⊥ Y |Z) ∧ (X ⊥⊥ W |Z ∪ Y) ==> X ⊥⊥ (Y ∪ W) |Z
     */
    public boolean contraction() {
        for (GraphoidIndFact fact1 : new HashSet<>(facts)) {
            Set<Node> X = fact1.getX();
            Set<Node> Y = fact1.getY();
            Set<Node> Z = fact1.getZ();

            Set<Node> ZY = new HashSet<>(Z);
            ZY.addAll(Y);

            for (GraphoidIndFact fact2 : new HashSet<>(facts)) {
                if (fact1.equals(fact2)) continue;

                if (fact2.getX().equals(X) && fact2.getZ().equals(ZY)) {
                    Set<Node> W = fact2.getY();
                    if (X.equals(W)) continue;
                    if (X.equals(ZY)) continue;
                    if (W.equals(ZY)) continue;

                    Set<Node> YW = new HashSet<>(Y);
                    YW.addAll(W);

                    GraphoidIndFact newFact = new GraphoidIndFact(X, YW, Z);

                    if (!facts.contains(newFact)) {
                        TetradLogger.getInstance().forceLogMessage("Contraction fails:" +
                                " Have " + fact1 + " and " + fact2 +
                                "; Missing " + newFact);
                        return false;
                    }
                }
            }
        }

        return true;
    }


    /**
     * (X ⊥⊥ Y | (Z ∪ W)) ∧ (X ⊥⊥ W | (Z ∪ Y)) ==> X ⊥⊥ (Y ∪ W) |Z
     */
    public boolean intersection() {
        for (GraphoidIndFact fact1 : facts) {
            final Set<Node> X = fact1.getX();
            Set<Node> Y = fact1.getY();
            Set<Node> ZW = fact1.getZ();

            List<Node> ZWList = new ArrayList<>(ZW);

            DepthChoiceGenerator gen = new DepthChoiceGenerator(ZWList.size(), ZWList.size());
            int[] choice;

            while ((choice = gen.next()) != null) {
                Set<Node> Z = GraphUtils.asSet(choice, ZWList);
                Set<Node> W = new HashSet<>(ZW);
                W.removeAll(Z);

                if (trivialtyAssumed && X.isEmpty() || W.isEmpty()) continue;

                Set<Node> ZY = new HashSet<>(Z);
                ZY.addAll(Y);

                GraphoidIndFact fact2 = new GraphoidIndFact(X, W, ZY);

                if ((facts.contains(fact2))) {
                    Set<Node> YW = new HashSet<>(Y);
                    YW.addAll(W);

                    if (YW.isEmpty()) continue;

                    GraphoidIndFact newFact = new GraphoidIndFact(X, YW, Z);

                    if (!facts.contains(newFact)) {
                        TetradLogger.getInstance().forceLogMessage("Intersection fails:" +
                                " Have " + fact1 + " and " + fact2 +
                                "; Missing " + newFact);
                        return false;
                    }
                }
            }
        }

        return true;
    }

    /**
     * (X ⊥⊥ Y | Z) ∧ (X ⊥⊥ W |Z) ==> X ⊥⊥ (Y ∪ W) |Z
     */
    public boolean composition() {
        for (GraphoidIndFact fact1 : facts) {
            Set<Node> X = fact1.getX();
            Set<Node> Y = fact1.getY();
            Set<Node> Z = fact1.getZ();

            for (GraphoidIndFact fact2 : new HashSet<>(facts)) {
                if (fact1.equals(fact2)) continue;

                if (fact2.getX().equals(X) && fact2.getZ().equals(Z)) {
                    Set<Node> W = fact2.Y;

                    Set<Node> YW = new HashSet<>(Y);
                    YW.addAll(W);

                    GraphoidIndFact newFact = new GraphoidIndFact(X, YW, Z);

                    if (!facts.contains(newFact)) {
                        TetradLogger.getInstance().forceLogMessage("Composition fails:" +
                                " Have " + fact1 + " and " + fact2 +
                                "; Missing " + newFact);
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public void setTrivialtyAssumed() {
        this.trivialtyAssumed = true;
    }

    /**
     * X ⊥⊥ Y | Z ==> Y ⊥⊥ X | Z
     */
    public void setSymmetryAssumed() {
        for (GraphoidIndFact IC : new HashSet<>(facts)) {
            Set<Node> X = IC.getX();
            Set<Node> Y = IC.getY();
            Set<Node> Z = IC.getZ();
            facts.add(new GraphoidIndFact(Y, X, Z));
        }
    }

    public static class GraphoidIndFact {
        private final Set<Node> X;
        private final Set<Node> Y;
        private final Set<Node> Z;

        public GraphoidIndFact(Set<Node> X, Set<Node> Y, Set<Node> Z) {
            if (X.isEmpty() || Y.isEmpty()) throw new IllegalArgumentException("X or Y is empty");
            if (!disjoint(X, Y, Z)) throw new IllegalArgumentException();

            this.X = new HashSet<>(X);
            this.Y = new HashSet<>(Y);
            this.Z = new HashSet<>(Z);
        }

        public Set<Node> getX() {
            return new HashSet<>(X);
        }

        public Set<Node> getY() {
            return new HashSet<>(Y);
        }

        public Set<Node> getZ() {
            return new HashSet<>(Z);
        }

        public int hashCode() {
            return 1;
        }

        public boolean equals(Object o) {
            if (!(o instanceof GraphoidIndFact)) return false;
            GraphoidIndFact _fact = (GraphoidIndFact) o;
            return X.equals(_fact.X) && Y.equals(_fact.Y) && Z.equals(_fact.Z);
        }

        public String toString() {
            return X + " : " + Y + " | " + Z;
        }

        private boolean disjoint(Set<Node> set1, Set<Node> set2, Set<Node> set3) {
            return intersection(set1, set2).isEmpty()
                    && intersection(set1, set3).isEmpty()
                    || !intersection(set2, set3).isEmpty();
        }

        private Set<Node> intersection(Set<Node> set1, Set<Node> set2) {
            Set<Node> W = new HashSet<>(set1);
            W.retainAll(set2);
            return W;
        }
    }
}
