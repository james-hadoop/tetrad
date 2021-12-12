package edu.cmu.tetrad.session;

import java.util.HashSet;
import java.util.Set;

public class AvailableModels {
    private static final AvailableModels unique = new AvailableModels();
    private final Set<Object> set = new HashSet<>();

    public static AvailableModels getUnique() {
        return unique;
    }

    public Set<Object> getSet() {
        return new HashSet<>(set);
    }

    public void add(Object model) {
        if (model == null) return;
        set.add(model);
    }

    public void remove(Object model) {
        if (model == null) return;
        set.remove(model);
    }
}
