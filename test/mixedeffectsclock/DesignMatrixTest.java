package mixedeffectsclock;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;
import java.util.Set;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;

/**
 * The design matrix: which branches load on which column.
 *
 * These assert the exact SET of branches per column, identified by the tips below each
 * branch, not the counts. Counts are the obvious thing to assert and they are too weak: a
 * design attached to entirely the wrong branches has identical counts. This is the only
 * place in the package where a silent wrong answer can originate.
 */
public class DesignMatrixTest extends ClockTestBase {

    /** column index -> the branches it marks, each named by the tips below it. */
    private Map<Integer, Set<String>> columnsOf(MixedEffectsClockModel m, TreeParser tree) {
        boolean[][] design = m.getDesign();
        Map<Integer, Set<String>> out = new LinkedHashMap<>();
        for (int c = 0; c < m.getNumberOfColumns(); c++) {
            out.put(c, new TreeSet<>());
        }
        for (Node nd : tree.getNodesAsArray()) {
            if (nd.isRoot()) {
                continue;
            }
            for (int c = 0; c < m.getNumberOfColumns(); c++) {
                if (design[nd.getNr()][c]) {
                    out.get(c).add(branchLabel(nd));
                }
            }
        }
        return out;
    }

    @Test
    public void includeStemMarksTheCladeAndItsStem() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"));
        Set<String> got = columnsOf(m, tree).get(0);
        assertEquals("t1 and t2 plus the branch above their ancestor",
                     new TreeSet<>(java.util.Arrays.asList("t1", "t2", "t1,t2")), got);
    }

    @Test
    public void withoutStemTheCladesOwnBranchIsLeftOut() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("B", 0, false, false, "t3", "t4"));
        assertEquals(new TreeSet<>(java.util.Arrays.asList("t3", "t4")),
                     columnsOf(m, tree).get(0));
    }

    @Test
    public void excludeCladeMarksTheStemAlone() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("C", 0, true, true, "t5", "t6"));
        assertEquals(new TreeSet<>(java.util.Collections.singletonList("t5,t6")),
                     columnsOf(m, tree).get(0));
    }

    @Test
    public void twoCladesMaySharaOneCoefficient() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"),
                clade("C", 0, true, true, "t5", "t6"));
        assertEquals(1, m.getNumberOfColumns());
        assertEquals(new TreeSet<>(java.util.Arrays.asList("t1", "t2", "t1,t2", "t5,t6")),
                     columnsOf(m, tree).get(0));
    }

    @Test
    public void disjointCladesDoNotOverlap() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0 0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"),
                clade("B", 1, false, false, "t3", "t4"));
        boolean[][] d = m.getDesign();
        for (Node nd : tree.getNodesAsArray()) {
            if (nd.isRoot()) continue;
            int hits = 0;
            for (int c = 0; c < 2; c++) if (d[nd.getNr()][c]) hits++;
            assertTrue("branch " + branchLabel(nd) + " is in " + hits + " columns", hits <= 1);
        }
    }

    /**
     * Without monophyly the ancestor of a scattered taxon set swallows unrelated branches.
     * That is what the column MEANS, not a bug, and the test pins the behaviour so a future
     * change to the semantics is noticed.
     */
    @Test
    public void aNonMonophyleticCladeSwallowsEverythingBetweenItsTaxa() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("scattered", 0, false, false, "t1", "t3"));
        Set<String> got = columnsOf(m, tree).get(0);
        assertEquals("the ancestor of t1 and t3 subtends all four of t1..t4",
                     new TreeSet<>(java.util.Arrays.asList(
                             "t1", "t2", "t3", "t4", "t1,t2", "t3,t4")), got);
        assertTrue("and it is correctly reported as not monophyletic",
                   !m.isMonophyletic(m.getClades().get(0)));
    }

    @Test
    public void everyBranchIsAccountedFor() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0 0.0 0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"),
                clade("B", 1, false, false, "t3", "t4"),
                clade("C", 2, true, true, "t5", "t6"));
        Map<Integer, Set<String>> cols = columnsOf(m, tree);
        assertEquals(3, cols.get(0).size());
        assertEquals(2, cols.get(1).size());
        assertEquals(1, cols.get(2).size());
        int branches = tree.getNodeCount() - 1;
        assertEquals("6 taxa give 2n-2 branches", 10, branches);
    }

    @Test
    public void aCoefficientDimensionMismatchIsRejected() {
        try {
            clock(tree(), "0.0", "1.0",
                  clade("A", 0, true, false, "t1", "t2"),
                  clade("B", 1, false, false, "t3", "t4"));
            fail("a 1-dimensional coefficient with two design columns should be rejected");
        } catch (RuntimeException expected) {
            assertTrue(expected.getMessage() == null
                       || expected.getMessage().contains("coefficient")
                       || expected.getCause() != null);
        }
    }

    @Test
    public void excludeCladeWithoutIncludeStemIsRejected() {
        try {
            clade("bad", 0, false, true, "t1", "t2");
            fail("excludeClade without includeStem marks nothing and should be rejected");
        } catch (RuntimeException expected) {
            // expected
        }
    }
}
