package mixedeffectsclock;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

/**
 * The rate a branch actually gets, and the design cache that feeds it.
 *
 * The cache is the highest-risk part of the class. The parent has its own tree-dirty check
 * commented out in the BEAST source, so the subclass reinstates it; if someone later tidies
 * that away, the model keeps running and is quietly wrong. These tests are what would notice.
 */
public class RateAndCacheTest extends ClockTestBase {

    private static final double EPS = 1e-12;

    private Node branchAbove(TreeParser tree, String tips) {
        for (Node nd : tree.getNodesAsArray()) {
            if (!nd.isRoot() && branchLabel(nd).equals(tips)) {
                return nd;
            }
        }
        throw new IllegalArgumentException("no branch above " + tips);
    }

    /** Swap two nodes between their parents: a topology change with no API for it. */
    private void swap(Node a, Node b) {
        Node pa = a.getParent(), pb = b.getParent();
        pa.removeChild(a);
        pb.removeChild(b);
        pa.addChild(b);
        pb.addChild(a);
    }

    @Test
    public void rateIsTheClockRateTimesTheProductOfFoldChanges() {
        TreeParser tree = tree();
        // log 4, log 0.2, log 10
        MixedEffectsClockModel m = clock(tree,
                "1.3862943611198906 -1.6094379124341003 2.302585092994046", "1.0E-3",
                clade("A", 0, true, false, "t1", "t2"),
                clade("B", 1, false, false, "t3", "t4"),
                clade("C", 2, true, true, "t5", "t6"));

        assertEquals("fold 4 on a 1e-3 background", 4.0E-3,
                     m.getRateForBranch(branchAbove(tree, "t1")), EPS);
        assertEquals(4.0E-3, m.getRateForBranch(branchAbove(tree, "t1,t2")), EPS);
        assertEquals("fold 0.2, a slowdown", 2.0E-4,
                     m.getRateForBranch(branchAbove(tree, "t3")), EPS);
        assertEquals("fold 10, on the stem alone", 1.0E-2,
                     m.getRateForBranch(branchAbove(tree, "t5,t6")), EPS);
        // background: the stem of B is excluded, and t5/t6 are inside an excluded clade
        assertEquals(1.0E-3, m.getRateForBranch(branchAbove(tree, "t3,t4")), EPS);
        assertEquals(1.0E-3, m.getRateForBranch(branchAbove(tree, "t5")), EPS);
        assertEquals("the root carries no branch", 1.0,
                     m.getRateForBranch(tree.getRoot()), EPS);
    }

    @Test
    public void aNegativeCoefficientSlowsTheBranchDown() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "-1.2039728043259361", "1.0E-3",  // log 0.3
                clade("slow", 0, true, false, "t1", "t2"));
        assertEquals(3.0E-4, m.getRateForBranch(branchAbove(tree, "t1")), EPS);
    }

    /** Overlapping columns add on the log scale, so the fold changes multiply. */
    @Test
    public void overlappingColumnsMultiply() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree,
                "0.6931471805599453 1.0986122886681098", "1.0E-3",   // log 2, log 3
                clade("outer", 0, false, false, "t1", "t3"),   // covers t1..t4
                clade("inner", 1, false, false, "t1", "t2"));
        assertEquals("t1 is in both columns: 2 x 3", 6.0E-3,
                     m.getRateForBranch(branchAbove(tree, "t1")), EPS);
        assertEquals("t3 is only in the outer one", 2.0E-3,
                     m.getRateForBranch(branchAbove(tree, "t3")), EPS);
    }

    @Test
    public void theDesignFollowsTheTopologyWhenRecomputed() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"));

        boolean[][] before = m.computeDesign();
        int markedBefore = count(before);

        swap(branchAbove(tree, "t1"), branchAbove(tree, "t5,t6"));

        boolean[][] after = m.computeDesign();
        assertFalse("moving t1 next to t5,t6 must change which branches the column marks",
                    same(before, after));
        assertTrue("and the clade is no longer monophyletic",
                   !m.isMonophyletic(m.getClades().get(0)));
        assertTrue("its ancestor now subtends more, so more branches are marked",
                   count(after) > markedBefore);
    }

    /**
     * The cache must not silently outlive the tree it was built for. getDesign() returns the
     * cached matrix; computeDesign() always rebuilds. After a topology change they disagree,
     * and that disagreement is exactly what DesignLogger's validate mode reports.
     */
    @Test
    public void theCachedDesignIsDetectablyStaleAfterATopologyChange() {
        TreeParser tree = tree();
        MixedEffectsClockModel m = clock(tree, "0.0", "1.0",
                clade("A", 0, true, false, "t1", "t2"));

        boolean[][] cached = deepCopy(m.getDesign());       // fills the cache
        swap(branchAbove(tree, "t1"), branchAbove(tree, "t5,t6"));

        assertFalse("a fresh computation must differ from the stale cache",
                    same(cached, m.computeDesign()));

        // requiresRecalculation reads the tree's dirty FLAG, not the tree itself, and
        // mutating nodes directly does not set it. BEAST sets it through the State when an
        // operator edits the tree; a test has to do it by hand.
        tree.setSomethingIsDirty(true);
        m.requiresRecalculation();                          // what BEAST calls
        assertTrue("after invalidation the cache must agree with a fresh computation",
                   same(m.getDesign(), m.computeDesign()));
    }

    @Test
    public void normalisationIsRejected() {
        TreeParser tree = tree();
        try {
            MixedEffectsClockModel m = new MixedEffectsClockModel();
            RealParameter rates = new RealParameter("1.0");
            rates.setDimension(tree.getNodeCount() - 1);
            m.initByName("tree", tree, "rates", rates,
                         "clock.rate", new RealParameter("1.0"),
                         "coefficient", new RealParameter("0.0"),
                         "distr", meanOneLogNormal(new RealParameter("0.3")),
                         "clade", java.util.Collections.singletonList(
                                 clade("A", 0, true, false, "t1", "t2")),
                         "normalize", true);
            fail("normalize=true would be undone by the fixed effect and must be rejected");
        } catch (RuntimeException expected) {
            // expected
        }
    }

    private static int count(boolean[][] d) {
        int n = 0;
        for (boolean[] row : d) for (boolean b : row) if (b) n++;
        return n;
    }

    private static boolean same(boolean[][] a, boolean[][] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] != b[i][j]) return false;
            }
        }
        return true;
    }

    private static boolean[][] deepCopy(boolean[][] a) {
        boolean[][] c = new boolean[a.length][];
        for (int i = 0; i < a.length; i++) c[i] = a[i].clone();
        return c;
    }
}
