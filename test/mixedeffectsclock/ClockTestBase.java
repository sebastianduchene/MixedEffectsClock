package mixedeffectsclock;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.parameter.RealParameter;

/**
 * Fixtures shared by the tests. A six-taxon tree whose design can be worked out by hand:
 *
 *     (((t1,t2)P1,(t3,t4)P2)X,(t5,t6)P3)root
 *
 * Ten branches, the root carrying none.
 */
public class ClockTestBase {

    static final String NEWICK = "(((t1:1,t2:1):1,(t3:1,t4:1):1):1,(t5:2,t6:2):1):0.0;";
    static final String[] TAXA = {"t1", "t2", "t3", "t4", "t5", "t6"};

    static TaxonSet taxonSet(String... names) {
        List<Taxon> t = new ArrayList<>();
        for (String n : names) {
            t.add(new Taxon(n));
        }
        TaxonSet ts = new TaxonSet();
        ts.initByName("taxon", t);
        ts.setID("taxa." + String.join("_", names));
        return ts;
    }

    static TreeParser tree() {
        TreeParser t = new TreeParser();
        t.initByName("newick", NEWICK, "IsLabelledNewick", true,
                     "adjustTipHeights", false, "taxonset", taxonSet(TAXA));
        return t;
    }

    static CladeDesign clade(String id, int category, boolean stem, boolean excl, String... taxa) {
        CladeDesign c = new CladeDesign();
        c.initByName("taxonset", taxonSet(taxa), "category", category,
                     "includeStem", stem, "excludeClade", excl);
        c.setID(id);
        return c;
    }

    static LogNormalDistributionModel meanOneLogNormal(RealParameter stdev) {
        LogNormalDistributionModel d = new LogNormalDistributionModel();
        d.initByName("S", stdev, "M", new RealParameter("1.0"), "meanInRealSpace", true);
        return d;
    }

    /**
     * Builds the clock. rates has one entry per branch; every value is 1 unless given.
     */
    static MixedEffectsClockModel clock(TreeParser tree, String coefficients,
                                        String clockRate, CladeDesign... clades) {
        int branches = tree.getNodeCount() - 1;
        RealParameter rates = new RealParameter("1.0");
        rates.setDimension(branches);
        RealParameter stdev = new RealParameter("0.3");
        MixedEffectsClockModel m = new MixedEffectsClockModel();
        m.initByName("tree", tree,
                     "rates", rates,
                     "clock.rate", new RealParameter(clockRate),
                     "coefficient", new RealParameter(coefficients),
                     "distr", meanOneLogNormal(stdev),
                     "clade", Arrays.asList(clades));
        return m;
    }

    /**
     * A label for the branch ABOVE a node that does not depend on node numbering: the
     * sorted set of tips below it. Node numbers shift as the topology moves, so asserting
     * on them would test the wrong thing.
     */
    static String branchLabel(Node node) {
        Set<String> tips = new TreeSet<>();
        collectTips(node, tips);
        return String.join(",", tips);
    }

    private static void collectTips(Node n, Set<String> out) {
        if (n.isLeaf()) {
            out.add(n.getID());
        } else {
            for (Node c : n.getChildren()) {
                collectTips(c, out);
            }
        }
    }
}
