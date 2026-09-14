package mixedeffectsclock;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.RealParameter;

/**
 * Mixed-effects molecular clock, Bletsa et al. 2019 Equation 1:
 *
 *     log r_i = beta_0 + sum_k X_ik beta_k + eps_i,   eps_i ~ Normal(0, sigma^2)
 *
 * The random-effect half is inherited from UCRelaxedClockModel: clock.rate carries
 * exp(beta_0) and the branch multipliers carry exp(eps_i), drawn from a lognormal
 * pinned to mean one in real space. The dispersion of that lognormal is the
 * LOG-SCALE standard deviation, which is the published sigma. Do NOT confuse it
 * with BEAST X's branchRates.scale, which is a real-space coefficient of
 * variation; see MIXED_EFFECTS_BEAST2_PLAN.md for the conversion.
 *
 * All this class adds is the fixed-effect factor exp(sum_k X_ik beta_k).
 *
 * STATUS 2026-09-13: compiles, installs, and passes stage 1 of the plan. The design
 * matrix reproduces the BEAST X reference counts and stays correct under a moving
 * topology. Prior-only and known-truth validation not yet done. See README.md.
 */
@Description("Mixed-effects relaxed clock: a background rate and per-branch random "
        + "effects from the uncorrelated relaxed clock, times a fixed-effect factor "
        + "built from a clade design matrix (Bletsa et al. 2019).")
public class MixedEffectsClockModel extends UCRelaxedClockModel {

    final public Input<RealParameter> coefficientInput = new Input<>("coefficient",
            "fixed-effect coefficients beta_k on the LOG scale, one per design column. "
          + "Unbounded: a negative value is a rate decrease. The intercept beta_0 is "
          + "NOT here, it is carried by clock.rate as exp(beta_0).",
            Input.Validate.REQUIRED);

    final public Input<List<CladeDesign>> cladeInput = new Input<>("clade",
            "clade design elements, one or more per coefficient",
            new ArrayList<>(), Input.Validate.REQUIRED);

    private RealParameter coefficient;
    private Tree tree;
    private List<CladeDesign> clades;

    /** design[nodeNr][k] = true when branch nodeNr loads on coefficient k. */
    private boolean[][] design;
    private boolean[][] storedDesign;
    // volatile: getRateForBranch uses double-checked locking on this flag, and the
    // threaded tree likelihood calls it from several threads. Without volatile a
    // thread can see designKnown=true before the design array is visible to it.
    private volatile boolean designKnown = false;
    private boolean monophylyChecked = false;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        coefficient = coefficientInput.get();
        tree        = treeInput.get();
        clades      = cladeInput.get();

        // Normalising the inherited rates to mean one would be undone by the
        // fixed-effect factor applied afterwards, so the combination is not meaningful.
        if (normalizeInput.get()) {
            throw new IllegalArgumentException(getID()
                    + ": normalize=\"true\" is not supported; the fixed effect is applied "
                    + "after normalisation and would break it.");
        }

        int maxCategory = -1;
        for (CladeDesign c : clades) {
            maxCategory = Math.max(maxCategory, c.getCategory());
        }
        if (coefficient.getDimension() != maxCategory + 1) {
            throw new IllegalArgumentException(getID() + ": coefficient has dimension "
                    + coefficient.getDimension() + " but the clade design uses "
                    + (maxCategory + 1) + " columns.");
        }

        designKnown = false;
    }

    /**
     * Rebuild the design matrix against the CURRENT tree. Necessary on every
     * topology change: monophyly fixes which taxa form a clade, not which branches
     * sit inside it, and node numbering is not stable across moves either.
     */
    private void buildDesign() {
        design = computeDesign();
        designKnown = true;
        if (!monophylyChecked) {
            monophylyChecked = true;
            warnIfNotMonophyletic();
        }
    }

    /**
     * Warn once, on the first design build, about any clade that is not monophyletic.
     *
     * This runs on the first build rather than in initAndValidate because at
     * initAndValidate the tree is still BEAST's default one; the starting tree supplied
     * by an init element has not been installed yet.
     *
     * It is a warning and not an exception on purpose. A proposal that transiently breaks
     * monophyly is evaluated BEFORE the prior rejects it, so throwing would kill
     * perfectly valid chains. A clade that is already broken at the start, though, is
     * almost always a misspecified taxon set or a bad starting tree.
     */
    private void warnIfNotMonophyletic() {
        for (CladeDesign c : clades) {
            if (isMonophyletic(c)) {
                continue;
            }
            int have = c.getTaxonSet().asStringList().size();
            Node mrca = TreeUtils.getCommonAncestorNode(
                    tree, new LinkedHashSet<>(c.getTaxonSet().asStringList()));
            Log.warning.println("WARNING: " + getID() + ": clade " + c.getID()
                    + " (column " + c.getCategory() + ") is NOT monophyletic on the "
                    + "starting tree. Its common ancestor subtends "
                    + (mrca == null ? 0 : mrca.getLeafNodeCount()) + " tips but the taxon "
                    + "set has " + have + ". The column will paint branches outside the "
                    + "clade, and columns may overlap, in which case a branch carries the "
                    + "SUM of their coefficients. Constrain it with an MRCAPrior "
                    + "monophyletic=\"true\".");
            if (c.includeStem()) {
                Log.warning.println("         includeStem=\"true\" makes this worse: the "
                        + "stem is then the branch above a larger, unintended group, and "
                        + "if the ancestor is the root the column can end up with no "
                        + "branches at all and its coefficient unidentified.");
            }
        }
    }

    /** Is this clade monophyletic on the current tree? */
    public boolean isMonophyletic(CladeDesign c) {
        Set<String> taxa = new LinkedHashSet<>(c.getTaxonSet().asStringList());
        Node mrca = TreeUtils.getCommonAncestorNode(tree, taxa);
        return mrca != null && mrca.getLeafNodeCount() == taxa.size();
    }

    /**
     * Computes the design against the current tree and returns it, without touching
     * the cache. Public so a validator can compare a freshly computed design against
     * the cached one; they must always agree.
     */
    public boolean[][] computeDesign() {
        final int nodeCount = tree.getNodeCount();
        final int k = coefficient.getDimension();
        final boolean[][] d = new boolean[nodeCount][k];

        for (CladeDesign c : clades) {
            Set<String> taxa = new LinkedHashSet<>(c.getTaxonSet().asStringList());
            Node mrca = TreeUtils.getCommonAncestorNode(tree, taxa);
            if (mrca == null) {
                throw new RuntimeException(getID() + ": no common ancestor found for "
                        + c.getTaxonSet().getID());
            }
            final int col = c.getCategory();

            if (!c.excludeClade()) {
                markSubtree(mrca, col, d);       // every branch strictly inside the clade
            }
            if (c.includeStem() && !mrca.isRoot()) {
                d[mrca.getNr()][col] = true;     // the branch subtending the MRCA
            }
        }
        return d;
    }

    /** Marks every branch below node, i.e. the branch subtending each descendant. */
    private void markSubtree(Node node, int col, boolean[][] d) {
        for (Node child : node.getChildren()) {
            d[child.getNr()][col] = true;
            markSubtree(child, col, d);
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot()) {
            return 1;   // the root has no branch; match the parent class
        }
        if (!designKnown) {
            synchronized (this) {
                if (!designKnown) {
                    buildDesign();
                }
            }
        }

        double linearPredictor = 0.0;
        final boolean[] row = design[node.getNr()];
        for (int k = 0; k < row.length; k++) {
            if (row[k]) {
                linearPredictor += coefficient.getValue(k);
            }
        }

        return super.getRateForBranch(node) * Math.exp(linearPredictor);
    }

    /**
     * The parent class has its tree-dirty check commented out, so it never
     * recomputes on a topology change. The design matrix does depend on the
     * topology, so that check has to be reinstated here.
     */
    @Override
    protected boolean requiresRecalculation() {
        boolean parentDirty = super.requiresRecalculation();

        if (tree.somethingIsDirty()) {
            designKnown = false;
            return true;
        }
        if (coefficient.somethingIsDirty()) {
            return true;
        }
        return parentDirty;
    }

    @Override
    public void store() {
        if (design != null) {
            storedDesign = new boolean[design.length][];
            for (int i = 0; i < design.length; i++) {
                storedDesign[i] = design[i].clone();
            }
        }
        super.store();
    }

    @Override
    public void restore() {
        if (storedDesign != null) {
            design = storedDesign;
            storedDesign = null;
            designKnown = true;
        } else {
            designKnown = false;
        }
        super.restore();
    }

    /** Exposed so a logger or a test can dump the realised design. */
    public boolean[][] getDesign() {
        if (!designKnown) {
            buildDesign();
        }
        return design;
    }

    /** Number of design columns, i.e. the dimension of the coefficient vector. */
    public int getNumberOfColumns() {
        return coefficient.getDimension();
    }

    /** The clade design elements, in the order they were given. */
    public List<CladeDesign> getClades() {
        return clades;
    }

    /** The tree this clock is attached to. */
    public Tree getTree() {
        return tree;
    }
}
