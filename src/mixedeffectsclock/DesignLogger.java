package mixedeffectsclock;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import java.util.List;

/**
 * Logs how many branches load on each column of a mixed-effects clock design
 * matrix, plus how many carry the intercept alone. Column order follows the
 * coefficient vector.
 *
 * This exists to check a design against a known reference. The BEAST X runs in
 * this project wrote the same quantity as n_painted_X1..XK, so the columns here
 * are named to match.
 */
@Description("Per-column branch counts of a mixed-effects clock design matrix.")
public class DesignLogger extends BEASTObject implements Loggable {

    final public Input<MixedEffectsClockModel> clockInput = new Input<>("clock",
            "the mixed-effects clock whose design should be logged",
            Input.Validate.REQUIRED);

    final public Input<Boolean> validateInput = new Input<>("validate",
            "also recompute the design from scratch at every logged state and report "
          + "how many entries disagree with the cached one. Always 0 if the cache is "
          + "invalidated correctly; anything else means a stale design. Default false.",
            false);

    private MixedEffectsClockModel clock;
    private int columns;
    private boolean validate;
    private List<CladeDesign> clades;

    @Override
    public void initAndValidate() {
        clock = clockInput.get();
        columns = clock.getNumberOfColumns();
        validate = validateInput.get();
        clades = clock.getClades();
    }

    @Override
    public void init(PrintStream out) {
        for (int c = 0; c < columns; c++) {
            out.print("nPainted.X" + (c + 1) + "\t");
        }
        out.print("nPainted.background\t");
        out.print("nPainted.multiple\t");
        // One indicator per clade, mirroring the monophyly(...) columns the BEAST X runs
        // in this project log and filter on. A column that is not monophyletic is not
        // painting what it was written to paint.
        for (CladeDesign c : clades) {
            out.print("monophyletic." + c.getID() + "\t");
        }
        if (validate) {
            out.print("design.staleEntries\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        final boolean[][] design = clock.getDesign();
        final int root = clock.getTree().getRoot().getNr();

        int[] counts = new int[columns];
        int background = 0;
        int multiple = 0;

        for (int n = 0; n < design.length; n++) {
            if (n == root) {
                continue;               // the root has no branch
            }
            int hits = 0;
            for (int c = 0; c < columns; c++) {
                if (design[n][c]) {
                    counts[c]++;
                    hits++;
                }
            }
            if (hits == 0) {
                background++;
            } else if (hits > 1) {
                multiple++;             // should normally be zero
            }
        }

        for (int c = 0; c < columns; c++) {
            out.print(counts[c] + "\t");
        }
        out.print(background + "\t");
        out.print(multiple + "\t");
        for (CladeDesign c : clades) {
            out.print((clock.isMonophyletic(c) ? 1 : 0) + "\t");
        }

        if (validate) {
            // Per-column branch counts cannot detect a stale design: if the cache is
            // not rebuilt after a topology move the same entries stay marked and the
            // counts are unchanged, merely attached to the wrong branches. Comparing
            // against a freshly computed design does detect it.
            final boolean[][] fresh = clock.computeDesign();
            int stale = 0;
            for (int n = 0; n < design.length; n++) {
                for (int c = 0; c < columns; c++) {
                    if (design[n][c] != fresh[n][c]) {
                        stale++;
                    }
                }
            }
            out.print(stale + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }
}
