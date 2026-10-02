package mixedeffectsclock;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;

@Description("One column entry of the mixed-effects clock design matrix: a clade, "
        + "optionally including its stem branch, optionally the stem only, optionally "
        + "only the terminal branches. Several CladeDesign elements sharing a category "
        + "load on the same coefficient.")
public class CladeDesign extends BEASTObject {

    final public Input<TaxonSet> taxonsetInput = new Input<>("taxonset",
            "taxa whose most recent common ancestor defines the clade",
            Input.Validate.REQUIRED);

    final public Input<Boolean> includeStemInput = new Input<>("includeStem",
            "also assign the branch subtending the clade's MRCA (default false)", false);

    final public Input<Boolean> includeTerminalInput = new Input<>("includeTerminal",
            "assign the terminal (tip) branches inside the clade (default true)", true);

    final public Input<Boolean> excludeCladeInput = new Input<>("excludeClade",
            "exclude the internal (non-terminal) branches inside the clade (default false)",
            false);

    final public Input<Integer> categoryInput = new Input<>("category",
            "index of the coefficient this clade loads on, counting from 0", 0);

    @Override
    public void initAndValidate() {
        if (categoryInput.get() < 0) {
            throw new IllegalArgumentException(getID() + ": category must be >= 0");
        }
        if (excludeCladeInput.get() && !includeStemInput.get() && !includeTerminalInput.get()) {
            throw new IllegalArgumentException(getID()
                    + ": excludeClade=\"true\" with includeStem=\"false\" and "
                    + "includeTerminal=\"false\" assigns no branches at all");
        }
    }

    public int getCategory()        { return categoryInput.get(); }
    public boolean includeStem()    { return includeStemInput.get(); }
    public boolean includeTerminal(){ return includeTerminalInput.get(); }
    public boolean excludeClade()   { return excludeCladeInput.get(); }
    public TaxonSet getTaxonSet()   { return taxonsetInput.get(); }
}
