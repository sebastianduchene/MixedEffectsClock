---
name: me-design-prior-check
description: End-to-end verifier of a BEAST 2 MixedEffectsClock design matrix - redact the sequences to Ns so the chain samples from the prior, run a short MCMC, then draw ten sampled trees with every branch coloured by the fixed-effect column (clade) it carries, into a PDF. Use when the user wants to CHECK or VERIFY the design matrix of a mixed-effects clock before launching a production run, wants to "sample from the prior" to confirm clades paint the branches they intended, asks for a prior-only / likelihood-flat version of a MixedEffectsClock analysis, wants to redact sequences in a BEAST 2 XML, wants to colour-code sampled trees by clade, or asks "does the mixed-effects design do what I think it does". BEAST 2 MixedEffectsClock stack only (not BEAST X <fixedEffects>). Automates exactly the all-N-sequence recipe that the sibling mixed-effects-clock-beast2 skill flags as point 1 of its "six things that will bite you".
metadata:
  type: tool
  engine: BEAST 2.7.7+, MixedEffectsClock v0.0.1+
  depends_on: ape (R), BEAST 2 on PATH as `beast`
---

# Prior-check for a mixed-effects clock design

Given a BEAST 2 XML that uses `mixedeffectsclock.MixedEffectsClockModel`, produce
a PDF that shows the **design matrix the model will actually use**, drawn on
trees sampled from the prior. The user eyeballs it to confirm each named clade
paints the branches they intended before committing a long run.

## When to use

- "I want to check the design matrix / clade colouring on this XML."
- "Run a prior-only version of this analysis."
- "Redact the sequences and sample the prior."
- "Does my ME-clock design do what I think?"
- Any time a new `<clade spec="mixedeffectsclock.CladeDesign" ...>` block was
  added, moved or had its `includeStem` flipped — the only honest way to confirm
  the change does what was intended is to look at it on a tree.

For colouring a single already-sampled tree, prefer the sibling skill
`colour-me-clock-tree`. For the full model + implementation reference, see the
sibling skill `mixed-effects-clock-beast2`.

Scope: **BEAST 2 only**. The BEAST X dialect uses `<fixedEffects>` with
`<category>` and `<clade>` elements and is a different parser; this skill will
refuse that XML and point elsewhere.

## Run it

```bash
python3 ~/.claude/skills/me-design-prior-check/scripts/sample_prior.py \
    <analysis.xml> [--states 1000000] [--trees 10] [--burnin 0.1] [--beast beast]
```

Produces, next to `<analysis.xml>`:

| File | What it is |
| --- | --- |
| `<stem>_prior.xml`   | redacted + chain-shortened copy |
| `<stem>_prior.log`   | BEAST trace log (prior sample) |
| `<stem>_prior.trees` | sampled trees |
| `<stem>_prior.pdf`   | the design-check PDF, one page per sampled tree + a summary page |

Defaults: 1,000,000 states, 10 trees drawn, first 10% of sampled trees
discarded as burn-in.

## What the script actually does

1. **Fails loud** if the XML has no `mixedeffectsclock.MixedEffectsClockModel`
   block, or if it looks like BEAST X (`<fixedEffects>` is present). Trying to
   check a design in the wrong dialect is strictly worse than erroring out.
2. **Replaces every `<sequence ... value="..."/>` content with `N`s** of
   matching length. Totalcount and all other attributes are left alone.
3. **Rewrites the top-level `<run>` `chainLength`** to `--states` and every
   `<logger ... mode="tree">` `logEvery` so the chain logs enough trees that,
   after burn-in, at least `--trees` remain. The trace logger is rewritten
   to the same interval to keep the log file small.
4. **Runs `beast -overwrite -seed 42 <stem>_prior.xml`** in the XML's directory.
5. **Calls the R script** `colour_prior_trees.R`, which parses the clade
   design from the XML the same way the `colour-me-clock-tree` skill does,
   reads the trees file (robust to a still-writing final line), picks the
   requested number of trees evenly across the post-burn-in portion, assigns
   each branch to its column by MRCA + strict-descendants (+ stem when
   `includeStem="true"`), and writes a multi-page PDF.

## Why all-N sequences, not `sampleFromPrior=true`

The `MixedEffectsClock` package subclasses the stock relaxed clock and its
design cache is invalidated through the `TreeLikelihood`. Removing the
likelihood breaks the invalidation chain, the design cache goes stale and
every number after that is wrong with no error. The sibling skill
`mixed-effects-clock-beast2` documents this as its first "will bite you".
The all-N trick keeps the wiring intact and leaves the likelihood flat so
topology moves freely under the prior.

## Reading the PDF

Page 1 is a **summary page**: mean branches-per-column across the sampled
trees, monophyly rate per clade, and the count of trees in which any two
columns overlapped.

Each following page is one sampled tree with branches coloured by column:
background = grey, each clade = its own colour. Confirm, going down the pages:

- Each clade paints the subtree you expected, and the colour palette is
  consistent across trees.
- `monophyly` is `T` for every clade on every tree. If it drifts to `F` on
  any tree, add or tighten the `MRCAPrior monophyletic="true"` constraint.
- A no-stem clade of `n` taxa paints `2n-2` branches; a `+stem` clade adds one.
- No **OVERLAP** warnings in the stdout / summary page. Any overlap means
  the design is not exclusive and two columns are loading on the same
  branches.

## Caveats

- The MCMC seed is pinned to 42 for reproducibility. Pass `--states` and rerun
  with a fresh output name if you want a second independent sample.
- The redaction is in-place on sequence values only. If the XML has a path to
  an external alignment file, this skill will not catch it; it only redacts
  sequences inlined in the XML. The project's convention is inlined `<sequence
  value="...">` blocks and the sibling skill's templates follow that.
- A design cache can still go stale at run time even if the one-off check
  looks fine. For an in-run guard, add the `mixedeffectsclock.DesignLogger`
  with `validate="true"` to the trace logger (see the sibling
  `mixed-effects-clock-beast2` skill, section "Checking the design before
  trusting anything"). The two checks are complementary.

## Install

The skill ships inside the `MixedEffectsClock` package. To make it available
globally the same way the sibling skill is installed:

```bash
ln -s "$PWD/skills/me-design-prior-check" ~/.claude/skills/
```

Only R package `ape` is required; the Python wrapper uses the standard
library only.
