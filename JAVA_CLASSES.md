# The Java in MixedEffectsClock, and what each class does

Three classes, 434 lines in total, in the package `mixedeffectsclock`. Written
2026-09-12 and 2026-09-13.

## Why the package exists

BEAST X already has the Bletsa et al. 2019 mixed-effects molecular clock. BEAST 2 does
not, and no combination of existing BEAST 2 components can express it: an exhaustive
survey of all 26 installed packages found eight branch-rate models, none of them keyed to
clades, and every one of them declares its per-branch rate vector as a *parameter* rather
than as a computed *function*, so a rate vector assembled from other parts cannot be
plugged in. The one thing BEAST 2 offers that BEAST X does not is the Mascot structured
coalescent together with the targeted tree operators. Running this clock with that tree
prior is the only reason to port it, and it needs Java.

The model is a log-link generalised linear model on branch rates:

    log r_i  =  beta_0  +  sum_k X_ik beta_k  +  eps_i,      eps_i ~ Normal(0, sigma^2)

`X_ik` is 1 when branch *i* belongs to clade *k*. Exponentiating, a branch rate is a
background rate, times a fold change for each clade the branch belongs to, times a
per-branch random deviation.

---

## 1. `MixedEffectsClockModel.java` — the clock (272 lines)

`extends beast.base.evolution.branchratemodel.UCRelaxedClockModel`

**The central decision is that it inherits almost everything.** Two of the three terms
above already exist in BEAST 2's stock uncorrelated relaxed clock: the clock rate is the
exponentiated intercept, and the branch multipliers are the exponentiated deviations,
drawn from a lognormal pinned to mean one in real space. Subclassing rather than starting
from scratch means the mean-one pinning that keeps the intercept identifiable comes for
free, the optimised relaxed clock operators keep working, and the dispersion parameter
keeps exactly the meaning it has in the existing BEAST 2 analyses.

So the entire model addition is one line: the rate for a branch is the parent class's
rate multiplied by the exponential of the sum of the coefficients whose columns that
branch loads on. Everything else in the file works out which columns those are, and makes
sure that answer is never stale.

**Inputs beyond the inherited ones**

| Input | Type | Meaning |
| --- | --- | --- |
| `coefficient` | RealParameter | the `beta_k`, on the log scale, one per design column, unbounded so a clade can slow down |
| `clade` | list of `CladeDesign` | the design matrix, one or more entries per column |

The intercept is **not** an input. It is carried by the inherited `clock.rate` as
`exp(beta_0)`, which is what keeps the existing up-and-down operators meaningful.

**What the methods do**

- `getRateForBranch(node)` returns the parent's rate times the clade factor. The root
  returns 1, matching the parent's convention that the root has no branch.
- `computeDesign()` builds the branch-by-column matrix against the current tree: for each
  clade it finds the most recent common ancestor of the taxon set, marks every branch in
  that subtree, and marks the stem when asked. It returns a fresh matrix without touching
  the cache, so a validator can compare the two.
- `requiresRecalculation()` invalidates the cached design whenever the tree changes.
  **The parent class has its own tree-dirty check commented out in the BEAST source**, so
  without this override the design would never be rebuilt after a topology move and the
  model would be quietly wrong.
- `store()` and `restore()` keep a copy of the design so a rejected proposal restores the
  previous matrix rather than recomputing it.
- `warnIfNotMonophyletic()` warns once, on the first design build, about any clade that is
  already broken on the starting tree, with a sharper message when the stem is included.
  It runs on the first build and not in `initAndValidate` because at that point the tree
  is still BEAST's default one and the starting tree has not been installed yet.
- `isMonophyletic(clade)` asks whether a clade is monophyletic on the current tree, by
  comparing the number of tips under its ancestor with the size of the taxon set.

---

## 2. `CladeDesign.java` — one entry of the design matrix (42 lines)

A small holder describing one clade and how it maps onto branches.

| Input | Default | Meaning |
| --- | --- | --- |
| `taxonset` | required | the taxa whose common ancestor defines the clade |
| `includeStem` | false | also assign the branch subtending that ancestor |
| `excludeClade` | false | assign **only** the stem, nothing inside |
| `category` | 0 | which coefficient this clade loads on |

Several entries may share a `category`, which is how two disjoint groups of branches share
one coefficient. It rejects the one combination that would assign no branches at all,
`excludeClade` without `includeStem`.

---

## 3. `DesignLogger.java` — the diagnostics (120 lines)

A `Loggable` that writes, at every sampled state:

| Column | Meaning |
| --- | --- |
| `nPainted.X1..XK` | branches loading on each column |
| `nPainted.background` | branches carrying the intercept alone |
| `nPainted.multiple` | branches in more than one column; must be 0 for disjoint clades |
| `monophyletic.<clade>` | 1 or 0 per clade, mirroring the `monophyly(...)` columns the BEAST X runs in this project already log and filter on |
| `design.staleEntries` | optional, see below |

**Why the last column exists.** Per-column counts cannot detect a stale design. If the
cache were never rebuilt after a topology move, the same entries would stay marked and the
totals would look perfectly normal while pointing at the wrong branches. With
`validate="true"` the logger recomputes the design from scratch at each sample and reports
how many entries disagree with the cache. It must always be 0. This column caught two real
wiring errors during development.

---

## How they fit together

```xml
<branchRateModel id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"
                 tree="@Tree" rates="@rates" clock.rate="@clockRate"
                 coefficient="@coefficient">
  <distr spec="...LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">
    <parameter estimate="false" name="M">1.0</parameter>
  </distr>
  <clade id="hyper" spec="mixedeffectsclock.CladeDesign" category="0" includeStem="true">
    <taxonset idref="taxa.hyper"/>
  </clade>
  <!-- one clade element per design column -->
</branchRateModel>
```

Three things about the surrounding XML matter and are not obvious.

**The clock must sit inside the posterior**, normally hanging off a tree likelihood. BEAST
invalidates a cache only for objects reachable from the posterior, so a clock declared
inside a logger is never told the tree moved and its design silently goes stale. For a
prior-only run, keep the tree likelihood but make the alignment entirely ambiguous: the
likelihood is then flat and the wiring stays correct.

**Constrain every clade that carries a column**, and point the constraint at the clock's
own `TaxonSet` by `idref` rather than redeclaring the taxa. Without monophyly a clade's
ancestor subsumes unrelated branches, columns overlap, and a branch in two columns carries
the sum of their coefficients. A second copy of the taxa would let the constraint and the
design drift apart with nothing to complain about.

**Every public class must be registered** in `version.xml` as a `beast.base.core.BEASTInterface`
service, or BEAST reports a missing class in a way that reads like a missing package.

---

## Decisions worth explaining

**No Hamiltonian Monte Carlo.** None exists in BEAST 2 or in any installed package, and
the project's own benchmark puts tuned gradient sampling at 7.4 seconds per effective
sample of the dispersion against 10.9 for plain Metropolis, while losing on the intercept,
the fixed effect and the mean rate. Little is given up.

**Overlapping columns are reported, never thrown.** A proposal that transiently breaks
monophyly is evaluated *before* the prior rejects it, so an exception would kill otherwise
valid chains.

**The dispersion is the log-scale standard deviation, which is the published sigma.** It is
not BEAST X's `branchRates.scale`, which is a real-space coefficient of variation. Convert
with `sigma = sqrt(log(1 + scale^2))`. The two agree to a couple of percent below about 0.3
and diverge badly above 1.

**Under a mixed-effects clock the dispersion is a residual**, left over after the named
clades have taken their share. In the BEAST X production posterior it is 0.26 against a
realised branch-rate coefficient of variation of 1.04. A small dispersion means the clade
assignments are doing work, not that the tree is clocklike.

---

## Status

Verified: the package builds and installs; the design matrix reproduces the BEAST X
reference counts in both stem configurations; the design stays correct while the topology
moves, with no stale cache entries; and branch rates are numerically exact, including a
slowdown and overlapping columns.

Not yet done: **the coefficient vector has never been sampled** by any operator, so the
moves are untested; the lazy rebuild uses double-checked locking on a non-volatile flag,
which matters only under the threaded likelihood; and prior-only recovery, the known-truth
simulation and agreement with BEAST X are all still ahead. No XML combining this clock with
Mascot exists yet, which is the actual objective.

The surrounding files are not Java: `examples/` holds small XMLs and two generators that
build design checks from the real data, `validation/` holds an independent rate checker
written with dendropy and an R script that draws the trees coloured by design column, and
`../lepromatosis/MIXED_EFFECTS_BEAST2_PLAN.md` is the specification and validation plan.
