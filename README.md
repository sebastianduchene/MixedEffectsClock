# MixedEffectsClock

A BEAST 2 package implementing the mixed-effects molecular clock of Bletsa et al. (2019),
so that it can be run together with the Mascot structured coalescent and the targeted tree
operators.

Version 0.0.1. Working and validated against known truth; not yet used for a published
analysis. See [Status](#status) for exactly what has and has not been tested.

## Why it exists

BEAST X already has this clock. BEAST 2 does not, and no combination of existing BEAST 2
components can express it. A survey of all 26 packages installed on the development machine
found eight concrete branch-rate models, none of them keyed to clades, and every one
declares its per-branch rate vector as a *parameter* rather than a computed *function*, so
a rate vector assembled from other parts cannot be connected. What BEAST 2 offers and
BEAST X does not is the Mascot structured coalescent alongside the targeted operators.
Running this clock with that tree prior is the reason the package exists.

## The model

Branch rates follow a log-link generalised linear model:

    log r_i  =  beta_0  +  sum_k X_ik beta_k  +  eps_i,      eps_i ~ Normal(0, sigma^2)

`X_ik` is 1 when branch *i* belongs to clade *k* and 0 otherwise. Exponentiating, a branch
rate is a background rate, times a fold change for each clade the branch belongs to, times
a per-branch random deviation:

    r_i  =  exp(beta_0) * prod_k exp(beta_k)^X_ik * exp(eps_i)

How that maps onto the XML:

| Model term | Where it lives |
| --- | --- |
| `exp(beta_0)`, the background rate | the inherited `clock.rate` |
| `beta_k`, log fold changes | the `coefficient` parameter, one entry per design column |
| `X_ik`, the design matrix | one `CladeDesign` element per column |
| `exp(eps_i)` | the inherited per-branch `rates`, from a mean-one lognormal |
| `sigma` | the log-scale standard deviation of that lognormal |

The coefficients are unbounded, so a clade can be slower than background as well as faster.
There is no separate intercept parameter: `clock.rate` carries it, which is what keeps the
standard up-and-down operators meaningful.

## The Java that was written

437 lines in three classes, in the package `mixedeffectsclock`. Everything else in this
repository is XML, generators or scoring scripts.

**The central decision was to inherit rather than start from scratch.**
`MixedEffectsClockModel` extends BEAST 2's stock `UCRelaxedClockModel` instead of
implementing `BranchRateModel` directly. Two of the three terms in the model already exist
there: the clock rate is the exponentiated intercept, and the branch multipliers are the
exponentiated deviations, drawn from a lognormal pinned to mean one in real space. Inheriting
means the mean-one pinning that keeps the intercept identifiable comes for free, the
optimised relaxed clock operators still apply, and the dispersion parameter keeps exactly the
meaning it has in every other BEAST 2 analysis.

So the entire model addition is one line: the rate for a branch is the parent's rate
multiplied by the exponential of the sum of the coefficients whose columns that branch loads
on. Everything else in the file works out which columns those are, and makes sure that answer
is never stale.

| Class | Lines | What it does |
| --- | --- | --- |
| `MixedEffectsClockModel` | 275 | the clock: the clade factor, the design matrix, and its cache |
| `CladeDesign` | 42 | one design-matrix entry: a taxon set, the stem and exclusion flags, and which coefficient it loads on |
| `DesignLogger` | 120 | diagnostics: branch counts per column, monophyly indicators, and a stale-cache check |

**Most of the complexity is cache invalidation, and that is where the risk lives.** The design
matrix depends on the topology, because monophyly fixes which taxa form a clade but not which
branches sit inside it, and internal node numbering shifts as subtrees move. So it cannot be
computed once. The class marks itself dirty when the tree changes, rebuilds lazily, and keeps
a copy so a rejected proposal restores the previous matrix rather than recomputing it.

One thing here was not obvious. The parent class has its own tree-dirty check **commented
out** in the BEAST source, so it never recomputes on a topology change. The subclass has to
reinstate it. Without that the model would look fine and be quietly wrong, which is why
`DesignLogger` can recompute the design from scratch and report disagreements: per-column
counts cannot detect a stale cache, since the same entries stay marked and only point at the
wrong branches.

**What is deliberately not ours.** No operators: the moves are BEAST's and ORC's, which is
what subclassing buys. No priors, no distribution: the mean-one lognormal is the stock one. No
gradients or Hamiltonian sampling, because BEAST 2 has none and the benchmark says little is
lost. `JAVA_CLASSES.md` covers all three classes method by method.

## Installation

**Prerequisites.**

| | | |
| --- | --- | --- |
| BEAST | 2.7.7 or newer | |
| A JDK | 17 or newer, to build | a Java 11 compiler fails with `class file has wrong version 61.0` |
| Apache Ant | any recent | |
| Python 3 | for the generators and checkers | `dendropy` only for `check_rates.py` |
| R | for the tree figures only | `ape` |

**BEAST packages.** Install these through BEAUti's package manager, or `packagemanager -add`:

| Package | Needed for |
| --- | --- |
| ORC | the dispersion operator, which is not optional; see below |
| feast | the dispersion conversion column |
| Mascot | the structured coalescent, only with `--mascot` |
| TargetedBeast | the targeted tree moves, only with `--targeted` |
| BEASTLabs | a metric used by ORC's exchange operators |

**Build and install this package.**

    ant install

That compiles with `--release 17`, jars the package and unpacks it into
`~/.beast/2.7/MixedEffectsClock`, where BEAST finds it beside the others. The build looks for
a JDK in the usual places; give it `-Djdk.home=/path/to/jdk` if it cannot find one. Other
targets: `build` for the jar alone, `addon` for an installable zip, `clean`.

Check it loaded by running any example: BEAST prints its package list at startup and
`MixedEffectsClock` should appear there.

## Running an analysis

Everything runs from this repository alone, with no external data.

    cd examples

    # 1. write an XML. This one is 50 taxa with known truth: a background rate of
    #    0.005, one clade at twice that, one at 0.3 times, and real per-branch noise.
    python3 gen_sim_re_slowdown.py > sim.xml

    # 2. run it
    beast -overwrite -seed 17 sim.xml

    # 3. score it against the truth: posterior median, 95% interval, whether it
    #    covers, effective sample size, and the design counts
    python3 ../validation/check_sim_re.py sim

For the full stack, the configuration this package exists for:

    python3 gen_sim_re_slowdown.py --targeted --mascot --boundaries "2.5" > stack.xml

To see which branches ended up with which rate, `validation/check_rates.py` recomputes every
branch rate independently with dendropy and writes the files that `plot_rate_trees.R` turns
into a tree per page, coloured by design column.

`validation/ANALYSES.md` lists every analysis run during development, the one command that
rebuilds each, and what it established. The XMLs themselves are not committed: each inlines
the alignment and they come to 9 MB, all of it regenerable.

**Reading a run.** Three columns are worth checking before anything else. `nPainted.multiple`
must be 0 whenever the clades are disjoint. Every `monophyletic.<clade>` must be 1. And if
`validate="true"` is set on the design logger, `design.staleEntries` must be 0 at every
sample. Any of those going wrong means the design is not what the XML says it is, and the
rest of the output is not worth reading.

`examples/design_check_6taxon.xml` is the smallest thing to read first: a six-taxon tree
whose design can be counted by hand.

### Generator options

The simulation generator has options that stack, so one script produces every configuration:

| Command | Tree moves | Tree prior |
| --- | --- | --- |
| `gen_sim_re_slowdown.py` | BEAUti's standard set | exponential-growth coalescent |
| `... --targeted` | targetedbeast | exponential-growth coalescent |
| `... --targeted --mascot` | targetedbeast | **Mascot**, two demes, one unsampled |

The Mascot skyline grid is set one of two ways, and never by hand:

| Option | Grid | Boundaries |
| --- | --- | --- |
| `--levels K` | fractions of the current root height | move with the tree, so an interval is not a fixed span |
| `--boundaries "2.5"` | times before the most recent tip | fixed calendar spans |

What makes the difference is whether the `rateShifts` element is given a `tree`. With one, the
values are read as fractions of the current root height; without one, as absolute times. An
absolute grid that does not reach the root is how the lepromatosis skyline was silently
collapsed to a constant, so the relative form is the safer default and the calendar form is
the one to reach for when the intervals need to mean fixed periods.

`--levels K` gives one more rate level than change points; the default is 2, so the rate
shifts once. **Do not set the grid by hand.**
The arithmetic is not what the XML suggests: `StructuredMigrationSkyline` caps its interval
index two below the number of shift values, so N values give N-1 levels, and `Skygrowth`
separately forces its parameter to dimension N+1 with the last entry never read by Mascot.
For K levels that is K+1 shift values and dimension K+2, of which K+1 are live. Verified by
perturbing each entry and watching the Mascot density.

A single shift value does not work and fails silently: the run exits 0, logs nothing, and
never produces a sample. The minimum is two values, which is one level.

The last is the configuration the package exists for: the mixed-effects clock, the targeted
tree proposals, ORC's rate moves and a structured coalescent with a ghost deme, all at once.
Note that the simulated data were not generated under population structure, so Mascot is
deliberately mis-specified there; the run tests that the stack holds together and that the
clock still recovers its coefficients, not that the demes mean anything.

## A minimal clock block

This is the clock element alone. For a whole analysis that runs, copy
`skills/mixed-effects-clock-beast2/templates/me_clock_minimal.xml` instead: it has the
priors, the constraints and the operators already wired, and the six traps below are what
it is wired to avoid.

```xml
<branchRateModel id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"
                 tree="@Tree" rates="@rates" clock.rate="@clockRate"
                 coefficient="@coefficient">
  <distr spec="beast.base.inference.distribution.LogNormalDistributionModel"
         S="@ucldStdev" meanInRealSpace="true">
    <parameter estimate="false" name="M">1.0</parameter>
  </distr>

  <clade id="hyper" spec="mixedeffectsclock.CladeDesign" category="0" includeStem="true">
    <taxonset idref="taxa.hyper"/>
  </clade>
  <!-- one clade element per design column; several may share a category -->
</branchRateModel>
```

`CladeDesign` takes `taxonset`, `category` (which coefficient it loads on), `includeStem`
(also assign the branch subtending the clade's ancestor) and `excludeClade` (assign only
that stem). Several elements sharing a category share one coefficient, which is how
disjoint groups of branches get a single shared effect.

## Six things that will bite you

**Anything with a cache must sit inside the posterior, or have caching turned off.** BEAST
invalidates a cache only for objects reachable from the posterior. A clock declared inside a
logger is never told the tree moved and its design silently goes stale, so hang it off a tree
likelihood; for a prior-only run, keep the likelihood but make the alignment entirely
ambiguous, which leaves it flat and the wiring correct. The same trap catches derived logger
columns: a feast `ExpCalculator` in a logger must be given `useCaching="false"` or it reports
its initial value for the whole run. Both failures are silent and produce plausible numbers.

**Constrain every clade that carries a column**, and point the constraint at the clock's own
`TaxonSet` by `idref` rather than redeclaring the taxa. Without monophyly a clade's ancestor
subsumes unrelated branches, columns overlap, and a branch in two columns carries the *sum*
of their coefficients. A second copy of the taxa lets the constraint and the design drift
apart with nothing to complain about. `DesignLogger` writes `nPainted.multiple`, which must
be 0 for disjoint clades, and a `monophyletic.<clade>` indicator per clade.

**The dispersion needs ORC's joint scaler**, `orc.consoperators.UcldScalerOperator`, not a
bare one. Changing the dispersion changes the prior density of every rate at once, so a move
touching it alone is nearly always rejected. Measured over the same 5,000,000 states:

| Operator on the dispersion | ESS | posterior sd, prior is 0.333 |
| --- | --- | --- |
| bare `BactrianScaleOperator` | 224 | 0.289, not converged |
| ORC `UcldScalerOperator` | 8875 | 0.336 |

**The branch rates need their own prior**, the same mean-one lognormal the clock declares:

```xml
<distribution id="ratesPrior" spec="beast.base.inference.distribution.Prior" x="@rates">
  <distr spec="beast.base.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">
    <parameter spec="beast.base.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>
  </distr>
</distribution>
```

It is easy to leave out precisely because the clock already names that distribution. Without
it the dispersion has nothing but its own prior to answer to, the joint scaler above has an
unbalanced Hastings ratio, and the dispersion climbs away with no error. On the six-taxon
template it reached a mean of 5.52 under a prior with mean 0.333, and took the tree height to
2.6e8 with it; with the prior in place the median is 0.218 and the tree height 0.95. Watch
also for an improper population-size prior, `OneOnX` under a flat likelihood, which lets the
tree run away by the same route.

**The coefficients need a random walk**, not a scale operator, since they are log-scale and
cross zero.

**Every parameter must start inside its prior.** An intercept left at 2e-3 under a Gamma with
mean 7e-9 makes the density underflow to NaN and BEAST cannot initialise at all.

## The dispersion is not the BEAST X quantity

`ucldStdev` here is a **log-scale standard deviation**, which is the sigma of the published
model. BEAST X's `branchRates.scale` is a **real-space coefficient of variation**. Convert:

    sigma = sqrt(log(1 + scale^2))        scale = sqrt(exp(sigma^2) - 1)

They agree to a couple of percent below about 0.3 and diverge badly above 1. Comparing the
two columns directly is wrong.

Note also that under a mixed-effects clock the dispersion is a **residual**, left after the
named clades have taken their share. A small value is evidence the design is working, not
evidence the tree is clocklike.

## Reporting the dispersion in both conventions

`examples/dispersion_conversion_6taxon.xml` logs the BEAST X quantity alongside ours, so
traces from the two codebases line up without hand conversion. The model is unchanged:
`ucldStdev` is still what is sampled and still carries its prior; the extra column is a
deterministic function of it, computed by feast.

```xml
<log id="branchRates.scale" spec="feast.expressions.ExpCalculator" useCaching="false"
     value="sqrt(exp(ucldStdev^2) - 1)">
  <arg idref="ucldStdev"/>
</log>
```

Verified against the closed form over 9001 samples, maximum absolute difference 1.4e-15.
The same column is now written by the simulation generators in `examples/`.

It also makes the prior mismatch visible rather than theoretical. Sampling the prior, the two
columns have almost the same median and very different tails:

| | median | 95% upper |
| --- | --- | --- |
| `ucldStdev`, ours | 0.227 | 1.01 |
| `branchRates.scale`, BEAST X convention | 0.230 | 1.33 |

BEAST X puts Exponential(mean 1/3) on the coefficient of variation; we put it on the
log-scale standard deviation. They agree closely where the posterior usually sits and differ
by about a third at the top, which is why this is documented rather than matched.

## What is in here

| Path | Contents |
| --- | --- |
| `JAVA_CLASSES.md` | what each class does, written for someone who has not followed the work |
| `src/mixedeffectsclock/MixedEffectsClockModel.java` | the clock; a subclass of the stock relaxed clock |
| `src/mixedeffectsclock/CladeDesign.java` | one design-matrix entry |
| `src/mixedeffectsclock/DesignLogger.java` | branch counts, monophyly indicators, stale-cache check |
| `version.xml` | package metadata **and service registration**; every public class must be listed or BEAST reports it as missing |
| `build.xml` | ant build |
| `examples/` | generators and small XMLs, one per validation stage |
| `validation/` | scoring scripts, written expectations, the figures, and `ANALYSES.md`, which records every run and the command that rebuilds it |
| `data/` | the small datasets the examples need; see `data/README.md` |
| `test/` | 15 JUnit tests over the Java, 10 over the generators |
| `skills/` | a Claude Code skill for setting the model up, and a minimal runnable template |

## The skill, and the shortest correct example

`skills/mixed-effects-clock-beast2/` is a [Claude Code](https://claude.com/claude-code)
skill covering the setup: the model and how each term maps to XML, the wiring that cannot be
skipped, the dispersion conversion, the Mascot grid arithmetic, and a triage table for the
silent failures. Install it with

    ln -s "$PWD/skills/mixed-effects-clock-beast2" ~/.claude/skills/

Its `templates/me_clock_minimal.xml` is useful with or without Claude: a complete six-taxon
analysis that runs on its own, samples the prior, and recovers each prior it declares. It is
the shortest correct example here, and the quickest check that an installation works.

## Status

Validated at the XML level against known truth, and covered by 33 unit tests; see
[Tests](#tests). The end-to-end checks came first and found every failure that actually
occurred here, all of them wiring or caching problems that a unit test on a six-taxon tree
would not have caught. The unit tests were added afterwards as regression protection.

| Check | Result |
| --- | --- |
| Design matrix vs the BEAST X implementation | reproduces 3, 25, 51, 3, 3 with 7 background, and the legacy stem configuration too |
| Design under a moving topology | correct at every state, no stale cache entries |
| Branch rates, including a slowdown and overlapping columns | exact |
| Prior recovery, all five parameters | passes, ESS near 9000 |
| 50-taxon local-clock simulation | background rate, coefficient and design all recovered |
| Random effect plus slowdown simulation | both coefficients recovered; see below |
| Compatibility with the full ORC operator suite | acceptance 0.40 to 0.48 |
| Compatibility with the targetedbeast tree moves | about 2x the mixing per state, same coverage |
| The full stack: clock, targeted moves, ORC and Mascot together | runs; no stale cache at any of 2138 samples |
| Clock estimates across three different tree priors | agree to within a percent or two |

**One known behaviour, and it is not a defect.** On a simulation with a genuine random
effect, estimating the tree and the rate together gave a background rate 21% high, a tree
15% short and a dispersion 35% low, while fitting the total substitutions per site to within
0.2%. Fixing the tree at truth recovers the rate and dispersion; fixing the rate at truth
recovers the tree and dispersion. So neither half is wrong: this is joint identifiability of
rate against time, a general property of relaxed clocks. It is one replicate, and a
ten-replicate coverage study is the outstanding piece of work.

## Tests

    ant test                      # 15 JUnit tests over the Java
    python3 test/test_generators.py   # 18 tests over the generators and the shipped template

Both suites were **mutation tested**: the code was deliberately broken and the tests were
checked to fail. That matters more than the pass count. Several times during development a
check passed while quietly testing nothing, so a green suite is only evidence once it has
been shown to go red.

| Mutation | Caught by |
| --- | --- |
| tree-dirty check removed from `requiresRecalculation` | 1 Java failure |
| clade factor dropped from `getRateForBranch` | 3 Java failures |
| skyline grid and parameter dimension drift apart | 2 generator failures |
| calendar grid keeps the `tree` attribute | 1 generator failure |
| the branch-rate prior removed from the template | 1 generator failure |
| a scale operator put on the coefficients | 1 generator failure |
| caching turned back on in the conversion column | 1 generator failure |

**Java, `test/mixedeffectsclock/`.** `DesignMatrixTest` asserts the exact *set* of branches in
each column, named by the tips below them rather than by node number. Counts are the obvious
thing to assert and are too weak: a design on entirely the wrong branches has identical
counts. It covers stems, stem-only columns, shared categories, disjointness, and pins the
behaviour of a clade that is not monophyletic. `RateAndCacheTest` checks the arithmetic
exactly, including a slowdown and a branch in two columns where the fold changes multiply,
and then the cache: that a fresh computation follows the topology, that the cached one is
detectably stale after a change, and that invalidation reconciles them.

One thing that test exposed and is worth knowing: `requiresRecalculation` reads the tree's
dirty *flag*, not the tree. Mutating nodes directly does not set it, so a test has to.

**Generators, `test/test_generators.py`.** The bugs that have actually happened in this
package are here, not in the Java. The tests assert that K levels produce K+1 grid values and
dimension K+2 together; that a one-value grid is never emitted, since it makes BEAST exit 0
and log nothing; that `--levels` passes a tree to the grid and `--boundaries` does not, which
is the entire difference between fractions and calendar times; that every flag combination
parses and carries exactly its own layers; that every clade carrying a column is constrained
by `idref` to the clock's own taxon set; that no comment contains a double hyphen; and that
the feast conversion has caching off.

A second class does the same for the template the skill ships, which is the first XML anyone
will copy: that it parses, that the clock sits inside the tree likelihood, that the branch
rates carry their own prior, that every clade is constrained, that the dispersion has ORC's
joint operator, and that the coefficients get a random walk. Each of those is a failure that
has actually happened here, and every one of them is silent.

**Still worth adding.** Store and restore across a rejected proposal. A golden test on the
50-taxon tree in `data/`, which is realistic size rather than a toy. Nothing needing an MCMC:
the operators are ORC's and BEAST's, and the end-to-end checks in `validation/` cover them
better than a unit test could.

## Not done

- The ten-replicate coverage study for the rate-against-time behaviour in [Status](#status)
  is not done. It is one replicate so far.
- Never run at real scale; every test is 5,000 or 10,000 sites against a 1.3 Mb target.
- The dispersion prior is documented rather than matched to BEAST X.
- No licence chosen. BEAST 2 and ORC are LGPL, which is the obvious candidate.
- No BEAUti integration; XMLs are written by the generators in `examples/`.
