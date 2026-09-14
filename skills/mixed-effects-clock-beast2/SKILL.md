---
name: mixed-effects-clock-beast2
description: Build a mixed-effects (fixed + random effects) molecular clock in BEAST 2 with the MixedEffectsClock package, where named clades get their own log fold change on top of a background rate and every branch also carries a lognormal random effect - the Bletsa et al. 2019 model, log r_i = beta_0 + sum_k X_ik beta_k + eps_i. Use when the user wants a local clock expressed as background plus per-clade fold changes in BEAST 2, wants to port a BEAST X arbitraryBranchRates / signTransform / fixedEffects analysis into BEAST 2, needs this clock combined with Mascot, ORC or the targeted tree operators, asks what ucldStdev means under a mixed-effects clock or how it relates to the BEAST X branchRates.scale, or hits a wrong, stale or silently expanding design matrix. Also covers the Mascot skyline grid arithmetic that this stack depends on.
metadata:
  type: reference
  engine: BEAST 2.7.7+, MixedEffectsClock v0.0.1
  repository: https://github.com/sebastianduchene/MixedEffectsClock
  template: templates/me_clock_minimal.xml
---

# A mixed-effects clock in BEAST 2

A **mixed-effects (ME) clock** gives named clades their own rate on top of a background
rate, and additionally lets every branch wander. BEAST X has had this for years; BEAST 2
gets it from the `MixedEffectsClock` package, which is a subclass of the stock relaxed
clock rather than a clock written from scratch.

```
log r_i  =  beta_0  +  sum_k X_ik beta_k  +  eps_i          eps_i ~ Normal(0, sigma^2)

r_i      =  exp(beta_0) * prod_k exp(beta_k)^X_ik * exp(eps_i)
```

- `exp(beta_0)` — the background rate. **It is the inherited `clock.rate`.** There is no
  separate intercept parameter, which is what keeps the up-down operators meaningful.
- `beta_k` — a **log fold change**, carried by design column *k*. `X_ik` is 1 when branch
  *i* is in that column. Unbounded, so a clade can slow down as well as speed up.
- `exp(eps_i)` — the inherited per-branch `rates`, drawn from a lognormal pinned to mean
  one in real space. Its log-scale SD is `ucldStdev`. Drive it to zero and the model
  collapses to a plain local clock.

This is the **log-scale (Bletsa) parameterisation only**. The package does not implement
BEAST X's real-space additive form `r_i = (beta_0 + sum_k beta_k z_ik) * exp(eps_i)`, the
one you get without `<signTransform/>`. If the analysis being ported is in the real-space
form, the coefficients do not transfer.

`templates/me_clock_minimal.xml` is a complete, runnable six-taxon skeleton with every
one of these pieces wired correctly. It samples the prior, so running it as-is checks the
installation and recovers each prior; the header lists what to replace to make it a real
analysis.

## Setup

Install BEAST 2.7.7 or newer, then these packages through BEAUti's package manager or
`packagemanager -add`:

| Package | Needed for |
| --- | --- |
| **ORC** | the dispersion operator. **Not optional**, see below |
| **feast** | the dispersion conversion column |
| **BEASTLabs** | a metric ORC's exchange operators use |
| Mascot | only if the tree prior is the structured coalescent |
| TargetedBeast | only if using the targeted tree moves |

Then the package itself, either from a release zip unpacked into
`~/.beast/2.7/MixedEffectsClock`, or from source:

```bash
git clone https://github.com/sebastianduchene/MixedEffectsClock
cd MixedEffectsClock
ant install          # needs JDK 17+; a Java 11 compiler fails with "class file has wrong version 61.0"
```

Confirm it loaded: BEAST prints its package list at startup and `MixedEffectsClock`
should be in it. There is **no BEAUti integration** — XMLs are written by hand or by the
generators in `examples/`.

Declare the dependency in the XML header:

```xml
<beast required="BEAST.base v2.7.7:MixedEffectsClock v0.0.1:ORC v1.2.1:feast v10.6.1"
       version="2.7">
```

## The clock block

```xml
<branchRateModel id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"
                 tree="@Tree" rates="@rates" clock.rate="@clockRate"
                 coefficient="@coefficient">

  <distr id="rateDistr" spec="beast.base.inference.distribution.LogNormalDistributionModel"
         S="@ucldStdev" meanInRealSpace="true">
    <parameter spec="beast.base.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>
  </distr>

  <clade id="cladeFast" spec="mixedeffectsclock.CladeDesign" category="0" includeStem="true">
    <taxonset id="taxa.fastClade" spec="beast.base.evolution.alignment.TaxonSet">
      <taxon idref="t1"/> <taxon idref="t2"/>
    </taxonset>
  </clade>

  <clade id="cladeSlow" spec="mixedeffectsclock.CladeDesign" category="1" includeStem="true">
    <taxonset id="taxa.slowClade" spec="beast.base.evolution.alignment.TaxonSet">
      <taxon idref="t5"/> <taxon idref="t6"/>
    </taxonset>
  </clade>
</branchRateModel>
```

`coefficient` is a `RealParameter` whose `dimension` is the number of **distinct
categories**, not the number of `<clade>` elements. `rates` has dimension 2n-2.

### CladeDesign

| Input | Default | Meaning |
| --- | --- | --- |
| `taxonset` | required | the taxa whose MRCA defines the clade |
| `category` | 0 | which coefficient this column loads on |
| `includeStem` | false | also paint the branch subtending the clade's ancestor |
| `excludeClade` | false | paint **only** that stem. Requires `includeStem` |

Several elements sharing a `category` share one coefficient. That is how a clade plus a
disjoint sister stem get a single shared effect, which is what BEAST X's `<category>`
grouping does. A branch that lands in two columns carries the **sum** of their
coefficients, so the fold changes multiply.

## Six things that will bite you

**1. The clock must sit inside the posterior.** BEAST invalidates a cache only for objects
reachable from the posterior, and the design matrix depends on the topology. A clock
declared inside a logger is never told the tree moved, so its design goes quietly stale
and every number after that is wrong with no error. Hang it off a `TreeLikelihood`. For a
prior-only run keep the likelihood and make the alignment entirely ambiguous (all `N`),
which leaves it flat and the wiring correct. Do not use random sequences for this: with
real random data the likelihood pins the topology and the tree never moves.

**2. Constrain every clade that carries a column, by `idref` to the clock's own TaxonSet.**

```xml
<distribution id="monophyly.fastClade" spec="beast.base.evolution.tree.MRCAPrior"
              monophyletic="true" tree="@Tree">
  <taxonset idref="taxa.fastClade"/>
</distribution>
```

Without monophyly the clade's MRCA subsumes unrelated branches, columns overlap and the
design silently expands. A *second copy* of the taxa lets the constraint and the design
drift apart with nothing to complain about, so reference the same object.

**3. The dispersion needs ORC's joint scaler**, not a bare one. Changing `ucldStdev`
changes the prior density of every branch rate at once, so a move touching it alone is
almost always rejected. Over the same 5,000,000 states:

| Operator on the dispersion | ESS |
| --- | --- |
| bare `BactrianScaleOperator` | 224, not converged |
| `orc.consoperators.UcldScalerOperator` | 8875 |

```xml
<operator id="ucldScalerJoint" spec="orc.consoperators.UcldScalerOperator"
          distr="@rateDistr" rates="@rates" stdev="@ucldStdev"
          scaleFactor="0.5" weight="3.0">
  <kernel spec="beast.base.inference.operator.kernel.KernelDistribution$Bactrian"/>
</operator>
```

**4. Coefficients need a random walk, never a scale operator.** They are log-scale and
cross zero. Use `beast.base.inference.operator.kernel.BactrianRandomWalkOperator`.

**5. The branch rates need their own prior**, the same mean-one lognormal the clock
uses. It is easy to leave out, because the clock already declares the distribution:

```xml
<distribution id="ratesPrior" spec="beast.base.inference.distribution.Prior" x="@rates">
  <distr spec="beast.base.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">
    <parameter spec="beast.base.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>
  </distr>
</distribution>
```

Without it `ucldStdev` has nothing but its own prior to answer to, ORC's joint move has an
unbalanced Hastings ratio, and the dispersion climbs away with no error. On the six-taxon
template this put `ucldStdev` at 5.5 under a prior with mean 0.333, and took the tree to
2.6e8 with it.

**6. Every parameter must start inside its prior.** A clock rate left at 2e-3 under a
Gamma with mean 7e-9 makes the density underflow to NaN and BEAST cannot initialise at
all. Start at the prior mean.

The Bactrian operator classes are split across two packages and this costs time
repeatedly:

| Class | Package |
| --- | --- |
| `BactrianScaleOperator`, `BactrianNodeOperator`, `BactrianSubtreeSlide` | `beast.base.evolution.operator.kernel` |
| `BactrianRandomWalkOperator`, `BactrianUpDownOperator`, `BactrianDeltaExchangeOperator`, `KernelDistribution$Bactrian` | `beast.base.inference.operator.kernel` |

## Checking the design before trusting anything

Add the design logger to the trace log:

```xml
<log id="designLog" spec="mixedeffectsclock.DesignLogger" validate="false" clock="@clock"/>
```

It writes `nPainted.X1..XK`, `nPainted.background`, `nPainted.multiple`, and a
`monophyletic.<clade>` indicator per clade. With `validate="true"` it also recomputes the
design from scratch at every sample and writes `design.staleEntries`.

Read these three columns first, before any parameter:

| Column | Must be |
| --- | --- |
| `nPainted.multiple` | 0, whenever the clades are disjoint |
| `monophyletic.<clade>` | 1, for every clade |
| `design.staleEntries` | 0 at every sample, if validation is on |

Any of those wrong means the design is not what the XML says, and the rest of the output
is not worth reading. Per-column counts alone **cannot** detect a stale cache: the same
number of entries stay marked, they just point at the wrong branches. That is what
`design.staleEntries` is for. Turn validation on for a short run, then off.

`examples/design_check_6taxon.xml` in the repository is a six-taxon tree whose design can
be counted by hand, and is the quickest way to confirm an installation works.

## The dispersion is not the BEAST X quantity

`ucldStdev` here is a **log-scale standard deviation**, the sigma of the published model.
BEAST X's `branchRates.scale` is a **real-space coefficient of variation**. Convert:

```
sigma = sqrt(log(1 + scale^2))          scale = sqrt(exp(sigma^2) - 1)
```

They agree to a couple of percent below about 0.3 and diverge badly above 1: sigma 1.0 is
a scale of 1.31, sigma 2.0 is 7.32. Comparing the two columns directly is wrong.

To log both, so traces from the two codebases line up without hand conversion:

```xml
<log id="branchRates.scale" spec="feast.expressions.ExpCalculator" useCaching="false"
     value="sqrt(exp(ucldStdev^2) - 1)">
  <arg idref="ucldStdev"/>
</log>
```

**`useCaching="false"` is mandatory.** A feast `ExpCalculator` in a logger is not
reachable from the posterior, so without it the column reports its initial value for the
entire run. The failure is silent and the number is plausible.

The priors are not matched either. BEAST X puts Exponential(mean 1/3) on the coefficient
of variation; this package puts it on the log-scale SD. Sampling the prior, the medians
agree (0.227 vs 0.230) and the 95% uppers do not (1.01 vs 1.33).

Interpretation differs too: under a mixed-effects clock the dispersion is a **residual**,
what is left after the named clades take their share. A small value is evidence the
design is working, not evidence the tree is clocklike.

## Combining with Mascot, ORC and the targeted operators

This is what the package exists for. ORC's moves and the targetedbeast moves both act on
the **raw** `rates` vector and node heights, exactly as the stock operators do, so the
per-branch clade factor is a constant they never touch. Nothing special is needed to
combine them beyond declaring the packages.

With targetedbeast, put the edge weights **inside the posterior** — they are a function of
the tree and follow the same invalidation rule as the clock:

```xml
<distribution id="consensusWeights" spec="targetedbeast.edgeweights.ParsimonyWeights"
              tree="@Tree" maxWeight="20">
  <data spec="targetedbeast.alignment.ConsensusAlignment" data="@alignment"/>
</distribution>
```

Use `ParsimonyWeights`, not `ParsimonyWeights2`, which throws
`ArrayIndexOutOfBoundsException` on real-sized data.

### The Mascot skyline grid arithmetic

**Do not set the grid by hand.** The arithmetic is not what the XML suggests.
`StructuredMigrationSkyline` caps its interval index two below the number of shift values,
so N values give N-1 levels; and `Skygrowth` separately forces its parameter to dimension
N+1, with the last entry never read by Mascot.

> For **K** rate levels: **K+1** shift values, and every `logNe` parameter at
> **dimension K+2**, of which K+1 are live.

A **single** shift value fails silently: the run exits 0, logs nothing, and never produces
a sample. The minimum is two values, which is one level.

Whether the values are fractions or calendar times is decided by one attribute:

| Form | `rateShifts` element | Values mean |
| --- | --- | --- |
| tree-relative | **with** `tree="@Tree"` | fractions of the current root height; boundaries move with the tree |
| calendar | **without** a `tree` | absolute times before the most recent tip |

An absolute grid that does not reach the root collapses the skyline to a constant with no
warning. Share **one** `RateShifts` object between the dynamics and every `Skygrowth`;
separate grids do the same thing.

For a ghost (unsampled) deme, name both types in `types=` — an unsampled type cannot be
read off the type trait — and leave `fromBeauti` at its default false, otherwise the
`types` string is ignored. Migration index order is fixed by the class: 0 is
outbreak-to-ghost, 1 is ghost-to-outbreak, both **forwards** in time.

```xml
<dynamics id="StructuredMigrationSkyline" spec="mascot.dynamics.StructuredMigrationSkyline"
          types="outbreak ghost">
  <NeDynamics spec="mascot.util.InitializedNeDynamicsList">
    <neDynamics spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineNe.outbreak" rateShifts="@rateShifts"/>
    <neDynamics spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineNe.ghost"    rateShifts="@rateShifts"/>
  </NeDynamics>
  <migrationDynamics spec="mascot.util.InitializedNeDynamicsList">
    <neDynamics spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineMig.outbreak_to_ghost" rateShifts="@rateShifts"/>
    <neDynamics spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineMig.ghost_to_outbreak" rateShifts="@rateShifts"/>
  </migrationDynamics>
  <rateShifts id="rateShifts" spec="mascot.dynamics.RateShifts">2.5 5 7.5</rateShifts>
  ...
</dynamics>
```

`NeDynamicsList` is not registered as a service; use `mascot.util.InitializedNeDynamicsList`.

## Generating an XML instead of writing one

The repository ships a generator that produces every configuration from one script:

```bash
cd examples
python3 gen_sim_re_slowdown.py                                  > plain.xml
python3 gen_sim_re_slowdown.py --targeted                       > targeted.xml
python3 gen_sim_re_slowdown.py --targeted --mascot --levels 2   > stack_relative.xml
python3 gen_sim_re_slowdown.py --targeted --mascot --boundaries "2.5" > stack_calendar.xml
```

`--levels K` gives a tree-relative grid with K levels; `--boundaries "a b c"` gives a
calendar grid at those times. Both derive the shift values and the parameter dimensions
together, which is the point of using them.

It writes a 50-taxon simulation with known truth: background rate 0.005, one clade at
twice that, one at 0.3 times, and a genuine per-branch random effect. Score a finished run
against the truth with `python3 ../validation/check_sim_re.py <basename>`.

## Interpreting the result

Exponentiate a coefficient to read it: `exp(beta_k)` is the fold change for that clade
relative to background, and the same transform applied to the 95% HPD bounds gives its
interval. A 95% interval on `beta_k` that excludes zero is the clade being distinguishable
from background.

**One known behaviour, and it is not a defect.** Estimating the tree and the background
rate together on a simulation with a real random effect gave a rate 21% high, a tree 15%
short and a dispersion 35% low, while fitting total substitutions per site to within 0.2%.
Fixing the tree at truth recovers the rate and dispersion; fixing the rate recovers the
tree and dispersion. Both fixed effects covered in every configuration. This is joint
identifiability of rate against time, a general property of relaxed clocks, not something
this model introduces. If the dates or the rate are weakly informed, expect it.

## Quick triage

| Symptom | Cause |
| --- | --- |
| Column reports one constant value for the whole run | something cacheable sits outside the posterior: the clock in a logger, or `ExpCalculator` without `useCaching="false"` |
| `nPainted.multiple` above 0 with disjoint clades | a clade is not constrained monophyletic, so its MRCA subsumes other branches |
| `design.staleEntries` above 0 | the design cache is not tracking the topology; per-column counts will still look right |
| Run exits 0, log file has a header and no samples | a one-value Mascot rate-shift grid |
| Skyline posterior flat, all intervals identical | absolute grid that does not reach the root, or separate `RateShifts` objects |
| Dispersion ESS in the low hundreds | a bare scale operator on `ucldStdev` instead of ORC's joint scaler |
| Dispersion climbs far above its prior, tree length with it | no `Prior` on `rates`, so the joint move is unbalanced |
| `class file has wrong version 61.0` | building with a Java 11 compiler; needs JDK 17+ |
| BEAST reports a class as missing | it is not registered in `version.xml` as a `beast.base.core.BEASTInterface` service |
| Cannot initialise, NaN density at state 0 | a parameter starts outside its prior |
| Element name not recognised, e.g. `<prior>` or `<Normal>` | no `<map>` declaration; write `spec="..."` in full instead |
| XML will not parse, no useful message | a comment contains a double hyphen, which is illegal |
