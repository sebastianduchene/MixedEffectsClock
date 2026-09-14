# Stage 3a, 50-taxon simulation: what we expect before looking

Written before the run, so the result can be judged against a stated prediction rather
than rationalised afterwards. Truth comes from `../../simulate_clocks.R`: four branches
were given a rate of 0.01 against a background of 0.005, with **no** per-branch random
variation at all, and the sequences were simulated under Jukes-Cantor.

## The three quantities with a known answer

| Parameter | Truth | Why |
| --- | --- | --- |
| `clockRate` | 0.005 | the background rate, which is `exp(beta_0)` |
| `coefficient` | 0.6931 | `log 2`, a two-fold increase on the four fast branches |
| `ucldStdev` | 0 | the generating process is a pure local clock, no random effect |

A pass means the 95% interval covers each of these. The coefficient is the one that
matters: it is the whole point of the model, and if it sits near 0 then either the design
is not reaching the fast branches or the effect is not identifiable from 5000 sites.

The dispersion cannot be "covered" in the usual sense, because 0 is on the boundary. What
we want is a posterior piled up against zero, with the bulk well below the prior's mean of
0.333. If it comes back near 0.3 it is simply following its prior and the data say nothing.

## The design, which is checkable exactly

The fast set is the two-taxon clade `t23_3.36`, `t9_2.75` *with* its stem, plus the stem
alone of its six-taxon sister. Both load on one shared column, as in BEAST X.

| Column | Expected branches |
| --- | --- |
| `nPainted.X1` | 4 |
| `nPainted.background` | 94 |
| `nPainted.multiple` | 0 |
| both monophyly indicators | 1 at every state |

98 branches for 50 taxa. If `nPainted.X1` is not 4, nothing else in the run is worth
reading.

## Deliberately mis-specified, and expected to show it

The data are Jukes-Cantor but the analysis fits HKY plus gamma, exactly as the BEAST X
run does. So:

- `kappa` should sit near **1**, its true value, against a prior whose median is
  `exp(1) = 2.7`. Seeing kappa near 2.7 would mean the data are not informing it.
- `gammaShape` has no true value, since there is no rate heterogeneity. Its prior,
  Exponential(mean 0.5), pulls hard toward small shapes, which is toward *more*
  heterogeneity. Expect a tug-of-war and a shape that is poorly determined. This is
  inherited mis-specification, not a fault of the port.
- `freqs` should be near 0.25 each.

## Interpretation traps to avoid

**A two-fold error in the rate would look like success on one parameter and failure on
another.** The intercept and the coefficient are not separately identifiable from branch
lengths alone: only the product of rate and time is. If the tree comes out half as deep
and the rate twice as high, `clockRate` would be 0.01 and the coefficient could still be
right. Check the root height too.

**The dispersion is a residual, not the rate variation.** Here the truth is zero because
the fixed effect explains everything. A small value is evidence the design is doing its
job.

## Mixing, and what is realistic in five minutes

The clock review measured this same model in BEAST X at 93 effective samples of the
dispersion per million states under Metropolis. So in a few minutes we should expect
effective sample sizes in the tens to low hundreds, not thousands. The question after five
minutes is whether the chain has reached stationarity and is heading to the right place,
not whether it is finished.
