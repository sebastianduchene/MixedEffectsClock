# Every analysis run here, and the command that rebuilds it

The XMLs themselves are not committed. Each is about half a megabyte because it inlines the
alignment, they come to 9 MB together, and every one is reproducible from a single command.
This file is the record; the generators in `../examples/` are the source of truth.

Run any of these from `examples/`, redirecting to a file, then score with the checker named
in the last column.

## The 50-taxon local-clock simulation, stage 3a

| Analysis | Command | What it tested |
| --- | --- | --- |
| `sim50_ME_beast2` | `gen_sim50_beast2.py` | port of the BEAST X build; recovered rate, coefficient and design |
| `sim50_ORC_test` | as above, plus the ORC operator block by hand | ORC's full suite is compatible; acceptance 0.40 to 0.48 |

## The random-effect and slowdown simulation

Truth: background 0.005, fast clade x2, slow clade x0.3, dispersion 0.2816 realised.

| Analysis | Command | What it tested |
| --- | --- | --- |
| `sim_re_slowdown` | `gen_sim_re_slowdown.py` | both coefficients recover, including the slowdown; rate and dispersion miss |
| `sim_re_JC` | as above with HKY+gamma replaced by JC | rules out substitution mis-specification as the cause |
| `sim_re_fixedtree` | as above, tree pinned at the simulating tree, tree operators removed | rate and dispersion both recover: the clock is sound |
| `sim_re_fixedrate` | as above, rate pinned at 0.005, its two operators removed | tree and dispersion recover: neither half is wrong |
| `sim_re_targeted` | `gen_sim_re_slowdown.py --targeted` | targetedbeast moves; about 2x the mixing per state, same coverage |

The fixed-tree and fixed-rate runs are hand-edited from the generated XML rather than
generated, since they exist to remove things. The edits are: delete the tree operators and
replace the starting tree with `../data/sim_re_slowdown_true.nwk`; or delete
`clockRateScaler` and `upDownRates` and set the rate to 0.005.

## The full stack

| Analysis | Command | What it tested |
| --- | --- | --- |
| `sim_full_stack` | `gen_sim_re_slowdown.py --targeted --mascot --levels 3` | all four layers compose; nothing stale |
| `sim_full_stack_1shift` | `... --targeted --mascot --levels 2` | one change point; dimensions follow the level count |
| `stack_cal_run` | `... --targeted --mascot --boundaries "2.5"` | calendar boundaries stay fixed while the root moves |

`stack_cal_run` is the one to look at first. 4.2M states, `design.staleEntries` 0 at all 2138
samples, boundaries pinned at 2.5, 5 and 7.5, design 4/19/75 throughout, and clock estimates
within a percent or two of the runs with a completely different tree prior.

## Scratch

`rs1`–`rs3`, `rs1chk`, `lv1`–`lv3`, `feast_conv_test` and `sim_re_default` were one-off
probes: the rate-shift grid arithmetic, the level-count option, and the feast caching
behaviour. Their findings are in the README and the lab notebook; the files are not kept.

## Checkers

| Script | For |
| --- | --- |
| `check_sim_re.py` | the random-effect and slowdown runs, and the full stack |
| `check_sim50.py` | the 50-taxon local clock |
| `check_prior_recovery.py` | prior-only runs, against analytic priors |
| `check_rates.py` | independent per-branch rate recomputation with dendropy |
| `plot_rate_tree.R`, `plot_rate_trees.R` | the design-coloured tree figures |
