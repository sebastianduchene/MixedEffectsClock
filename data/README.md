# Data shipped with the package

Small datasets the example generators need, so the repository is self-contained.

| File | What it is |
| --- | --- |
| `sim_local_MEclock_Bletsa.xml` | The BEAST X analysis of the 50-taxon local-clock simulation. `../examples/gen_sim50_beast2.py` reads the taxa, dates, sequences, clade sets and priors out of it, so the port is like for like rather than retyped. |
| `simulate_re_slowdown.R` | Generates the random-effect plus slowdown simulation. Needs R with NELSI, ape and phangorn. |
| `sim_re_slowdown.fasta` | Its output: 50 taxa, 10,000 sites, Jukes-Cantor. |
| `sim_re_slowdown_true.nwk` | The simulating tree, used by the fixed-tree diagnostic. |
| `sim_re_slowdown_truth.csv` | Per-column truth: branches, fold change, coefficient, expected substitutions. |
| `sim_re_slowdown_taxa.tsv` | The three clade taxon sets. |

**Not here: the lepromatosis data.** `gen_design_check_lepromatosis.py` and
`gen_design_check_moving.py` read the real analysis from `../../lepromatosis`, which is
unpublished and whose alignment is 61 MB. Those two generators fail with a clear message
when it is absent; every other example runs from this folder alone.
