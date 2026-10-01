#!/usr/bin/env python3
"""
sample_prior.py -- end-to-end verifier of a BEAST 2 MixedEffectsClock design.

Takes a BEAST 2 XML that uses mixedeffectsclock.MixedEffectsClockModel, redacts
every sequence to Ns, runs a short MCMC that samples from the prior, and then
calls the companion R script to colour N evenly-spaced sampled trees by their
fixed-effect column (clade).

Usage:
    sample_prior.py <analysis.xml>
                    [--states 1000000]
                    [--trees 10]
                    [--burnin 0.1]
                    [--beast beast]
                    [--out-pdf PATH]
                    [--keep]

Writes next to <analysis.xml>:
    <stem>_prior.xml      the redacted, chain-shortened XML
    <stem>_prior.log      BEAST trace log
    <stem>_prior.trees    sampled trees
    <stem>_prior.pdf      the design-check PDF (one page per sampled tree)

Why all-N sequences, not sampleFromPrior=true: with the MixedEffectsClock stack
the clock sits inside the posterior through the TreeLikelihood, and ripping the
likelihood out breaks the invalidation chain that keeps the design cache fresh.
The sibling skill mixed-effects-clock-beast2 documents this as point 1 of its
"Six things that will bite you". The all-N trick keeps the wiring intact and
leaves the likelihood flat so topology moves freely.
"""
from __future__ import annotations
import argparse, os, re, shutil, subprocess, sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
R_SCRIPT   = SCRIPT_DIR / "colour_prior_trees.R"


def redact_and_rechain(xml_text: str, chain_length: int, tree_log_every: int) -> str:
    """Replace every <sequence ... value="..."> with Ns, and rewrite the MCMC
    chain length and tree-logger logEvery."""

    # 1. redact sequences. Preserve the original length so totalcount/weight stay
    #    consistent. BEAST 2 syntax: <sequence ... value="ACGT..." .../>
    def _redact(m: re.Match) -> str:
        return 'value="' + "N" * len(m.group(1)) + '"'
    xml_text, n_seq = re.subn(r'value="([ACGTUMRWSYKVHDBN\-\?acgtumrwsykvhdbn\-\?]+)"',
                              _redact, xml_text)
    if n_seq == 0:
        print("WARN: no <sequence ... value='...'> blocks matched; nothing redacted.",
              file=sys.stderr)
    else:
        print(f"redacted {n_seq} sequences to Ns")

    # 2. chainLength on the top-level MCMC run element
    xml_text, n_run = re.subn(
        r'(<run[^>]*spec="beast\.base\.inference\.MCMC"[^>]*chainLength=")\d+(")',
        lambda m: f"{m.group(1)}{chain_length}{m.group(2)}", xml_text, count=1)
    if n_run == 0:
        # chainLength may come before the spec attribute; try the other order
        xml_text, n_run = re.subn(
            r'(<run[^>]*chainLength=")\d+("[^>]*spec="beast\.base\.inference\.MCMC")',
            lambda m: f"{m.group(1)}{chain_length}{m.group(2)}", xml_text, count=1)
    if n_run == 0:
        sys.exit("ERROR: could not find <run ... spec='...MCMC' chainLength='...'>")
    print(f"set chainLength = {chain_length}")

    # 3. logEvery on every tree logger (mode="tree"). Rewrite them all.
    def _retree(m: re.Match) -> str:
        head, rest = m.group(1), m.group(2)
        # swap logEvery="..."; mode="tree" was already matched in the enclosing group
        rest2 = re.sub(r'logEvery="\d+"', f'logEvery="{tree_log_every}"', rest)
        return head + rest2
    xml_text, n_tree = re.subn(r'(<logger\b)([^>]*mode="tree"[^>]*>)',
                                _retree, xml_text)
    if n_tree == 0:
        sys.exit("ERROR: no tree logger (<logger ... mode='tree'>) found")
    print(f"set treelog logEvery = {tree_log_every} ({n_tree} tree logger(s) rewritten)")

    # 4. also keep the trace logger from drowning the disk: logEvery = tree_log_every.
    def _retrace(m: re.Match) -> str:
        head, rest = m.group(1), m.group(2)
        if 'mode="tree"' in rest:
            return head + rest
        rest2 = re.sub(r'logEvery="\d+"', f'logEvery="{tree_log_every}"', rest)
        return head + rest2
    xml_text, _ = re.subn(r'(<logger\b)([^>]*>)', _retrace, xml_text)

    return xml_text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xml", help="BEAST 2 analysis XML using MixedEffectsClockModel")
    ap.add_argument("--states", type=int, default=1_000_000,
                    help="MCMC chain length (default 1,000,000)")
    ap.add_argument("--trees",  type=int, default=10,
                    help="number of trees to draw in the PDF (default 10)")
    ap.add_argument("--burnin", type=float, default=0.1,
                    help="fraction of trees discarded as burn-in before picking (default 0.1)")
    ap.add_argument("--beast",  default="beast",
                    help="BEAST 2 executable on PATH (default 'beast')")
    ap.add_argument("--out-pdf", default=None,
                    help="override the output PDF path")
    ap.add_argument("--keep",   action="store_true",
                    help="keep _prior.log and _prior.trees (default: keep them anyway; this is a no-op kept for symmetry)")
    args = ap.parse_args()

    xml_in = Path(args.xml).resolve()
    if not xml_in.is_file():
        sys.exit(f"not a file: {xml_in}")
    stem   = xml_in.with_suffix("").name
    outdir = xml_in.parent
    prior_xml   = outdir / f"{stem}_prior.xml"
    prior_trees = outdir / f"{stem}_prior.trees"
    prior_pdf   = Path(args.out_pdf).resolve() if args.out_pdf else outdir / f"{stem}_prior.pdf"

    src = xml_in.read_text()

    # Fail loud if this is not a MixedEffectsClock XML.
    if "mixedeffectsclock.MixedEffectsClockModel" not in src:
        sys.exit("ERROR: no mixedeffectsclock.MixedEffectsClockModel block in "
                 f"{xml_in}. This skill is for BEAST 2 MixedEffectsClock XMLs only; "
                 "BEAST X <fixedEffects> is a different stack.")
    if "<fixedEffects" in src:
        sys.exit("ERROR: this looks like a BEAST X XML (<fixedEffects> block). "
                 "This skill is BEAST 2 only.")

    # tree_log_every so that chainLength / logEvery = args.trees / (1 - burnin).
    # i.e. we log enough trees that after discarding burnin we still have >= args.trees.
    keep_fraction = max(0.0, 1.0 - args.burnin)
    n_logged = max(args.trees + 1, int(round(args.trees / keep_fraction)) + 1)
    tree_log_every = max(1, args.states // n_logged)

    prior_xml.write_text(redact_and_rechain(src, args.states, tree_log_every))
    print(f"wrote {prior_xml}")

    # Run BEAST 2. Use -overwrite so a stale .log/.trees from a previous attempt
    # does not stop it; -D sets $(filebase) to our stem.
    cmd = [args.beast, "-overwrite", "-seed", "42", prior_xml.name]
    print(f"\n>>> {' '.join(cmd)}   (cwd={outdir})\n")
    r = subprocess.run(cmd, cwd=outdir)
    if r.returncode != 0:
        sys.exit(f"BEAST exited {r.returncode}")

    if not prior_trees.is_file():
        sys.exit(f"BEAST ran but did not produce {prior_trees}")

    # Hand off to the R script.
    rcmd = ["Rscript", str(R_SCRIPT), str(prior_xml), str(prior_trees),
            str(prior_pdf), str(args.trees), str(args.burnin)]
    print(f"\n>>> {' '.join(rcmd)}\n")
    r = subprocess.run(rcmd)
    if r.returncode != 0:
        sys.exit(f"R script exited {r.returncode}")
    print(f"\nwrote {prior_pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
