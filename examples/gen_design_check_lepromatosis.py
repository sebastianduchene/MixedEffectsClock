#!/usr/bin/env python3
"""
Stage 1 of MIXED_EFFECTS_BEAST2_PLAN.md on the real data.

Emits a BEAST 2 XML that builds nothing but the tree, the clade taxon sets and the
mixed-effects design, then logs the per-column branch counts. No alignment is needed:
the design depends only on the topology and the taxon sets.

The counts are checked against the BEAST X implementation, which wrote the same
quantity as n_painted_X1..X5. Two references are known for the 47-taxon tree:

  stems  hyper chile Argentina only : 3, 24, 50, 3, 3   background 9
  stems  + pddc and squirrels       : 3, 25, 51, 3, 3   background 7

The second is the corrected design of 2026-09-10.

    python3 gen_design_check_lepromatosis.py [--legacy-stems] > out.xml
"""
import re, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
LEP  = os.path.join(HERE, "..", "..", "lepromatosis")

if not os.path.isdir(LEP):
    sys.exit("This example needs the lepromatosis analysis, which is not part of this\n"
             "repository: the data are unpublished and the alignment is 61 MB. Expected\n"
             "to find it at:\n    %s\n"
             "Every other example in this folder runs from ../data alone." % LEP)
ORIG = os.path.join(LEP, "MSC_fakeOutgroup_posterior_fakeAr_FINAL_v7_expGrowth_"
                         "pddc_squirrels_ME_bletsa_hmc_b0eqIngroup_sd2_prior.xml")
TREE = os.path.join(LEP, "starting_tree_dated.nwk")

legacy = "--legacy-stems" in sys.argv

txt = open(ORIG, encoding="utf-8", errors="replace").read()

def taxset(name):
    m = re.search(r'<taxa id="%s">(.*?)</taxa>' % re.escape(name), txt, re.S)
    return re.findall(r'<taxon idref="([^"]+)"', m.group(1))

CLADES = ["hyper", "pddc", "squirrels", "chile", "Argentina"]
sets = {c: taxset(c) for c in CLADES}
all_taxa = taxset("ingroup") + ["fake_outgroup"]
all_taxa = sorted(set(all_taxa) | {t for c in CLADES for t in sets[c]})

# includeStem as the BEAST X file has it. pddc and squirrels gained their stems on
# 2026-09-10; --legacy-stems reproduces the design from before that.
STEM = {"hyper": True, "chile": True, "Argentina": True,
        "pddc": not legacy, "squirrels": not legacy}

newick = open(TREE, encoding="utf-8").read().strip()
nbranch = 2 * len(all_taxa) - 2
sys.stderr.write("taxa=%d  branches=%d  stems=%s\n"
                 % (len(all_taxa), nbranch, {k: v for k, v in STEM.items() if v}))

P = "beast.base"
L = []; w = L.append
w('<beast namespace="%s.evolution.alignment:%s.evolution.tree:%s.inference"' % (P, P, P))
w('       required="BEAST.base v2.7.7:MixedEffectsClock v0.0.1" version="2.7">')
w('  <run id="mcmc" spec="%s.inference.MCMC" chainLength="0">' % P)
w('    <state id="state" spec="%s.inference.State">' % P)
w('      <tree id="Tree" spec="%s.evolution.tree.Tree" name="stateNode">' % P)
w('        <taxonset id="TaxonSet" spec="%s.evolution.alignment.TaxonSet">' % P)
for t in all_taxa:
    w('          <taxon id="%s" spec="%s.evolution.alignment.Taxon"/>' % (t, P))
w('        </taxonset>')
w('      </tree>')
w('      <parameter id="coefficient" spec="%s.inference.parameter.RealParameter" dimension="%d" name="stateNode">0.0</parameter>' % (P, len(CLADES)))
w('      <parameter id="rates" spec="%s.inference.parameter.RealParameter" dimension="%d" lower="0.0" name="stateNode">1.0</parameter>' % (P, nbranch))
w('      <parameter id="ucldStdev" spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">0.3</parameter>' % P)
w('      <parameter id="clockRate" spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">7.0E-9</parameter>' % P)
w('    </state>')
w('    <init id="startingTree" spec="%s.evolution.tree.TreeParser" initial="@Tree"' % P)
w('          taxonset="@TaxonSet" IsLabelledNewick="true" adjustTipHeights="false"')
w('          newick="%s"/>' % newick)
w('    <distribution id="posterior" spec="%s.inference.CompoundDistribution">' % P)
w('      <distribution id="coefficientPrior" spec="%s.inference.distribution.Prior" x="@coefficient">' % P)
w('        <distr spec="%s.inference.distribution.Normal">' % P)
w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">0.0</parameter>' % P)
w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="sigma">2.0</parameter>' % P)
w('        </distr>')
w('      </distribution>')
w('    </distribution>')
w('    <operator id="dummy" spec="%s.inference.operator.kernel.BactrianRandomWalkOperator" parameter="@coefficient" scaleFactor="0.1" weight="1.0"/>' % P)
w('    <logger id="tracelog" spec="%s.inference.Logger" fileName="$(filebase).log" logEvery="1">' % P)
w('      <log idref="posterior"/>')
w('      <log id="designLog" spec="mixedeffectsclock.DesignLogger">')
w('        <clock id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"')
w('               tree="@Tree" rates="@rates" clock.rate="@clockRate" coefficient="@coefficient">')
w('          <distr id="rateDistr" spec="%s.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">' % P)
w('            <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>' % P)
w('          </distr>')
for i, c in enumerate(CLADES):
    w('          <clade id="clade.%s" spec="mixedeffectsclock.CladeDesign" category="%d" includeStem="%s">'
      % (c, i, "true" if STEM[c] else "false"))
    w('            <taxonset id="taxa.%s" spec="%s.evolution.alignment.TaxonSet">' % (c, P))
    for t in sets[c]:
        w('              <taxon idref="%s"/>' % t)
    w('            </taxonset>')
    w('          </clade>')
w('        </clock>')
w('      </log>')
w('    </logger>')
w('  </run>')
w('</beast>')
print("\n".join(L))
