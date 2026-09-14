#!/usr/bin/env python3
"""
Design matrix under a MOVING topology, the check that state-0 tests cannot do.

A constant-size coalescent over the real 47-taxon lepromatosis tree with topology
operators active, plus a SMALL DUMMY ALIGNMENT.

The alignment is not optional and its content is irrelevant. BEAST only calls
requiresRecalculation() on objects reachable from the posterior, so a clock declared
inside a logger is never told the tree changed and its cache is never invalidated. The
first version of this script did exactly that and produced a stale design. The clock
must hang off a tree likelihood, as it does in production.

Why per-column counts are not enough. If the cached design were never rebuilt after a
topology move, the same entries would stay marked and the counts would not change; they
would simply be attached to the wrong branches. The logger's validate mode therefore also
recomputes the design from scratch at every logged state and reports how many entries
disagree with the cache. That number must be 0 at every state.

Two runs, and they test opposite things:

  --constrained    monophyly enforced on all five clades. The clades are then fixed, so
                   the counts must stay at 3, 25, 51, 3, 3 with 7 background. Catches a
                   design that drifts when it should not.
  --unconstrained  no monophyly. The clades break up, their MRCAs subsume unrelated
                   branches and the counts MUST move. A positive control: if these stay
                   constant, the design is not tracking the tree at all.

    python3 gen_design_check_moving.py --constrained   > moving_constrained.xml
    python3 gen_design_check_moving.py --unconstrained > moving_unconstrained.xml
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

if "--unconstrained" in sys.argv:
    constrained = False
elif "--constrained" in sys.argv:
    constrained = True
else:
    sys.exit("give --constrained or --unconstrained")

CHAIN, LOG_EVERY = 200000, 1000

txt = open(ORIG, encoding="utf-8", errors="replace").read()

def taxset(name):
    m = re.search(r'<taxa id="%s">(.*?)</taxa>' % re.escape(name), txt, re.S)
    return re.findall(r'<taxon idref="([^"]+)"', m.group(1))

CLADES = ["hyper", "pddc", "squirrels", "chile", "Argentina"]
sets = {c: taxset(c) for c in CLADES}
all_taxa = sorted(set(taxset("ingroup")) | {"fake_outgroup"} | {t for c in CLADES for t in sets[c]})
STEM = {c: True for c in CLADES}

NSITE = 200
import random
_rng = random.Random(20260913)
align = {t: "".join(_rng.choice("ACGT") for _ in range(NSITE)) for t in all_taxa}

newick = open(TREE, encoding="utf-8").read().strip()
nbranch = 2 * len(all_taxa) - 2
sys.stderr.write("taxa=%d branches=%d constrained=%s chain=%d\n"
                 % (len(all_taxa), nbranch, constrained, CHAIN))

P = "beast.base"
L = []; w = L.append
w('<beast namespace="%s.evolution.alignment:%s.evolution.tree:%s.inference"' % (P, P, P))
w('       required="BEAST.base v2.7.7:MixedEffectsClock v0.0.1" version="2.7">')
w('  <!-- Dummy alignment: content irrelevant. It exists only so the clock sits inside')
w('       the posterior, which is what makes BEAST invalidate its design cache. -->')
w('  <data id="alignment" spec="%s.evolution.alignment.Alignment" name="alignment" dataType="nucleotide">' % P)
for t in all_taxa:
    w('    <sequence spec="%s.evolution.alignment.Sequence" taxon="%s" totalcount="4" value="%s"/>' % (P, t, align[t]))
w('  </data>')
w('  <run id="mcmc" spec="%s.inference.MCMC" chainLength="%d">' % (P, CHAIN))
w('    <state id="state" spec="%s.inference.State" storeEvery="10000">' % P)
w('      <tree id="Tree" spec="%s.evolution.tree.Tree" name="stateNode">' % P)
w('        <taxonset id="TaxonSet" spec="%s.evolution.alignment.TaxonSet">' % P)
for t in all_taxa:
    w('          <taxon id="%s" spec="%s.evolution.alignment.Taxon"/>' % (t, P))
w('        </taxonset>')
w('      </tree>')
w('      <parameter id="popSize" spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">5000.0</parameter>' % P)
w('      <parameter id="coefficient" spec="%s.inference.parameter.RealParameter" dimension="%d" name="stateNode">0.0</parameter>' % (P, len(CLADES)))
w('      <parameter id="rates" spec="%s.inference.parameter.RealParameter" dimension="%d" lower="0.0" name="stateNode">1.0</parameter>' % (P, nbranch))
w('      <parameter id="ucldStdev" spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">0.3</parameter>' % P)
w('      <parameter id="clockRate" spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">7.0E-9</parameter>' % P)
w('    </state>')
w('    <init id="startingTree" spec="%s.evolution.tree.TreeParser" initial="@Tree"' % P)
w('          taxonset="@TaxonSet" IsLabelledNewick="true" adjustTipHeights="false"')
w('          newick="%s"/>' % newick)
w('    <distribution id="posterior" spec="%s.inference.CompoundDistribution">' % P)
w('      <distribution id="coalescent" spec="%s.evolution.tree.coalescent.Coalescent">' % P)
w('        <treeIntervals id="TreeIntervals" spec="%s.evolution.tree.TreeIntervals" tree="@Tree"/>' % P)
w('        <populationModel id="constant" spec="%s.evolution.tree.coalescent.ConstantPopulation" popSize="@popSize"/>' % P)
w('      </distribution>')
w('      <distribution id="likelihood" spec="%s.evolution.likelihood.TreeLikelihood" data="@alignment" tree="@Tree">' % P)
w('        <siteModel id="siteModel" spec="%s.evolution.sitemodel.SiteModel">' % P)
w('          <substModel id="jc" spec="%s.evolution.substitutionmodel.JukesCantor"/>' % P)
w('        </siteModel>')
w('        <branchRateModel id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"')
w('                         tree="@Tree" rates="@rates" clock.rate="@clockRate" coefficient="@coefficient">')
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
w('        </branchRateModel>')
w('      </distribution>')
if constrained:
    for c in CLADES:
        w('      <distribution id="monophyly.%s" spec="%s.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree">' % (c, P))
        w('        <taxonset id="mono.%s" spec="%s.evolution.alignment.TaxonSet">' % (c, P))
        for t in sets[c]:
            w('          <taxon idref="%s"/>' % t)
        w('        </taxonset>')
        w('      </distribution>')
w('    </distribution>')
# topology and height moves
w('    <operator id="narrowExchange" spec="%s.evolution.operator.Exchange" tree="@Tree" isNarrow="true" weight="15.0"/>' % P)
w('    <operator id="wideExchange" spec="%s.evolution.operator.Exchange" tree="@Tree" isNarrow="false" weight="3.0"/>' % P)
w('    <operator id="subtreeSlide" spec="%s.evolution.operator.SubtreeSlide" tree="@Tree" size="200.0" weight="15.0"/>' % P)
w('    <operator id="wilsonBalding" spec="%s.evolution.operator.WilsonBalding" tree="@Tree" weight="3.0"/>' % P)
w('    <operator id="uniformHeights" spec="%s.evolution.operator.Uniform" tree="@Tree" weight="30.0"/>' % P)
w('    <operator id="popSizeScaler" spec="%s.evolution.operator.kernel.BactrianScaleOperator" parameter="@popSize" scaleFactor="0.5" weight="1.0"/>' % P)
w('    <logger id="tracelog" spec="%s.inference.Logger" fileName="$(filebase).log" logEvery="%d">' % (P, LOG_EVERY))
w('      <log idref="posterior"/>')
w('      <log id="treeHeight" spec="%s.evolution.tree.TreeStatLogger" tree="@Tree"/>' % P)
w('      <log id="designLog" spec="mixedeffectsclock.DesignLogger" validate="true" clock="@clock"/>')
w('    </logger>')
w('    <logger id="screenlog" spec="%s.inference.Logger" logEvery="50000">' % P)
w('      <log idref="posterior"/>')
w('    </logger>')
w('  </run>')
w('</beast>')
print("\n".join(L))
