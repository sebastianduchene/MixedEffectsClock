#!/usr/bin/env python3
"""
A 50-taxon simulation with a genuine random effect AND a genuine slowdown.

Reads everything out of the BEAST X file so nothing is retyped:
  ../../sim_local_MEclock_Bletsa.xml
That file is the log-scale (Bletsa) build, which is the parameterisation this BEAST 2
package implements, so the port is like for like.

Settings carried over unchanged
  data            50 taxa, 5000 sites, tip-dated forwards in years
  substitution    HKY + 4 gamma categories, estimated frequencies
                  kappa ~ LogNormal(1.0, 1.25) on the log scale
                  gamma shape ~ Exponential(mean 0.5)
  clock           beta_0 ~ Normal(-5.3, 2) on the LOG scale. In BEAST 2 the intercept is
                  clock.rate = exp(beta_0), so the same prior is written as a LogNormal
                  with M=-5.3, S=2 in log space. Identical distribution.
                  beta_1 ~ Normal(0, 1); dispersion ~ Exponential(mean 1/3)
                  per-branch multipliers ~ mean-one lognormal
  tree prior      exponential growth; popSize ~ Exponential(mean 1e5),
                  growthRate ~ Laplace(0, 100)
  constraints     monophyly on fastClade and fastSisterStem

One deliberate difference. BEAST X's branchRates.scale is a real-space coefficient of
variation; BEAST 2's ucldStdev is a log-scale standard deviation. The Exponential(mean
1/3) is carried across as written. They agree to about 1% at these values and the truth
is 0, so it does not matter here, but the two columns are not the same quantity.

    python3 gen_sim_re_slowdown.py > sim50_ME_beast2.xml
"""
import re, sys, os

# --targeted swaps BEAUti's tree moves for the targetedbeast set the lepromatosis
# BEAST 2 builds use, and adds the edge weights they need.
TARGETED = "--targeted" in sys.argv
# --mascot swaps the exponential-growth coalescent for a Mascot structured coalescent
# with two demes, "outbreak" holding every sample and "ghost" unsampled. With --targeted
# this is the full stack: mixed-effects clock, targeted tree moves, ORC, and Mascot.
MASCOT   = "--mascot" in sys.argv

# How many distinct rate LEVELS the Mascot skyline should have, i.e. one more than the
# number of change points. --levels 2 means the rate shifts once.
#
# The arithmetic is not what the XML suggests, so it is derived here rather than typed.
# StructuredMigrationSkyline caps its interval index two below the number of shift values,
# so N shift values give N-1 levels; and Skygrowth forces its parameter to dimension N+1,
# with the last entry never read by Mascot. For K levels:
#       shift values  = K+1      (K-1 interior boundaries, then the root, then past it)
#       parameter dim = K+2      (of which K+1 are live)
# ONE shift value does not work. It exits cleanly, logs nothing, and the tip types come
# back NaN, so K must be at least 1 and the grid at least two values.
MASCOT_LEVELS = 2
if "--levels" in sys.argv:
    MASCOT_LEVELS = int(sys.argv[sys.argv.index("--levels") + 1])
assert MASCOT_LEVELS >= 1, "need at least one rate level"
_frac = ["%.10g" % ((i + 1) / MASCOT_LEVELS) for i in range(MASCOT_LEVELS - 1)] + ["1", "1.5"]
MASCOT_GRID = " ".join(_frac)            # K+1 values
MASCOT_DIM  = MASCOT_LEVELS + 2          # K+2 entries, last one inert

HERE = os.path.dirname(os.path.abspath(__file__))
SIM  = os.path.join(HERE, "..", "data")

seqs, name = [], None
for line in open(os.path.join(SIM, "sim_re_slowdown.fasta")):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs.append([name, ""])
    elif line:
        seqs[-1][1] += line.upper()
seqs = [(n, sq) for n, sq in seqs]
assert len(seqs) == 50, len(seqs)
NCHAR = len(seqs[0][1])
assert all(len(sq) == NCHAR for _, sq in seqs)

# the tip label carries its own sampling time: t9_2.75 was sampled at 2.75
dates = {n: n.rsplit("_", 1)[1] for n, _ in seqs}

sets = {}
for line in open(os.path.join(SIM, "sim_re_slowdown_taxa.tsv")):
    k, v = line.rstrip("\n").split("\t"); sets[k] = v.split(",")
fast_clade = sets["fastClade"]
fast_stem  = sets["fastSisterStem"]
slow_clade = sets["slowClade"]
order = [n for n, _ in seqs]

# A Taxon is declared with id= the first time it is used and referenced with idref=
# thereafter. The alignment supplies sequence labels, not Taxon objects, so a bare
# idref into a clade TaxonSet fails with "Could not find object associated with idref".
_seen = set()
def taxon_line(t, indent):
    if t in _seen:
        return '%s<taxon idref="%s"/>' % (indent, t)
    _seen.add(t)
    return '%s<taxon id="%s" spec="beast.base.evolution.alignment.Taxon"/>' % (indent, t)
NBRANCH = 2 * len(order) - 2

sys.stderr.write("taxa=%d sites=%d branches=%d fastClade=%d fastSisterStem=%d\n"
                 % (len(order), NCHAR, NBRANCH, len(fast_clade), len(fast_stem)))
sys.stderr.write("slowClade=%d taxa\n" % len(slow_clade))

P  = "beast.base"
EV = "%s.evolution.operator" % P
IN = "%s.inference.operator" % P
L = []; w = L.append

w('<!--')
w('  Stage 3a: 50-taxon local-clock simulation with KNOWN TRUTH, ported from')
w('  sim_local_MEclock_Bletsa.xml. Generated by gen_sim_re_slowdown.py; do not hand-edit.')
w('')
w('  Truth: 4 of the %d branches evolve at 0.01, the rest at 0.005, with NO' % NBRANCH)
w('  per-branch random variation. In this parameterisation that is')
w('      clock.rate   = 0.005          (exp of the intercept)')
w('      coefficient  = log 2 = 0.6931 (a two-fold increase)')
w('      ucldStdev    = 0')
w('  The 4 fast branches are the clade {%s} with its stem,' % ", ".join(fast_clade))
w('  and the stem alone of its six-tip sister. Both load on the SAME column.')
w('')
w('  The sequences were simulated under Jukes-Cantor, but HKY + gamma is fitted, as in')
w('  the BEAST X run. So kappa should sit near 1 and the gamma shape should be large.')
w('-->')
w('<beast namespace="%s.evolution.alignment:%s.evolution.tree:%s.inference"' % (P, P, P))
w('       required="BEAST.base v2.7.7:MixedEffectsClock v0.0.1:ORC v1.2.1:feast v10.6.1%s%s" version="2.7">' % (":TargetedBeast v2.0.0" if TARGETED else "", ":Mascot v3.0.7" if MASCOT else ""))
w('')
w('  <data id="alignment" spec="%s.evolution.alignment.Alignment" name="alignment" dataType="nucleotide">' % P)
for n, s in seqs:
    w('    <sequence spec="%s.evolution.alignment.Sequence" taxon="%s" totalcount="4" value="%s"/>' % (P, n, s))
w('  </data>')
w('')
w('  <run id="mcmc" spec="%s.inference.MCMC" chainLength="20000000">' % P)
w('    <state id="state" spec="%s.inference.State" storeEvery="5000">' % P)
w('      <tree id="Tree" spec="%s.evolution.tree.Tree" name="stateNode">' % P)
w('        <trait id="dateTrait" spec="%s.evolution.tree.TraitSet" traitname="date-forward"' % P)
w('               value="%s">' % ",".join("%s=%s" % (t, dates[t]) for t in order))
w('          <taxa id="TaxonSet" spec="%s.evolution.alignment.TaxonSet" alignment="@alignment"/>' % P)
w('        </trait>')
w('        <taxonset idref="TaxonSet"/>')
w('      </tree>')
w('      <parameter id="clockRate"   spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">0.005</parameter>' % P)
w('      <parameter id="coefficient" spec="%s.inference.parameter.RealParameter" dimension="2" name="stateNode">0.0</parameter>' % P)
w('      <parameter id="rates"       spec="%s.inference.parameter.RealParameter" dimension="%d" lower="0.0" name="stateNode">1.0</parameter>' % (P, NBRANCH))
w('      <parameter id="ucldStdev"   spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">0.15</parameter>' % P)
w('      <parameter id="kappa"       spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>' % P)
w('      <parameter id="gammaShape"  spec="%s.inference.parameter.RealParameter" lower="0.1" name="stateNode">1.0</parameter>' % P)
w('      <parameter id="freqs"       spec="%s.inference.parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25</parameter>' % P)
if MASCOT:
    # Ne and migration are scaled to THIS tree, not copied from lepromatosis. That tree is
    # thousands of years deep with logNe near 7; this one is ~5.5 years, so logNe near 0
    # (Ne of order 1) and migration near log(0.1) are the sensible starting points.
    w('      <!-- Mascot skyline: %d rate level(s), from %d shift values, parameter' % (MASCOT_LEVELS, MASCOT_LEVELS + 1))
    w('           dimension %d of which the last entry is never read by Mascot. -->' % MASCOT_DIM)
    for d in ['outbreak','ghost']:
        w('      <parameter id="SkylineNe.%s" spec="%s.inference.parameter.RealParameter" dimension="%d" name="stateNode">0.0</parameter>' % (d, P, MASCOT_DIM))
    for m in ['outbreak_to_ghost','ghost_to_outbreak']:
        w('      <parameter id="SkylineMig.%s" spec="%s.inference.parameter.RealParameter" dimension="%d" name="stateNode">-2.302585</parameter>' % (m, P, MASCOT_DIM))
else:
    w('      <parameter id="popSize"     spec="%s.inference.parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>' % P)
    w('      <parameter id="growthRate"  spec="%s.inference.parameter.RealParameter" name="stateNode">0.0</parameter>' % P)
w('    </state>')
w('')
w('    <!-- Random coalescent starting tree honouring both monophyly constraints. -->')
w('    <init id="startingTree" spec="%s.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree" taxa="@alignment">' % P)
w('      <populationModel spec="%s.evolution.tree.coalescent.ConstantPopulation">' % P)
w('        <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="popSize">1.0</parameter>' % P)
w('      </populationModel>')
w('      <constraint idref="monophyly.fastClade"/>')
w('      <constraint idref="monophyly.fastSisterStem"/>')
w('      <constraint idref="monophyly.slowClade"/>')
w('    </init>')
w('')
w('    <distribution id="posterior" spec="%s.inference.CompoundDistribution">' % P)
if MASCOT:
    w('      <!-- Two demes. Every sample is in "outbreak"; "ghost" is unsampled, so it cannot')
    w('           be read off the type trait and must be named in types=. That string is only')
    w('           honoured when fromBeauti is false, which is its default, so it is NOT written.')
    w('           StructuredMigrationSkyline, not StructuredSkyline: migration varies over time.')
    w('           One shared rateShifts object for the dynamics and all four Skygrowths; separate')
    w('           grids silently collapse the skyline to a constant. -->')
    w('      <distribution id="Mascot" spec="mascot.distribution.Mascot" tree="@Tree">')
    w('        <dynamics id="StructuredMigrationSkyline" spec="mascot.dynamics.StructuredMigrationSkyline" types="outbreak ghost">')
    w('          <NeDynamics id="NeDynamicsList" spec="mascot.util.InitializedNeDynamicsList">')
    for d in ['outbreak','ghost']:
        w('            <neDynamics id="NeDynamics.%s" spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineNe.%s" rateShifts="@rateShifts"/>' % (d,d))
    w('          </NeDynamics>')
    w('          <!-- Index order is fixed by the class: 0 is outbreak to ghost, 1 is ghost to')
    w('               outbreak, both FORWARDS in time. -->')
    w('          <migrationDynamics id="MigDynamicsList" spec="mascot.util.InitializedNeDynamicsList">')
    for m in ['outbreak_to_ghost','ghost_to_outbreak']:
        w('            <neDynamics id="MigDynamics.%s" spec="mascot.parameterdynamics.Skygrowth" logNe="@SkylineMig.%s" rateShifts="@rateShifts"/>' % (m,m))
    w('          </migrationDynamics>')
    w('          <rateShifts id="rateShifts" spec="mascot.dynamics.RateShifts" tree="@Tree">%s</rateShifts>' % MASCOT_GRID)
    w('          <indicators id="indicators" spec="%s.inference.parameter.BooleanParameter" dimension="2" estimate="false">true</indicators>' % P)
    w('          <typeTrait id="typeTraitSet" spec="mascot.util.InitializedTraitSet" traitname="type"')
    w('                     value="%s">' % ",".join("%s=outbreak" % t for t in order))
    w('            <taxa idref="TaxonSet"/>')
    w('          </typeTrait>')
    w('        </dynamics>')
    w('        <structuredTreeIntervals id="StructuredTreeIntervals" spec="mascot.distribution.StructuredTreeIntervals" tree="@Tree"/>')
    w('      </distribution>')
    for nm, ctr in [('SkylineNe.outbreak','0.0'),('SkylineNe.ghost','0.0'),
                    ('SkylineMig.outbreak_to_ghost','-2.302585'),('SkylineMig.ghost_to_outbreak','-2.302585')]:
        w('      <distribution id="%s.FirstPrior" spec="%s.inference.distribution.Prior">' % (nm, P))
        w('        <x id="first.%s" spec="mascot.util.First" arg="@%s"/>' % (nm, nm))
        w('        <distr spec="%s.inference.distribution.Normal">' % P)
        w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">%s</parameter>' % (P, ctr))
        w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="sigma">2.0</parameter>' % P)
        w('        </distr>')
        w('      </distribution>')
        w('      <distribution id="%s.SmoothPrior" spec="%s.inference.distribution.Prior">' % (nm, P))
        w('        <x id="diff.%s" spec="mascot.util.Difference" arg="@%s"/>' % (nm, nm))
        w('        <distr spec="%s.inference.distribution.Normal">' % P)
        w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">0.0</parameter>' % P)
        w('          <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="sigma">1.0</parameter>' % P)
        w('        </distr>')
        w('      </distribution>')
else:
    w('      <distribution id="coalescent" spec="%s.evolution.tree.coalescent.Coalescent">' % P)
    w('        <treeIntervals id="TreeIntervals" spec="%s.evolution.tree.TreeIntervals" tree="@Tree"/>' % P)
    w('        <populationModel id="expGrowth" spec="%s.evolution.tree.coalescent.ExponentialGrowth" popSize="@popSize" growthRate="@growthRate"/>' % P)
    w('      </distribution>')
PRIORS = [
    ("clockRatePrior", "@clockRate",
     ['<distr spec="%s.inference.distribution.LogNormalDistributionModel" meanInRealSpace="false">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="M">-5.3</parameter>' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="S">2.0</parameter>' % P,
      '</distr>']),
    ("coefficientPrior", "@coefficient",
     ['<distr spec="%s.inference.distribution.Normal">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">0.0</parameter>' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="sigma">1.0</parameter>' % P,
      '</distr>']),
    ("ucldStdevPrior", "@ucldStdev",
     ['<distr spec="%s.inference.distribution.Exponential">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">0.3333333333333333</parameter>' % P,
      '</distr>']),
    ("ratesPrior", "@rates",
     ['<distr spec="%s.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>' % P,
      '</distr>']),
    ("kappaPrior", "@kappa",
     ['<distr spec="%s.inference.distribution.LogNormalDistributionModel" meanInRealSpace="false">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="S">1.25</parameter>' % P,
      '</distr>']),
    ("gammaShapePrior", "@gammaShape",
     ['<distr spec="%s.inference.distribution.Exponential">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">0.5</parameter>' % P,
      '</distr>']),
    ("freqsPrior", "@freqs",
     ['<distr spec="%s.inference.distribution.Dirichlet">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" dimension="4" estimate="false" name="alpha">1.0 1.0 1.0 1.0</parameter>' % P,
      '</distr>']),
]
if not MASCOT:
    PRIORS += [
    ("popSizePrior", "@popSize",
     ['<distr spec="%s.inference.distribution.Exponential">' % P,
      '  <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="mean">1.0E5</parameter>' % P,
      '</distr>']),
    ("growthRatePrior", "@growthRate",
     ['<distr spec="%s.inference.distribution.LaplaceDistribution" mu="0.0" scale="100.0"/>' % P]),
]
for pid, spec, body in PRIORS:
    w('      <distribution id="%s" spec="%s.inference.distribution.Prior" x="%s">' % (pid, P, spec))
    for b in body:
        w('        ' + b)
    w('      </distribution>')
w('')
if TARGETED:
    w('      <!-- Edge weights for the targeted moves. Inside the posterior on purpose: they')
    w('           are a function of the tree, and BEAST invalidates a cache only for objects')
    w('           reachable from the posterior. ParsimonyWeights, not ParsimonyWeights2: the')
    w('           latter throws ArrayIndexOutOfBoundsException on the lepromatosis data. -->')
    w('      <distribution id="consensusWeights" spec="targetedbeast.edgeweights.ParsimonyWeights" tree="@Tree" maxWeight="20">')
    w('        <data spec="targetedbeast.alignment.ConsensusAlignment" data="@alignment"/>')
    w('      </distribution>')
w('      <distribution id="likelihood" spec="%s.evolution.likelihood.TreeLikelihood" data="@alignment" tree="@Tree">' % P)
w('        <siteModel id="siteModel" spec="%s.evolution.sitemodel.SiteModel" gammaCategoryCount="4" shape="@gammaShape">' % P)
w('          <substModel id="hky" spec="%s.evolution.substitutionmodel.HKY" kappa="@kappa">' % P)
w('            <frequencies id="freqModel" spec="%s.evolution.substitutionmodel.Frequencies" frequencies="@freqs"/>' % P)
w('          </substModel>')
w('        </siteModel>')
w('        <branchRateModel id="clock" spec="mixedeffectsclock.MixedEffectsClockModel"')
w('                         tree="@Tree" rates="@rates" clock.rate="@clockRate" coefficient="@coefficient">')
w('          <distr id="rateDistr" spec="%s.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">' % P)
w('            <parameter spec="%s.inference.parameter.RealParameter" estimate="false" name="M">1.0</parameter>' % P)
w('          </distr>')
w('          <!-- Both sets load on column 0: a single shared coefficient, as in BEAST X. -->')
w('          <clade id="cladeFast" spec="mixedeffectsclock.CladeDesign" category="0" includeStem="true">')
w('            <taxonset id="taxa.fastClade" spec="%s.evolution.alignment.TaxonSet">' % P)
for t in fast_clade:
    w(taxon_line(t, '              '))
w('            </taxonset>')
w('          </clade>')
w('          <clade id="cladeFastStem" spec="mixedeffectsclock.CladeDesign" category="0" includeStem="true" excludeClade="true">')
w('            <taxonset id="taxa.fastSisterStem" spec="%s.evolution.alignment.TaxonSet">' % P)
for t in fast_stem:
    w(taxon_line(t, '              '))
w('            </taxonset>')
w('          </clade>')
w('          <!-- Column 1: the SLOW clade, with its stem. Truth is a fold change of 0.3. -->')
w('          <clade id="cladeSlow" spec="mixedeffectsclock.CladeDesign" category="1" includeStem="true">')
w('            <taxonset id="taxa.slowClade" spec="%s.evolution.alignment.TaxonSet">' % P)
for t in slow_clade:
    w(taxon_line(t, '              '))
w('            </taxonset>')
w('          </clade>')
w('        </branchRateModel>')
w('      </distribution>')
w('')
w('      <!-- Constraints reference the clock\'s own TaxonSets, so the two cannot drift apart. -->')
for nm, ts in [("fastClade", "taxa.fastClade"), ("fastSisterStem", "taxa.fastSisterStem"), ("slowClade", "taxa.slowClade")]:
    w('      <distribution id="monophyly.%s" spec="%s.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree">' % (nm, P))
    w('        <taxonset idref="%s"/>' % ts)
    w('      </distribution>')
w('    </distribution>')
w('')
w('    <!-- clock and substitution parameters -->')
w('    <operator id="clockRateScaler" spec="%s.kernel.BactrianScaleOperator" parameter="@clockRate" scaleFactor="0.5" weight="3.0"/>' % EV)
w('    <operator id="coefficientWalk" spec="%s.kernel.BactrianRandomWalkOperator" parameter="@coefficient" scaleFactor="0.5" weight="3.0"/>' % IN)
w('    <operator id="ratesScaler"     spec="%s.kernel.BactrianScaleOperator" parameter="@rates" scaleFactor="0.5" weight="15.0"/>' % EV)
w('    <!-- ORC joint move: without it the dispersion is about 40x less efficient. -->')
w('    <operator id="ucldScalerJoint" spec="orc.consoperators.UcldScalerOperator" distr="@rateDistr" rates="@rates" stdev="@ucldStdev" scaleFactor="0.5" weight="3.0">')
w('      <kernel spec="%s.kernel.KernelDistribution$Bactrian"/>' % IN)
w('    </operator>')
w('    <operator id="kappaScaler"      spec="%s.kernel.BactrianScaleOperator" parameter="@kappa" scaleFactor="0.5" weight="1.0"/>' % EV)
w('    <operator id="gammaShapeScaler" spec="%s.kernel.BactrianScaleOperator" parameter="@gammaShape" scaleFactor="0.5" weight="1.0"/>' % EV)
w('    <operator id="freqsExchange"    spec="%s.kernel.BactrianDeltaExchangeOperator" parameter="@freqs" delta="0.01" weight="1.0"/>' % IN)
if MASCOT:
    for nm in ['SkylineNe.outbreak','SkylineNe.ghost','SkylineMig.outbreak_to_ghost','SkylineMig.ghost_to_outbreak']:
        w('    <operator id="%s.Walker" spec="%s.kernel.BactrianRandomWalkOperator" parameter="@%s" scaleFactor="0.5" weight="1.0"/>' % (nm, IN, nm))
else:
    w('    <operator id="popSizeScaler"    spec="%s.kernel.BactrianScaleOperator" parameter="@popSize" scaleFactor="0.5" weight="3.0"/>' % EV)
    w('    <!-- Small window: the tree is only a few years deep, so a wide random walk on the -->')
    w('    <!-- growth rate is rejected almost always and floods the log with -Infinity.      -->')
    w('    <operator id="growthRateWalk"   spec="%s.kernel.BactrianRandomWalkOperator" parameter="@growthRate" scaleFactor="0.1" weight="3.0"/>' % IN)
w('')
TB = "targetedbeast.operators"
if TARGETED:
    w('    <!-- targetedbeast tree moves, the set the lepromatosis BEAST 2 builds use. -->')
    w('    <!-- They act on the RAW rates and node heights, exactly as ORC does, so the')
    w('         clade factor is a per-branch constant the moves never touch. -->')
    w('    <operator id="WeightBasedNodeRandomizer" spec="%s.WeightBasedNodeRandomizer" percentage="0.01" optimise="true" tree="@Tree" weight="0.5" edgeWeights="@consensusWeights"/>' % TB)
    w('    <operator id="HeightBasedNodeRandomizer" spec="%s.HeightBasedNodeRandomizer" percentage="0.01" optimise="true" tree="@Tree" weight="0.5"/>' % TB)
    w('    <operator id="RangeByLength" spec="%s.RangeSlide" tree="@Tree" weight="10.0" size="0.1" edgeWeights="@consensusWeights" weightByBranchLength="true"/>' % TB)
    w('    <operator id="Range" spec="%s.RangeSlide" tree="@Tree" weight="30.0" size="0.1" edgeWeights="@consensusWeights" sqrtWeights="true"/>' % TB)
    w('    <operator id="RangeUniform" spec="%s.RangeSlide" tree="@Tree" weight="10.0" size="0.1" edgeWeights="@consensusWeights" uniform="true"/>' % TB)
    w('    <operator id="UnTargetedWide" spec="%s.WeightedWideOperator" tree="@Tree" weight="5.0"/>' % TB)
    w('    <operator id="TargetedWide" spec="%s.WeightedWideOperator" tree="@Tree" weight="15.0" edgeWeights="@consensusWeights"/>' % TB)
    w('    <operator id="TargetedWilsonBalding" spec="%s.TargetedWilsonBaldingRates" rates="@rates" tree="@Tree" weight="5" edgeWeights="@consensusWeights"/>' % TB)
    w('    <operator id="TargetedWilsonBalding2" spec="%s.TargetedWilsonBaldingRates" rates="@rates" tree="@Tree" weight="1" edgeWeights="@consensusWeights" useEdgeLength="true"/>' % TB)
    w('    <operator id="ScaleAll" spec="%s.IntervalScaleOperator" tree="@Tree" weight="0.5" scaleFactor="0.9">' % TB)
    w('      <down idref="clockRate"/>')
    w('    </operator>')
    w('    <operator id="ScaleIntervalsRandomly" spec="%s.IntervalScaleOperator" tree="@Tree" weight="0.5" scaleFactor="0.1" scaleAllNodesIndependently="true"/>' % TB)
    w('    <operator id="ScaleIntervalsWithRates" spec="%s.IntervalRateCoScaler" scaleFactor="0.25" weight="3.0">' % TB)
    w('      <branchRates idref="rates"/>')
    w('      <tree idref="Tree"/>')
    w('    </operator>')
    w('    <operator id="ScaleIntervalsWithRatesAndStdev" spec="%s.IntervalRateCoScaler" scaleFactor="0.25" weight="3.0">' % TB)
    w('      <branchRates idref="rates"/>')
    w('      <stdev idref="ucldStdev"/>')
    w('      <tree idref="Tree"/>')
    w('    </operator>')
else:
    w('    <!-- BEAUti\'s standard tree move set -->')
    w('    <operator id="treeScaler"     spec="%s.kernel.BactrianScaleOperator" tree="@Tree" scaleFactor="0.5" upper="10.0" weight="3.0"/>' % EV)
    w('    <operator id="treeRootScaler" spec="%s.kernel.BactrianScaleOperator" tree="@Tree" rootOnly="true" scaleFactor="0.5" upper="10.0" weight="3.0"/>' % EV)
    w('    <operator id="uniformHeights" spec="%s.kernel.BactrianNodeOperator" tree="@Tree" weight="30.0"/>' % EV)
    w('    <operator id="subtreeSlide"   spec="%s.kernel.BactrianSubtreeSlide" tree="@Tree" weight="15.0"/>' % EV)
    w('    <operator id="narrowExchange" spec="%s.Exchange" tree="@Tree" weight="15.0"/>' % EV)
    w('    <operator id="wideExchange"   spec="%s.Exchange" tree="@Tree" isNarrow="false" weight="3.0"/>' % EV)
    w('    <operator id="wilsonBalding"  spec="%s.WilsonBalding" tree="@Tree" weight="3.0"/>' % EV)
    w('    <operator id="upDownRates"    spec="%s.kernel.BactrianUpDownOperator" scaleFactor="0.75" weight="3.0">' % IN)
    w('      <up idref="clockRate"/>')
    w('      <down idref="Tree"/>')
    w('    </operator>')
    w('')

w('    <logger id="tracelog" spec="%s.inference.Logger" fileName="$(filebase).log" logEvery="2000">' % P)
_tree_prior_logs = ["Mascot"] if MASCOT else ["coalescent", "popSize", "growthRate"]
for r in ["posterior", "likelihood", "clockRate", "coefficient",
          "ucldStdev", "kappa", "gammaShape", "freqs"] + _tree_prior_logs:
    w('      <log idref="%s"/>' % r)
if MASCOT:
    w('      <!-- Migration counts between the two demes, from a stochastic map declared')
    w('           once here and reused by the events tree logger. -->')
    w('      <log id="nrEventsLogger" spec="mascot.logger.MigrationCountLogger">')
    w('        <mappedMascot id="mascotEventsTreelog" spec="mascot.distribution.MappedMascot"')
    w('                      dynamics="@StructuredMigrationSkyline"')
    w('                      structuredTreeIntervals="@StructuredTreeIntervals" tree="@Tree"/>')
    w('      </log>')
    w('      <!-- the realised interval boundaries in years: the grid is a fraction of the root -->')
    w('      <log idref="rateShifts"/>')
    for nm in ['SkylineNe.outbreak','SkylineNe.ghost','SkylineMig.outbreak_to_ghost','SkylineMig.ghost_to_outbreak']:
        w('      <log idref="%s"/>' % nm)
w('      <!-- the BEAST X convention for the same quantity: scale = sqrt(exp(sigma^2)-1).')
w('           useCaching="false" is required: a calculator in a logger is never told its')
w('           argument moved, and with caching on it logs its initial value forever. -->')
w('      <log id="branchRates.scale" spec="feast.expressions.ExpCalculator" useCaching="false" value="sqrt(exp(ucldStdev^2) - 1)">')
w('        <arg idref="ucldStdev"/>')
w('      </log>')
w('      <log id="treeHeight" spec="%s.evolution.tree.TreeStatLogger" tree="@Tree"/>' % P)
w('      <log id="rateStat" spec="%s.evolution.RateStatistic" branchratemodel="@clock" tree="@Tree"/>' % P)
w('      <log id="designLog" spec="mixedeffectsclock.DesignLogger" validate="false" clock="@clock"/>')
w('    </logger>')
w('    <logger id="screenlog" spec="%s.inference.Logger" logEvery="100000">' % P)
w('      <log idref="posterior"/>')
w('      <log idref="clockRate"/>')
w('      <log idref="coefficient"/>')
w('    </logger>')
w('    <logger id="treelog" spec="%s.inference.Logger" fileName="$(filebase).trees" logEvery="20000" mode="tree">' % P)
w('      <log id="treeLogger" spec="%s.evolution.TreeWithMetaDataLogger" tree="@Tree" branchratemodel="@clock"/>' % P)
w('    </logger>')
w('  </run>')
w('</beast>')
print("\n".join(L))
