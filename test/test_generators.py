#!/usr/bin/env python3
"""
Tests for the XML generators.

Everything in test/mixedeffectsclock is Java. No Java bug has been found in this package;
the bugs that have actually happened are in the generators, and they are arithmetic and
wiring. These are cheap and catch that class.

    python3 test/test_generators.py
"""
import os, re, subprocess, sys, unittest
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
GEN  = os.path.join(HERE, "..", "examples", "gen_sim_re_slowdown.py")


def generate(*flags):
    r = subprocess.run([sys.executable, GEN] + list(flags),
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("generator failed: " + r.stderr.strip()[-400:])
    return r.stdout


class GridArithmetic(unittest.TestCase):
    """
    K rate levels need K+1 shift values and parameter dimension K+2.

    The class caps its interval index two below the number of shift values, so N values give
    N-1 levels, and Skygrowth separately forces dimension N+1 with the last entry never read.
    Letting the grid and the dimensions drift apart is the single easiest way to produce a
    silently wrong Mascot model, which is why the grid is derived rather than typed.
    """

    def _grid_and_dims(self, xml):
        grid = re.search(r'<rateShifts[^>]*>([^<]*)<', xml).group(1).split()
        dims = set(re.findall(r'id="Skyline[^"]*"[^>]*dimension="(\d+)"', xml))
        self.assertEqual(1, len(dims), "all four skyline vectors must share one dimension")
        return grid, int(dims.pop())

    def test_levels_set_grid_and_dimension_together(self):
        for k in (1, 2, 3, 4):
            xml = generate("--targeted", "--mascot", "--levels", str(k))
            grid, dim = self._grid_and_dims(xml)
            self.assertEqual(k + 1, len(grid), "%d levels need %d shift values" % (k, k + 1))
            self.assertEqual(k + 2, dim, "%d levels need dimension %d" % (k, k + 2))

    def test_boundaries_set_grid_and_dimension_together(self):
        xml = generate("--targeted", "--mascot", "--boundaries", "2.5")
        grid, dim = self._grid_and_dims(xml)
        self.assertEqual(["2.5", "5", "7.5"], grid)
        self.assertEqual(4, dim, "one boundary is two levels, so dimension 4")

        xml = generate("--targeted", "--mascot", "--boundaries", "2 4")
        grid, dim = self._grid_and_dims(xml)
        self.assertEqual(3, len(grid) - 1, "two boundaries are three levels")
        self.assertEqual(5, dim)

    def test_a_single_shift_value_is_never_emitted(self):
        """One shift value makes BEAST exit 0, log nothing and report NaN tip types."""
        for flags in (("--levels", "1"), ("--boundaries", "2.5")):
            xml = generate("--targeted", "--mascot", *flags)
            grid, _ = self._grid_and_dims(xml)
            self.assertGreaterEqual(len(grid), 2,
                                    "a one-value grid fails silently and must never be written")


class GridKind(unittest.TestCase):
    """
    Relative and calendar grids differ by ONE attribute and are otherwise identical on
    inspection. An absolute grid that does not reach the root is how the lepromatosis
    skyline was silently collapsed to a constant.
    """

    def _rateshifts_element(self, xml):
        return re.search(r'<rateShifts[^>]*>', xml).group(0)

    def test_levels_gives_a_tree_relative_grid(self):
        el = self._rateshifts_element(generate("--targeted", "--mascot", "--levels", "2"))
        self.assertIn('tree=', el, "--levels must pass a tree, making the values fractions")

    def test_boundaries_gives_an_absolute_grid(self):
        el = self._rateshifts_element(generate("--targeted", "--mascot", "--boundaries", "2.5"))
        self.assertNotIn('tree=', el, "--boundaries must NOT pass a tree, or the times are"
                                      " read as fractions of the root height")


class Layers(unittest.TestCase):
    """Each flag must add exactly its own layer, and no flag may drop another."""

    CASES = {
        (): {"targetedbeast": False, "mascot.dynamics": False},
        ("--targeted",): {"targetedbeast": True, "mascot.dynamics": False},
        ("--targeted", "--mascot"): {"targetedbeast": True, "mascot.dynamics": True},
        ("--mascot",): {"targetedbeast": False, "mascot.dynamics": True},
    }

    def test_every_combination_parses_and_has_the_right_layers(self):
        for flags, want in self.CASES.items():
            xml = generate(*flags)
            ET.fromstring(xml)                       # must be well-formed
            for token, expected in want.items():
                self.assertEqual(expected, token in xml,
                                 "%s with flags %s" % (token, flags or "(none)"))
            # the clock and ORC are in every configuration
            self.assertIn("mixedeffectsclock.MixedEffectsClockModel", xml)
            self.assertIn("orc.consoperators.UcldScalerOperator", xml)

    def test_the_design_has_every_clade_the_simulation_defines(self):
        xml = generate()
        clades = re.findall(r'<clade id="(\w+)"[^>]*category="(\d+)"', xml)
        self.assertEqual(3, len(clades), "two fast elements and one slow")
        self.assertEqual({"0", "1"}, {c for _, c in clades})
        dim = re.search(r'id="coefficient"[^>]*dimension="(\d+)"', xml).group(1)
        self.assertEqual("2", dim, "two columns, so two coefficients")

    def test_every_clade_is_constrained_monophyletic(self):
        xml = generate()
        clade_sets = set(re.findall(r'<clade id="\w+"[^>]*>\s*<taxonset id="([\w.]+)"', xml))
        constrained = set(re.findall(r'monophyletic="true"[^>]*>\s*<taxonset idref="([\w.]+)"', xml))
        self.assertEqual(clade_sets, constrained,
                         "every clade carrying a column must be constrained, by idref to the"
                         " clock's OWN taxon set so the two cannot drift apart")


class XmlHygiene(unittest.TestCase):
    def test_no_double_hyphen_inside_a_comment(self):
        """Illegal in XML, and it has broken generated files three times."""
        for flags in ((), ("--targeted",), ("--targeted", "--mascot")):
            xml = generate(*flags)
            bad = [m.group(0)[:70] for m in re.finditer(r'<!--.*?-->', xml, re.S)
                   if "--" in m.group(0)[4:-3]]
            self.assertEqual([], bad, "illegal comment with flags %s" % (flags or "(none)",))

    def test_the_feast_conversion_never_caches(self):
        """A calculator in a logger is never told its argument moved; with caching on it
        logs its initial value for the whole run."""
        xml = generate()
        el = re.search(r'<log[^>]*ExpCalculator[^>]*>', xml).group(0)
        self.assertIn('useCaching="false"', el)


class ShippedTemplate(unittest.TestCase):
    """
    The skill's template is the first XML anyone will copy, so the non-negotiables have to
    be in it. Each assertion here is a failure that has actually happened, and every one of
    them is silent: the run finishes and the numbers look plausible.
    """

    def setUp(self):
        path = os.path.join(HERE, "..", "skills", "mixed-effects-clock-beast2",
                            "templates", "me_clock_minimal.xml")
        with open(path) as fh:
            self.xml = fh.read()
        self.root = ET.fromstring(self.xml)

    def test_it_parses(self):
        self.assertEqual("beast", self.root.tag)

    def test_the_clock_is_inside_the_tree_likelihood(self):
        """In a logger instead, it is never told the tree moved and its design goes stale."""
        lik = self.root.find(".//distribution[@id='likelihood']")
        self.assertIsNotNone(lik)
        self.assertIsNotNone(lik.find("branchRateModel"))

    def test_the_branch_rates_carry_a_prior(self):
        """Without it ucldStdev answers only to its own prior, the joint move is unbalanced,
        and the dispersion climbs away."""
        pri = self.root.find(".//distribution[@id='ratesPrior']")
        self.assertIsNotNone(pri, "no Prior on @rates")
        self.assertEqual("@rates", pri.get("x"))

    def test_every_clade_is_constrained_by_idref_to_the_clock_taxonset(self):
        clades = self.root.findall(".//clade")
        self.assertTrue(clades)
        declared = {c.find("taxonset").get("id") for c in clades}
        constrained = {t.get("idref")
                       for m in self.root.findall(".//distribution")
                       if "MRCAPrior" in (m.get("spec") or "")
                       and m.get("monophyletic") == "true"
                       for t in m.findall("taxonset")}
        self.assertEqual(declared, declared & constrained,
                         "unconstrained design columns: %s" % (declared - constrained))

    def test_the_dispersion_has_the_orc_joint_operator(self):
        specs = [o.get("spec") for o in self.root.findall(".//operator")]
        self.assertIn("orc.consoperators.UcldScalerOperator", specs)

    def test_the_coefficients_get_a_random_walk_not_a_scaler(self):
        """They are log scale and cross zero, so a scale operator cannot move them past it."""
        ops = [o.get("spec") for o in self.root.findall(".//operator")
               if o.get("parameter") == "@coefficient"]
        self.assertTrue(ops, "nothing operates on @coefficient")
        for spec in ops:
            self.assertIn("RandomWalk", spec)

    def test_the_feast_conversion_never_caches(self):
        el = re.search(r'<log[^>]*ExpCalculator[^>]*>', self.xml).group(0)
        self.assertIn('useCaching="false"', el)

    def test_no_double_hyphen_inside_a_comment(self):
        bad = [m.group(0)[:70] for m in re.finditer(r'<!--.*?-->', self.xml, re.S)
               if "--" in m.group(0)[4:-3]]
        self.assertEqual([], bad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
