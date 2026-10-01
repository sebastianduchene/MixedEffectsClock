# Claude Code skills

Two [Claude Code](https://claude.com/claude-code) skills ship with this package:

`mixed-effects-clock-beast2/` is a setup reference: the model and how each term maps to XML,
the wiring that cannot be skipped, the dispersion conversion to the BEAST X convention, the
Mascot skyline grid arithmetic, and a triage table for the silent failures.
`templates/me_clock_minimal.xml` is a complete six-taxon skeleton that runs on its own.

`me-design-prior-check/` is an end-to-end verifier of the design matrix: it redacts the
sequences in a MixedEffectsClock XML to Ns, runs a short MCMC that samples from the prior,
then draws ten sampled trees with every branch coloured by the clade column it carries.
Use it on any analysis XML before launching a production run.

To install either, copy or link the directory into `~/.claude/skills/`:

    ln -s "$PWD/mixed-effects-clock-beast2" ~/.claude/skills/
    ln -s "$PWD/me-design-prior-check"     ~/.claude/skills/

Then ask Claude for a mixed-effects clock or a design check, or invoke a skill by name.
