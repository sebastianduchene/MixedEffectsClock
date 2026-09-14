# Claude Code skill

`mixed-effects-clock-beast2/` is a [Claude Code](https://claude.com/claude-code) skill for
setting up this package: the model and how each term maps to XML, the wiring that cannot be
skipped, the dispersion conversion to the BEAST X convention, the Mascot skyline grid
arithmetic, and a triage table for the silent failures.

To install it, copy or link the directory into `~/.claude/skills/`:

    ln -s "$PWD/mixed-effects-clock-beast2" ~/.claude/skills/

Then ask Claude for a mixed-effects clock, or invoke it by name.

`templates/me_clock_minimal.xml` is a complete six-taxon skeleton that runs on its own. It
samples the prior, so running it checks the installation end to end and should recover every
prior it declares. It is also useful without Claude, as the shortest correct example in the
repository.
