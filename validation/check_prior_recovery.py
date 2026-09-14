#!/usr/bin/env python3
"""
Stage 2: did every sampled parameter recover its prior?

With an all-ambiguous alignment the likelihood is flat, so the posterior IS the prior and
every prior here has a closed form. The script compares the posterior mean and standard
deviation with the analytic ones, and judges the mean against its Monte Carlo standard
error rather than by eye: |posterior mean - prior mean| must be within 3 MCSE, where
MCSE = sd / sqrt(ESS).

    python3 check_prior_recovery.py <stem> [burninFraction]
"""
import math, sys

stem = sys.argv[1] if len(sys.argv) > 1 else "prior_recovery_6taxon"
burnin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10

rows, hdr = [], None
for line in open(stem + ".log"):
    if line.startswith("#") or not line.strip():
        continue
    p = line.rstrip("\n").split("\t")
    if hdr is None:
        hdr = p
        continue
    if len(p) == len(hdr):
        rows.append(p)
rows = rows[int(len(rows) * burnin):]
col = {n: i for i, n in enumerate(hdr)}
def v(name):
    return [float(r[col[name]]) for r in rows]

def ess(x):
    """Standard autocorrelation-time ESS, truncated on the first negative pair sum."""
    n = len(x)
    m = sum(x) / n
    d = [xi - m for xi in x]
    var = sum(t * t for t in d) / n
    if var <= 0:
        return float(n)
    rho, s, k = [], 0.0, 1
    while k < n // 2:
        r = sum(d[i] * d[i + k] for i in range(n - k)) / (n * var)
        rho.append(r)
        if k % 2 == 0 and len(rho) >= 2 and rho[-1] + rho[-2] < 0:
            break
        k += 1
    act = 1 + 2 * sum(rho)
    return n / act if act > 0 else float(n)

# name -> (analytic mean, analytic sd, label)
TARGETS = {
    "coefficient.1": (0.0, 2.0, "Normal(0, 2)"),
    "coefficient.2": (0.0, 2.0, "Normal(0, 2)"),
    "coefficient.3": (0.0, 2.0, "Normal(0, 2)"),
    "clockRate":     (70 * 1e-10, math.sqrt(70) * 1e-10, "Gamma(70, 1e-10)"),
    "ucldStdev":     (1 / 3.0, 1 / 3.0, "Exponential(mean 1/3)"),
}

print("%s: %d samples after %.0f%% burn-in\n" % (stem, len(rows), burnin * 100))
print("%-14s %-22s %11s %11s %9s %9s %7s  %s"
      % ("parameter", "prior", "post.mean", "prior.mean", "post.sd", "prior.sd", "ESS", "mean within 3 MCSE"))
bad = 0
for name, (pm, psd, label) in TARGETS.items():
    if name not in col:
        print("  %-12s MISSING from the log" % name); bad += 1; continue
    x = v(name)
    n = len(x)
    mean = sum(x) / n
    sd = math.sqrt(sum((xi - mean) ** 2 for xi in x) / (n - 1))
    e = ess(x)
    mcse = sd / math.sqrt(max(e, 1.0))
    ok = abs(mean - pm) <= 3 * mcse
    if not ok:
        bad += 1
    print("%-14s %-22s %11.4g %11.4g %9.4g %9.4g %9.0f  %s"
          % (name, label, mean, pm, sd, psd, e, "yes" if ok else "NO"))

print()
# the design must have stayed correct throughout
for name, want in [("nPainted.X1", 3), ("nPainted.X2", 2), ("nPainted.X3", 1),
                   ("nPainted.background", 4), ("nPainted.multiple", 0),
                   ("design.staleEntries", 0)]:
    if name in col:
        u = sorted(set(int(x) for x in v(name)))
        ok = u == [want]
        if not ok:
            bad += 1
        print("  %-22s %-14s %s" % (name, u, "ok" if ok else "EXPECTED %d" % want))
for name in [c for c in hdr if c.startswith("monophyletic.")]:
    u = sorted(set(int(x) for x in v(name)))
    ok = u == [1]
    if not ok:
        bad += 1
    print("  %-22s %-14s %s" % (name, u, "ok" if ok else "CONSTRAINT VIOLATED"))

print("\n%s" % ("PASS" if bad == 0 else "%d PROBLEM(S)" % bad))
sys.exit(1 if bad else 0)
