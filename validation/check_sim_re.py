#!/usr/bin/env python3
"""Score the random-effect + slowdown simulation. Usage: check_sim_re.py <stem> [burnin]"""
import math, sys

stem   = sys.argv[1] if len(sys.argv) > 1 else "sim_re_slowdown"
burnin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.20

hdr, rows = None, []
for line in open(stem + ".log"):
    if line.startswith("#") or not line.strip(): continue
    p = line.rstrip("\n").split("\t")
    if hdr is None: hdr = p; continue
    if len(p) == len(hdr): rows.append(p)
n_all = len(rows)
rows = rows[int(n_all * burnin):]
col = {c: i for i, c in enumerate(hdr)}
def v(c): return [float(r[col[c]]) for r in rows]

def ess(x):
    n = len(x); m = sum(x)/n; d=[t-m for t in x]
    var = sum(t*t for t in d)/n
    if var <= 0: return float(n)
    rho, k = [], 1
    while k < n//2:
        rho.append(sum(d[i]*d[i+k] for i in range(n-k))/(n*var))
        if k % 2 == 0 and len(rho) >= 2 and rho[-1]+rho[-2] < 0: break
        k += 1
    act = 1 + 2*sum(rho)
    return n/act if act > 0 else float(n)

def hpd(x, p=0.95):
    s = sorted(x); n = len(s); k = max(1, int(math.floor(p*n)))
    best = min(range(n-k), key=lambda i: s[i+k]-s[i])
    return s[best], s[best+k]

TRUTH = {"clockRate": 0.005,
         "coefficient.1": math.log(2.0),     # fast clade, 4 branches
         "coefficient.2": math.log(0.3)}     # SLOW clade, 19 branches
# ucldStdev is handled separately: its truth is 0, which sits ON THE BOUNDARY of a
# strictly positive parameter, so no interval can ever contain it and an HPD coverage
# test is meaningless. The right question is whether the posterior has been pulled well
# below its prior (Exponential with mean 1/3) and piles up against zero.
TRUE_DISPERSION = 0.2816   # realised sd of log eps in this dataset (nominal 0.300)
EXTRA = ["kappa", "gammaShape", "Tree.height", "rateStat.mean", "rateStat.coefficientOfVariation"]

print("%s: %d samples logged, %d after %.0f%% burn-in (state %s)\n"
      % (stem, n_all, len(rows), burnin*100, rows[-1][0]))
print("%-14s %>10s %11s %24s %8s  %s".replace(">","") % ("parameter","truth","median","95% HPD","ESS","covers truth"))
bad = 0
for name, t in TRUTH.items():
    if name not in col: print("  %s MISSING" % name); bad += 1; continue
    x = v(name); x.sort()
    med = x[len(x)//2]; lo, hi = hpd(x); e = ess(v(name))
    ok = lo <= t <= hi
    if not ok: bad += 1
    print("%-14s %10.5g %11.5g  [%9.5g,%9.5g] %8.0f  %s"
          % (name, t, med, lo, hi, e, "yes" if ok else "NO"))
if "ucldStdev" in col:
    x = sorted(v("ucldStdev")); med = x[len(x)//2]; lo, hi = hpd(x); e = ess(v("ucldStdev"))
    ok = lo <= TRUE_DISPERSION <= hi
    if not ok: bad += 1
    print("%-14s %10.5g %11.5g  [%9.5g,%9.5g] %8.0f  %s"
          % ("ucldStdev", TRUE_DISPERSION, med, lo, hi, e, "yes" if ok else "NO"))
print()
for name in EXTRA:
    if name not in col: continue
    x = sorted(v(name)); med = x[len(x)//2]; lo, hi = hpd(x); e = ess(v(name))
    print("%-14s %10s %11.5g  [%9.5g,%9.5g] %8.0f" % (name, "-", med, lo, hi, e))
print()
for name, want in [("nPainted.X1",4),("nPainted.X2",19),("nPainted.background",75),("nPainted.multiple",0)]:
    if name in col:
        u = sorted(set(int(x) for x in v(name)))
        ok = u == [want]
        if not ok: bad += 1
        print("  %-24s %-12s %s" % (name, u, "ok" if ok else "EXPECTED %d" % want))
for name in [c for c in hdr if c.startswith("monophyletic.")]:
    u = sorted(set(int(x) for x in v(name)))
    ok = u == [1]
    if not ok: bad += 1
    print("  %-24s %-12s %s" % (name, u, "ok" if ok else "CONSTRAINT VIOLATED"))
print("\n%s" % ("all known quantities covered" if bad == 0 else "%d problem(s)" % bad))
