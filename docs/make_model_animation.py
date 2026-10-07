#!/usr/bin/env python3
"""Animated explainer for the mixed-effects clock, for the README.

Draws, side by side, a clade-coloured tree and the per-clade lognormal
distributions each branch rate is drawn from, then loops by resampling the
per-branch random effect e_i so the dots slide and the branches re-thicken --
the model's claim that every branch rate is a *draw* from its clade's
distribution, made visible.

    log r_i = beta_0 + sum_k X_ik beta_k + e_i ,   e_i ~ Normal(0, sigma^2)

Regenerate:  python3 docs/make_model_animation.py
Writes:      docs/mixed_effects_clock.gif
No data or BEAST needed. Mirrors examples/gen_sim_re_slowdown.py's known truth:
background 0.005, clade A x2 (speed-up), clade B x0.3 (slow-down), clade C x1.4.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation, PillowWriter

OUT = os.path.join(os.path.dirname(__file__), "mixed_effects_clock.gif")

# ---- palette (reads on white; matches the interactive illustration) ----------
COL = {"bg": "#8c98a7", "A": "#dd5430", "B": "#3b7cc9", "C": "#2f9c6a"}
FG, MUT, PANEL, LINE = "#1a2330", "#586675", "#ffffff", "#d6dde4"
EPS = "#6f52c4"   # the per-branch random effect, kept distinct from any clade hue

BG_RATE = 0.005
FOLD = {"bg": 1.0, "A": 2.0, "B": 0.3, "C": 1.4}   # exp(beta_k)
SIGMA = 0.33                                        # log-scale dispersion
STEM_INCL = {"A": True, "B": False, "C": False}     # includeStem per clade

# ---- tree topology (nested); leaves tagged by clade --------------------------
_uid = [0]
def _id():
    _uid[0] += 1
    return _uid[0]

def L(clade):
    return {"id": _id(), "clade": clade, "children": None}

def Nn(*ch, mrca=None):
    return {"id": _id(), "children": list(ch), "mrca": mrca}

cladeA = Nn(Nn(L("A"), L("A")), Nn(L("A"), L("A")), mrca="A")
cladeB = Nn(Nn(L("B"), L("B")), Nn(L("B"), L("B")), mrca="B")
cladeC = Nn(Nn(L("C"), L("C")), L("C"), mrca="C")
root = Nn(cladeA, Nn(L("bg"), Nn(Nn(cladeB, cladeC), L("bg"))))

def tag(n):
    if n["children"] is None:
        n["cl"] = n["clade"]
    else:
        cs = [tag(c) for c in n["children"]]
        n["cl"] = cs[0] if (all(c == cs[0] for c in cs) and cs[0] != "bg") else "bg"
    return n["cl"]
tag(root)

_maxd = [0]
def depth(n):
    n["d"] = 0 if n["children"] is None else 1 + max(depth(c) for c in n["children"])
    _maxd[0] = max(_maxd[0], n["d"])
    return n["d"]
depth(root)
MAXD = _maxd[0]

_order = [0]
def place(n):
    if n["children"] is None:
        n["y"] = _order[0] + 0.5
        _order[0] += 1
    else:
        for c in n["children"]:
            place(c)
        n["y"] = (n["children"][0]["y"] + n["children"][-1]["y"]) / 2.0
    n["x"] = (MAXD - n["d"]) / MAXD
place(root)
NLEAF = _order[0]

edges = []   # every node except root
def collect(n, parent):
    if parent is not None:
        c = n["cl"]
        if n.get("mrca") and not STEM_INCL[n["mrca"]]:
            c = "bg"
        edges.append({"node": n, "parent": parent, "clade": c})
    if n["children"] is not None:
        for ch in n["children"]:
            collect(ch, n)
collect(root, None)

# ---- rate axis (log) ---------------------------------------------------------
rmin, rmax = BG_RATE * 0.14, BG_RATE * 5.2
lmin, lmax = np.log(rmin), np.log(rmax)
def centre(cl):
    return BG_RATE * FOLD[cl]
def width_of(r):
    t = np.clip((np.log(r) - lmin) / (lmax - lmin), 0, 1)
    return 1.3 + t * (7.0 - 1.3)

# ---- resample keyframes, interpolated for a seamless loop --------------------
rng = np.random.default_rng(7)
ne = len(edges)
KEYS = 5                       # distinct draws
STEP = 16                      # frames between draws
zkeys = rng.standard_normal((KEYS, ne))
zkeys[-1] = zkeys[0]           # close the loop
FRAMES = (KEYS - 1) * STEP

def z_at(frame):
    seg = frame // STEP
    f = (frame % STEP) / STEP
    f = f * f * (3 - 2 * f)    # smoothstep ease
    return (1 - f) * zkeys[seg] + f * zkeys[seg + 1]

# ---- figure ------------------------------------------------------------------
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
fig = plt.figure(figsize=(9.2, 4.5), dpi=112)
fig.patch.set_facecolor(PANEL)
gs = fig.add_gridspec(1, 2, width_ratios=[0.92, 1.08], left=0.02, right=0.985,
                      top=0.78, bottom=0.115, wspace=0.08)
axT = fig.add_subplot(gs[0]); axR = fig.add_subplot(gs[1])
for ax in (axT, axR):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

fig.text(0.02, 0.955, "The mixed-effects clock: every branch rate is a draw",
         ha="left", va="center", fontsize=13, color=FG, weight="bold")

# colour-coded equation — each parameter in the colour of its role in the figure:
#   beta_0 grey (background), clade effect in the clade hue, epsilon_i / sigma violet.
def colored_eq(x, y, frags, fs):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for s, c, w in frags:
        t = fig.text(x, y, s, color=c, fontsize=fs, ha="left", va="center",
                     weight=w)
        bb = t.get_window_extent(renderer=rend)
        x += inv.transform((bb.width, 0))[0] - inv.transform((0, 0))[0]

colored_eq(0.02, 0.895, [
    (r"$\log r_i = $", FG, "normal"),
    (r"$\beta_0$", COL["bg"], "bold"),
    (r"$\, + \,$", MUT, "normal"),
    (r"$\sum_k X_{ik}\,\beta_k$", COL["A"], "bold"),
    (r"$\, + \,$", MUT, "normal"),
    (r"$\varepsilon_i$", EPS, "bold"),
    (r"$,\ \ \ \varepsilon_i \sim \mathrm{Normal}(0,\, $", MUT, "normal"),
    (r"$\sigma^2$", EPS, "bold"),
    (r"$\,)$", MUT, "normal"),
], 11.5)

axT.set_title("branches coloured by clade", fontsize=10, color=MUT, pad=5)
axR.set_title("drawn from these distributions", fontsize=10, color=MUT, pad=5)

# --- tree: static connectors + tips, dynamic branch lines ---
def TX(x):
    return 0.09 + x * 0.86
def TY(y):
    return 0.93 - (y / NLEAF) * 0.86
def draw_conn(n):
    if n["children"] is not None:
        ys = [TY(c["y"]) for c in n["children"]]
        axT.plot([TX(n["x"])] * 2, [min(ys), max(ys)], color=LINE, lw=1.1, zorder=1)
        for c in n["children"]:
            draw_conn(c)
draw_conn(root)
axT.plot([0.03, TX(root["x"])], [TY(root["y"])] * 2, color=LINE, lw=1.1, zorder=1)
branch_lines = []
for e in edges:
    y = TY(e["node"]["y"])
    ln = Line2D([TX(e["parent"]["x"]), TX(e["node"]["x"])], [y, y],
                color=COL[e["clade"]], lw=2, solid_capstyle="round", zorder=2)
    axT.add_line(ln); branch_lines.append(ln)
for e in edges:
    if e["node"]["children"] is None:
        axT.plot(TX(e["node"]["x"]) + 0.004, TY(e["node"]["y"]), "o",
                 ms=3, color=COL[e["clade"]], zorder=3)
axT.set_xlim(0, 1); axT.set_ylim(0, 1)

# --- ridgeline: static bells + axis, dynamic dots ---
ROWS = ["A", "C", "bg", "B"]        # top -> bottom, descending centre
def RX(lr):
    return 0.17 + (lr - lmin) / (lmax - lmin) * 0.80
def RY(i, h=0.0):
    return 0.86 - i * 0.205 - h
ROWH = 0.185
k = ROWH * 1.0 * SIGMA * np.sqrt(2 * np.pi)
CAP = ROWH * 1.35
def pdf(lr, mu):
    return np.exp(-0.5 * ((lr - mu) / SIGMA) ** 2) / (SIGMA * np.sqrt(2 * np.pi))
# fold-change grid
for fac in (0.25, 0.5, 1, 2, 4):
    x = RX(np.log(BG_RATE * fac))
    axR.plot([x, x], [0.08, 0.90], color=LINE, lw=1.3 if fac == 1 else 0.8,
             ls="-" if fac == 1 else (0, (2, 3)), zorder=1)
    axR.text(x, 0.045, "×1" if fac == 1 else f"×{fac:g}", ha="center", va="center",
             fontsize=8.5, color=MUT)
axR.text(0.57, 0.005, "branch rate  (relative to background, log scale)",
         ha="center", va="center", fontsize=8.5, color=MUT)
lrs = np.linspace(lmin, lmax, 160)
for i, cl in enumerate(ROWS):
    mu = np.log(centre(cl)); base = RY(i)
    h = np.minimum(pdf(lrs, mu) * k, CAP)
    axR.fill_between([RX(l) for l in lrs], base, base + h, color=COL[cl],
                     alpha=0.16, zorder=2)
    axR.plot([RX(l) for l in lrs], base + h, color=COL[cl], lw=1.5, alpha=0.9, zorder=3)
    axR.plot([RX(mu)] * 2, [base, base + min(pdf(mu, mu) * k, CAP)],
             color=COL[cl], lw=1.0, ls=(0, (1, 2)), alpha=0.7, zorder=3)
    name = {"A": "Clade A", "C": "Clade C", "bg": "Background", "B": "Clade B"}[cl]
    axR.text(0.155, base + 0.028, name, ha="right", va="center", fontsize=9.5, color=FG)
    fac = FOLD[cl]
    axR.text(0.155, base - 0.005, "×1" if cl == "bg" else f"×{fac:g}",
             ha="right", va="center", fontsize=8, color=MUT)
axR.set_xlim(0, 1); axR.set_ylim(0, 1)

# dynamic dots, grouped by clade row
dot_artists = []   # (PathCollection, [edge idx], base, jitter array)
for i, cl in enumerate(ROWS):
    idx = [j for j, e in enumerate(edges) if e["clade"] == cl]
    base = RY(i)
    jit = np.array([((m % 3) - 1) * 0.011 for m in range(len(idx))])
    sc = axR.scatter(np.zeros(len(idx)), np.full(len(idx), base + 0.02) + jit,
                     s=34, color=COL[cl], edgecolor=PANEL, linewidth=0.6, zorder=5)
    dot_artists.append((sc, idx, base, jit))

def update(frame):
    z = z_at(frame)
    rates = np.array([centre(e["clade"]) * np.exp(SIGMA * z[j]) for j, e in enumerate(edges)])
    for ln, r in zip(branch_lines, rates):
        ln.set_linewidth(width_of(r))
    for sc, idx, base, jit in dot_artists:
        xs = np.array([RX(np.log(rates[j])) for j in idx])
        sc.set_offsets(np.column_stack([xs, np.full(len(idx), base + 0.02) + jit]))
    return branch_lines + [sc for sc, *_ in dot_artists]

anim = FuncAnimation(fig, update, frames=FRAMES, interval=55, blit=False)
anim.save(OUT, writer=PillowWriter(fps=18))
print(f"wrote {OUT}  ({FRAMES} frames, {os.path.getsize(OUT)/1024:.0f} KB)")
