#!/usr/bin/env python3
"""
Independent check that BEAST put the right rate on every branch.

Reads the clade design, the coefficients and the intercept straight out of the XML, so
nothing here restates a constant by hand, then recomputes the expected rate for every
branch of every sampled tree and compares it with the rate BEAST annotated. Also writes
the files the plotting script needs.

    python3 check_rates.py <run stem>

Expects <stem>.xml and <stem>.trees; writes <stem>_plot.newick, <stem>_nodes.csv and
<stem>_cols.csv.
"""
import csv, math, sys, os
import xml.etree.ElementTree as ET
import dendropy

stem = sys.argv[1] if len(sys.argv) > 1 else "branch_rate_map_6taxon_run"
TOL = 1e-9

# ---------------------------------------------------------------- read the model
root = ET.parse(stem + ".xml").getroot()
byid = {e.get("id"): e for e in root.iter() if e.get("id")}

def param(pid):
    return [float(x) for x in byid[pid].text.split()]

clock_rate = param("clockRate")[0]
beta       = param("coefficient")
multiplier = param("rates")

clades = []
for e in root.iter("clade"):
    taxa = [t.get("idref") or t.get("id") for t in e.iter("taxon")]
    clades.append(dict(
        name=e.get("id"),
        col=int(e.get("category")),
        stem=e.get("includeStem") == "true",
        excl=e.get("excludeClade") == "true",
        taxa=[t for t in taxa if t]))

print("from %s.xml" % stem)
print("  intercept (clock.rate) = %g" % clock_rate)
for i, b in enumerate(beta):
    print("  beta[%d] = %+.6f   fold change %g" % (i, b, math.exp(b)))
print("  branch multipliers all 1: %s" % all(abs(m - 1.0) < TOL for m in multiplier))
for c in clades:
    print("  column %d  %-8s taxa=%s stem=%s excludeClade=%s"
          % (c["col"], c["name"], ",".join(c["taxa"]), c["stem"], c["excl"]))
print()

# ---------------------------------------------------------------- read the trees
# rooting="force-rooted" is essential. Without it dendropy treats the tree as unrooted
# and collapses the basal node on read, which silently removes one branch: the first
# version of this script "passed" because the stem-only column was simply not there.
trees = dendropy.TreeList.get(path=stem + ".trees", schema="nexus",
                              extract_comment_metadata=True,
                              preserve_underscores=True,
                              rooting="force-rooted")
ntax = len(trees.taxon_namespace)
expected_branches = 2 * ntax - 2
print("read %d sampled trees, %d taxa, so %d branches expected per tree"
      % (len(trees), ntax, expected_branches))

def design_of(tree):
    """node -> set of design columns, recomputed independently of the Java code."""
    cols = {nd: set() for nd in tree.preorder_node_iter()}
    for c in clades:
        mrca = tree.mrca(taxon_labels=c["taxa"])
        if mrca is None:
            raise SystemExit("no MRCA for " + str(c["taxa"]))
        if not c["excl"]:
            for nd in mrca.preorder_iter():
                if nd is not mrca:
                    cols[nd].add(c["col"])
        if c["stem"] and mrca.parent_node is not None:
            cols[mrca].add(c["col"])
    return cols

bad = 0
empty_columns = []
checked = 0        # real branches
root_checked = 0   # BEAST also annotates the root, with the placeholder value 1
for ti, tree in enumerate(trees):
    nodes = list(tree.preorder_node_iter())
    if len(nodes) != expected_branches + 1:
        raise SystemExit("tree %d has %d nodes, expected %d. The tree was probably "
                         "derooted on read, which would hide a branch."
                         % (ti, len(nodes), expected_branches + 1))
    cols = design_of(tree)
    if not any(cols[nd] for nd in nodes):
        raise SystemExit("tree %d: the design is empty" % ti)
    for c in clades:
        if not any(c["col"] in cols[nd] for nd in nodes):
            # Legitimate without monophyly constraints: if a clade's ancestor becomes
            # the root, a stem-only column has no branch left to mark. Worth reporting,
            # not an error. The hard check is the branch count above.
            empty_columns.append((ti, c["name"]))
    for nd in tree.preorder_node_iter():
        got = nd.annotations.get_value("rate")
        if got is None:
            continue
        got = float(got)
        if nd.parent_node is None:
            want = 1.0                       # BEAST convention: the root has no branch
            root_checked += 1
        else:
            want = clock_rate * math.exp(sum(beta[k] for k in cols[nd]))
            checked += 1
        if abs(got - want) > TOL * max(1.0, abs(want)):
            bad += 1
            label = nd.taxon.label if nd.taxon else "internal"
            print("  MISMATCH tree %d %-10s cols=%s got %.10g want %.10g"
                  % (ti, label, sorted(cols[nd]), got, want))

if empty_columns:
    print("  note: a column marked no branch in %d tree/column cases, e.g. %s. That is "
          "expected without monophyly constraints." % (len(empty_columns), empty_columns[:4]))

per_tree = checked / len(trees)
print("checked %d branch rates across %d trees (%g per tree), plus %d root placeholders: "
      "%d mismatches" % (checked, len(trees), per_tree, root_checked, bad))
if per_tree != expected_branches:
    raise SystemExit("only %g branches per tree were checked, expected %d. A branch is "
                     "missing, so this run proves nothing." % (per_tree, expected_branches))

# ---------------------------------------------------------------- files for plotting
# Every sampled tree is written, one per line, so a multi-page figure can show how the
# design follows the topology. Internal nodes are labelled n1..nK within each tree.
names = {c["col"]: c["name"] for c in clades}

def col_label(cset):
    return "background" if not cset else "+".join(names[k] for k in sorted(cset))

rows = []
newicks = []
for ti, tree in enumerate(trees):
    cols = design_of(tree)
    counter = [0]
    for nd in tree.preorder_node_iter():
        if nd.taxon:
            nid = nd.taxon.label
        else:
            counter[0] += 1
            nid = "n%d" % counter[0]
            nd.label = nid
        rate = nd.annotations.get_value("rate")
        rows.append((ti, nid,
                     "root" if nd.parent_node is None else col_label(cols[nd]),
                     float(rate) if rate is not None else float("nan")))
    for nd in tree.preorder_node_iter():
        nd.annotations.clear()
    newicks.append(tree.as_string(schema="newick", suppress_rooting=True,
                                  unquoted_underscores=True,
                                  suppress_internal_node_labels=False).strip())

with open(stem + "_plot_all.newick", "w") as fh:
    fh.write("\n".join(newicks) + "\n")
with open(stem + "_nodes_all.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["tree", "label", "col", "rate"]); w.writerows(rows)

# legend covers every column combination that actually occurred
order, seen = [], set()
for _, _, c, _ in rows:
    if c not in seen and c != "root":
        seen.add(c); order.append(c)
order.sort(key=lambda c: (c != "background", c.count("+"), c))
with open(stem + "_cols.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["col", "label"])
    for c in order:
        if c == "background":
            w.writerow([c, "background (intercept only)"])
        else:
            fold = 1.0
            for nm in c.split("+"):
                k = [d["col"] for d in clades if d["name"] == nm][0]
                fold *= math.exp(beta[k])
            w.writerow([c, "%s  x%.4g" % (c, fold)])

print("wrote %s_plot_all.newick (%d trees), %s_nodes_all.csv, %s_cols.csv"
      % (stem, len(newicks), stem, stem))
sys.exit(1 if bad else 0)
