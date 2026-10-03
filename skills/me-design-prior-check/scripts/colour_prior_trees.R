#!/usr/bin/env Rscript
# ============================================================================
# colour_prior_trees.R  --  colour N sampled trees from a prior-only run by
# the fixed-effect DESIGN COLUMN (clade) that each branch carries.
#
# Usage:
#   Rscript colour_prior_trees.R <xml> <trees> <out.pdf> [n_trees] [burnin_frac]
#
# <xml>          BEAST2 XML with a mixedeffectsclock.MixedEffectsClockModel
#                block (CladeDesign children: category, includeStem, taxonset
#                idref). The design is parsed from the XML, not hardcoded.
# <trees>        BEAST .trees file (robust to a still-writing final line).
# <out.pdf>      output PDF; one page per sampled tree, plus a summary page.
# [n_trees]      number of trees to draw (default 10).
# [burnin_frac]  fraction of trees to discard before picking (default 0.1).
#
# Needs only ape. Design rule (same as the sibling colour-me-clock-tree skill):
# a branch is in column k if its CHILD node is a strict descendant of clade k's
# MRCA; the stem branch (child = MRCA) is added only if includeStem. Branches
# in no clade carry the intercept (background).
# ============================================================================
suppressPackageStartupMessages(library(ape))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("usage: colour_prior_trees.R <xml> <trees> <out.pdf> [n_trees] [burnin_frac]")
xml_path   <- a[1]; trees_path <- a[2]; out_pdf <- a[3]
n_trees    <- if (length(a) >= 4) as.integer(a[4])   else 10L
burnin_frac<- if (length(a) >= 5) as.numeric(a[5])   else 0.1

# ---- 1. read the design from the XML ---------------------------------------
# Handles both idioms: a <clade> whose <taxonset> is referenced by idref (as in
# the lepromatosis analyses), and a <clade> that defines its <taxonset>
# inline (as in the sibling skill's six-taxon template).
xml <- paste(readLines(xml_path, warn = FALSE), collapse = "\n")
clades <- list()

# Grab each full <clade ...>...</clade> block that declares CladeDesign.
clade_blocks <- regmatches(
  xml,
  gregexpr('(?s)<clade\\b[^>]*spec="mixedeffectsclock\\.CladeDesign"[^>]*>.*?</clade>',
           xml, perl = TRUE))[[1]]

# Fallback: self-closing <clade .../> is rare but legal.
if (!length(clade_blocks))
  clade_blocks <- regmatches(
    xml,
    gregexpr('(?s)<clade\\b[^>]*spec="mixedeffectsclock\\.CladeDesign"[^>]*/>',
             xml, perl = TRUE))[[1]]

taxa_from_taxonset_block <- function(tsblk) {
  m <- regmatches(tsblk,
                  gregexpr('<taxon\\b[^>]*(?:idref|id)="([^"]+)"', tsblk, perl = TRUE))[[1]]
  tx <- sub('.*"([^"]+)"', "\\1", m)
  tx[!grepl("^taxa\\.", tx)]
}

for (blk in clade_blocks) {
  nm   <- sub('.*<clade\\b[^>]*\\bid="([^"]+)".*', "\\1", blk)
  nm   <- sub("^design\\.?", "", nm)
  cat_m <- regmatches(blk, regexpr('category="(\\d+)"', blk, perl = TRUE))
  cat_  <- as.integer(sub('.*"(\\d+)".*', "\\1", cat_m))
  stem    <- grepl('includeStem="true"', blk)
  incTerm <- !grepl('includeTerminal="false"', blk)   # default true
  exclCl  <- grepl('excludeClade="true"', blk)        # default false

  # taxon set: inline <taxonset id="..." ...>...</taxonset>, or <taxonset idref="...">
  tsblk <- regmatches(blk,
                      regexpr('(?s)<taxonset\\b[^>]*>.*?</taxonset>', blk, perl = TRUE))
  if (length(tsblk) && nzchar(tsblk)) {
    taxa <- taxa_from_taxonset_block(tsblk)
  } else {
    ref <- sub('.*<taxonset\\b[^>]*idref="([^"]+)".*', "\\1", blk)
    tsblk2 <- regmatches(
      xml,
      regexpr(sprintf('(?s)<taxonset\\b[^>]*\\bid="%s"[^>]*>.*?</taxonset>',
                      gsub("([.])", "\\\\\\1", ref)), xml, perl = TRUE))
    taxa <- if (length(tsblk2) && nzchar(tsblk2))
              taxa_from_taxonset_block(tsblk2) else character(0)
  }
  clades[[nm]] <- list(cat = cat_, stem = stem, includeTerminal = incTerm,
                       excludeClade = exclCl, taxa = taxa)
}
if (!length(clades))
  stop("no CladeDesign blocks found in ", xml_path,
       " -- is this a BEAST 2 MixedEffectsClock XML?")
clades <- clades[order(sapply(clades, function(c) c$cat))]
cat(sprintf("design: %d columns\n", length(clades)))
for (nm in names(clades))
  cat(sprintf("  cat %d  %-14s stem=%s terminal=%s excludeInternal=%s  (%d taxa)\n",
              clades[[nm]]$cat, nm, clades[[nm]]$stem, clades[[nm]]$includeTerminal,
              clades[[nm]]$excludeClade, length(clades[[nm]]$taxa)))

# ---- 2. read the sampled trees ---------------------------------------------
raw <- readLines(trees_path, warn = FALSE)
ti  <- grep("^\\s*[Tt]ranslate", raw)
t1  <- grep("^\\s*tree ", raw)[1]
tmap <- list()
if (length(ti) && !is.na(t1)) {
  for (ln in raw[(ti + 1):(t1 - 1)]) {
    m <- regmatches(ln, regexec("^\\s*(\\d+)\\s+([^,;]+)", ln))[[1]]
    if (length(m) == 3) tmap[[m[2]]] <- gsub("['\" ]", "", m[3])
  }
}
treelines <- raw[grepl("^\\s*tree ", raw) & grepl(";\\s*$", raw)]
if (!length(treelines)) stop("no complete tree lines in ", trees_path)

# discard burnin, then pick n_trees evenly spaced
n_total  <- length(treelines)
first    <- max(1L, floor(n_total * burnin_frac) + 1L)
usable   <- treelines[first:n_total]
if (length(usable) < 1) stop("burn-in left no trees")
n_draw   <- min(n_trees, length(usable))
picks    <- round(seq(1, length(usable), length.out = n_draw))
states   <- as.numeric(sub(".*STATE_(\\d+).*", "\\1", usable[picks]))
cat(sprintf("read %d trees, discarding first %d as burn-in; drawing %d\n",
            n_total, first - 1L, n_draw))

read_one <- function(line) {
  nwk <- sub("^[^(]*", "", line)                 # drop everything before first (
  nwk <- gsub("\\[&[^]]*\\]", "", nwk)           # strip [&rate=...] metadata
  t   <- read.tree(text = nwk)
  if (length(tmap)) {
    nm <- unlist(tmap[t$tip.label])
    if (!any(is.na(nm))) t$tip.label <- nm
  }
  t
}
trees <- lapply(usable[picks], read_one)

# ---- 3. assign the design column per edge ----------------------------------
# Steiner subtree: only branches on paths from MRCA to taxon tips are
# candidates. For monophyletic groups this equals the full subtree.
descendants <- function(tree, m) {
  kids <- function(n) tree$edge[tree$edge[, 1] == n, 2]
  out <- c(); st <- m
  while (length(st)) { n <- st[1]; st <- st[-1]; out <- c(out, n); st <- c(st, kids(n)) }
  out
}

steiner_set <- function(tree, mrca_node, tip_idx) {
  n_tips <- Ntip(tree)
  max_node <- max(tree$edge)
  tip_in_set <- logical(max_node); tip_in_set[tip_idx] <- TRUE
  has_taxon  <- logical(max_node); has_taxon[tip_idx]  <- TRUE

  children <- vector("list", max_node)
  for (i in seq_len(nrow(tree$edge)))
    children[[tree$edge[i, 1]]] <- c(children[[tree$edge[i, 1]]], tree$edge[i, 2])

  post_order <- function(node) {
    if (node <= n_tips) return(node)
    result <- c()
    for (ch in children[[node]]) result <- c(result, post_order(ch))
    c(result, node)
  }
  nodes <- post_order(mrca_node)

  for (nd in nodes) {
    if (nd <= n_tips) next
    for (ch in children[[nd]])
      if (has_taxon[ch]) { has_taxon[nd] <- TRUE; break }
  }

  on_steiner <- logical(max_node)
  for (nd in nodes) if (has_taxon[nd]) on_steiner[nd] <- TRUE
  on_steiner[mrca_node] <- FALSE
  on_steiner
}

assign_design <- function(tree) {
  ec   <- rep(NA_integer_, nrow(tree$edge))
  mono <- setNames(logical(length(clades)), names(clades))
  overlap <- c()
  n_tips <- Ntip(tree)
  for (nm in names(clades)) {
    cl   <- clades[[nm]]
    tp   <- match(cl$taxa, tree$tip.label)
    if (any(is.na(tp))) {
      cat(sprintf("  ** %d taxa of clade %s not in tree\n", sum(is.na(tp)), nm))
      next
    }
    mnode <- if (length(tp) == 1) tp else getMRCA(tree, cl$taxa)
    on_st <- steiner_set(tree, mnode, tp)

    ie <- c()
    for (i in seq_len(nrow(tree$edge))) {
      child <- tree$edge[i, 2]
      if (!on_st[child]) next
      is_leaf <- child <= n_tips
      if (is_leaf) {
        if (cl$includeTerminal && tree$tip.label[child] %in% cl$taxa)
          ie <- c(ie, i)
      } else {
        if (!cl$excludeClade)
          ie <- c(ie, i)
      }
    }
    if (cl$stem) ie <- c(ie, which(tree$edge[, 2] == mnode))
    if (any(!is.na(ec[ie]))) overlap <- c(overlap, nm)
    ec[ie] <- cl$cat
    desc <- descendants(tree, mnode)
    mono[nm] <- setequal(tree$tip.label[desc[desc <= n_tips]], cl$taxa)
  }
  list(edge_cat = ec, mono = mono, overlap = overlap)
}

# ---- 4. draw ---------------------------------------------------------------
pal  <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
          "#A65628", "#F781BF", "#999999")
cats <- sapply(clades, function(c) c$cat)
names_by_cat <- names(clades)[order(cats)]
sorted_cats  <- sort(cats)
col_by_cat   <- setNames(pal[seq_along(sorted_cats)], as.character(sorted_cats))

draw_one <- function(tree, res, idx_shown, idx_total, state) {
  ecol <- rep("grey70", nrow(tree$edge))
  for (k in sorted_cats) ecol[which(res$edge_cat == k)] <- col_by_cat[as.character(k)]
  tipcol <- rep("grey60", Ntip(tree))
  for (nm in names(clades))
    tipcol[match(clades[[nm]]$taxa, tree$tip.label)] <-
      col_by_cat[as.character(clades[[nm]]$cat)]

  nb <- sapply(sorted_cats, function(k) sum(res$edge_cat == k, na.rm = TRUE))
  bg <- sum(is.na(res$edge_cat))
  mono_str <- paste(sprintf("%s=%s", names(res$mono), ifelse(res$mono, "T", "F")),
                    collapse = "  ")

  plot(tree, edge.color = ecol, edge.width = 2.2, tip.color = tipcol,
       cex = 0.5, no.margin = FALSE,
       main = sprintf("sample %d/%d (STATE_%s)  bg=%d  %s",
                      idx_shown, idx_total, state, bg,
                      paste(sprintf("%s=%d", names_by_cat, nb), collapse = "  ")))
  mtext(mono_str, side = 1, line = 2, cex = 0.7)
  legend("bottomleft", bty = "n", cex = 0.7, lwd = 2.4,
         legend = c("background (intercept only)",
                    sprintf("cat %d: %s%s", sorted_cats, names_by_cat,
                            sapply(clades[names_by_cat], function(c) {
                              flags <- c()
                              if (c$stem) flags <- c(flags, "+stem")
                              if (!c$includeTerminal) flags <- c(flags, "-term")
                              if (c$excludeClade) flags <- c(flags, "-internal")
                              if (length(flags)) paste0(" (", paste(flags, collapse=","), ")")
                              else ""
                            }))),
         col = c("grey70", col_by_cat[as.character(sorted_cats)]))
}

pdf(out_pdf, width = 9, height = 11, pointsize = 10)
# summary page: branches-per-column across the sampled trees, monophyly rate
res_all   <- lapply(trees, assign_design)
nb_matrix <- sapply(res_all, function(r)
  c(bg = sum(is.na(r$edge_cat)),
    sapply(sorted_cats, function(k) sum(r$edge_cat == k, na.rm = TRUE))))
rownames(nb_matrix)[-1] <- names_by_cat
mono_rate <- rowMeans(sapply(res_all, function(r) r$mono))
overlap_any <- sum(sapply(res_all, function(r) length(r$overlap) > 0))

plot.new()
title("Design check across sampled prior trees")
y <- 0.95
mtext_lines <- function(lines, start_y) {
  for (i in seq_along(lines))
    text(0.02, start_y - 0.045 * (i - 1), lines[i], adj = 0, cex = 0.9, family = "mono")
}
summary_lines <- c(
  sprintf("XML:    %s", basename(xml_path)),
  sprintf("Trees:  %s   (%d sampled, after burn-in)",
          basename(trees_path), length(trees)),
  "",
  "Branches per column (mean across sampled trees):",
  sprintf("  %-14s %6.2f", "background", mean(nb_matrix["bg", ])),
  sapply(names_by_cat, function(nm) {
    cl <- clades[[nm]]
    flags <- c()
    if (cl$stem) flags <- c(flags, "+stem")
    if (!cl$includeTerminal) flags <- c(flags, "-term")
    if (cl$excludeClade) flags <- c(flags, "-internal")
    fstr <- if (length(flags)) paste0(", ", paste(flags, collapse=",")) else ""
    sprintf("  %-14s %6.2f   (cat %d%s)", nm, mean(nb_matrix[nm, ]), cl$cat, fstr)
  }),
  "",
  "Monophyly rate (fraction of sampled trees where the clade's MRCA equals its taxa set):",
  sapply(names(mono_rate), function(nm)
    sprintf("  %-14s %6.2f", nm, mono_rate[nm])),
  "",
  sprintf("Trees with any OVERLAP (branches in >1 column): %d / %d",
          overlap_any, length(trees))
)
mtext_lines(summary_lines, 0.92)
mtext_lines(c("Each following page shows one sampled tree, branches coloured by",
              "the column that carries them. Monophyly bools printed under each tree."),
            0.08)

for (i in seq_along(trees))
  draw_one(trees[[i]], res_all[[i]], i, length(trees), states[i])

invisible(dev.off())

# ---- 5. echo a short per-column report to stdout ---------------------------
cat("\n=== design check: per-column branch counts across sampled trees ===\n")
for (rn in rownames(nb_matrix))
  cat(sprintf("  %-14s mean=%6.2f  min=%d  max=%d\n", rn,
              mean(nb_matrix[rn, ]), min(nb_matrix[rn, ]), max(nb_matrix[rn, ])))
cat("\n=== monophyly rate across sampled trees ===\n")
for (nm in names(mono_rate))
  cat(sprintf("  %-14s %5.1f%%\n", nm, 100 * mono_rate[nm]))
if (overlap_any > 0) {
  cat(sprintf("\n** OVERLAP in %d/%d trees -- the design is not exclusive\n",
              overlap_any, length(trees)))
}
cat("\nwrote ", out_pdf, "\n", sep = "")
