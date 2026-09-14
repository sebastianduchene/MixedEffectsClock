#!/usr/bin/env Rscript
# One page per sampled tree, every branch coloured by its mixed-effects design column
# and labelled with the rate BEAST gave it.
#
#   Rscript plot_rate_trees.R <stem> ["title"]
#
# Inputs come from check_rates.py, which recomputes the design independently of the Java
# code and verifies every rate before writing these files:
#   <stem>_plot_all.newick  one tree per line, internal nodes labelled n1..nK
#   <stem>_nodes_all.csv    tree, label, design column, rate
#   <stem>_cols.csv         column -> legend label, including overlapping combinations
#   <stem>.trees            read only to recover the MCMC state number of each sample
#
# Where columns overlap, a branch carries the SUM of their coefficients, so the fold
# change is the product. Those combinations appear in the legend in their own right.

library(ape)
a    <- commandArgs(trailingOnly = TRUE)
stem <- if (length(a) >= 1) a[1] else "topology_rate_map_6taxon"
main <- if (length(a) >= 2) a[2] else "Mixed-effects design under a moving topology"

trees <- read.tree(paste0(stem, "_plot_all.newick"))
if (inherits(trees, "phylo")) trees <- structure(list(trees), class = "multiPhylo")
nodes <- read.csv(paste0(stem, "_nodes_all.csv"), stringsAsFactors = FALSE)
legd  <- read.csv(paste0(stem, "_cols.csv"),      stringsAsFactors = FALSE)

states <- grep("^tree STATE_", readLines(paste0(stem, ".trees")), value = TRUE)
states <- sub("^tree STATE_([0-9]+).*$", "\\1", states)

lev <- legd$col
pal <- c("grey60", "#0072B2", "#D55E00", "#009E73",
         "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999")[seq_along(lev)]
names(pal) <- lev

pdf(paste0(stem, "_rates.pdf"), width = 12, height = 7)
for (i in seq_along(trees)) {
  phy <- trees[[i]]
  nd  <- nodes[nodes$tree == i - 1, ]
  lab <- c(phy$tip.label, phy$node.label)
  idx <- match(lab, nd$label)
  stopifnot(!any(is.na(idx)))

  phy   <- ladderize(phy, right = FALSE)
  child <- phy$edge[, 2]
  ecol  <- pal[nd$col[idx][child]]
  erate <- nd$rate[idx][child]
  n     <- table(factor(nd$col[idx][child], levels = lev))

  par(mar = c(4.5, 1, 4.5, 16), xpd = NA)
  plot(phy, edge.color = ecol, edge.width = 5, cex = 1.0,
       label.offset = max(node.depth.edgelength(phy)) * 0.02)
  edgelabels(formatC(erate, format = "g", digits = 3),
             frame = "none", adj = c(0.5, -0.6), cex = 0.7, col = "grey15")
  axisPhylo(backward = TRUE)
  mtext("time before present", side = 1, line = 2.6, cex = 0.9)
  title(sprintf("%s\nsample %d of %d, MCMC state %s", main, i, length(trees),
                if (i <= length(states)) states[i] else "?"),
        cex.main = 0.95)
  # legend goes in the right margin, never over the tree
  shown <- lev[as.integer(n[lev]) > 0]
  u <- par("usr")
  legend(x = u[2] + 0.06 * diff(u[1:2]), y = u[4], bty = "n", cex = 0.75,
         lwd = 5, seg.len = 1.6, col = pal[shown], title = "design column",
         title.adj = 0, xpd = NA,
         legend = sprintf("%s  [%d branch%s]", legd$label[match(shown, legd$col)],
                          as.integer(n[shown]),
                          ifelse(as.integer(n[shown]) == 1, "", "es")))
}
invisible(dev.off())

cat("wrote", paste0(stem, "_rates.pdf"), "with", length(trees), "pages\n")
