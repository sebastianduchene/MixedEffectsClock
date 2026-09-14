#!/usr/bin/env Rscript
# Draw the tree with every branch coloured by its mixed-effects design column and
# labelled with the rate BEAST gave it.
#
#   Rscript plot_rate_tree.R <stem> ["title"]
#
# Inputs come from check_rates.py, which recomputes the design independently of the
# Java code and verifies the rates before writing them:
#   <stem>_plot.newick   internal nodes labelled n1..nK
#   <stem>_nodes.csv     label, design column, rate for every node
#   <stem>_cols.csv      column -> legend label
#
# Colour is the fixed effect. Under the Bletsa form the branch rate is
#   r_i = exp(beta_0 + sum_k X_ik beta_k) * L_i
# so a colour says which beta_k is added to the background log-rate. Here every L_i
# is 1, so each branch shows its class exactly.

library(ape)
a    <- commandArgs(trailingOnly = TRUE)
stem <- if (length(a) >= 1) a[1] else "branch_rate_map_6taxon_run"
main <- if (length(a) >= 2) a[2] else "Mixed-effects clock: branch rates by design column"

phy   <- read.tree(paste0(stem, "_plot.newick"))
nodes <- read.csv(paste0(stem, "_nodes.csv"), stringsAsFactors = FALSE)
legd  <- read.csv(paste0(stem, "_cols.csv"),  stringsAsFactors = FALSE)

ntip <- Ntip(phy)
lab  <- c(phy$tip.label, phy$node.label)
idx  <- match(lab, nodes$label)
stopifnot(!any(is.na(idx)))
node_col  <- nodes$col[idx]
node_rate <- nodes$rate[idx]

lev <- legd$col
pal <- c("grey60", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")[seq_along(lev)]
names(pal) <- lev

phy   <- ladderize(phy, right = FALSE)
child <- phy$edge[, 2]
ecol  <- pal[node_col[child]]
erate <- node_rate[child]
n     <- table(factor(node_col[child], levels = lev))

draw <- function() {
  par(mar = c(4.5, 1, 4, 12), xpd = NA)
  plot(phy, edge.color = ecol, edge.width = 5, cex = 1.0,
       label.offset = max(node.depth.edgelength(phy)) * 0.02)
  edgelabels(formatC(erate, format = "g", digits = 3),
             frame = "none", adj = c(0.5, -0.6), cex = 0.72, col = "grey15")
  axisPhylo(backward = TRUE)
  mtext("time before present", side = 1, line = 2.6, cex = 0.9)
  title(paste0(main,
               "\ncolour = design column, number on branch = rate from BEAST"),
        cex.main = 0.95)
  legend("topleft", inset = c(0, 0.04), bty = "n", cex = 0.8, lwd = 5,
         seg.len = 1.6, col = pal[lev], title = "design column", title.adj = 0,
         legend = sprintf("%s  [%d branch%s]", legd$label, as.integer(n[lev]),
                          ifelse(as.integer(n[lev]) == 1, "", "es")))
}

pdf(paste0(stem, "_rates.pdf"), width = 10, height = 7); draw(); invisible(dev.off())

cat("branches per design column:\n"); print(n)
cat(sprintf("total non-root branches: %d (expected %d)\n", sum(n), 2 * ntip - 2))
cat("wrote", paste0(stem, "_rates.pdf"), "\n")
