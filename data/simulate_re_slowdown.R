#!/usr/bin/env Rscript
# A 50-taxon simulation with BOTH a genuine random effect and a genuine SLOWDOWN.
#
# Neither existing simulation tests either. The 50-taxon local clock and the 47-taxon
# lepromatosis simulation both have a true dispersion of essentially zero, so the random
# effect has only ever been checked against its prior; and neither has a clade that is
# really slower than background, which is the case the log-scale parameterisation exists
# to represent.
#
# Same tree as sim_local_MEclock_Bletsa.xml (rtree(50) under seed 22, NELSI tip dates), so
# the fast clade and its taxon names carry straight over from that build.
#
#   log r_i = log(0.005) + log(2)*fast_i + log(0.3)*slow_i + eps_i,  eps_i ~ N(0, 0.3^2)

library(NELSI); library(ape); library(phangorn)

BACKGROUND <- 0.005
FAST_FOLD  <- 2.0      # as in the existing 50-taxon build
SLOW_FOLD  <- 0.3      # the new case
SD_LOG     <- 0.3      # the new case: a real random effect
NSITE      <- 10000    # more sites than the 5000 of the older build: three quantities now
SLOW_NODE  <- 72       # 10 tips, disjoint from the fast set, ~10.4 time units of branch

set.seed(22)
tr <- rtree(50)
tr$tip.label <- paste0(tr$tip.label, "_", round(all.node.times(tr, tipsonly = TRUE), 2))

## edges of a clade, optionally including its stem -----------------------------
clade_edges <- function(phy, node, include_stem = TRUE) {
  desc <- node
  repeat {
    nw <- phy$edge[phy$edge[, 1] %in% desc, 2]
    if (all(nw %in% desc)) break
    desc <- unique(c(desc, nw))
  }
  inside <- which(phy$edge[, 2] %in% setdiff(desc, node))
  stem   <- which(phy$edge[, 2] == node)
  if (include_stem) c(stem, inside) else inside
}

fast_edges <- c(4, 5, 6, 7)                     # the BEAST X fast set, unchanged
slow_edges <- clade_edges(tr, SLOW_NODE, TRUE)
stopifnot(length(intersect(fast_edges, slow_edges)) == 0)

## true rates ------------------------------------------------------------------
fixed <- rep(BACKGROUND, nrow(tr$edge))
fixed[fast_edges] <- BACKGROUND * FAST_FOLD
fixed[slow_edges] <- BACKGROUND * SLOW_FOLD

set.seed(101)
eps  <- rlnorm(nrow(tr$edge), meanlog = -SD_LOG^2 / 2, sdlog = SD_LOG)  # mean one
rate <- fixed * eps

## simulate under Jukes-Cantor --------------------------------------------------
phylo <- tr
phylo$edge.length <- tr$edge.length * rate
set.seed(7)
seqs <- as.DNAbin(simSeq(phylo, l = NSITE))
write.dna(seqs, "sim_re_slowdown.fasta", format = "fasta", nbcol = -1, colsep = "")
write.tree(tr, "sim_re_slowdown_true.nwk")

## taxon sets the analysis will need --------------------------------------------
tips_of <- function(node) extract.clade(tr, node)$tip.label
fast_clade <- tips_of(tr$edge[4, 2])          # the 2-taxon cherry, edges 5 and 6
fast_sister <- tips_of(tr$edge[7, 2])         # 6-taxon sister, STEM ONLY
slow_clade <- tips_of(SLOW_NODE)
writeLines(c(paste0("fastClade\t", paste(fast_clade, collapse = ",")),
             paste0("fastSisterStem\t", paste(fast_sister, collapse = ",")),
             paste0("slowClade\t", paste(slow_clade, collapse = ","))),
           "sim_re_slowdown_taxa.tsv")

## truth table -------------------------------------------------------------------
realised_sd   <- sd(log(eps))
realised_mean <- mean(eps)
bg_edges <- setdiff(seq_len(nrow(tr$edge)), c(fast_edges, slow_edges))
tab <- data.frame(
  column       = c("background", "fast", "slow"),
  branches     = c(length(bg_edges), length(fast_edges), length(slow_edges)),
  fold_change  = c(1, FAST_FOLD, SLOW_FOLD),
  true_beta    = c(NA, log(FAST_FOLD), log(SLOW_FOLD)),
  mean_rate    = c(mean(rate[bg_edges]), mean(rate[fast_edges]), mean(rate[slow_edges])),
  expected_subs= c(sum(tr$edge.length[bg_edges]   * rate[bg_edges]),
                   sum(tr$edge.length[fast_edges] * rate[fast_edges]),
                   sum(tr$edge.length[slow_edges] * rate[slow_edges])) * NSITE)
write.csv(tab, "sim_re_slowdown_truth.csv", row.names = FALSE)

cat("\n=== TRUTH ===\n")
cat(sprintf("background rate  %.6g   (clock.rate, = exp(beta_0))\n", BACKGROUND))
cat(sprintf("fast  beta = log %.3g = %+.4f   on %d branches\n", FAST_FOLD, log(FAST_FOLD), length(fast_edges)))
cat(sprintf("slow  beta = log %.3g = %+.4f   on %d branches\n", SLOW_FOLD, log(SLOW_FOLD), length(slow_edges)))
cat(sprintf("dispersion sd(log eps)  nominal %.3f, realised %.4f (mean of eps %.4f)\n",
            SD_LOG, realised_sd, realised_mean))
cat(sprintf("background branches %d\n", length(bg_edges)))
cat(sprintf("variable sites %d of %d\n", length(seg.sites(seqs)), NSITE))
print(tab, row.names = FALSE)
