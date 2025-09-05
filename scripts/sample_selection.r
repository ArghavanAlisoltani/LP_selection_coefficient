############################################################
# Sample selection for ARGweaver algorithms
############################################################
# This script chooses a single representative child from each mother.
# For every child, similarities to all *outsiders* (individuals not
# sharing the same mother) are computed. We then minimise the following
# metrics to pick the most distinct child per mother:
#   * `max_sim`  – highest similarity to any outsider (primary metric)
#   * `q95_sim`  – 95th percentile similarity to outsiders (tie‑breaker)
#   * `mean_sim` – average similarity to outsiders (secondary tie‑breaker)
# The child with the smallest `max_sim` (followed by `q95_sim` and
# `mean_sim`) is selected for each mother.

############################################################
# run Algorithm
############################################################

setwd("~")  # ensure relative paths resolve correctly

# Read data
GRM <- as.matrix(read.csv("Cappa_GRM.csv", check.names = FALSE, row.names = 1))
meta <- read.csv("metadata.matchesCappa_GRM.csv", stringsAsFactors = FALSE)
meta <- meta[, c(2, 4)]                             # keep sample and mum columns
mum_freq <- data.frame(table(meta$mum))             # counts per mum (diagnostic)
id2mom <- setNames(meta$mum, meta$name)             # map sample ID → mum
target_n <- 67                                     # desired number of picks

# get stats in GRM
# Compute similarity statistics between each child and outsiders.
# Returns a data.frame with one row per sample.
rel_stats_to_other_mums_G <- function(G, id2mom) {
  ids <- rownames(G); id2mom <- id2mom[ids]  # reorder id2mom to match matrix
  G <- as.matrix(G); diag(G) <- NA           # ignore self-similarity
  do.call(rbind, lapply(ids, function(i){
    # outsiders are individuals with different mothers
    others <- ids[id2mom != id2mom[i]]
    if (!length(others))
      return(data.frame(sample_id = i, mum = id2mom[i],
                        n_others = 0L, max_sim = NA,
                        q95_sim = NA, mean_sim = NA))
    s <- G[i, others]; n_ok <- sum(!is.na(s))
    data.frame(sample_id = i, mum = id2mom[i], n_others = as.integer(n_ok),
               max_sim = max(s, na.rm = TRUE),
               q95_sim = as.numeric(quantile(s, 0.95, na.rm = TRUE)),
               mean_sim = mean(s, na.rm = TRUE))
  }))
}

#per-mum pick: minimize closest outsider (max_sim)
# Select the best representative per mother based on similarity stats.
pick_one_per_mum_maxsim_safe <- function(stats_df, min_n_others = 5) {
  out <- lapply(split(stats_df, stats_df$mum), function(df) {
    df <- df[df$n_others >= min_n_others, , drop = FALSE]  # require sufficient outsiders
    if (!nrow(df)) return(NULL)
    ord <- order(df$max_sim, -df$n_others, df$sample_id,
                 decreasing = FALSE, na.last = TRUE)
    df[ord[1], c("sample_id","mum","max_sim","q95_sim","mean_sim","n_others")]
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

# usage/run
stats_rel <- rel_stats_to_other_mums_G(GRM, id2mom)             # compute metrics
picks <- pick_one_per_mum_maxsim_safe(stats_rel, min_n_others = 5)  # choose samples

#visualize selected samples on heatmap
library(ComplexHeatmap)
library(circlize)
library(mixOmics)
meta$selected <- ifelse(meta$name %in% picks$sample_id,
                        "Selected", "Not_selected")

## ==== HEATMAP → PNG ====
library(ComplexHeatmap)
library(circlize)
library(grid)

# build the heatmap object (same settings you used)
ht <- Heatmap(
  GRM, name = " ",                           # main heatmap of GRM
  column_km = 1, row_km = 1,
  col = colorRamp2(c(0, 0.2, 0.5, 0.8, 1),
                   c("black","green","yellow1","orange","red")),
  show_row_names = FALSE, show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 5.5),
  cluster_rows = TRUE, cluster_columns = TRUE
) + Heatmap(
  meta$selected, name = "Selection",         # annotation showing picks
  width = unit(8, "mm"),
  col = c("gray","red4")
)

# write PNG (increase width/height if your matrix is huge)
png("selected_heatmap.png", width = 2400, height = 2200,
    res = 300, type = "cairo")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

  

# PCA plot of selected samples
snp <- read.table("snp_1490x25099_CharlesSEP2025.txt")
snp[1:5,1:5]                          # peek at raw SNP data
snpm <- as.matrix(snp)
snpm[1:5,1:5]                         # confirm conversion
var_snp <- apply(snpm, 2, var, na.rm = TRUE)
X_filtered <- snpm[, var_snp > 0]     # remove constant SNPs
X_scaled <- scale(X_filtered, center = TRUE, scale = TRUE)
pca.res <- pca(X_scaled, center = F, ncomp = 3, scale = F)
df <- data.frame(IDs = row.names(snp), order = c(1:1490))
merged <- merge(df, meta, by.x = "IDs", by.y = "name", all.x = T)
merged <- merged[order(merged$order),]

png("pca_selected.png", width = 1800, height = 1400,
    res = 300, type = "cairo")
plotIndiv(
  pca.res, group = merged$selected,
  title = " ", col = c("gray","red4"),
  centroid = FALSE, label = FALSE, star = FALSE, ind.names = FALSE,
  ellipse = TRUE, legend = TRUE, legend.title = "Selection"
)
dev.off()
