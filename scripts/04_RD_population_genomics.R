# ============================================================
# RD set1 — Population Genomics Pipeline
# Species: Rhynchobatus djiddensis (RD set1: 47 individuals, 27 366 SNPs)
# Populations: SA (SB + ZK + DU), MO, AR (Arabian region:
#   Arabian Gulf + Strait of Hormuz + Arabian Sea; n=1–2 per
#   sub-location; merged due to insufficient per-site n)
# Dataset: full SNP dataset (no outlier splits)
# Sections:
#   1. Setup and data loading
#   2. Pairwise FST + heatmap
#   3. DAPC (including the 5-PC DAPC chosen by the a-score)
#   4. sNMF
#   5. Diversity statistics (HO, HE, FIS, AR, PA)
#   6. IBD (five populations)
#   7. Unmerged FST and divMigrate (excluding the Arabian region)
#   8. Sampling map
# ============================================================

# ── Libraries ────────────────────────────────────────────────
library(vcfR)
library(adegenet)
library(hierfstat)
library(poppr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(LEA)
library(pophelper)
library(RColorBrewer)
library(nlme)
library(vegan)
library(gridExtra)
library(grid)

# ============================================================
# SECTION 1: Setup and data loading
# ============================================================

# ── Popmap ───────────────────────────────────────────────────
# Expected columns: ID, pop (original fine-scale pop codes)
popmap <- read.table("RDpopmap.txt", header = FALSE,
                     col.names = c("ID", "pop"))

# Merged population assignments
# SA = SB + ZK + DU (South African coastal locations)
# AR = RS sub-locations (Arabian Gulf, Strait of Hormuz,
#      Arabian Sea; n=1-2 per sub-location; merged due to
#      insufficient per-site sample sizes)
# MO = Mozambique (unchanged)
popmap_merged <- popmap
popmap_merged$pop[popmap_merged$pop %in% c("SB", "ZK", "DU")] <- "SA"
popmap_merged$pop[popmap_merged$pop == "RS"] <- "AR"
pop_levels <- c("SA", "MO", "AR")

# ── Plot aesthetics ───────────────────────────────────────────
pop_colours <- c("SA" = "#7F8A52",
                 "MO" = "#F0E442",
                 "AR" = "#7A1F2B")

pop_labels  <- c("SA" = "South Africa",
                 "MO" = "Mozambique",
                 "AR" = "Arabian region")

# ── Load VCF and build genotype matrix ───────────────────────
obj.vcfR <- read.vcfR("RDpop_set1.snps.vcf", verbose = FALSE)
geno_raw  <- extract.gt(obj.vcfR)
position  <- getPOS(obj.vcfR)
chrom     <- getCHROM(obj.vcfR)

G <- matrix(NA, nrow = nrow(geno_raw), ncol = ncol(geno_raw))
G[geno_raw %in% c("0/0", "0|0")]             <- 0
G[geno_raw %in% c("0/1","1/0","1|0","0|1")]  <- 1
G[geno_raw %in% c("1/1", "1|1")]             <- 2

# Mean imputation (for sNMF / DAPC; FST uses non-imputed via genind)
G_imp <- apply(G, 1, function(x) {
  x[is.na(x)] <- as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
  x
})
if (nrow(G_imp) != ncol(geno_raw)) G_imp <- t(G_imp)

locus_names <- paste(chrom, position, sep = "_")
colnames(G_imp) <- locus_names
rownames(G_imp) <- colnames(geno_raw)
cat("G_imp dimensions:", dim(G_imp), "\n")

# ── Population factor (in VCF sample order) ──────────────────
pop_factor <- popmap_merged$pop[match(rownames(G_imp),
                                      popmap_merged$ID)]
pop_full   <- factor(pop_factor, levels = pop_levels)
cat("Population counts:\n"); print(table(pop_full))

cat("Section 1 complete\n")

# ============================================================
# SECTION 2: Pairwise FST + Heatmap
# ============================================================

# ── Compute pairwise FST ─────────────────────────────────────
hf_full    <- cbind(pop = as.numeric(pop_full), as.data.frame(G_imp))
global_fst <- wc(hf_full)
cat("Global FST:", round(global_fst$FST, 4), "\n")

pairwise_fst <- pairwise.WCfst(hf_full)
rownames(pairwise_fst) <- pop_levels
colnames(pairwise_fst) <- pop_levels
print(round(pairwise_fst, 4))

# ── Bootstrap CIs ────────────────────────────────────────────
cat("Bootstrapping pairwise FST (1000 replicates)...\n")
boot_fst <- boot.ppfst(hf_full, nboot = 1000, quant = c(0.025, 0.975))

# Significance: CI excludes zero → p < 0.05 (BH-corrected)
n_pops    <- length(pop_levels)
pvals_mat <- matrix(NA, n_pops, n_pops,
                    dimnames = list(pop_levels, pop_levels))
for (i in 1:(n_pops - 1)) {
  for (j in (i + 1):n_pops) {
    pvals_mat[i, j] <- ifelse(!is.na(boot_fst$ll[i, j]) &
                                boot_fst$ll[i, j] > 0, 0.001, 1)
  }
}
pvals_vec <- pvals_mat[upper.tri(pvals_mat)]
pvals_bh  <- p.adjust(pvals_vec, method = "BH")
pvals_mat[upper.tri(pvals_mat)] <- pvals_bh
pvals_mat[lower.tri(pvals_mat)] <- t(pvals_mat)[lower.tri(pvals_mat)]

write.csv(as.data.frame(pairwise_fst),
          "RD_set1_pairwiseFST.csv")
write.csv(as.data.frame(pvals_mat),
          "RD_set1_pairwiseFST_pvals_BH.csv")
save(pairwise_fst, pvals_mat, global_fst,
     file = "RD_set1_FST_results.RData")

# ============================================================
# SECTION 2 (updated): Pairwise FST heatmap
# Upper triangle: FST values
# Lower triangle: significance symbols (BH-corrected p-values)
# Diagonal: white / blank
# ============================================================

# ── Significance symbol function ─────────────────────────────
sig_symbol <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("***")
  if (p < 0.01)    return("**")
  if (p < 0.05)    return("*")
  return("ns")
}

# ── Build long-format data frame ─────────────────────────────
n_pops   <- length(pop_levels)
fst_mat  <- as.matrix(pairwise_fst)
fst_mat[fst_mat < 0] <- 0  # floor negative FST at 0 for display

df_heat <- data.frame()

for (i in 1:n_pops) {
  for (j in 1:n_pops) {
    if (i == j) {
      # Diagonal — blank
      df_heat <- rbind(df_heat, data.frame(
        Pop1      = pop_levels[i],
        Pop2      = pop_levels[j],
        fill_val  = NA,
        label     = "",
        is_diag   = TRUE,
        stringsAsFactors = FALSE))
    } else if (i < j) {
      # Upper triangle — FST value
      df_heat <- rbind(df_heat, data.frame(
        Pop1      = pop_levels[i],
        Pop2      = pop_levels[j],
        fill_val  = fst_mat[i, j],
        label     = as.character(round(fst_mat[i, j], 3)),
        is_diag   = FALSE,
        stringsAsFactors = FALSE))
    } else {
      # Lower triangle — significance symbol
      df_heat <- rbind(df_heat, data.frame(
        Pop1      = pop_levels[i],
        Pop2      = pop_levels[j],
        fill_val  = fst_mat[i, j],  # still fill by FST for colour continuity
        label     = sig_symbol(pvals_mat[i, j]),
        is_diag   = FALSE,
        stringsAsFactors = FALSE))
    }
  }
}

df_heat$Pop1 <- factor(df_heat$Pop1, levels = pop_levels)
df_heat$Pop2 <- factor(df_heat$Pop2, levels = rev(pop_levels))

fst_max <- max(fst_mat, na.rm = TRUE)

p_fst <- ggplot(df_heat, aes(x = Pop1, y = Pop2, fill = fill_val)) +
  geom_tile(data = subset(df_heat, !is_diag), color = NA) +
  geom_tile(data = subset(df_heat, is_diag),
            fill = "white", color = NA) +
  geom_text(aes(label = label), size = 3, color = "black") +
  scale_fill_gradientn(
    colors   = c("#FFFFFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    na.value = "white",
    name     = expression(F[ST]),
    limits   = c(0, fst_max)) +
  scale_x_discrete(labels = pop_labels) +
  scale_y_discrete(labels = pop_labels) +
  theme_bw(base_size = 11) +
  theme(
    axis.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "right")

pdf("RD_set1_FST_heatmap.pdf", width = 5.5, height = 4.5)
print(p_fst); dev.off()
png("RD_set1_FST_heatmap.png", width = 5.5, height = 4.5,
    units = "in", res = 300)
print(p_fst); dev.off()

cat("FST heatmap saved\n")

# ============================================================
# SECTION 3: DAPC
# ============================================================

gl_RD <- vcfR2genlight(obj.vcfR)

# Fix locus names: colon → underscore, 0-based → 1-based position
loc_fixed <- gsub(":", "_", locNames(gl_RD))
loc_split <- strsplit(loc_fixed, "_")
loc_fixed <- sapply(loc_split, function(x)
  paste0(x[1], "_", as.integer(x[2]) + 1))
locNames(gl_RD) <- loc_fixed
pop(gl_RD)      <- pop_full

geno_mat <- tab(gl_RD, NA.method = "mean")

# ── PCA axes explaining ≥80% variance ────────────────────────
pca_full <- glPca(gl_RD, nf = min(nInd(gl_RD) - 1, nLoc(gl_RD)))
var_exp  <- cumsum(pca_full$eig / sum(pca_full$eig))
n_pc_80  <- which(var_exp >= 0.80)[1]
cat("PCs explaining 80% variance:", n_pc_80, "\n")

# ── Optimal K via find.clusters (10 replicates) ───────────────
best_k_vec <- integer(10)
for (rep in 1:10) {
  fc <- find.clusters(gl_RD, max.n.clust = 10, n.pca = n_pc_80,
                      choose.n.clust = FALSE, criterion = "diffNgroup",
                      n.iter = 1e6, n.start = 10)
  best_k_vec[rep] <- length(levels(fc$grp))
}
kmeans_k <- as.integer(names(sort(table(best_k_vec),
                                   decreasing = TRUE)[1]))
cat("Optimal K (k-means):", kmeans_k, "\n")

# ── Optimal PC retention via cross-validation ─────────────────
set.seed(42)
xval <- xvalDapc(geno_mat, pop(gl_RD),
                 n.pca        = 1:min(n_pc_80 * 2, 50),
                 training.set = 0.9,
                 result       = "groupMean",
                 center       = TRUE, scale = FALSE,
                 n.rep        = 100, xval.plot = FALSE)
n_pc_xval <- as.integer(xval$`Number of PCs Achieving Lowest MSE`)
cat("Cross-validated PCs:", n_pc_xval, "\n")

# ── A-score optimisation (for reporting) ─────────────────────
dapc_tmp    <- dapc(gl_RD, pop(gl_RD), n.pca = n_pc_80,
                    n.da = nlevels(pop(gl_RD)) - 1)
ascore      <- optim.a.score(dapc_tmp, n.sim = 100, plot = FALSE)
n_pc_ascore <- ascore$best
cat("A-score optimal PCs:", n_pc_ascore, "\n")

# ── DAPC with population prior ────────────────────────────────
dapc_prior <- dapc(gl_RD, pop(gl_RD),
                   n.pca = n_pc_xval,
                   n.da  = nlevels(pop(gl_RD)) - 1)

# ── DAPC without population prior ────────────────────────────
fc_final    <- find.clusters(gl_RD, n.clust = kmeans_k,
                              n.pca = n_pc_80,
                              choose.n.clust = FALSE,
                              n.iter = 1e6, n.start = 10)
dapc_noprior <- dapc(gl_RD, fc_final$grp,
                     n.pca = n_pc_xval,
                     n.da  = kmeans_k - 1)

# ── DAPC plot function ────────────────────────────────────────
make_dapc_plot <- function(dapc_obj, filename_pdf) {
  scores_df     <- as.data.frame(dapc_obj$ind.coord)
  scores_df$pop <- dapc_obj$grp
  eig_df        <- data.frame(axis = seq_along(dapc_obj$eig),
                               eig  = dapc_obj$eig)
  xlab <- paste0("DA1 (",
                  round(dapc_obj$eig[1] / sum(dapc_obj$eig) * 100, 1), "%)")
  ylab <- if (length(dapc_obj$eig) > 1)
    paste0("DA2 (",
           round(dapc_obj$eig[2] / sum(dapc_obj$eig) * 100, 1), "%)")
  else "DA2"

  has_ld2 <- "LD2" %in% colnames(scores_df)

  p_main <- ggplot(scores_df,
                   aes(x = LD1,
                       y = if (has_ld2) LD2 else 0,
                       fill = pop)) +
    geom_point(shape = 21, size = 3.5, stroke = 0.4, color = "grey30") +
    scale_fill_manual(values = pop_colours, labels = pop_labels,
                      name = "Population") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    xlab(xlab) + ylab(ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "right")

  p_eig <- ggplot(eig_df, aes(x = factor(axis), y = eig)) +
    geom_bar(stat = "identity", fill = "grey40", width = 0.6) +
    ggtitle("DA eigenvalues") +
    theme_bw(base_size = 8) +
    theme(panel.grid  = element_blank(),
          axis.title  = element_blank(),
          plot.title  = element_text(size = 8, hjust = 0.5))

  p_noleg   <- p_main + theme(legend.position = "none")
  leg       <- cowplot::get_legend(p_main)
  right_col <- plot_grid(leg, p_eig, ncol = 1, rel_heights = c(0.6, 0.4))
  p_final   <- plot_grid(p_noleg, right_col, ncol = 2,
                          rel_widths = c(0.78, 0.22))

  pdf(filename_pdf, width = 8, height = 6)
  print(p_final); dev.off()
  return(p_final)
}

p_dapc_prior   <- make_dapc_plot(dapc_prior,
                                  "RD_set1_dapc_prior.pdf")

pdf("RD_set1_dapc_noprior.pdf", width = 8, height = 6)
scatter(dapc_noprior, bg = "white", solid = 0.7, cex = 1.5,
        clab = 0, leg = TRUE)
dev.off()

# ── Summary table ─────────────────────────────────────────────
dapc_summary <- data.frame(
  N_SNPs          = nLoc(gl_RD),
  PCA_axes_80pct  = n_pc_80,
  Kmeans_optimal_K = kmeans_k,
  Crossval_PCs    = n_pc_xval,
  Ascore_PCs      = n_pc_ascore)
write.csv(dapc_summary, "RD_set1_DAPC_summary.csv", row.names = FALSE)

save(dapc_prior, dapc_noprior, fc_final, xval, ascore,
     file = "RD_set1_dapc_results.RData")

cat("Section 3 complete\n")


# ======================================================================
# 3b. DAPC with 5 PCs (a-score)
# ======================================================================

# ============================================================
# DAPC with PC = 5 as per a-score
# ============================================================

#obtain genlight object

pop_colours <- c("SA" = "#7F8A52",
                 "MO" = "#F0E442",
                 "AR" = "#7A1F2B")

pop_labels  <- c("SA" = "South Africa",
                 "MO" = "Mozambique",
                 "AR" = "Arabian region")

# ── DAPC with 5 PCs ──────────────────────────────────────────
dapc_5pc <- dapc(gl_RD, n.pca = 5, n.da = 2)

# ── Build plot dataframes ─────────────────────────────────────
scores_df     <- as.data.frame(dapc_5pc$ind.coord)
scores_df$pop <- dapc_5pc$grp

eig_df <- data.frame(axis = seq_along(dapc_5pc$eig),
                     eig  = dapc_5pc$eig)

xlab <- paste0("DA1 (",
               round(dapc_5pc$eig[1] / sum(dapc_5pc$eig) * 100, 1), "%)")
ylab <- paste0("DA2 (",
               round(dapc_5pc$eig[2] / sum(dapc_5pc$eig) * 100, 1), "%)")

# ── Main scatter plot ─────────────────────────────────────────
p_main <- ggplot(scores_df,
                 aes(x = LD1, y = LD2, fill = pop)) +
    geom_point(shape = 21, size = 3.5, stroke = 0.4, color = "grey30") +
    scale_fill_manual(values = pop_colours, labels = pop_labels,
                      name = "Population") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    xlab(xlab) + ylab(ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "right")

# ── Eigenvalue barplot ────────────────────────────────────────
p_eig <- ggplot(eig_df, aes(x = factor(axis), y = eig)) +
    geom_bar(stat = "identity", fill = "grey40", width = 0.6) +
    ggtitle("DA eigenvalues") +
    theme_bw(base_size = 8) +
    theme(panel.grid  = element_blank(),
          axis.title  = element_blank(),
          plot.title  = element_text(size = 8, hjust = 0.5))

# ── Combine ───────────────────────────────────────────────────
p_noleg   <- p_main + theme(legend.position = "none")
leg       <- cowplot::get_legend(p_main)
right_col <- plot_grid(leg, p_eig, ncol = 1, rel_heights = c(0.6, 0.4))
p_5pc     <- plot_grid(p_noleg, right_col, ncol = 2,
                        rel_widths = c(0.78, 0.22))

# ── Save ──────────────────────────────────────────────────────
cairo_pdf("RD_set1_dapc_5pc.pdf", width = 8, height = 6)
print(p_5pc)
dev.off()

png("RD_set1_dapc_5pc.png", width = 8, height = 6,
    units = "in", res = 300)
print(p_5pc)
dev.off()

cat("DAPC 5PC plot saved\n")
# ============================================================
# SECTION 4: sNMF
# ============================================================

# Order individuals by population for barplot
ind_meta <- data.frame(ID  = rownames(G_imp),
                        pop = pop_factor,
                        stringsAsFactors = FALSE)
ind_meta$pop <- factor(ind_meta$pop, levels = pop_levels)
ind_meta     <- ind_meta[order(ind_meta$pop), ]
grp_labels   <- data.frame(pop = as.character(ind_meta$pop),
                             stringsAsFactors = FALSE)

# Subset and reorder G_imp to population order
G_ordered <- G_imp[ind_meta$ID, ]
G_t       <- t(G_ordered)
G_t[is.na(G_t)] <- 9  # sNMF missing data code

geno_file <- "RD_set1_full.geno"
write.table(G_t, file = geno_file,
            row.names = FALSE, col.names = FALSE,
            sep = "", quote = FALSE)

# ── Run sNMF ─────────────────────────────────────────────────
proj_file <- "RD_set1_full.snmfProject"
if (file.exists(proj_file)) remove.snmfProject(proj_file)

snmf_proj <- snmf(geno_file, K = 1:10, rep = 10,
                  alpha = 100, entropy = TRUE, ploidy = 2,
                  project = "new")

# ── Cross-entropy plot ────────────────────────────────────────
ce_vals <- sapply(1:10, function(k)
  mean(cross.entropy(snmf_proj, K = k)))
best_k  <- which.min(ce_vals)
cat("Optimal K (cross-entropy):", best_k, "\n")

pdf("RD_set1_snmf_crossentropy.pdf", width = 6, height = 5)
plot(1:10, ce_vals, type = "b", pch = 19,
     xlab = "K", ylab = "Cross-entropy", bty = "l", las = 1)
abline(v = best_k, lty = 2, col = "red")
dev.off()

# ── NOTE: inspect cross-entropy plot and set best_k below ────
# Default set to best_k from cross-entropy; override if needed
# e.g. best_k_plot <- 2
best_k_plot <- best_k

best_run <- which.min(cross.entropy(snmf_proj, K = best_k_plot))

qmat <- as.data.frame(Q(snmf_proj, K = best_k_plot, run = best_run))
colnames(qmat) <- paste0("Cluster", 1:best_k_plot)
rownames(qmat) <- ind_meta$ID
qmat           <- qmat[ind_meta$ID, ]  # ensure population order

snmf_cols <- c("#7a1f2b", "#7F8A52")  # WIO cluster, Arabian region cluster

# Labelled barplot
plotQ(as.qlist(list(qmat)),
      clustercol      = snmf_cols,
      grplab          = grp_labels, selgrp = "pop",
      ordergrp        = FALSE,   # population order already set
      showgrplab      = TRUE,
      grplabsize      = 4, grplabangle = 0, grplabjust = 0.5,
      grplabcol       = "grey30", showindlab = FALSE,
      sortind         = "all",
      barbordercolour = "black", barbordersize = 0.1,
      divtype = 1, divcol = "white", divsize = 0.5,
      panelratio      = c(8, 1), grplabheight = 0.8,
      splab           = paste0("K = ", best_k_plot),
      splabsize = 4, showsp = TRUE,
      outputfilename  = paste0("RD_set1_snmf_barplot_K",
                                best_k_plot, "_labelled"),
      imgtype = "pdf", exportplot = TRUE, exportpath = getwd(),
      height = 6, width = 30, units = "cm", dpi = 300)

# Unlabelled barplot (for publication figures)
plotQ(as.qlist(list(qmat)),
      clustercol      = snmf_cols,
      grplab          = grp_labels, selgrp = "pop",
      ordergrp        = FALSE,
      showgrplab      = FALSE,
      showindlab      = FALSE,
      sortind         = "all",
      barbordercolour = "black", barbordersize = 0.1,
      divtype = 1, divcol = "white", divsize = 0.5,
      panelratio      = c(8, 1), grplabheight = 0.8,
      splab           = paste0("K = ", best_k_plot),
      splabsize = 4, showsp = TRUE,
      outputfilename  = paste0("RD_set1_snmf_barplot_K",
                                best_k_plot, "_unlabelled"),
      imgtype = "pdf", exportplot = TRUE, exportpath = getwd(),
      height = 6, width = 30, units = "cm", dpi = 300)

save(snmf_proj, qmat, best_k_plot, best_run,
     file = "RD_set1_snmf_results.RData")

cat("Section 4 complete\n")

# ============================================================
# SECTION 5: Diversity statistics (HO, HE, FIS, AR, PA)
# ============================================================

gind_full <- vcfR2genind(obj.vcfR)

# Fix locus names (colon → underscore, 0-based → 1-based)
loc_fixed2 <- gsub(":", "_", locNames(gind_full))
loc_split2 <- strsplit(loc_fixed2, "_")
loc_fixed2 <- sapply(loc_split2, function(x)
  paste0(x[1], "_", as.integer(x[2]) + 1))
locNames(gind_full) <- loc_fixed2
pop(gind_full)      <- factor(pop_factor, levels = pop_levels)

hf_div  <- genind2hierfstat(gind_full)
bs      <- basic.stats(hf_div)

ho_pop  <- colMeans(bs$Ho,  na.rm = TRUE)
he_pop  <- colMeans(bs$Hs,  na.rm = TRUE)
fis_pop <- colMeans(bs$Fis, na.rm = TRUE)

# ── Bootstrap FIS CIs (locus resampling) ─────────────────────
# Note: boot.ppfst computes FST CIs only; custom bootstrap
# required for per-population FIS CIs
cat("Bootstrapping FIS CIs (500 replicates)...\n")
set.seed(42)
hf_mat      <- hf_div[, -1]
pop_vec     <- hf_div[, 1]
fis_boot    <- matrix(NA, 500, length(pop_levels),
                       dimnames = list(NULL, pop_levels))

for (b in 1:500) {
  boot_idx  <- sample(ncol(hf_mat), ncol(hf_mat), replace = TRUE)
  hf_boot   <- cbind(pop = as.numeric(pop_vec),
                      as.data.frame(hf_mat[, boot_idx]))
  bs_boot   <- basic.stats(hf_boot)
  fis_boot[b, ] <- colMeans(bs_boot$Fis, na.rm = TRUE)
}
fis_ll <- apply(fis_boot, 2, quantile, 0.025, na.rm = TRUE)
fis_ul <- apply(fis_boot, 2, quantile, 0.975, na.rm = TRUE)

# ── Allelic richness (rarefied to min n) ─────────────────────
min_n  <- min(table(pop(gind_full)))
ar     <- allelic.richness(hf_div, min.n = min_n)
ar_pop <- colMeans(ar$Ar, na.rm = TRUE)

# ── Private alleles ───────────────────────────────────────────
pa_mat <- private_alleles(gind_full, count.alleles = FALSE, drop = TRUE)
pa_pop <- rowSums(pa_mat)

# ── Diversity table ───────────────────────────────────────────
div_table <- data.frame(
  Population = pop_levels,
  n          = as.integer(table(pop(gind_full))[pop_levels]),
  HO         = round(ho_pop[pop_levels],  4),
  HE         = round(he_pop[pop_levels],  4),
  FIS        = round(fis_pop[pop_levels], 4),
  FIS_lower  = round(fis_ll[pop_levels],  4),
  FIS_upper  = round(fis_ul[pop_levels],  4),
  AR         = round(ar_pop[pop_levels],  4),
  PA         = as.integer(pa_pop[pop_levels]),
  row.names  = NULL)

print(div_table)
write.csv(div_table, "RD_set1_diversity.csv", row.names = FALSE)
save(div_table, bs, ar, file = "RD_set1_diversity_results.RData")

cat("Section 5 complete\n")


# ======================================================================
# 6. Isolation by distance: five populations (SB, ZK, DU, MO, AR)
# ======================================================================


library(vcfR)
library(adegenet)
library(hierfstat)
library(nlme)
library(vegan)
library(ggplot2)

# ============================================================
# SECTION 1: Load VCF and assign unmerged populations
# ============================================================

vcf_full <- read.vcfR("RDpop_set1.snps.vcf", verbose = FALSE)

popmap <- read.table("RDpopmap.txt", header = FALSE,
                     col.names = c("SampleID", "Population"))

# Keep original unmerged population codes
pop_levels <- c("SB", "ZK", "DU", "MO", "AR")

# Recode RS → AR
popmap$Population[popmap$Population == "RS"] <- "AR"

pop_factor <- factor(popmap$Population[match(
    colnames(extract.gt(vcf_full)), popmap$SampleID)],
    levels = pop_levels)

cat("Population counts:\n")
print(table(pop_factor))

# ============================================================
# SECTION 2: Pairwise FST — all five populations
# ============================================================

genind_obj <- vcfR2genind(vcf_full)
pop(genind_obj) <- pop_factor

X     <- tab(genind_obj, freq = TRUE, NA.method = "mean")
G_imp <- X

hf_full <- cbind(pop = as.numeric(pop_factor),
                 as.data.frame(G_imp))

cat("Calculating pairwise FST (5 populations)...\n")
pairwise_fst_5pop <- pairwise.WCfst(hf_full)
rownames(pairwise_fst_5pop) <- pop_levels
colnames(pairwise_fst_5pop) <- pop_levels

cat("Global FST:\n")
print(round(wc(hf_full)$FST, 4))

cat("\nPairwise FST:\n")
print(round(pairwise_fst_5pop, 4))

write.csv(as.data.frame(pairwise_fst_5pop),
          "RD_set1_5pop_pairwiseFST.csv")

# ============================================================
# SECTION 3: Marine geographic distances — 5 populations
# ============================================================

load("RD_set1_spatial_data.RData")  # provides marine_dist_RD

n_pops       <- length(pop_levels)
geo_dist_pop <- matrix(NA, n_pops, n_pops,
                        dimnames = list(pop_levels, pop_levels))

for (i in 1:n_pops) {
    for (j in 1:n_pops) {
        idx_i <- which(pop_factor == pop_levels[i])
        idx_j <- which(pop_factor == pop_levels[j])
        geo_dist_pop[i, j] <- mean(marine_dist_RD[idx_i, idx_j])
    }
}
diag(geo_dist_pop) <- 0

cat("Population geographic distance matrix (km):\n")
print(round(geo_dist_pop, 0))

# ============================================================
# SECTION 4: Linearised FST
# ============================================================

linearize_fst <- function(m) {
    m[m < 0] <- 0
    m / (1 - m)
}

fst_lin <- linearize_fst(as.matrix(pairwise_fst_5pop))

# ============================================================
# SECTION 5: Extract lower triangle pairs
# ============================================================

pairs_idx <- which(lower.tri(geo_dist_pop), arr.ind = TRUE)
geo_vec   <- geo_dist_pop[pairs_idx]
gen_vec   <- fst_lin[pairs_idx]

cat("Number of pairwise comparisons:", length(geo_vec), "\n")

# ============================================================
# SECTION 6: Permutation test function
# ============================================================

perm_r2 <- function(y, x, n_perm = 10000) {
    obs   <- cor(y, x, use = "complete.obs")^2
    perms <- replicate(n_perm,
                       cor(sample(y), x, use = "complete.obs")^2)
    p_val <- mean(perms >= obs)
    list(r2 = obs, p = p_val)
}

# ============================================================
# SECTION 7: GLS IBD
# ============================================================

ibd_gls  <- gls(gen_vec ~ geo_vec, method = "ML")
ibd_perm <- perm_r2(gen_vec, geo_vec)

cat("IBD GLS — R²:", round(ibd_perm$r2, 3),
    " p =", round(ibd_perm$p, 4), "\n")
cat("IBD GLS — AIC:", round(AIC(ibd_gls), 2), "\n")

ibd_results <- data.frame(
    IBD_R2    = round(ibd_perm$r2, 3),
    IBD_p     = round(ibd_perm$p,  4),
    Slope     = round(coef(ibd_gls)[2], 6),
    Intercept = round(coef(ibd_gls)[1], 6),
    AIC       = round(AIC(ibd_gls), 2))

print(ibd_results)
write.csv(ibd_results, "RD_set1_5pop_IBD_results.csv",
          row.names = FALSE)

# ============================================================
# SECTION 8: IBD plot
# ============================================================

ibd_df <- data.frame(Geo = geo_vec, FST = gen_vec)

p_ibd <- ggplot(ibd_df, aes(x = Geo, y = FST)) +
    geom_point(size = 3, colour = "#2171B5") +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#2171B5", fill = "#C6DBEF",
                linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("R² = ", round(ibd_perm$r2, 3),
                            "; p = ", round(ibd_perm$p, 4)),
             hjust = 1.05, vjust = -0.5,
             size = 3, colour = "grey30") +
    xlab("Geographic distance (km)") +
    ylab(expression("Linearised genetic distance (F"[ST]*
                        "/(1\u2212F"[ST]*"))")) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

cairo_pdf("RD_set1_5pop_IBD.pdf", width = 7, height = 6)
print(p_ibd)
dev.off()

png("RD_set1_5pop_IBD.png", width = 7, height = 6,
    units = "in", res = 300)
print(p_ibd)
dev.off()

# ============================================================
# SECTION 9: Save
# ============================================================

save(pairwise_fst_5pop, geo_dist_pop, fst_lin,
     ibd_results, ibd_df, p_ibd,
     file = "RD_set1_5pop_IBD_results.RData")


# ======================================================================
# 7. Pairwise FST and divMigrate for the unmerged South African and Mozambican sites (excluding the Arabian region)
# ======================================================================

# Additional
# RD set1 — Unmerged pairwise FST + divMigrate (excl. RS)
# Populations: SB, ZK, DU, MO
# Conventions: hierfstat pairwise.WCfst, 1000 boot, BH correction
# ============================================================

library(vcfR)
library(adegenet)
library(hierfstat)
library(diveRsity)
library(poppr)
library(Cairo)

# ── Load popmap (excl. RS) ───────────────────────────────────
popmap <- read.table("RDpopmap.txt",
                     header = FALSE,
                     col.names = c("SampleID", "Population"),
                     stringsAsFactors = FALSE)

popmap_noRS <- popmap[!popmap$Population == "RS", ]
pop_levels  <- c("SB", "ZK", "DU", "MO")
popmap_noRS$Population <- factor(popmap_noRS$Population,
                                  levels = pop_levels)

cat("Samples per population (excl. RS):\n")
print(table(popmap_noRS$Population))

# ── Load and filter VCF ──────────────────────────────────────
vcf_full <- read.vcfR("RDpop_set1.snps.vcf")

# Retain only non-RS individuals
samples_noRS <- popmap_noRS$SampleID
vcf_noRS     <- vcf_full[, c("FORMAT", samples_noRS)]

cat("VCF dimensions (excl. RS):",
    nrow(vcf_noRS@gt), "loci x",
    ncol(vcf_noRS@gt) - 1, "individuals\n")

# ── Convert to genind ────────────────────────────────────────
genind_noRS <- vcfR2genind(vcf_noRS)
pop(genind_noRS) <- popmap_noRS$Population

# Impute missing data with mean allele frequency per locus
# (consistent with previous RD pipeline)
X <- tab(genind_noRS, freq = TRUE, NA.method = "mean")
G_imp <- X

# ── Build hierfstat dataframe ────────────────────────────────
pop_full   <- pop(genind_noRS)
hf_full    <- cbind(pop = as.numeric(pop_full),
                    as.data.frame(G_imp))

# ── Global FST ───────────────────────────────────────────────
global_fst <- wc(hf_full)
cat("Global FST (SB+ZK+DU+MO):", round(global_fst$FST, 4), "\n")

# ── Pairwise FST ─────────────────────────────────────────────
cat("Calculating pairwise FST...\n")
pairwise_fst <- pairwise.WCfst(hf_full)
rownames(pairwise_fst) <- pop_levels
colnames(pairwise_fst) <- pop_levels
print(round(pairwise_fst, 4))

# ── Bootstrap CIs + BH-corrected significance ────────────────
cat("Bootstrapping pairwise FST (1000 replicates)...\n")
boot_fst <- boot.ppfst(hf_full,
                        nboot = 1000,
                        quant = c(0.025, 0.975))

n_pops    <- length(pop_levels)
pvals_mat <- matrix(NA, n_pops, n_pops,
                    dimnames = list(pop_levels, pop_levels))

for (i in 1:(n_pops - 1)) {
  for (j in (i + 1):n_pops) {
    pvals_mat[i, j] <- ifelse(!is.na(boot_fst$ll[i, j]) &
                                boot_fst$ll[i, j] > 0, 0.001, 1)
  }
}

pvals_vec        <- pvals_mat[upper.tri(pvals_mat)]
pvals_bh         <- p.adjust(pvals_vec, method = "BH")
pvals_mat[upper.tri(pvals_mat)] <- pvals_bh
pvals_mat[lower.tri(pvals_mat)] <- t(pvals_mat)[lower.tri(pvals_mat)]

cat("BH-corrected significance matrix:\n")
print(round(pvals_mat, 4))

# ── Save FST results ─────────────────────────────────────────
write.csv(as.data.frame(pairwise_fst),
          "RD_set1_noRS_pairwiseFST.csv")
write.csv(as.data.frame(pvals_mat),
          "RD_set1_noRS_pairwiseFST_pvals_BH.csv")
save(pairwise_fst, pvals_mat, global_fst,
     file = "RD_set1_noRS_FST_results.RData")

# heatmap
# Fix significance function
sig_symbol <- function(p) {
    if (is.na(p))    return("")
    if (p <= 0.001)  return("***")
    if (p <= 0.01)   return("**")
    if (p <= 0.05)   return("*")
    return("ns")
}

# Fix label formatting — no scientific notation
df_heat <- data.frame()

for (i in 1:n_pops) {
    for (j in 1:n_pops) {
        if (i == j) {
            df_heat <- rbind(df_heat, data.frame(
                Pop1 = pop_levels[i], Pop2 = pop_levels[j],
                fill_val = NA, label = "", is_diag = TRUE,
                stringsAsFactors = FALSE))
        } else if (i < j) {
            df_heat <- rbind(df_heat, data.frame(
                Pop1 = pop_levels[i], Pop2 = pop_levels[j],
                fill_val = fst_mat[i, j],
                label = formatC(pairwise_fst[i, j], 
                               format = "f", digits = 4),
                is_diag = FALSE, stringsAsFactors = FALSE))
        } else {
            df_heat <- rbind(df_heat, data.frame(
                Pop1 = pop_levels[i], Pop2 = pop_levels[j],
                fill_val = fst_mat[i, j],
                label = sig_symbol(pvals_mat[i, j]),
                is_diag = FALSE, stringsAsFactors = FALSE))
        }
    }
}

df_heat$Pop1 <- factor(df_heat$Pop1, levels = pop_levels)
df_heat$Pop2 <- factor(df_heat$Pop2, levels = rev(pop_levels))
fst_max <- max(fst_mat, na.rm = TRUE)

p_fst_noRS <- ggplot(df_heat, aes(x = Pop1, y = Pop2, fill = fill_val)) +
    geom_tile(data = subset(df_heat, !is_diag), color = NA) +
    geom_tile(data = subset(df_heat, is_diag),
              fill = "white", color = NA) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_gradientn(
        colors   = c("#FFFFFF","#C6DBEF","#6BAED6","#2171B5","#08306B"),
        na.value = "white",
        name     = expression(F[ST]),
        limits   = c(0, fst_max)) +
    scale_x_discrete(labels = pop_labels) +
    scale_y_discrete(labels = pop_labels) +
    theme_bw(base_size = 11) +
    theme(
        axis.title      = element_blank(),
        axis.text.x     = element_text(angle = 45, hjust = 1),
        panel.grid      = element_blank(),
        panel.border    = element_blank(),
        axis.ticks      = element_blank(),
        legend.position = "right")

cairo_pdf("RD_set1_noRS_pairwiseFST_heatmap_labelled.pdf",
          width = 5, height = 4)
print(p_fst_noRS)
dev.off()

png("RD_set1_noRS_pairwiseFST_heatmap_labelled.png",
    width = 5, height = 4,
    units = "in", res = 300)
print(p_fst_noRS)
dev.off()

cat("Heatmap saved\n")

# divMigrate (diveRsity v1.9.90, which was removed from CRAN):
# install.packages("https://cran.r-project.org/src/contrib/Archive/diveRsity/diveRsity_1.9.90.tar.gz",
#                  repos = NULL, type = "source")
#
# Workflow (see docs/04_RD_population_genomics.md)
#   1. On the cluster: remove the Arabian region (RS) individuals and draw a random subset of 5 000 SNPs
#      (vcftools); download RD_set1_noRS.recode.vcf (full dataset, 27 366 SNPs, for relative migration) and
#      RD_set1_noRS_5000.recode.vcf (subset, for significance testing and the network plot)
#   2. Convert each VCF to GENEPOP on Galaxy Europe (usegalaxy.eu): genepop.txt (full) and genepop_5000.txt
#   3. Restructure both files in R (below), because Galaxy puts every individual in a single Pop section
#   4. Relative migration on the full dataset in R, without bootstraps (bootstrapping fails silently)
#   5. Significance (500 bootstraps, alpha = 0.05) on the 5 000-SNP subset in the divMigrate online application
#   6. Network plot in R from the 5 000-SNP subset

library(diveRsity)
library(Cairo)

# ── Restructure the GENEPOP files ────────────────────────────
popmap <- read.table("RDpopmap.txt", header = FALSE,
                     col.names = c("SampleID", "Population"))
popmap_noRS <- popmap[popmap$Population != "RS", ]

genepop_files <- c("genepop_5000.txt" = "RD_set1_noRS_5000_divMigrate.gen",
                   "genepop.txt"      = "RD_set1_noRS_full_divMigrate.gen")

for (genepop_in in names(genepop_files)) {
    genepop_out <- genepop_files[[genepop_in]]

    lines       <- readLines(genepop_in)
    pop_line    <- grep("^Pop$", lines)
    locus_lines <- lines[2:(pop_line[1] - 1)]
    ind_lines   <- lines[(pop_line[1] + 1):length(lines)]
    ind_lines   <- ind_lines[!grepl("^Pop$", ind_lines)]
    ind_lines   <- ind_lines[ind_lines != ""]

    cat(genepop_in, "- loci:", length(locus_lines),
        "| individuals:", length(ind_lines), "\n")

    # Write restructured genepop with correct population sections
    pops <- c("SB", "ZK", "DU", "MO")
    con  <- file(genepop_out, "w")
    writeLines(sub("_divMigrate\\.gen$", "", genepop_out), con)
    writeLines(locus_lines, con)

    for (p in pops) {
        writeLines("Pop", con)
        pop_inds <- popmap_noRS$SampleID[popmap_noRS$Population == p]
        for (ind in pop_inds) {
            ind_data <- ind_lines[grep(paste0("^", ind, " ,"), ind_lines)]
            if (length(ind_data) > 0) writeLines(ind_data, con)
        }
    }
    close(con)

    # Verify
    new_lines     <- readLines(genepop_out)
    new_pop_lines <- grep("^Pop$", new_lines)
    cat("Pop sections:", length(new_pop_lines), "\n")
    for (i in new_pop_lines) {
        cat("Pop", which(new_pop_lines == i),
            "- first individual:", substr(new_lines[i+1], 1, 10), "\n")
    }
}
# Expected: Pop 1 = SB, Pop 2 = ZK, Pop 3 = DU (starts with SB9, which is coded DU in the popmap), Pop 4 = MO

# ── Relative migration on the full dataset (27,366 SNPs), no bootstraps ──────
dm_d <- divMigrate(
    infile           = "RD_set1_noRS_full_divMigrate.gen",
    boots            = 0,
    stat             = "d",
    filter_threshold = 0,
    para             = TRUE)

dm_gst <- divMigrate(
    infile           = "RD_set1_noRS_full_divMigrate.gen",
    boots            = 0,
    stat             = "gst",
    filter_threshold = 0,
    para             = TRUE)

dm_nm <- divMigrate(
    infile           = "RD_set1_noRS_full_divMigrate.gen",
    boots            = 0,
    stat             = "Nm",
    filter_threshold = 0,
    para             = TRUE)

# Assign population names to the matrices
pops_ordered <- c("SB", "ZK", "DU", "MO")
rownames(dm_d$dRelMig)   <- colnames(dm_d$dRelMig)   <- pops_ordered
rownames(dm_gst$gRelMig) <- colnames(dm_gst$gRelMig) <- pops_ordered
rownames(dm_nm$nmRelMig) <- colnames(dm_nm$nmRelMig) <- pops_ordered

print(round(dm_d$dRelMig, 3))
print(round(dm_gst$gRelMig, 3))
print(round(dm_nm$nmRelMig, 3))

save(dm_d, dm_gst, dm_nm,
     file = "RD_set1_noRS_divMigrate_results.RData")

# ── Network plot from the 5,000-SNP subset (Jost's D, following Manuzzi et al. 2019) ──
# Node labels show the first individual of each population (SB1t, ZK1, SB9, MO1);
# they were corrected to population names manually after export.
cairo_pdf("RD_set1_noRS_divMigrate_network_labelled.pdf",
          width = 7, height = 7)
divMigrate(
    infile           = "RD_set1_noRS_5000_divMigrate.gen",
    outfile          = NULL,
    boots            = 0,
    stat             = "d",
    filter_threshold = 0,
    plot_network     = TRUE,
    plot_col         = "darkblue",
    para             = TRUE
)
dev.off()

png("RD_set1_noRS_divMigrate_network_labelled.png",
    width = 7, height = 7, units = "in", res = 300)
divMigrate(
    infile           = "RD_set1_noRS_5000_divMigrate.gen",
    outfile          = NULL,
    boots            = 0,
    stat             = "d",
    filter_threshold = 0,
    plot_network     = TRUE,
    plot_col         = "darkblue",
    para             = TRUE
)
dev.off()

# ======================================================================
# 8. Sampling map
# ======================================================================

### figure 1 map

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# ── Data ─────────────────────────────────────────────────────
world <- ne_countries(scale = "medium", returnclass = "sf")

# ── Main map ──────────────────────────────────────────────────
main_map <- ggplot() +
    geom_sf(data = world,
            fill = "gray95", color = "gray70", linewidth = 0.2) +
    coord_sf(xlim = c(25, 65), ylim = c(-35, 30), expand = FALSE) +
    annotation_scale(
        location   = "br",
        width_hint = 0.15,
        line_width = 0.3,
        text_cex   = 0.7
    ) +
    scale_x_continuous(
        breaks = seq(25, 65, by = 10),
        labels = function(x) paste0(abs(x), "\u00b0E")
    ) +
    scale_y_continuous(
        breaks = seq(-30, 30, by = 10),
        labels = function(y) ifelse(y < 0,
                                    paste0(abs(y), "\u00b0S"),
                                    paste0(y, "\u00b0N"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        axis.title       = element_blank(),
        axis.text        = element_text(size = 9, color = "black"),
        text             = element_text(color = "black"),
        panel.border     = element_rect(color = "gray40", fill = NA,
                                        linewidth = 0.5)
    )

# ── Inset ─────────────────────────────────────────────────────
main_bbox <- st_as_sfc(st_bbox(c(
    xmin = 25, xmax = 65,
    ymin = -35, ymax = 30
), crs = 4326))

inset_map <- ggplot() +
    geom_sf(data = world,
            fill = "gray95", color = "gray60", linewidth = 0.15) +
    geom_sf(data = main_bbox,
            fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(xlim = c(-20, 120), ylim = c(-50, 40), expand = FALSE) +
    theme_void() +
    theme(
        panel.border     = element_rect(color = "gray40", fill = NA,
                                        linewidth = 0.3),
        panel.background = element_rect(fill = "white", color = NA)
    )

# ── Export main map ───────────────────────────────────────────
ggsave("RD_main_basemap.png",
       plot   = main_map,
       width  = 10,
       height = 8,
       dpi    = 600)

ggsave("RD_main_basemap.pdf",
       plot   = main_map,
       width  = 10,
       height = 8)

# ── Export inset ──────────────────────────────────────────────
ggsave("RD_inset_map.png",
       plot   = inset_map,
       width  = 4,
       height = 3,
       dpi    = 600)

ggsave("RD_inset_map.pdf",
       plot   = inset_map,
       width  = 4,
       height = 3)

message("Done — main and inset maps saved separately")
