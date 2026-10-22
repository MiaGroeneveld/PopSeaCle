# =============================================================================
# Hybridisation analyses (R. djiddensis): PCA, FST-based diagnostic loci, sNMF admixture, NewHybrids
# Part 1: Mozambique (R. djiddensis x R. australiae; individual MO3, with MZ16 as a suspected hybrid)
# Part 2: Arabian region (R. djiddensis x R. laevis; individual RS2)
# NewHybrids is run outside R (Windows executable) through the parallelnewhybrid package.
# =============================================================================

# =============================================================================
# Following Nevatte et al. 2023 methods
# Inputs:
#   hyb_noMO3.vcf    — noMO3 dataset (35 ind: RA n=28, RD n=7)
#   hyb_withMO3.vcf  — withMO3 dataset (36 ind: RA n=28, RD n=8)

library(vcfR)
library(pcadapt)
library(ggplot2)
library(ggrepel)

# setwd("path/to/working/directory")

# Colours
col_RA      <- "#F5A623"   # lighter orange
col_RD      <- "#2166AC"   # darker blue
col_hybrid  <- "#7B2D8B"   # purple for MO3 and MZ16

# Individuals to highlight as suspected/detected hybrids
hybrid_inds <- c("MO3", "MZ16")

# =============================================================================
# Function to run PCA and build plot dataframe
# =============================================================================
run_pca <- function(vcf_file, K = 2) {
  bed <- read.pcadapt(vcf_file, type = "vcf")
  pca <- pcadapt(bed, K = K)
  vcf <- read.vcfR(vcf_file, verbose = FALSE)

  ind_names  <- colnames(vcf@gt)[-1]
  pop_assign <- ifelse(grepl("^MZ", ind_names), "RA", "RD")

  # Override population for hybrid individuals
  pop_assign[ind_names %in% hybrid_inds] <- "Hybrid"

  var_exp <- round(pca$singular.values^2 /
                     sum(pca$singular.values^2) * 100, 2)

  df <- data.frame(
    ind  = ind_names,
    pop  = factor(pop_assign, levels = c("RA", "RD", "Hybrid")),
    PC1  = pca$scores[, 1],
    PC2  = pca$scores[, 2]
  )

  list(df = df, var_exp = var_exp)
}

# =============================================================================
# Function to build final plot
# =============================================================================
make_pca_plot <- function(df, var_exp, labelled = FALSE) {
  p <- ggplot(df, aes(x = PC1, y = PC2,
                      fill = pop, colour = pop)) +
    geom_point(size = 4, shape = 21,
               colour = "black",      # border
               stroke = 0.6,
               aes(fill = pop)) +
    scale_fill_manual(
      values = c("RA" = col_RA, "RD" = col_RD, "Hybrid" = col_hybrid),
      labels = c("RA" = expression(italic("R. australiae")),
                 "RD" = expression(italic("R. djiddensis")),
                 "Hybrid" = "Putative hybrid")
    ) +
    scale_colour_manual(
      values = c("RA" = col_RA, "RD" = col_RD, "Hybrid" = col_hybrid),
      guide  = "none"
    ) +
    labs(
      x    = paste0("PC1 (", var_exp[1], "%)"),
      y    = paste0("PC2 (", var_exp[2], "%)"),
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "right",
      legend.text.align = 0,
      legend.key        = element_blank(),
      axis.line         = element_line(colour = "black"),
      panel.border      = element_blank()
    )

  if (labelled) {
    p <- p + geom_text_repel(
      aes(label = ind),
      size         = 2.5,
      colour       = "black",
      max.overlaps = 30,
      box.padding  = 0.3
    )
  }
  p
}

# =============================================================================
# withMO3 PCA
# =============================================================================
res_MO3 <- run_pca("hyb_withMO3.vcf")

p_MO3_unlab <- make_pca_plot(res_MO3$df, res_MO3$var_exp, labelled = FALSE)
p_MO3_lab   <- make_pca_plot(res_MO3$df, res_MO3$var_exp, labelled = TRUE)

ggsave("hyb_PCA_withMO3_final.pdf",     p_MO3_unlab, width = 7, height = 5)
ggsave("hyb_PCA_withMO3_labelled.pdf",  p_MO3_lab,   width = 7, height = 5)
cat("withMO3 PCA saved\n")

# =============================================================================
# noMO3 PCA
# =============================================================================
res_noMO3 <- run_pca("hyb_noMO3.vcf")

p_noMO3_unlab <- make_pca_plot(res_noMO3$df, res_noMO3$var_exp, labelled = FALSE)
p_noMO3_lab   <- make_pca_plot(res_noMO3$df, res_noMO3$var_exp, labelled = TRUE)

ggsave("hyb_PCA_noMO3_final.pdf",    p_noMO3_unlab, width = 7, height = 5)
ggsave("hyb_PCA_noMO3_labelled.pdf", p_noMO3_lab,   width = 7, height = 5)
cat("noMO3 PCA saved\n")

cat("\nDone. Final figures:\n")
cat("  hyb_PCA_withMO3_final.pdf   — publication figure (no labels)\n")
cat("  hyb_PCA_withMO3_labelled.pdf — for reference\n")
cat("  hyb_PCA_noMO3_final.pdf     — species separation (no labels)\n")
cat("  hyb_PCA_noMO3_labelled.pdf  — for reference\n")


# =============================================================================
# SECTION 2: FST-based diagnostic loci filtering
# Following Nevatte et al. — use noMO3 dataset to identify diagnostic loci,
# then apply to withMO3 dataset
# =============================================================================

# Read noMO3 VCF
vcf_noMO3 <- read.vcfR("hyb_noMO3.vcf")
cat("noMO3 dataset:", nrow(vcf_noMO3@gt), "loci,",
    ncol(vcf_noMO3@gt) - 1, "individuals\n")

# Extract genotypes
geno_noMO3 <- extract.gt(vcf_noMO3, element = "GT")
fix_noMO3  <- getFIX(vcf_noMO3)

****
# Correct hierfstat input format
# hierfstat needs two-digit allele coding: 11=homozygous ref, 12=het, 22=homozygous alt
# G_noMO3 is 0/1/2 — convert to 11/12/22

hf_geno <- G_noMO3
hf_geno[G_noMO3 == 0] <- 11
hf_geno[G_noMO3 == 1] <- 12
hf_geno[G_noMO3 == 2] <- 22
hf_geno[is.na(G_noMO3)] <- NA

# Population vector (1=RA, 2=RD)
pop_num <- ifelse(pop_noMO3 == "RA", 1, 2)

# Build hierfstat dataframe — individuals as rows, loci as columns
hf_input <- as.data.frame(cbind(pop_num, t(hf_geno)))

# Calculate per-locus FST
bs <- basic.stats(hf_input)
fst_vals <- bs$perloc$Fst

cat("FST distribution:\n")
print(summary(fst_vals))
cat("Loci with FST = 1:", sum(fst_vals == 1, na.rm = TRUE), "\n")
cat("Loci with FST > 0.9:", sum(fst_vals > 0.9, na.rm = TRUE), "\n")

****

# Plot FST distribution
pdf("hyb_FST_distribution.pdf", width = 6, height = 4)
hist(fst_vals, breaks = 50, main = "Per-locus FST (RA vs RD)",
     xlab = "FST", col = "steelblue", border = "white")
abline(v = 1, col = "red", lty = 2)
dev.off()

# Select diagnostic loci: FST = 1 (following Nevatte et al.)
# If too few loci, consider FST > 0.9 as per Nevatte's note
diag_idx  <- which(fst_vals == 1)
cat("Diagnostic loci (FST=1):", length(diag_idx), "\n")

# If FST=1 gives fewer than 200 loci, expand to FST > 0.9
if (length(diag_idx) < 200) {
  cat("Fewer than 200 loci at FST=1 — expanding to FST>0.9\n")
  diag_idx <- which(fst_vals > 0.9)
  cat("Diagnostic loci (FST>0.9):", length(diag_idx), "\n")
}

# Get locus IDs for diagnostic loci
diag_loci <- fix_noMO3[diag_idx, 1]  # CHROM column = STACKS locus IDs
cat("Diagnostic locus IDs extracted:", length(diag_loci), "\n")

# =============================================================================
# SECTION 3: Apply diagnostic loci to withMO3 dataset + remove missing data
# =============================================================================

# Read the withMO3 VCF (also used for the NewHybrids input below)
vcf_MO3 <- read.vcfR("hyb_withMO3.vcf", verbose = FALSE)

# Filter withMO3 VCF to diagnostic loci only
geno_MO3 <- extract.gt(vcf_MO3, element = "GT")
fix_MO3   <- getFIX(vcf_MO3)

# Match loci by CHROM field
diag_rows_MO3 <- which(fix_MO3[, 1] %in% diag_loci)
cat("Diagnostic loci found in withMO3 dataset:", length(diag_rows_MO3), "\n")

geno_diag <- geno_MO3[diag_rows_MO3, ]
fix_diag  <- fix_MO3[diag_rows_MO3, ]

# Remove loci with ANY missing data (following Nevatte et al.)
G_diag <- matrix(NA, nrow = nrow(geno_diag), ncol = ncol(geno_diag))
G_diag[geno_diag %in% c("0/0","0|0")] <- 0
G_diag[geno_diag %in% c("0/1","1/0","1|0","0|1")] <- 1
G_diag[geno_diag %in% c("1/1","1|1")] <- 2

miss_per_locus <- rowMeans(is.na(G_diag))
no_miss_idx    <- which(miss_per_locus == 0)
cat("Loci with zero missing data:", length(no_miss_idx), "\n")

geno_final <- geno_diag[no_miss_idx, ]
fix_final  <- fix_diag[no_miss_idx, ]
G_final    <- G_diag[no_miss_idx, ]

cat("\n=== FINAL DIAGNOSTIC DATASET ===\n")
cat("SNPs:", nrow(geno_final), "\n")
cat("Individuals:", ncol(geno_final), "\n")

# save
saveRDS(list(geno = geno_final,
             fix  = fix_final,
             G    = G_final,
             inds = colnames(geno_final),
             pops = ifelse(grepl("^MZ", colnames(geno_final)), "RA", "RD")),
        "hyb_diagnostic_dataset.rds")

cat("Saved: hyb_diagnostic_dataset.rds\n")
cat("Individuals:", colnames(geno_final), "\n")

# =============================================================================
# Hybridisation Analysis: sNMF admixture
# Following Nevatte et al. 2023 — LEA package
# Input: hyb_diagnostic_dataset.rds (from hyb_analysis_1_PCA_FST.R)
# =============================================================================

library(LEA)
library(ggplot2)
library(dplyr)
library(pophelper)  # for consistent sNMF bar plots with the main analyses

## setwd("path/to/working/directory")

# Population colours — consistent with main datasets
col_RA <- "#F5A623"   # update to match RA colour
col_RD <- "#2166AC"   # update to match RD colour

# Load final diagnostic dataset
dat <- readRDS("hyb_diagnostic_dataset.rds")
G   <- dat$G
ind_names <- dat$inds
pop_assign <- dat$pops

cat("Individuals:", ncol(G), "\n")
cat("Loci:", nrow(G), "\n")

# =============================================================================
# Convert genotype matrix to .geno format for LEA
# .geno format: one row per locus, one column per individual, values 0/1/2/9
# =============================================================================

geno_lea <- t(G)  # LEA needs individuals as rows, loci as columns
geno_lea[is.na(geno_lea)] <- 9  # missing data = 9 in LEA

# Write .geno file
write.geno(geno_lea, "hyb_diagnostic.geno")
cat("Geno file written: hyb_diagnostic.geno\n")

# =============================================================================
# sNMF: cross-entropy to select K
# K=1-5, 100 repetitions — following Nevatte et al.
# =============================================================================

cat("Running sNMF K=1-5 (100 reps each)...\n")
snmf_obj <- snmf("hyb_diagnostic.geno",
                 K          = 1:5,
                 entropy    = TRUE,
                 repetitions = 100,
                 project    = "new",
                 seed       = 123)

# Cross-entropy plot
pdf("hyb_sNMF_crossentropy.pdf", width = 6, height = 4)
plot(snmf_obj, col = "steelblue", cex = 1.2, pch = 19,
     main = "sNMF cross-entropy — hybridisation dataset")
dev.off()

# Print cross-entropy values
cat("\n--- Cross-entropy values ---\n")
for (k in 1:5) {
  ce <- cross.entropy(snmf_obj, K = k)
  cat(sprintf("K=%d: min=%.4f  mean=%.4f\n", k, min(ce), mean(ce)))
}

# =============================================================================
# Best K — based on cross-entropy knee (expect K=2 as per Nevatte)
# Update best_K if the cross-entropy plot suggests otherwise
# =============================================================================

best_K <- 2  # set from the cross-entropy plot

# Best run for best_K
ce_best_K <- cross.entropy(snmf_obj, K = best_K)
best_run   <- which.min(ce_best_K)
cat(sprintf("\nBest run for K=%d: run %d (CE=%.4f)\n",
            best_K, best_run, min(ce_best_K)))

# Extract Q matrix
Q_matrix <- Q(snmf_obj, K = best_K, run = best_run)
colnames(Q_matrix) <- paste0("K", 1:best_K)

# Build plot dataframe
q_df <- data.frame(
  ind  = ind_names,
  pop  = pop_assign,
  Q_matrix
)

# Order individuals: RA first, then RD (MO3 last within RD)
q_df$pop <- factor(q_df$pop, levels = c("RA", "Hybrid", "RD"))
q_df <- q_df %>% arrange(pop, ind)

# Identify which K column corresponds to RA vs RD
# (whichever has higher mean Q in RA individuals)
ra_idx <- which(q_df$pop == "RA")
k1_mean_ra <- mean(q_df$K1[ra_idx])
if (k1_mean_ra > 0.5) {
  q_df$Q_RA <- q_df$K1
  q_df$Q_RD <- q_df$K2
} else {
  q_df$Q_RA <- q_df$K2
  q_df$Q_RD <- q_df$K1
}

# plot

# Extract Q matrix for best run
best_run <- which.min(cross.entropy(snmf_obj, K = best_K))
Q_matrix <- Q(snmf_obj, K = best_K, run = best_run)
colnames(Q_matrix) <- paste0("K", 1:best_K)

# Individual order: RA, then MO3 (hybrid), then RD
ind_order <- c(which(pop_assign == "RA"),
               which(ind_names == "MO3"),
               which(pop_assign == "RD" & ind_names != "MO3"))

qmat <- as.data.frame(Q_matrix[ind_order, ])
rownames(qmat) <- ind_names[ind_order]

# Ensure RA column is first (higher Q in RA individuals)
ra_idx <- which(pop_assign[ind_order] == "RA")
if (mean(qmat[ra_idx, 1]) < 0.5) qmat <- qmat[, c(2, 1)]

# Colours — RA orange, RD blue (matching PCA)
hyb_cols <- c("#F5A623", "#2166AC")

plotQ(as.qlist(list(qmat)),
      clustercol      = hyb_cols,
      showindlab      = FALSE,
      showgrplab      = FALSE,
      barbordercolour = "black",
      barbordersize   = 0.1,
      divtype         = 1,
      divcol          = "white",
      divsize         = 0.5,
      panelratio      = c(8, 1),
      splab           = "K = 2",
      splabsize       = 4,
      showsp          = TRUE,
      outputfilename  = "hyb_withMO3_snmf_barplot_K2",
      imgtype         = "pdf",
      exportplot      = TRUE,
      exportpath      = getwd(),
      height          = 6,
      width           = 30,
      units           = "cm",
      dpi             = 300)

# Print Q values for key individuals
cat("\n--- Q values for MO3 ---\n")
print(qmat[rownames(qmat) == "MO3", ])
cat("\n--- Q values for MZ16 ---\n")
print(qmat[rownames(qmat) == "MZ16", ])
cat("\n--- RD individuals Q summary ---\n")
rd_rows <- which(pop_assign[ind_order] == "RD" & ind_names[ind_order] != "MO3")
print(qmat[rd_rows, ])

# NewHybrids prep
library(vcfR)
library(hierfstat)

# =============================================================================
# Step 1: Load data and recalculate FST
# =============================================================================

# Load withMO3 VCF (for NewHybrids input)
vcf_MO3   <- read.vcfR("hyb_withMO3.vcf", verbose = FALSE)
geno_MO3  <- extract.gt(vcf_MO3, element = "GT")
fix_MO3   <- getFIX(vcf_MO3)
ind_names <- colnames(geno_MO3)
pop_assign <- ifelse(grepl("^MZ", ind_names), "RA", "RD")

cat("withMO3 dataset:", ncol(geno_MO3), "individuals |",
    nrow(geno_MO3), "loci\n")

# Convert to numeric
G_MO3 <- matrix(NA, nrow = nrow(geno_MO3), ncol = ncol(geno_MO3))
G_MO3[geno_MO3 %in% c("0/0","0|0")] <- 0
G_MO3[geno_MO3 %in% c("0/1","1/0","1|0","0|1")] <- 1
G_MO3[geno_MO3 %in% c("1/1","1|1")] <- 2

# Load noMO3 VCF for FST calculation
vcf_noMO3  <- read.vcfR("hyb_noMO3.vcf", verbose = FALSE)
geno_noMO3 <- extract.gt(vcf_noMO3, element = "GT")
ind_noMO3  <- colnames(geno_noMO3)
pop_noMO3  <- ifelse(grepl("^MZ", ind_noMO3), "RA", "RD")

G_noMO3 <- matrix(NA, nrow = nrow(geno_noMO3), ncol = ncol(geno_noMO3))
G_noMO3[geno_noMO3 %in% c("0/0","0|0")] <- 0
G_noMO3[geno_noMO3 %in% c("0/1","1/0","1|0","0|1")] <- 1
G_noMO3[geno_noMO3 %in% c("1/1","1|1")] <- 2

# =============================================================================
# Step 2: Calculate FST and identify diagnostic loci (FST = 1)
# =============================================================================

cat("Calculating FST...\n")
hf_geno <- G_noMO3
hf_geno[G_noMO3 == 0] <- 11
hf_geno[G_noMO3 == 1] <- 12
hf_geno[G_noMO3 == 2] <- 22
hf_geno[is.na(G_noMO3)] <- NA

pop_num  <- ifelse(pop_noMO3 == "RA", 1, 2)
hf_input <- as.data.frame(cbind(pop_num, t(hf_geno)))

bs       <- basic.stats(hf_input)
fst_vals <- bs$perloc$Fst

cat("FST=1 loci:", sum(fst_vals == 1, na.rm = TRUE), "\n")

# Select FST=1 diagnostic loci
diag_idx  <- which(fst_vals == 1)
fix_noMO3 <- getFIX(vcf_noMO3)
diag_loci <- fix_noMO3[diag_idx, 1]

# Apply to withMO3 dataset
diag_rows_MO3 <- which(fix_MO3[, 1] %in% diag_loci)
cat("Diagnostic loci in withMO3:", length(diag_rows_MO3), "\n")

G_diag   <- G_MO3[diag_rows_MO3, ]
fix_diag <- fix_MO3[diag_rows_MO3, ]

# Remove loci with any missing data
miss_per_locus <- rowMeans(is.na(G_diag))
no_miss_idx    <- which(miss_per_locus == 0)
G_final        <- G_diag[no_miss_idx, ]
fix_final      <- fix_diag[no_miss_idx, ]

cat("Final diagnostic loci (zero missing):", nrow(G_final), "\n")

# =============================================================================
# Step 3: Convert to NewHybrids allele coding (11/12/22)
# =============================================================================

geno_nh <- matrix(NA, nrow = ncol(G_final), ncol = nrow(G_final))
for (j in 1:nrow(G_final)) {
  for (i in 1:ncol(G_final)) {
    g <- G_final[j, i]
    if (is.na(g)) {
      geno_nh[i, j] <- 0
    } else if (g == 0) {
      geno_nh[i, j] <- 11
    } else if (g == 1) {
      geno_nh[i, j] <- 12
    } else if (g == 2) {
      geno_nh[i, j] <- 22
    }
  }
}

cat("Genotype matrix ready:", nrow(geno_nh), "individuals x",
    ncol(geno_nh), "loci\n")

# =============================================================================
# Step 4: Write individuals file (same for all seeds)
# =============================================================================

ind_file <- "hyb_individuals.txt"
con <- file(ind_file, "w")
for (i in 1:length(ind_names)) {
  writeLines(paste(ind_names[i], 0), con)
}
close(con)
cat("Individuals file written:", ind_file, "\n")

# =============================================================================
# Step 5: Write 200-loci input files for multiple seeds
# =============================================================================

write_nh_input <- function(seed_val) {
  set.seed(seed_val)
  sub_idx  <- sample(1:nrow(G_final), 200)
  geno_sub <- geno_nh[, sub_idx]
  fix_sub  <- fix_final[sub_idx, ]
  locus_names <- paste(fix_sub[, 1], fix_sub[, 2], sep = "_")

  outfile <- paste0("hyb_newhybrids_200loci_seed", seed_val, ".txt")
  con <- file(outfile, "w")
  writeLines(paste("NumIndivs", length(ind_names)), con)
  writeLines("NumLoci 200", con)
  writeLines("Digits 1", con)
  writeLines("Format Lumped", con)
  writeLines("", con)
  writeLines(paste(locus_names, collapse = " "), con)
  writeLines("", con)
  for (i in 1:length(ind_names)) {
    geno_row <- c(ind_names[i], geno_sub[i, ])
    writeLines(paste(geno_row, collapse = " "), con)
  }
  close(con)

  # Report MO3 genotype pattern
  mo3_idx <- which(ind_names == "MO3")
  mo3_geno <- geno_sub[mo3_idx, ]
  cat("\nSeed", seed_val, "—", outfile, "\n")
  cat("  MO3: 11 count:", sum(mo3_geno == 11),
      "| 12 count:", sum(mo3_geno == 12),
      "| 22 count:", sum(mo3_geno == 22), "\n")

  return(invisible(outfile))
}

# Generate input files for three seeds
write_nh_input(321)
write_nh_input(456)
write_nh_input(789)

cat("\nAll input files written.\n")
cat("Copy each file 5 times into separate newhybrids_data folders\n")
cat("and run parallelnh_WIN for each seed.\n")

library(parallelnewhybrid)

path.data <- "path/to/newhybrids_data/"
path.NH   <- "path/to/newhybrids/"

# Full dataset (no subsetting). NOTE: this run produced floating-point underflow and was not used;
# NewHybrids was run on the 200-locus subsets above.
locus_names_full <- paste(fix_final[, 1], fix_final[, 2], sep = "_")

outfile <- "hyb_newhybrids_full.txt"
con <- file(outfile, "w")
writeLines(paste("NumIndivs", length(ind_names)), con)
writeLines(paste("NumLoci", nrow(G_final)), con)
writeLines("Digits 1", con)
writeLines("Format Lumped", con)
writeLines("", con)
writeLines(paste(locus_names_full, collapse = " "), con)
writeLines("", con)
for (i in 1:length(ind_names)) {
  geno_row <- c(ind_names[i], geno_nh[i, ])
  writeLines(paste(geno_row, collapse = " "), con)
}
close(con)

cat("Full dataset written:", outfile, "\n")
cat("NumLoci:", nrow(G_final), "\n")

library(parallelnewhybrid)

path.data <- "path/to/newhybrids_data/"
path.NH   <- "path/to/newhybrids/"

# Confirm files are in place
list.files(path.data)
list.files(path.NH)

parallelnh_WIN(
  folder.data = path.data,
  where.NH    = path.NH,
  burnin      = 10000,
  sweeps      = 100000
)

# =============================================================================
# RS2 outlier investigation
# Pipeline: PCA → FST=1 diagnostic loci → sNMF → NewHybrids prep
# Input files:
#   rs_noRS2.vcf  — without RS2 (for FST calculation: RD vs RL)
#   rs_withRS2.vcf — with RS2 (for PCA, sNMF, NewHybrids)
# De novo RS-only assembly
# =============================================================================

library(vcfR)
library(pcadapt)
library(hierfstat)
library(LEA)
library(ggplot2)
library(ggrepel)
library(pophelper)

# --- Population assignments (manual — all RS prefix) ---
# noRS2 dataset
rd_inds      <- c("RS1","RS3","RS5","RS6","RS7","RS9","RS10","RS11","RS12","RS13","RS14")
rl_inds      <- c("RS4","RS8","RS10")

# withRS2 dataset — RS2 treated as "Unknown" for PCA/sNMF
# Population colours — consistent with combined assembly plots
col_RD      <- "#2166AC"   # blue — R. djiddensis
col_RL      <- "#1B7837"   # green — R. laevis (cryptic)
col_RS2     <- "#7B2D8B"   # purple — RS2 putative intermediate

# =============================================================================
# SECTION 1: PCA
# =============================================================================

# --- Helper: run pcadapt PCA and return plot dataframe ---
run_pca <- function(vcf_file, pop_map, K = 2) {
    bed       <- read.pcadapt(vcf_file, type = "vcf")
    pca       <- pcadapt(bed, K = K)
    vcf_obj   <- read.vcfR(vcf_file, verbose = FALSE)
    ind_names <- colnames(vcf_obj@gt)[-1]

    # Assign populations from pop_map
    pop_assign <- pop_map$pop[match(ind_names, pop_map$ind)]
    pop_assign[is.na(pop_assign)] <- "Unknown"

    var_exp <- round(pca$singular.values^2 /
                         sum(pca$singular.values^2) * 100, 2)

    df <- data.frame(
        ind = ind_names,
        pop = factor(pop_assign, levels = c("RD", "RL", "RS2")),
        PC1 = pca$scores[, 1],
        PC2 = pca$scores[, 2]
    )
    list(df = df, var_exp = var_exp)
}

# --- Helper: build PCA plot ---
make_pca_plot <- function(df, var_exp, labelled = FALSE) {
    p <- ggplot(df, aes(x = PC1, y = PC2, fill = pop)) +
        geom_point(size = 4, shape = 21,
                   colour = "black", stroke = 0.6) +
        scale_fill_manual(
            values = c("RD"  = col_RD,
                       "RL"  = col_RL,
                       "RS2" = col_RS2),
            labels = c("RD"  = expression(italic("R. djiddensis")),
                       "RL"  = expression(italic("R. laevis")),
                       "RS2" = "RS2 (putative intermediate)"),
            na.value = "grey60"
        ) +
        labs(
            x    = paste0("PC1 (", var_exp[1], "%)"),
            y    = paste0("PC2 (", var_exp[2], "%)"),
            fill = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
            legend.position   = "right",
            legend.text.align = 0,
            legend.key        = element_blank()
        )

    if (labelled) {
        p <- p + geom_text_repel(
            aes(label = ind),
            size         = 2.5,
            colour       = "black",
            max.overlaps = 30,
            box.padding  = 0.3
        )
    }
    p
}

# --- Pop maps ---
# noRS2
popmap_noRS2 <- data.frame(
    ind = c("RS1","RS3","RS4","RS5","RS6","RS7",
            "RS8","RS9","RS10","RS11","RS12","RS13","RS14"),
    pop = c("RD","RD","RL","RD","RD","RD",
            "RL","RD","RL","RD","RD","RD","RD")
)

# withRS2 — RS2 labelled separately
popmap_withRS2 <- rbind(
    data.frame(ind = "RS2", pop = "RS2"),
    popmap_noRS2
)

# --- Run PCAs ---
res_noRS2   <- run_pca("rs_noRS2.vcf",   popmap_noRS2)
res_withRS2 <- run_pca("rs_withRS2.vcf", popmap_withRS2)

# noRS2 plots
p_noRS2_lab   <- make_pca_plot(res_noRS2$df,   res_noRS2$var_exp,   labelled = TRUE)
p_noRS2_unlab <- make_pca_plot(res_noRS2$df,   res_noRS2$var_exp,   labelled = FALSE)

# withRS2 plots
p_withRS2_lab   <- make_pca_plot(res_withRS2$df, res_withRS2$var_exp, labelled = TRUE)
p_withRS2_unlab <- make_pca_plot(res_withRS2$df, res_withRS2$var_exp, labelled = FALSE)

# Save
ggsave("RS2_PCA_noRS2_labelled.pdf",    p_noRS2_lab,     width = 7, height = 5, device = cairo_pdf)
ggsave("RS2_PCA_noRS2_final.pdf",       p_noRS2_unlab,   width = 7, height = 5, device = cairo_pdf)
ggsave("RS2_PCA_withRS2_labelled.pdf",  p_withRS2_lab,   width = 7, height = 5, device = cairo_pdf)
ggsave("RS2_PCA_withRS2_final.pdf",     p_withRS2_unlab, width = 7, height = 5, device = cairo_pdf)
cat("PCA figures saved.\n")

# =============================================================================
# SECTION 2: FST=1 diagnostic loci
# Calculated on noRS2 dataset (RD vs RL only)
# Applied to withRS2 dataset
# =============================================================================

# --- Load noRS2 VCF ---
vcf_noRS2  <- read.vcfR("rs_noRS2.vcf", verbose = FALSE)
geno_noRS2 <- extract.gt(vcf_noRS2, element = "GT")
fix_noRS2  <- getFIX(vcf_noRS2)
ind_noRS2  <- colnames(geno_noRS2)

cat("noRS2 dataset:", ncol(geno_noRS2), "individuals |",
    nrow(geno_noRS2), "loci\n")

# Convert to numeric
G_noRS2 <- matrix(NA, nrow = nrow(geno_noRS2), ncol = ncol(geno_noRS2))
G_noRS2[geno_noRS2 %in% c("0/0","0|0")] <- 0
G_noRS2[geno_noRS2 %in% c("0/1","1/0","0|1","1|0")] <- 1
G_noRS2[geno_noRS2 %in% c("1/1","1|1")] <- 2

# Population vector from popmap
pop_noRS2 <- popmap_noRS2$pop[match(ind_noRS2, popmap_noRS2$ind)]
pop_num   <- ifelse(pop_noRS2 == "RD", 1, 2)

# hierfstat format
hf_geno <- G_noRS2
hf_geno[G_noRS2 == 0] <- 11
hf_geno[G_noRS2 == 1] <- 12
hf_geno[G_noRS2 == 2] <- 22
hf_geno[is.na(G_noRS2)] <- NA

hf_input <- as.data.frame(cbind(pop_num, t(hf_geno)))
bs        <- basic.stats(hf_input)
fst_vals  <- bs$perloc$Fst

cat("FST summary:\n"); print(summary(fst_vals))
cat("Loci with FST = 1:", sum(fst_vals == 1, na.rm = TRUE), "\n")
cat("Loci with FST > 0.9:", sum(fst_vals > 0.9, na.rm = TRUE), "\n")

# Diagnostic loci selection
diag_idx <- which(fst_vals == 1)
if (length(diag_idx) < 50) {
    cat("Fewer than 50 loci at FST=1 — expanding to FST > 0.9\n")
    diag_idx <- which(fst_vals > 0.9)
}
diag_loci <- fix_noRS2[diag_idx, 1]
cat("Diagnostic loci selected:", length(diag_loci), "\n")

# --- Apply to withRS2 dataset ---
vcf_withRS2  <- read.vcfR("rs_withRS2.vcf", verbose = FALSE)
geno_withRS2 <- extract.gt(vcf_withRS2, element = "GT")
fix_withRS2  <- getFIX(vcf_withRS2)
ind_withRS2  <- colnames(geno_withRS2)

diag_rows <- which(fix_withRS2[, 1] %in% diag_loci)
cat("Diagnostic loci found in withRS2 dataset:", length(diag_rows), "\n")

G_withRS2 <- matrix(NA, nrow = nrow(geno_withRS2), ncol = ncol(geno_withRS2))
G_withRS2[geno_withRS2 %in% c("0/0","0|0")] <- 0
G_withRS2[geno_withRS2 %in% c("0/1","1/0","0|1","1|0")] <- 1
G_withRS2[geno_withRS2 %in% c("1/1","1|1")] <- 2

G_diag   <- G_withRS2[diag_rows, ]
fix_diag <- fix_withRS2[diag_rows, ]

# Remove loci with any missing data
no_miss_idx <- which(rowMeans(is.na(G_diag)) == 0)
G_final     <- G_diag[no_miss_idx, ]
fix_final   <- fix_diag[no_miss_idx, ]

cat("Final diagnostic loci (zero missing):", nrow(G_final), "\n")
cat("Individuals:", ncol(G_final), "\n")

# Save diagnostic dataset
saveRDS(list(G    = G_final,
             fix  = fix_final,
             inds = ind_withRS2,
             pops = popmap_withRS2$pop[match(ind_withRS2, popmap_withRS2$ind)]),
        "RS2_diagnostic_dataset.rds")
cat("Saved: RS2_diagnostic_dataset.rds\n")

# =============================================================================
# SECTION 3: sNMF (K=2: RD vs RL)
# =============================================================================

dat       <- readRDS("RS2_diagnostic_dataset.rds")
G         <- dat$G
ind_names <- dat$inds
pop_assign <- dat$pops

# LEA geno format
geno_lea <- t(G)
geno_lea[is.na(geno_lea)] <- 9
write.geno(geno_lea, "RS2_diagnostic.geno")

cat("Running sNMF K=1-4...\n")
snmf_obj <- snmf("RS2_diagnostic.geno",
                 K           = 1:4,
                 entropy     = TRUE,
                 repetitions = 100,
                 project     = "new",
                 seed        = 123)

# Cross-entropy plot
pdf("RS2_sNMF_crossentropy.pdf", width = 6, height = 4, onefile = FALSE)
plot(snmf_obj, col = "steelblue", cex = 1.2, pch = 19,
     main = "sNMF cross-entropy — RS2 dataset")
dev.off()

for (k in 1:4) {
    ce <- cross.entropy(snmf_obj, K = k)
    cat(sprintf("K=%d: min=%.4f  mean=%.4f\n", k, min(ce), mean(ce)))
}

# Best K — expect K=2 but confirm from cross-entropy plot
best_K   <- 2
best_run <- which.min(cross.entropy(snmf_obj, K = best_K))
Q_matrix <- Q(snmf_obj, K = best_K, run = best_run)
colnames(Q_matrix) <- paste0("K", 1:best_K)

# Individual order: RD first, RS2, then RL
ind_order <- c(which(pop_assign == "RD"),
               which(ind_names  == "RS2"),
               which(pop_assign == "RL"))

qmat <- as.data.frame(Q_matrix[ind_order, ])
rownames(qmat) <- ind_names[ind_order]

# Ensure RD = col 1 (blue)
rd_idx <- which(pop_assign[ind_order] == "RD")
if (mean(qmat[rd_idx, 1]) < 0.5) qmat <- qmat[, c(2, 1)]

# Colours: RD blue, RL green
snmf_cols <- c(col_RD, col_RL)

plotQ(as.qlist(list(qmat)),
      clustercol      = snmf_cols,
      showindlab      = TRUE,
      indlabsize      = 4,
      showgrplab      = FALSE,
      barbordercolour = "black",
      barbordersize   = 0.1,
      divtype         = 1,
      divcol          = "white",
      divsize         = 0.5,
      panelratio      = c(8, 1),
      splab           = "K = 2",
      splabsize       = 4,
      showsp          = TRUE,
      outputfilename  = "RS2_snmf_barplot_K2",
      imgtype         = "pdf",
      exportplot      = TRUE,
      exportpath      = getwd(),
      height          = 6,
      width           = 20,
      units           = "cm",
      dpi             = 300)

# Print Q values for RS2 specifically
cat("\n--- Q values for RS2 ---\n")
print(qmat["RS2", ])
cat("\n--- RD individuals Q summary ---\n")
print(qmat[pop_assign[ind_order] == "RD", ])
cat("\n--- RL individuals Q summary ---\n")
print(qmat[pop_assign[ind_order] == "RL", ])

# =============================================================================
# SECTION 4: NewHybrids input preparation
# Maximum 200 loci — subset across multiple seeds
# =============================================================================

# NewHybrids allele coding
geno_nh <- matrix(NA, nrow = ncol(G_final), ncol = nrow(G_final))
for (j in seq_len(nrow(G_final))) {
    for (i in seq_len(ncol(G_final))) {
        g <- G_final[j, i]
        geno_nh[i, j] <- if (is.na(g)) 0 else if (g == 0) 11 else if (g == 1) 12 else 22
    }
}

# Individuals file — no population code column (same format as for MO3)
con <- file("RS2_individuals.txt", "w")
for (nm in ind_withRS2) writeLines(paste(nm, 0), con)
close(con)
cat("Individuals file written.\n")

# Write 200-loci subsets for multiple seeds
write_nh_input <- function(seed_val) {
    set.seed(seed_val)

    # If fewer than 200 diagnostic loci, use all
    n_loci  <- min(200, nrow(G_final))
    sub_idx <- sample(seq_len(nrow(G_final)), n_loci)

    geno_sub    <- geno_nh[, sub_idx]
    fix_sub     <- fix_final[sub_idx, ]
    locus_names <- paste(fix_sub[, 1], fix_sub[, 2], sep = "_")

    outfile <- paste0("RS2_newhybrids_", n_loci, "loci_seed", seed_val, ".txt")
    con <- file(outfile, "w")
    writeLines(paste("NumIndivs", length(ind_withRS2)), con)
    writeLines(paste("NumLoci", n_loci), con)
    writeLines("Digits 1", con)
    writeLines("Format Lumped", con)
    writeLines("", con)
    writeLines(paste(locus_names, collapse = " "), con)
    writeLines("", con)
    for (i in seq_along(ind_withRS2)) {
        writeLines(paste(c(ind_withRS2[i], geno_sub[i, ]), collapse = " "), con)
    }
    close(con)

    # RS2 genotype pattern check
    rs2_idx  <- which(ind_withRS2 == "RS2")
    rs2_geno <- geno_sub[rs2_idx, ]
    cat(sprintf("Seed %d — %s\n", seed_val, outfile))
    cat(sprintf("  RS2: 11=%d | 12=%d | 22=%d\n",
                sum(rs2_geno == 11),
                sum(rs2_geno == 12),
                sum(rs2_geno == 22)))
    invisible(outfile)
}

write_nh_input(321)
write_nh_input(456)
write_nh_input(789)
write_nh_input(101)

cat("\nAll NewHybrids input files written.\n")
cat("Run parallelnh_WIN for each seed independently.\n")
cat("Expect: RD individuals → pure RD (posterior ~1)\n")
cat("        RL individuals → pure RL (posterior ~1)\n")
cat("        RS2 → examine posterior across F1/backcross/pure classes\n")

library(parallelnewhybrid)

path.data <- "path/to/RS_RS_only/"
path.NH   <- "path/to/newhybrids/"

# Confirm files are in place
list.files(path.data)
list.files(path.NH)

parallelnh_WIN(
  folder.data = path.data,
  where.NH    = path.NH,
  burnin      = 10000,
  sweeps      = 100000
)
