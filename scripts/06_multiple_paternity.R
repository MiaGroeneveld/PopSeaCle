# =============================================================================
# COLONY Input File Preparation — RD Family Multiple Paternity Analysis
# Following Nevatte et al. 2023 parameters
# Family: SC12 (mother) + SC1-SC11 (11 pups)
# Dataset: 1378 SNPs of dataset D (family extracted); the STRUCTURE file is written by Stacks populations
# =============================================================================
#
# WORKFLOW:
# 1. Run this R script to generate all input files
# 2. Open the COLONY GUI (Colony2.exe)
# 3. File -> New Project -> Empirical Data Analysis
# 4. Fill in each tab as described below
# 5. Save Data -> Run -> Start Running
#
# COLONY PARAMETERS (following Nevatte et al. 2023):
#   - Female Polygamy (sharks can have multiple mates)
#   - Male Polygamy
#   - Without Inbreeding
#   - Without Clone
#   - Dioecious
#   - Diploid
#   - Full-Likelihood (FL) method
#   - Medium run length
#   - 5 runs
#   - Random seed: 1234
#   - No sibship prior
#   - Sibship scaling: Yes
#   - Update allele frequency: No (estimated during analysis)
#   - Genotype error rate: 0.01
#   - Allelic dropout rate: 0
#   - Mother SC12 specified as known female with prob = 1
#   - All pups specified as known maternal sibship group
# =============================================================================

library(vcfR)

# The STRUCTURE file is the output of Stacks populations for the 12 family individuals
# (whitelist of the 1378 SNPs of dataset D), renamed MP_family_1378.structure
structure_file <- "MP_family_1378.structure"

# Verify: 2 header lines, then 2 lines for each of the 12 individuals
lines <- readLines(structure_file)[-(1:2)]
cat("Total lines:", length(lines), "\n")  # should be 24 (12 x 2)
cols  <- length(strsplit(lines[1], "\t")[[1]])
cat("Columns:", cols, "\n")  # number of loci + ID + pop
cat("Loci:", cols - 2, "\n")  # SNPs in the STRUCTURE file

# Check individual IDs
ind_ids <- sapply(seq(1, length(lines), by=2), function(i)
  strsplit(lines[i], "\t")[[1]][1])
cat("Individual IDs:\n")
print(ind_ids)


# =============================================================================
# Step 1: Read STACKS STRUCTURE file
# Format: 2 header rows, 2 rows per individual, col 2 = population (remove)
# Alleles coded as 1/2/3/4 (A/C/G/T) — correct multi-allele coding for related
# skip = 2 if structure downloaded directly from hpc
# =============================================================================

struct <- read.table("MP_family_1378.structure",
                     skip = 2,
                     header = FALSE,
                     stringsAsFactors = FALSE)

# Remove population column (col 2)
struct <- struct[, -2]

n_rows <- nrow(struct)
n_ind  <- n_rows / 2
n_loci <- ncol(struct) - 1

cat("Individuals:", n_ind, "\n")
cat("Loci:", n_loci, "\n")

# Individual names from odd rows
ind_names <- struct[seq(1, n_rows, by = 2), 1]
cat("Individuals found:", ind_names, "\n")

# =============================================================================
# Step 2: Build genotype matrix
# COLONY format: IndivID  a1_L1  a2_L1  a1_L2  a2_L2 ...
# =============================================================================

geno_mat <- matrix(NA, nrow = n_ind, ncol = 2 * n_loci)

for (i in 1:n_ind) {
  row1 <- as.integer(struct[2*i - 1, -1])  # allele 1 for all loci
  row2 <- as.integer(struct[2*i,     -1])  # allele 2 for all loci
  for (j in 1:n_loci) {
    geno_mat[i, 2*j - 1] <- row1[j]
    geno_mat[i, 2*j]     <- row2[j]
  }
}

rownames(geno_mat) <- ind_names

# Replace NA with 0 (COLONY missing data code)
geno_mat[is.na(geno_mat)] <- 0

cat("Genotype matrix:", nrow(geno_mat), "individuals x",
    ncol(geno_mat), "allele columns\n")

# =============================================================================
# Step 3: Separate mother and pups
# Mother = SC12 (confirmed from field collection)
# Pups   = SC1-SC11
# =============================================================================

mother_id  <- "SC12"
pup_ids    <- ind_names[ind_names != mother_id]

cat("Mother:", mother_id, "\n")
cat("Pups:", paste(pup_ids, collapse = ", "), "\n")
cat("Number of pups:", length(pup_ids), "\n")

mother_geno <- geno_mat[mother_id, , drop = FALSE]
pup_geno    <- geno_mat[pup_ids,   , drop = FALSE]

# =============================================================================
# Step 4: Write offspring genotype file
# COLONY format: IndivID  a1_L1  a2_L1  a1_L2  a2_L2 ...
# =============================================================================

write.table(pup_geno,
            "OffspringGenotype.txt",
            row.names = TRUE,
            col.names = FALSE,
            sep       = " ",
            quote     = FALSE)

cat("Offspring genotype file written: OffspringGenotype.txt\n")
cat("  Individuals:", nrow(pup_geno), "\n")
cat("  Allele columns:", ncol(pup_geno), "(=", n_loci, "loci x 2)\n")

# =============================================================================
# Step 5: Write mother (female candidate) genotype file
# COLONY format: IndivID  a1_L1  a2_L1  a1_L2  a2_L2 ...
# SC12 is the known mother — probability of being the mother = 1
# =============================================================================

write.table(mother_geno,
            "FemaleGenotype.txt",
            row.names = TRUE,
            col.names = FALSE,
            sep       = " ",
            quote     = FALSE)

cat("Female genotype file written: FemaleGenotype.txt\n")

# =============================================================================
# Step 6: Write marker types file
# COLONY format: 4 lines, all values space-separated
#   Line 1: Locus names (L1 L2 ... L1923)
#   Line 2: Marker types (0 = codominant for all)
#   Line 3: Allelic dropout rates (0 for all)
#   Line 4: Other error rates (0.01 for all, following Nevatte et al.)
#
# NOTE: Requires Windows decimal separator to be "." not ","
# Check: Control Panel -> Region -> Additional Settings -> Decimal symbol = "."
# =============================================================================

con <- file("MarkerTypeErrorRate.txt", "w")
writeLines(paste(paste0("L", 1:n_loci), collapse = " "), con)  # locus names
writeLines(paste(rep(0,    n_loci), collapse = " "), con)       # codominant
writeLines(paste(rep(0,    n_loci), collapse = " "), con)       # dropout = 0
writeLines(paste(rep(0.01, n_loci), collapse = " "), con)       # error = 0.01
close(con)

cat("Marker file written: MarkerTypeErrorRate.txt\n")
cat("  Loci:", n_loci, "\n")
cat("  Marker type: codominant (0)\n")
cat("  Dropout rate: 0\n")
cat("  Error rate: 0.01\n")

# =============================================================================
# Step 7: Write known maternal sibship file
# All 11 pups share SC12 as their mother
# COLONY format: MotherID  Offspring1  Offspring2 ... OffspringN
# This tells COLONY all pups are maternal siblings — critical constraint
# =============================================================================

con <- file("KnownMaternalSibs.txt", "w")
writeLines(paste(c(mother_id, pup_ids), collapse = " "), con)
close(con)

cat("Known maternal sibship file written: KnownMaternalSibs.txt\n")
cat("  Mother:", mother_id, "\n")
cat("  Pups:", paste(pup_ids, collapse = " "), "\n")

# ── Sensitivity analysis: exclude SC3 (putative MZ twin) ─────────────────────
# SC2 and SC3 show near-identical genotypes (TrioML = 1.003; bitwise dist = 0.004-0.014)
# consistent with MZ twinning. Including both may violate COLONY's independence
# assumption and drive spurious sibship assignment.
# This sensitivity analysis tests whether the two-father result persists
# when one putative twin is excluded.
pup_ids_noSC3    <- pup_ids[pup_ids != "SC3"]
pup_geno_noSC3   <- geno_mat[pup_ids_noSC3, , drop = FALSE]

write.table(pup_geno_noSC3, "OffspringGenotype_noSC3.txt",
            row.names = TRUE, col.names = FALSE,
            sep = " ", quote = FALSE)
cat("OffspringGenotype_noSC3.txt written —",
    nrow(pup_geno_noSC3), "pups (SC3 excluded)\n")

# Update known maternal sibship file for sensitivity analysis
con <- file("KnownMaternalSibs_noSC3.txt", "w")
writeLines(paste(c(mother_id, pup_ids_noSC3), collapse = " "), con)
close(con)
cat("KnownMaternalSibs_noSC3.txt written\n")


# =============================================================================
# COLONY GUI SETUP INSTRUCTIONS
# =============================================================================
#
# 1. Open COLONY -> File -> New Project -> Empirical Data Analysis -> OK
#
# 2. PARAMETERS TAB:
#    - Female Polygamy
#    - Male Polygamy
#    - Without Inbreeding
#    - Without Clone
#    - Dioecious
#    - Diploid
#    - Full-Likelihood (FL)
#    - Medium run length
#    - Number of Run: 5
#    - Random Number Seed: 1234
#    - Update Allele Frequency: No
#    - Sibship Scaling: Yes
#    - Sibship Prior: No Prior
#
# 3. MARKERS TAB:
#    - Number of Loci: 1923
#    - Click Load -> select MarkerTypeErrorRate.txt
#    - Allele Frequency: Unknown
#    - Click Check Data -> OK
#
# 4. OFFSPRING GENOTYPES TAB:
#    - Number of Offspring: 11
#    - Click Load Genotypes -> select OffspringGenotype.txt
#    - Click Check Data -> OK
#
# 5. MALE GENOTYPES TAB:
#    - Number of Males: 0 (no candidate fathers — COLONY will reconstruct)
#
# 6. FEMALE GENOTYPES TAB:
#    - Number of Females: 1 (SC12 is the known mother)
#    - Click Load Genotypes -> select FemaleGenotype.txt
#    - Prob of a Mum in Candidates: 1
#    - Click Check Data -> OK
#
# 7. KNOWN PATERNAL SIBS TAB:
#    - Number: 0
#
# 8. KNOWN MATERNAL SIBS TAB:
#    - Number: 1 (one maternal sibship group = all pups share SC12)
#    - Maternity Exc Threshold: 50
#      (allows up to 50 mismatches between mother-offspring dyads;
#       with 1923 loci at 0.01 error rate, expect ~19 mismatches by chance)
#    - Click Load -> select KnownMaternalSibs.txt
#    - Click Check Data -> OK
#
# 9. EXC. PATERNITY TAB: skip
#
# 10. EXC. MATERNITY TAB:
#     - Number: 0
#
# 11. EXC. PATERNAL SIBS TAB:
#     - Number: 0
#
# 12. EXC. MATERNAL SIBS TAB:
#     - Number: 0
#
# 13. SAVE DATA TAB:
#     - Click Save Data -> OK
#
# 14. Run -> Start Running -> threads: 4 -> OK
#
# =============================================================================
# KEY OUTPUT FILES (for the COLONY project folder)
# =============================================================================
#
# MP_RD.BestFSFamily     — Full-sibling family assignments (MAIN RESULT)
#                          Number of families = minimum number of sires
# MP_RD.FullSibDyad      — Pairwise full-sibling probabilities
# MP_RD.HalfSibDyad      — Pairwise half-sibling probabilities
# MP_RD.DadGenotype      — Reconstructed paternal genotypes
# MP_RD.Maternity        — Maternity assignments (confirm SC12)
# MP_RD.PairwisePaternity — Paternity probabilities
#
# =============================================================================

cat("\nAll input files generated successfully.\n")
cat("Copy all files to the COLONY project folder before opening COLONY.\n")
cat("Follow the GUI instructions above to set up and run the analysis.\n")

# pairwise dist and KING coefficients

library(vcfR)
library(adegenet)
library(poppr)
library(SNPRelate)
library(ggplot2)
library(reshape2)
library(cowplot)

# ── Individual order — pups only, natural order ───────────────
pup_order <- paste0("SC", 1:11)


# ============================================================
# DATASET D (family + SWIO reference panel, 47 individuals, 1378 SNPs)
# ============================================================

# ── Load combined VCF ─────────────────────────────────────────
vcf_comb <- read.vcfR("MP_combined.vcf", verbose = FALSE)
cat("Combined VCF — individuals:", ncol(extract.gt(vcf_comb)), "\n")
cat("Combined VCF — loci:", nrow(extract.gt(vcf_comb)), "\n")

# ── Bitwise distance — pups only ──────────────────────────────
gl_comb      <- vcfR2genlight(vcf_comb)
gl_pups_comb <- gl_comb[indNames(gl_comb) %in% pup_order, ]
bit_comb     <- as.matrix(bitwise.dist(gl_pups_comb, percent = TRUE,
                                        mat = TRUE, missing_match = FALSE))
bit_comb     <- bit_comb[pup_order, pup_order]

# ── KING kinship — all individuals, subset to pups ───────────
snpgdsVCF2GDS("MP_combined.vcf", "MP_combined_noPD1.gds",
               method = "biallelic.only", verbose = FALSE)
gds_comb     <- snpgdsOpen("MP_combined_noPD1.gds")
king_comb_obj <- snpgdsIBDKING(gds_comb, type = "KING-robust",
                                 num.thread = 2, verbose = FALSE)
snpgdsClose(gds_comb)

king_comb <- king_comb_obj$kinship
rownames(king_comb) <- colnames(king_comb) <- king_comb_obj$sample.id
king_comb <- king_comb[pup_order, pup_order]

# ============================================================
# HEATMAP FUNCTION
# ============================================================

make_lower_heatmap <- function(mat, legend_name,
                                low_col, high_col,
                                mid_col  = "white",
                                midpoint = NULL,
                                title    = "") {
  # Upper triangle + diagonal to NA, flip so lower triangle shows
  mat_t <- t(mat)
  mat_t[lower.tri(mat_t, diag = TRUE)] <- NA

  df <- melt(mat_t, varnames = c("Ind1", "Ind2"),
             value.name = "value")
  df <- df[!is.na(df$value), ]

  df$Ind1 <- factor(df$Ind1, levels = pup_order)
  df$Ind2 <- factor(df$Ind2, levels = rev(pup_order))

  if (is.null(midpoint)) {
    p <- ggplot(df, aes(x = Ind1, y = Ind2, fill = value)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient(low = low_col, high = high_col,
                          name = legend_name)
  } else {
    p <- ggplot(df, aes(x = Ind1, y = Ind2, fill = value)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(low = low_col, mid = mid_col,
                            high     = high_col,
                            midpoint = midpoint,
                            name     = legend_name)
  }

  p + ggtitle(title) +
    theme_bw(base_size = 11) +
    theme(
      axis.title      = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.ticks      = element_blank(),
      plot.title      = element_text(size = 11),
      legend.position = "right")
}

# ============================================================
# DATASET D HEATMAPS
# ============================================================

p_bit_comb <- make_lower_heatmap(
  bit_comb,
  legend_name = "Genetic\ndistance",
  low_col     = "#ECF0F1",
  high_col    = "#2980B9",
  title       = "Pairwise genetic distance - dataset D")

p_king_comb <- make_lower_heatmap(
  king_comb,
  legend_name = "KING\nkinship \u03c6",
  low_col     = "#4B0082",
  mid_col     = "white",
  high_col    = "#2980B9",
  midpoint    = 0,
  title       = "Pairwise KING kinship - dataset D")

combined_comb <- plot_grid(p_bit_comb, p_king_comb, ncol = 2)

cairo_pdf("MP_heatmaps_combined.pdf", width = 12, height = 6)
print(combined_comb); dev.off()
png("MP_heatmaps_combined.png", width = 12, height = 6,
    units = "in", res = 300)
print(combined_comb); dev.off()
cat("Combined dataset heatmaps saved\n")

# Empirical relatedness (LynchRD, TrioML, DyadML) and the relatedness simulation were run on the cluster:
#   scripts/hpc/MP_related_empirical.R and scripts/hpc/MP_related_simulation.R
#   (inputs made by 06a_MP_related_input.R). The results are read back in below.

# process in R

library(ggplot2)
library(reshape2)
library(cowplot)

## lynchrd input
emp_lynchrd <- readRDS("MP_emp_lynchrd.rds")
emprel      <- emp_lynchrd$relatedness
family_ids <- paste0("SC", 1:12)
fam_pairs  <- emprel[emprel$ind1.id %in% family_ids &
                      emprel$ind2.id %in% family_ids, ]


# Build pairwise matrix — all family including mother SC12
all_fam_ids <- c(paste0("SC", 1:12))
n_fam       <- length(all_fam_ids)

rel_mat_full <- matrix(NA, nrow = n_fam, ncol = n_fam,
                        dimnames = list(all_fam_ids, all_fam_ids))

for (i in 1:nrow(fam_pairs)) {
  id1 <- fam_pairs$ind1.id[i]
  id2 <- fam_pairs$ind2.id[i]
  val <- fam_pairs$lynchrd[i]
  if (id1 %in% all_fam_ids & id2 %in% all_fam_ids) {
    rel_mat_full[id1, id2] <- val
    rel_mat_full[id2, id1] <- val
  }
}

cat("Matrix dimensions:", dim(rel_mat_full), "\n")
cat("Range of values:", round(range(rel_mat_full, na.rm=TRUE), 3), "\n")
print(round(rel_mat_full, 3))

# ── Build lower triangle dataframe ───────────────────────────
rel_mat_plot <- rel_mat_full
rel_mat_plot[lower.tri(rel_mat_plot, diag = TRUE)] <- NA

df <- melt(rel_mat_plot, varnames = c("Ind1", "Ind2"),
           value.name = "value")
df <- df[!is.na(df$value), ]

df$Ind1 <- factor(df$Ind1, levels = all_fam_ids)
df$Ind2 <- factor(df$Ind2, levels = rev(all_fam_ids))

p_lynchrd <- ggplot(df, aes(x = Ind1, y = Ind2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3)),
            size = 2.8, color = "black") +
  scale_fill_gradient(low  = "white",
                      high = "#1B7837", name = "Lynch-Ritland\nrelatedness") +
  theme_bw(base_size = 11) +
  theme(
    axis.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.ticks      = element_blank(),
    plot.title      = element_text(size = 11),
    legend.position = "right")

cairo_pdf("MP_relatedness_lynchrd.pdf", width = 7, height = 6)
print(p_lynchrd); dev.off()

png("MP_relatedness_lynchrd.png", width = 7, height = 6,
    units = "in", res = 300)
print(p_lynchrd); dev.off()

cat("Saved\n")

## trioml input
# change input here
emp_ml <- readRDS("MP_emp_ML.rds")          # change file
emprel      <- emp_ml$relatedness

# Family pairs only
family_ids <- paste0("SC", 1:12)
fam_pairs  <- emprel[emprel$ind1.id %in% family_ids &
                      emprel$ind2.id %in% family_ids, ]

# Build pairwise matrix
all_fam_ids  <- c(paste0("SC", 1:12))
n_fam        <- length(all_fam_ids)
rel_mat_full <- matrix(NA, nrow = n_fam, ncol = n_fam,
                        dimnames = list(all_fam_ids, all_fam_ids))

for (i in 1:nrow(fam_pairs)) {
  id1 <- fam_pairs$ind1.id[i]
  id2 <- fam_pairs$ind2.id[i]
  val <- fam_pairs$trioml[i]                     # change estimator
  if (id1 %in% all_fam_ids & id2 %in% all_fam_ids) {
    rel_mat_full[id1, id2] <- val
    rel_mat_full[id2, id1] <- val
  }
}

# ── Build lower triangle dataframe ───────────────────────────
rel_mat_plot <- rel_mat_full
rel_mat_plot[lower.tri(rel_mat_plot, diag = TRUE)] <- NA

df <- melt(rel_mat_plot, varnames = c("Ind1", "Ind2"),
           value.name = "value")
df <- df[!is.na(df$value), ]

df$Ind1 <- factor(df$Ind1, levels = all_fam_ids)
df$Ind2 <- factor(df$Ind2, levels = rev(all_fam_ids))

p_trioml <- ggplot(df, aes(x = Ind1, y = Ind2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3)),
            size = 2.8, color = "black") +
  scale_fill_gradient(low  = "white",
                      high = "#4B0082", name = "Triadic likelihood\nrelatedness") +
  theme_bw(base_size = 11) +
  theme(
    axis.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.ticks      = element_blank(),
    plot.title      = element_text(size = 11),
    legend.position = "right")

cairo_pdf("MP_relatedness_trioml.pdf", width = 7, height = 6)
print(p_trioml); dev.off()

png("MP_relatedness_trioml.png", width = 7, height = 6,
    units = "in", res = 300)
print(p_trioml); dev.off()

cat("Saved\n")

# both on one panel

coancestry_all <- plot_grid(p_lynchrd, p_trioml, 
                             ncol = 2)

cairo_pdf("MP_coancestry_all.pdf", width = 12, height = 6)
print(coancestry_all); dev.off()

png("MP_coancestry_all.png", width = 12, height = 6,
    units = "in", res = 300)
print(coancestry_all); dev.off()

cat("Saved\n")

#PISibs

calc_pisibs <- function(vcf_file, label) {
  vcf  <- read.vcfR(vcf_file, verbose = FALSE)
  geno <- extract.gt(vcf, element = "GT")
  
  pisibs_loci <- apply(geno, 1, function(gt) {
    gt  <- gt[!is.na(gt)]
    a1  <- ifelse(grepl("^0", gt), 1, 2)
    a2  <- ifelse(grepl("0$", gt), 1, 2)
    alleles <- c(a1, a2)
    n   <- length(alleles)
    p   <- sum(alleles == 1) / n
    q   <- sum(alleles == 2) / n
    # Waits et al. 2001 equation 3
    sum_pi2 <- p^2 + q^2
    sum_pi4 <- p^4 + q^4
    0.25 + (0.5 * sum_pi2) + (0.5 * sum_pi2^2) - (0.25 * sum_pi4)
  })
  
  pisibs_cumulative <- prod(pisibs_loci, na.rm = TRUE)
  
  cat(label, "\n")
  cat("  N loci:", length(pisibs_loci), "\n")
  cat("  Mean PISibs per locus:", round(mean(pisibs_loci, na.rm=TRUE), 6), "\n")
  cat("  Cumulative PISibs:", formatC(pisibs_cumulative, format="e", digits=3), "\n\n")
  
  return(list(per_locus = pisibs_loci, cumulative = pisibs_cumulative))
}

# Dataset D: the 12 family individuals extracted (1378 SNPs)
pisibs_D <- calc_pisibs("fam_only.vcf",
                         "Dataset D (family individuals extracted, 1378 SNPs)")

#Identity by descent
library(SNPRelate)

# Define the VCF file and the GDS file to create
vcf_input <- "MP_combined.vcf"        
gds_output <- "temporary_data.gds"    # R will create this file from scratch right now

# 3. This step automatically creates the GDS file for you
snpgdsVCF2GDS(vcf_input, gds_output, method = "biallelic.only")

# 4. Open the file that R just built
genofile <- snpgdsOpen(gds_output)

# 5. Prune for Linkage Disequilibrium (highly recommended for VCF data)
set.seed(123) 
ld_pruned <- snpgdsLDpruning(genofile, ld.threshold = 0.2, autosome.only = TRUE)
snp_ids <- unlist(ld_pruned)

# Calculate the IBD metrics
ibd_output <- snpgdsIBDMoM(genofile, snp.id = snp_ids, num.thread = 2)

# 7. Turn the results into a readable R dataframe
ibd_dataframe <- snpgdsIBDSelection(ibd_output)

# 8. Close the file connection safely
snpgdsClose(genofile)

# View the relationship matrix
# View(ibd_dataframe)

# Vector of the target sample IDs
target_ids <- c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", 
                "SC7", "SC8", "SC9", "SC10", "SC11", "SC12")

# Subset the dataframe where BOTH individuals are in the target list
# (all pairwise combinations within the SC1-12 group)
sc_subset <- subset(ibd_dataframe, ID1 %in% target_ids & ID2 %in% target_ids)

# Export this group to a CSV file
write.csv(sc_subset, file = "SC1_to_SC12_IBD_matrix.csv", row.names = FALSE)

# ML
genofile <- snpgdsOpen(gds_output)

# 2. Run the Maximum Likelihood Estimation (MLE)
# Note: num.thread utilizes multiple CPU cores to speed up the intense MLE math
ibd_mle_output <- snpgdsIBDMLE(genofile, 
                               snp.id = snp_ids, 
                               num.thread = 4)

# 3. Format the raw MLE output into a data frame
ibd_mle_df <- snpgdsIBDSelection(ibd_mle_output)

# 4. Always close the file connection safely
snpgdsClose(genofile)

# Filter for the SC1-12 individuals
target_ids <- c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", 
                "SC7", "SC8", "SC9", "SC10", "SC11", "SC12")

sc_mle_subset <- subset(ibd_mle_df, ID1 %in% target_ids & ID2 %in% target_ids)

# 7. Export the MLE results to a new CSV file
write.csv(sc_mle_subset, file = "SC1_to_SC12_IBD_MLE_results.csv", row.names = FALSE)