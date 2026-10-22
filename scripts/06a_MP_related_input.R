# Prepare the relatedness inputs from dataset D (MP_combined.vcf: 47 individuals, 1378 SNPs)
# Outputs, uploaded to the cluster: MP_combined_related_input.txt and MP_WIO_frequencies.csv

library(vcfR)
library(related)

# ── Load combined VCF (all individuals, zero missingness) ─────
combined_vcf  <- read.vcfR("MP_combined.vcf",
                             verbose = FALSE)
geno_combined <- extract.gt(combined_vcf, element = "GT")
cat("Individuals:", ncol(geno_combined), "\n")
cat("Loci:", nrow(geno_combined), "\n")

# ── Convert to related input format ──────────────────────────
n_loci <- nrow(geno_combined)
n_ind  <- ncol(geno_combined)
out    <- matrix(0, nrow = n_ind, ncol = 1 + 2 * n_loci)
out[, 1] <- colnames(geno_combined)

for (i in 1:n_loci) {
  gt <- geno_combined[i, ]
  a1 <- ifelse(grepl("^0", gt), 1, 2)
  a2 <- ifelse(grepl("0$", gt), 1, 2)
  a1[is.na(gt)] <- 0
  a2[is.na(gt)] <- 0
  out[, 2 * i]     <- a1
  out[, 2 * i + 1] <- a2
}

out_df <- as.data.frame(out, stringsAsFactors = FALSE)
out_df[, 2:ncol(out_df)] <- lapply(out_df[, 2:ncol(out_df)], as.numeric)

write.table(out_df, "MP_combined_related_input.txt",
            row.names = FALSE, col.names = FALSE,
            sep = "\t", quote = FALSE)
cat("Written MP_combined_related_input.txt\n")

# ── Extract WIO allele frequencies ───────────────────────────
# WIO individuals = all except SC1-SC12
family_ids <- paste0("SC", 1:12)
wio_idx    <- !colnames(geno_combined) %in% family_ids
geno_wio   <- geno_combined[, wio_idx]
cat("WIO individuals for frequency estimation:", ncol(geno_wio), "\n")

wio_freqs <- apply(geno_wio, 1, function(gt) {
  gt <- gt[!is.na(gt)]
  a1 <- ifelse(grepl("^0", gt), 1, 2)
  a2 <- ifelse(grepl("0$", gt), 1, 2)
  sum(c(a1, a2) == 2) / length(c(a1, a2))
})

wio_freq_df <- data.frame(
  locus    = names(wio_freqs),
  freq_al1 = 1 - wio_freqs,
  freq_al2 = wio_freqs)
write.csv(wio_freq_df, "MP_WIO_frequencies.csv", row.names = FALSE)
cat("Written MP_WIO_frequencies.csv —", nrow(wio_freq_df), "loci\n")

cat("\nUpload MP_combined_related_input.txt and MP_WIO_frequencies.csv to HPC\n")
