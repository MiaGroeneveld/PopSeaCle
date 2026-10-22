# HPC step: relatedness simulation on dataset D (related package)
# Allele frequencies of the SWIO reference panel replace the frequencies in the combined input,
# 400 loci are sampled at random (Fortran backend memory limit), and 150 dyads per relationship
# class are simulated and scored with seven estimators.
# Inputs: MP_combined_related_input.txt, MP_WIO_frequencies.csv (from 06a_MP_related_input.R)
# Outputs: MP_simrel_150.rds, MP_sim_output_150.rds

# .libPaths("<path to R library>")
library(related)

# setwd("<working directory containing the input files>")

# Read combined genotypes natively (all individuals, zero missingness)
input_combined <- readgenotypedata("MP_combined_related_input.txt")
cat("Individuals:", input_combined$ninds, "\n")
cat("Loci:", input_combined$nloci, "\n")

# Substitute WIO frequencies
wio_freqs_csv   <- read.csv("MP_WIO_frequencies.csv")
locus_names     <- names(input_combined$freqs)
match_idx       <- match(wio_freqs_csv$locus, locus_names)

input_wio <- input_combined
for (i in which(!is.na(match_idx))) {
  idx <- match_idx[i]
  input_wio$freqs[[idx]] <- data.frame(
    allele    = as.integer(c(1, 2)),
    frequency = c(wio_freqs_csv$freq_al1[i], wio_freqs_csv$freq_al2[i])
  )
}
cat("WIO frequencies substituted:", sum(!is.na(match_idx)), "loci\n")

# Thin to 400 loci — Fortran backend memory limit
set.seed(123)
thin_idx   <- sample(1:length(input_wio$freqs), 400)
freqs_thin <- input_wio$freqs[thin_idx]
cat("Loci after thinning:", length(freqs_thin), "\n")

# Simulate 150 dyads per relationship class
cat("Simulating dyads...\n")
sim <- familysim(freqs_thin, 150)
cat("familysim complete\n")

# All 7 estimators
cat("Running coancestry (trioml and dyadml will be slow)...\n")
sim_output <- coancestry(sim,
                          lynchli              = 1,
                          lynchrd              = 1,
                          quellergt            = 1,
                          ritland              = 1,
                          wang                 = 1,
                          trioml               = 1,
                          dyadml               = 1,
                          error.rates          = 0.01,
                          trioml.num.reference = 100)

simrel <- cleanuprvals(sim_output$relatedness, 150)
saveRDS(simrel,     "MP_simrel_150.rds")
saveRDS(sim_output, "MP_sim_output_150.rds")
cat("Saved\n")

# Pearson correlations
relvals    <- c(rep(0.5, 150), rep(0.5, 150), rep(0.25, 150), rep(0, 150))
estimators <- c("trioml","wang","lynchli","lynchrd","ritland","quellergt","dyadml")
est_cols   <- c(5, 6, 7, 8, 9, 10, 11)

cat("\n--- Pearson correlations ---\n")
for (i in seq_along(estimators)) {
  cat(sprintf("%-12s %.4f\n", estimators[i],
              cor(relvals, simrel[, est_cols[i]], use = "complete.obs")))
}

# Means and SDs
classes <- list(PO = 1:150, FS = 151:300, HS = 301:450, UR = 451:600)
cat("\n--- Means and SDs per class ---\n")
for (i in seq_along(estimators)) {
  cat("\n", estimators[i], "\n")
  for (cls in names(classes)) {
    vals <- simrel[classes[[cls]], est_cols[i]]
    cat(sprintf("  %s: mean=%.4f  SD=%.4f\n",
                cls, mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE)))
  }
}
cat("\nDone — download MP_simrel_150.rds\n")