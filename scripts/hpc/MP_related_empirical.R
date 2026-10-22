# HPC step: empirical relatedness on dataset D (related package)
# Inputs: MP_combined_related_input.txt (from 06a_MP_related_input.R)
# Outputs: MP_emp_lynchrd.rds, MP_emp_ML.rds

# .libPaths("<path to R library>")
library(related)

input_combined <- readgenotypedata("MP_combined_related_input.txt")

# lynchrd best
# run empirically
emp_lynchrd <- coancestry(input_combined$gdata,
                           lynchrd     = 1,
                           error.rates = 0.01)

emprel_lynchrd <- emp_lynchrd$relatedness

# Family pairs only
family_ids <- paste0("SC", 1:12)
fam_pairs  <- emprel_lynchrd[
  emprel_lynchrd$ind1.id %in% family_ids &
  emprel_lynchrd$ind2.id %in% family_ids, ]

cat("Mother-pup pairs (lynchrd):\n")
mother_pairs <- fam_pairs[fam_pairs$ind1.id == "SC12" |
                            fam_pairs$ind2.id == "SC12", ]
print(mother_pairs[, c("ind1.id","ind2.id","lynchrd")])

saveRDS(emp_lynchrd, "MP_emp_lynchrd.rds")
cat("LynchRD empirical complete\n")

# likelihood estimates
# TrioML and DyadML empirical
emp_ML <- coancestry(input_combined$gdata,
                      trioml               = 1,
                      dyadml               = 1,
                      error.rates          = 0.01,
                      trioml.num.reference = 100)

emprel_ML <- emp_ML$relatedness

fam_pairs_ML <- emprel_ML[
  emprel_ML$ind1.id %in% family_ids &
  emprel_ML$ind2.id %in% family_ids, ]

cat("Mother-pup pairs (trioml and dyadml):\n")
mother_pairs_ML <- fam_pairs_ML[
  fam_pairs_ML$ind1.id == "SC12" |
  fam_pairs_ML$ind2.id == "SC12", ]
print(mother_pairs_ML[, c("ind1.id","ind2.id","trioml","dyadml")])

saveRDS(emp_ML, "MP_emp_ML.rds")
cat("ML empirical complete\n")
