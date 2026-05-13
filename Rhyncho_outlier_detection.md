# Seascape Genomics Pipeline
## *Rhynchobatus australiae* and *Rhynchobatus djiddensis*
### Outlier Detection and GEA Pipeline

---

## Table of Contents
1. [Data Preparation — Genotype Matrix](#1-data-preparation--genotype-matrix)
2. [OutFLANK — Differentiation-based Outlier Detection](#2-outflank--differentiation-based-outlier-detection)
3. [pcadapt — Structure-based Outlier Detection](#3-pcadapt--structure-based-outlier-detection)
4. [Environmental Variable Extraction](#4-environmental-variable-extraction)
5. [Environmental Variable Correlation Heatmap](#5-environmental-variable-correlation-heatmap)
6. [Marine Distances and dbMEMs](#6-marine-distances-and-dbmems)
7. [RDA Pipeline — Variable Selection](#7-rda-pipeline--variable-selection)
8. [Partial RDA and Variance Partitioning](#8-partial-rda-and-variance-partitioning)
9. [Outlier Detection from Partial RDA](#9-outlier-detection-from-partial-rda)
10. [SNP Cluster Plot](#10-snp-cluster-plot)
11. [Candidate-Environment Correlations](#11-candidate-environment-correlations)
12. [Top 10 Candidate Loci — Pearson Correlation Plots](#12-top-10-candidate-loci--pearson-correlation-plots)
13. [LFMM — Latent Factor Mixed Models](#13-lfmm--latent-factor-mixed-models)
14. [Four SNP Dataset Compilation](#14-four-snp-dataset-compilation)
15. [Venn Diagram](#15-venn-diagram)
16. [RDA Triplot](#16-rda-triplot)
17. [Environmental Variable Distribution Maps](#17-environmental-variable-distribution-maps)

---

## Notes

- **Primary datasets:** RA_set1 (MAC>=3, 27,562 SNPs, n=71), RD_set1 (MAC>=3, 27,366 SNPs, n=47)
- **Robustness datasets:** RA_set2 (MAC>=2, 30,380 SNPs, n=71), RD_set2 (DP8, 18,942 SNPs, n=47)
- **Change file paths and species prefixes (RA/RD) as needed for each species/dataset**
- **K values:** RA = 4, RD = 2 (from pcadapt GIF analysis)
- **Bio-ORACLE** temperature and oxygen variables extracted using **surface layers** (depthsurf) to avoid deep-water extraction artefacts at sites with steep bathymetric gradients
- **VIF thresholds:** 10 for env-only, 15 for combined env+spatial

---

## 1. Data Preparation — Genotype Matrix

```r
library(vcfR)
library(readxl)

obj.vcfR   <- read.vcfR("RApop_set1.snps.vcf")
geno       <- extract.gt(obj.vcfR)
position   <- getPOS(obj.vcfR)
chromosome <- getCHROM(obj.vcfR)

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

G_imp <- apply(G, 1, function(x){
  x[is.na(x)] <- as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
  x
})
G_imp <- t(G_imp)

locus_names_full <- paste(chromosome, position, sep = "_")
colnames(G_imp)  <- locus_names_full
rownames(G_imp)  <- colnames(geno)

popmap <- read.table("RApopmap.txt", header = FALSE)
colnames(popmap) <- c("ID", "pop")

coords <- read_excel("RA_coordinates.xlsx", col_types = "text")
colnames(coords) <- c("ID", "Location", "lat", "lon")
coords$lat <- as.numeric(gsub(",", ".", trimws(coords$lat)))
coords$lon <- as.numeric(gsub(",", ".", trimws(coords$lon)))

all(rownames(G_imp) == coords$ID)  # must be TRUE
dim(G_imp)                          # 71 x 27562

save(G, G_imp, locus_names_full, chromosome, position,
     file = "RA_set1_genotype_matrix.RData")
```

---

## 2. OutFLANK — Differentiation-based Outlier Detection

```r
library(OutFLANK)

G_of <- G
G_of[is.na(G_of)] <- 9

my_fst <- MakeDiploidFSTMat(
  SNPmat     = t(G_of),
  locusNames = locus_names_full,
  popNames   = popmap$pop
)

out <- OutFLANK(
  my_fst,
  LeftTrimFraction  = 0.05,
  RightTrimFraction = 0.05,
  Hmin              = 0.1,
  NumberOfSamples   = 7,    # RA: 7 populations; change for RD
  qthreshold        = 0.05
)

P1 <- pOutlierFinderChiSqNoCorr(
  my_fst,
  Fstbar     = out$FSTNoCorrbar,
  dfInferred = out$dfInferred,
  qthreshold = 0.05,
  Hmin       = 0.1
)

outliers_outflank <- P1$LocusName[which(P1$OutlierFlag == TRUE)]
cat("OutFLANK outliers:", length(outliers_outflank), "\n")

pdf("outflank_fst_distribution.pdf", width = 7, height = 6)
  OutFLANKResultsPlotter(out, withOutliers = TRUE, NoCorr = TRUE,
                         Hmin = 0.1, binwidth = 0.005, titletext = "")
dev.off()

save(G_of, my_fst, out, P1, outliers_outflank,
     file = "RA_set1_outflank_results.RData")
```

---

## 3. pcadapt — Structure-based Outlier Detection

```r
library(pcadapt)
library(qvalue)

ra_pcadapt <- read.pcadapt("RApop_set1.snps.vcf", type = "vcf")

# Determine optimal K via GIF
gif_values <- data.frame(K = 1:20, GIF = NA)
for(k in 1:20){
  x_temp <- pcadapt(input = ra_pcadapt, K = k)
  gif_values$GIF[k] <- x_temp$gif
}
print(gif_values, digits = 6)

x_K20 <- pcadapt(input = ra_pcadapt, K = 20)
plot(x_K20, option = "screeplot")

# Run with optimal K
x <- pcadapt(input = ra_pcadapt, K = 4)  # change K for RD
plot(x, option = "scores", pop = popmap$pop)

# Extract outliers — Storey q-value < 0.05
qval                   <- qvalue(x$pvalues)$qvalues
outliers_pcadapt       <- which(qval < 0.05)
outliers_pcadapt_names <- locus_names_full[outliers_pcadapt]
cat("pcadapt outliers (q < 0.05):", length(outliers_pcadapt_names), "\n")

# Manhattan plot
fdr_threshold <- max(x$pvalues[qval < 0.05], na.rm = TRUE)
fdr_line      <- -log10(fdr_threshold)

pdf("manhattan_pcadapt.pdf", width = 10, height = 5)
  plot(seq_along(x$pvalues), -log10(x$pvalues),
       pch = 19, cex = 0.3, col = "#050404",
       xlab = "SNP index",
       ylab = expression(-log[10](italic(p)~value)),
       cex.lab = 1.2, cex.axis = 1.1, las = 1, bty = "l")
  abline(h = fdr_line, col = "#E69F00", lty = 2, lwd = 1.5)
  legend("topright", legend = "FDR (q < 0.05)",
         col = "#E69F00", lty = 2, lwd = 1.5, bty = "n", cex = 0.9)
dev.off()

save(x, gif_values, qval, outliers_pcadapt, outliers_pcadapt_names,
     file = "RA_set1_pcadapt_results.RData")
```

---

## 4. Environmental Variable Extraction

```r
library(sdmpredictors)
library(biooracler)
library(raster)

# MARSPEC variables
marspec_layers  <- c("MS_biogeo13_sst_mean_5m", "MS_biogeo16_sst_range_5m",
                     "MS_biogeo08_sss_mean_5m", "MS_biogeo11_sss_range_5m",
                     "MS_bathy_5m")
marspec_rasters <- load_layers(marspec_layers)
marspec_values  <- raster::extract(marspec_rasters, coords[, c("lon", "lat")])
marspec_df      <- data.frame(ID = coords$ID, marspec_values)

# Bio-ORACLE point extraction function
extract_biooracle_points <- function(dataset_id, variables, coords){
  results <- data.frame(ID = coords$ID)
  for(var in variables){
    vals <- numeric(nrow(coords))
    for(i in 1:nrow(coords)){
      tryCatch({
        tiny <- list(latitude  = c(coords$lat[i] - 0.05, coords$lat[i] + 0.05),
                     longitude = c(coords$lon[i] - 0.05, coords$lon[i] + 0.05))
        d <- download_layers(dataset_id = dataset_id, variables = var,
                             constraints = tiny, fmt = "csv")
        vals[i] <- mean(d[[var]], na.rm = TRUE)
      }, error = function(e){ vals[i] <<- NA })
    }
    results[[var]] <- vals
  }
  return(results)
}

# Extract all Bio-ORACLE variables — surface layers only
thetao_df <- extract_biooracle_points("thetao_baseline_2000_2019_depthsurf",
               c("thetao_mean", "thetao_ltmin", "thetao_ltmax", "thetao_range"), coords)
chl_df    <- extract_biooracle_points("chl_baseline_2000_2018_depthsurf",
               c("chl_mean", "chl_range"), coords)
o2_df     <- extract_biooracle_points("o2_baseline_2000_2018_depthsurf",
               c("o2_mean", "o2_ltmin"), coords)
phyc_df   <- extract_biooracle_points("phyc_baseline_2000_2020_depthsurf",
               c("phyc_mean"), coords)
sws_df    <- extract_biooracle_points("sws_baseline_2000_2019_depthsurf",
               c("sws_mean", "sws_range"), coords)
ph_df     <- extract_biooracle_points("ph_baseline_2000_2018_depthsurf",
               c("ph_mean"), coords)

# Combine all variables
env_data <- Reduce(function(x, y) merge(x, y, by = "ID"),
                   list(marspec_df, thetao_df, chl_df, o2_df, phyc_df, sws_df, ph_df))
env_data <- merge(env_data, popmap, by = "ID")

colSums(is.na(env_data))
dim(env_data)  # 71 x 19 (ID + 17 vars + pop)

save(env_data, coords, file = "RA_environmental_data.RData")
```

---

## 5. Environmental Variable Correlation Heatmap

```r
library(corrplot)

# Clean labels for plotting — does not rename env_data columns
clean_names <- c(
  "MS_biogeo13_sst_mean_5m"  = "SST mean",
  "MS_biogeo16_sst_range_5m" = "SST range",
  "MS_biogeo08_sss_mean_5m"  = "SSS mean",
  "MS_biogeo11_sss_range_5m" = "SSS range",
  "MS_bathy_5m"              = "Bathymetry",
  "thetao_mean"              = "Temp mean",
  "thetao_ltmin"             = "Temp min",
  "thetao_ltmax"             = "Temp max",
  "thetao_range"             = "Temp range",
  "chl_mean"                 = "Chl mean",
  "chl_range"                = "Chl range",
  "o2_mean"                  = "DO mean",
  "o2_ltmin"                 = "DO min",
  "phyc_mean"                = "Phytoplankton",
  "sws_mean"                 = "Current mean",
  "sws_range"                = "Current range",
  "ph_mean"                  = "pH"
)

env_vars_all <- env_data[, !colnames(env_data) %in% c("ID", "pop")]
rownames(env_vars_all) <- coords$ID

env_matrix <- as.matrix(env_vars_all)
colnames(env_matrix) <- clean_names[colnames(env_matrix)]

# Reorder so pH (last variable) appears as a row — avoids clipping
new_order  <- rev(clean_names)
cor_matrix <- cor(env_matrix, method = "pearson", use = "complete.obs")
cor_matrix <- cor_matrix[new_order, new_order]

pdf("env_correlation_heatmap_all17.pdf", width = 10, height = 10)
  corrplot(cor_matrix,
           method      = "color",
           type        = "upper",
           tl.cex      = 0.85,
           tl.col      = "black",
           tl.srt      = 45,
           addCoef.col = "black",
           number.cex  = 0.55,
           col         = colorRampPalette(c("#0072B2", "white", "#D55E00"))(200),
           diag        = FALSE,
           mar         = c(0, 0, 0, 2))
dev.off()
```

---

## 6. Marine Distances and dbMEMs

```r
library(adespatial)

write.table(data.frame(Name = coords$ID, Lon = coords$lon, Lat = coords$lat),
            file = "RA_individuals.txt", sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run marineDistances (requires GSHHS_h_L1.shp; runtime ~8 hours on laptop)
source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")
contour(global.polygon = "GSHHS_h_L1.shp", file = "RA_individuals.txt",
        file.sep = " ", file.dec = ".", file.strucutre = 1,
        file.header = FALSE, resolution = 0.5, buffer = c(5, 5, 5, 5),
        export.file = TRUE)

marine_dist <- as.matrix(read.table("Contour _ Pairwise Marine Distances.txt",
                                     header = TRUE, row.names = 1,
                                     sep = " ", dec = "."))

# PCoA to Euclidean distances
marine_pcoa   <- cmdscale(as.dist(marine_dist), k = nrow(marine_dist) - 1, eig = TRUE)
n_pos         <- sum(marine_pcoa$eig > 0)
marine_euclid <- marine_pcoa$points[, 1:n_pos]
rownames(marine_euclid) <- coords$ID

# dbMEMs — positive eigenvalues only
dbmems <- dbmem(dist(marine_euclid), MEM.autocor = "positive")
cat("dbMEMs computed:", ncol(dbmems), "\n")

save(marine_dist, marine_euclid, dbmems, file = "RA_spatial_data.RData")
```

---

## 7. RDA Pipeline — Variable Selection

```r
library(vegan)
library(usdm)

env_data     <- env_data[match(coords$ID, env_data$ID), ]
rownames(env_data) <- env_data$ID

env_vars_all <- env_data[, !colnames(env_data) %in% c("ID", "pop")]
rownames(env_vars_all) <- coords$ID

mem_df <- data.frame(MEM1 = as.numeric(dbmems[, 1]),
                     MEM2 = as.numeric(dbmems[, 2]),
                     row.names = coords$ID)

# VIF pre-filter (threshold 10) on all 17 variables
vif_result <- vifstep(env_vars_all, th = 10)
env_vars10  <- env_vars_all[, vif_result@results$Variables]
cat("VIF-retained variables:", colnames(env_vars10), "\n")

# ordistep ENV only — use 5000 SNP subset if LAPACK errors occur
rda_null_env <- rda(G_imp ~ 1, data = as.data.frame(env_vars10))
rda_full_env <- rda(G_imp ~ ., data = as.data.frame(env_vars10))
rda_step_env <- ordistep(rda_null_env, scope = formula(rda_full_env),
                          direction = "both", permutations = 500, trace = TRUE)
selected_env <- attr(rda_step_env$terms, "term.labels")
cat("Selected env variables:", selected_env, "\n")

# ordistep MEMs only
rda_null_sp <- rda(G_imp ~ 1, data = mem_df)
rda_full_sp <- rda(G_imp ~ ., data = mem_df)
rda_step_sp <- ordistep(rda_null_sp, scope = formula(rda_full_sp),
                         direction = "both", permutations = 500, trace = TRUE)
selected_sp <- attr(rda_step_sp$terms, "term.labels")
cat("Selected spatial variables:", selected_sp, "\n")

# VIF on combined set (threshold 15)
env_retained <- as.data.frame(env_vars10)[, selected_env, drop = FALSE]
mem_retained <- mem_df[, selected_sp, drop = FALSE]
env_combined <- data.frame(env_retained, mem_retained)
rownames(env_combined) <- coords$ID
print(vif(env_combined))

env_only_final <- env_combined[, !grepl("^MEM", colnames(env_combined)), drop = FALSE]
mem_only_final <- env_combined[,  grepl("^MEM", colnames(env_combined)), drop = FALSE]

cat("Final env variables:", colnames(env_only_final), "\n")
cat("Final spatial variables:", colnames(mem_only_final), "\n")
```

---

## 8. Partial RDA and Variance Partitioning

```r
rda_partial <- rda(G_imp ~ . + Condition(as.matrix(mem_only_final)),
                   data = env_only_final)

anova_partial      <- anova(rda_partial, permutations = 999)
anova_partial_axes <- anova(rda_partial, by = "axis", permutations = 999)
print(anova_partial)
print(anova_partial_axes)

vp <- varpart(G_imp, env_only_final, mem_only_final)
print(vp)

# Test individual fractions
rda_sp_partial   <- rda(G_imp ~ . + Condition(as.matrix(env_only_final)),
                        data = mem_only_final)
anova_sp_partial <- anova(rda_sp_partial, permutations = 999)
print(anova_sp_partial)

rda_env_only <- rda(G_imp ~ ., data = env_only_final)
anova_env    <- anova(rda_env_only, permutations = 999)
cat("Env-only adj R2:", RsquareAdj(rda_env_only)$adj.r.squared, "\n")

rda_full   <- rda(G_imp ~ ., data = data.frame(env_only_final, mem_only_final))
anova_full <- anova(rda_full, permutations = 999)
cat("Full model adj R2:", RsquareAdj(rda_full)$adj.r.squared, "\n")

save(rda_partial, anova_partial, anova_partial_axes,
     vp, env_combined, env_only_final, mem_only_final,
     file = "RA_set1_RDA_corrected_results.RData")
```

---

## 9. Outlier Detection from Partial RDA

```r
sig_axes <- which(anova_partial_axes[, "Pr(>F)"] < 0.05)
sig_axes <- sig_axes[sig_axes <= rda_partial$CCA$rank]
cat("Significant axes:", sig_axes, "\n")

locus_scores <- scores(rda_partial, display = "species", choices = sig_axes)

outliers_rda <- c()
for(ax in 1:ncol(locus_scores)){
  ax_scores  <- locus_scores[, ax]
  thresh_pos <- mean(ax_scores) + 3 * sd(ax_scores)
  thresh_neg <- mean(ax_scores) - 3 * sd(ax_scores)
  outliers_ax <- rownames(locus_scores)[ax_scores > thresh_pos | ax_scores < thresh_neg]
  cat("Axis", ax, ":", length(outliers_ax), "outliers\n")
  outliers_rda <- unique(c(outliers_rda, outliers_ax))
}
cat("Total RDA environmental candidates:", length(outliers_rda), "\n")

write.csv(data.frame(locus = outliers_rda),
          "outliers_rda_RA_set1.csv", row.names = FALSE)
```

---

## 10. SNP Cluster Plot

```r
# Arrow label order (top to bottom of clustered arrows):
# SSS range (RDA2 = -0.0142), Phytoplankton (-0.0159), DO min (-0.0200)
# Labels added manually in Canva

all_scores   <- scores(rda_partial, display = "species", choices = c(1, 2))
is_candidate <- rownames(all_scores) %in% outliers_rda

ax1_var <- round(rda_partial$CCA$eig[1] / sum(rda_partial$CCA$eig) * 100, 1)
ax2_var <- round(rda_partial$CCA$eig[2] / sum(rda_partial$CCA$eig) * 100, 1)

env_arrows <- scores(rda_partial, display = "bp")

scale <- 0.25
xlim  <- c(-0.18, 0.18)
ylim  <- c(-0.18, 0.18)

pdf("partial_RDA_australiae.pdf", width = 7, height = 7)

  plot(all_scores[!is_candidate, 1], all_scores[!is_candidate, 2],
       pch = 19, col = adjustcolor("#85C0F9", alpha.f = 0.25),
       cex = 0.3, xlim = xlim, ylim = ylim,
       xlab = paste0("RDA1 (", ax1_var, "%)"),
       ylab = paste0("RDA2 (", ax2_var, "%)"),
       cex.lab = 1.2, cex.axis = 1.1, las = 1)

  points(all_scores[is_candidate, 1], all_scores[is_candidate, 2],
         pch = 19, col = adjustcolor("#F5793A", alpha.f = 0.75), cex = 0.3)

  arrows(0, 0, env_arrows[, 1] * scale, env_arrows[, 2] * scale,
         length = 0.08, col = "black", lwd = 1.8)

  legend("topright",
         legend = c(paste0("Neutral SNPs (n=", formatC(sum(!is_candidate), big.mark = ",")),
                    paste0("Candidate SNPs (n=", formatC(sum(is_candidate), big.mark = ","), ")")),
         col = c("#85C0F9", "#F5793A"), pch = 19, pt.cex = 1.2, cex = 0.9, bty = "n")

dev.off()
```

---

## 11. Candidate-Environment Correlations

```r
outlier_indices <- which(locus_names_full %in% outliers_rda)
G_candidates    <- G_imp[, outlier_indices]

cor_results <- data.frame()

for(var in colnames(env_only_final)){
  env_var <- env_only_final[, var]
  r_vals  <- apply(G_candidates, 2, function(g) cor(g, env_var, method = "pearson"))
  p_vals  <- apply(G_candidates, 2, function(g) cor.test(g, env_var)$p.value)

  cor_results <- rbind(cor_results, data.frame(
    Environmental_variable = var,
    n_significant          = sum(p_vals < 0.05, na.rm = TRUE),
    pct_significant        = round(sum(p_vals < 0.05, na.rm = TRUE) / ncol(G_candidates) * 100, 1),
    mean_abs_r             = round(mean(abs(r_vals), na.rm = TRUE), 3),
    max_abs_r              = round(max(abs(r_vals), na.rm = TRUE), 3),
    n_strong               = sum(abs(r_vals) > 0.5, na.rm = TRUE)
  ))
}

cor_results <- cor_results[order(cor_results$pct_significant, decreasing = TRUE), ]
print(cor_results)
write.csv(cor_results, "RA_set1_candidate_env_correlations.csv", row.names = FALSE)
```

---

## 12. Top 10 Candidate Loci — Pearson Correlation Plots

```r
r_matrix <- matrix(NA, nrow = ncol(G_candidates), ncol = ncol(env_only_final))
rownames(r_matrix) <- colnames(G_candidates)
colnames(r_matrix) <- colnames(env_only_final)

for(v in 1:ncol(env_only_final)){
  env_var <- env_only_final[, v]
  r_matrix[, v] <- apply(G_candidates, 2, function(g)
    cor(g, env_var, method = "pearson"))
}

max_r_per_locus <- apply(abs(r_matrix), 1, max)
top10_idx       <- order(max_r_per_locus, decreasing = TRUE)[1:10]
top10_loci      <- rownames(r_matrix)[top10_idx]

cat("Top 10 candidate loci:\n")
print(data.frame(
  locus     = top10_loci,
  max_abs_r = round(max_r_per_locus[top10_idx], 3),
  best_var  = colnames(r_matrix)[apply(abs(r_matrix[top10_idx, ]), 1, which.max)]
))

var_labels <- c("MS_biogeo11_sss_range_5m" = "SSS range",
                "MS_bathy_5m"              = "Bathymetry",
                "phyc_mean"                = "Phytoplankton",
                "ph_mean"                  = "pH",
                "o2_ltmin"                 = "DO min")

pop_levels <- c("TZ", "MD", "SY", "AU", "RI", "MZ", "MA")
pop_cols   <- c("#E69F00", "#56B4E9", "#009E73", "#8B008B", "#D55E00", "#0072B2", "#FF69B4")
names(pop_cols) <- pop_levels
pop_factor <- popmap$pop[match(coords$ID, popmap$ID)]
pt_cols    <- pop_cols[pop_factor]

pdf("top10_candidate_correlations.pdf", width = 7, height = 14)
par(mfrow = c(5, 2), mar = c(4, 4, 3, 1))

for(i in 1:10){
  locus   <- top10_loci[i]
  best_v  <- colnames(r_matrix)[which.max(abs(r_matrix[locus, ]))]
  env_var <- env_only_final[, best_v]
  geno    <- G_candidates[, locus]
  r_val   <- round(r_matrix[locus, best_v], 3)
  p_val   <- cor.test(geno, env_var)$p.value

  plot(env_var, geno,
       pch = 19, col = pt_cols, cex = 1.0,
       xlab = var_labels[best_v], ylab = "Genotype frequency",
       main = paste0("Locus ", i, "\nr = ", r_val,
                     ", p = ", formatC(p_val, format = "e", digits = 2)),
       cex.main = 0.85, cex.lab = 0.9, las = 1, bty = "l")
  abline(lm(geno ~ env_var), col = "red", lwd = 1.5)
}

dev.off()
```

---

## 13. LFMM — Latent Factor Mixed Models

```r
library(LEA)
library(qvalue)
library(vegan)

# PCA of all 17 environmental variables
env_all_lfmm <- env_data[, !colnames(env_data) %in% c("ID", "pop")]
rownames(env_all_lfmm) <- coords$ID

env_pca     <- prcomp(env_all_lfmm, scale = TRUE, center = TRUE)
eigenvalues <- env_pca$sdev^2 / sum(env_pca$sdev^2)
bstick_vals <- bstick(length(env_pca$sdev))

n_pcs <- sum(eigenvalues > bstick_vals)
cat("PCs retained:", n_pcs, "\n")
cat("Cumulative variance explained:", round(sum(eigenvalues[1:n_pcs]) * 100, 1), "%\n")
if(n_pcs == 0){ n_pcs <- 1; cat("Warning: defaulting to 1 PC\n") }

env_pcs <- env_pca$x[, 1:n_pcs, drop = FALSE]

# Broken stick plot
pdf("lfmm_brokenstack_set1.pdf", width = 6, height = 5)
  plot(eigenvalues, type = "b", pch = 19, xlab = "PC",
       ylab = "Proportion of variance",
       ylim = c(0, max(c(eigenvalues, bstick_vals)) + 0.05),
       cex.lab = 1.2, cex.axis = 1.1, las = 1, bty = "l")
  lines(bstick_vals, col = "red", lty = 2)
  legend("topright", legend = c("Observed", "Broken stick"),
         col = c("black", "red"), lty = c(1, 2), pch = c(19, NA),
         bty = "n", cex = 0.9)
dev.off()

# Run LFMM2
mod_lfmm <- lfmm2(input = G_imp, env = env_pcs, K = 4)

pv  <- lfmm2.test(object = mod_lfmm, input = G_imp, env = env_pcs, full = TRUE)
gif <- pv$gif
cat("GIF:", round(gif, 3), "\n")

# FDR < 1%
qv                  <- qvalue(pv$pvalues)$qvalues
outliers_lfmm       <- which(qv < 0.01)
outliers_lfmm_names <- locus_names_full[outliers_lfmm]
cat("LFMM outliers (FDR < 1%):", length(outliers_lfmm_names), "\n")

# Manhattan plot
fdr_threshold <- max(pv$pvalues[qv < 0.01], na.rm = TRUE)
fdr_line      <- -log10(fdr_threshold)

pdf("lfmm_manhattan_set1.pdf", width = 10, height = 5)
  plot(seq_along(pv$pvalues), -log10(pv$pvalues),
       pch = 19, cex = 0.3, col = "#050404",
       xlab = "SNP index",
       ylab = expression(-log[10](italic(p)~value)),
       cex.lab = 1.2, cex.axis = 1.1, las = 1, bty = "l")
  abline(h = fdr_line, col = "#87c1f7", lty = 2, lwd = 1.5)
  legend("topright", legend = "FDR (q < 0.01)",
         col = "#87c1f7", lty = 2, lwd = 1.5, bty = "n", cex = 0.9)
dev.off()

save(mod_lfmm, pv, gif, qv, outliers_lfmm, outliers_lfmm_names,
     env_pca, env_pcs, n_pcs,
     file = "RA_set1_LFMM_corrected_results.RData")
```

---

## 14. Four SNP Dataset Compilation

```r
all_outliers   <- union(union(outliers_rda, outliers_pcadapt_names), outliers_lfmm_names)
env_candidates <- outliers_rda
high_conf      <- union(union(
                    intersect(outliers_rda, outliers_pcadapt_names),
                    intersect(outliers_rda, outliers_lfmm_names)),
                    intersect(outliers_pcadapt_names, outliers_lfmm_names))
neutral        <- locus_names_full[!locus_names_full %in% all_outliers]

cat("All putative outliers:   ", length(all_outliers),
    "(", round(length(all_outliers) / length(locus_names_full) * 100, 1), "%)\n")
cat("Environmental candidates:", length(env_candidates),
    "(", round(length(env_candidates) / length(locus_names_full) * 100, 1), "%)\n")
cat("High confidence (>=2):   ", length(high_conf),
    "(", round(length(high_conf) / length(locus_names_full) * 100, 1), "%)\n")
cat("Neutral:                 ", length(neutral),
    "(", round(length(neutral) / length(locus_names_full) * 100, 1), "%)\n")

save(all_outliers, env_candidates, high_conf, neutral,
     outliers_rda, outliers_pcadapt_names, outliers_lfmm_names,
     file = "RA_set1_four_datasets.RData")
```

---

## 15. Venn Diagram

```r
library(ggVennDiagram)
library(ggplot2)

venn_list <- list(RDA     = outliers_rda,
                  pcadapt = outliers_pcadapt_names,
                  LFMM    = outliers_lfmm_names)

p_venn <- ggVennDiagram(venn_list,
                         label_alpha = 0,
                         label       = "count",
                         set_size    = 4.5) +
  scale_fill_gradientn(colours = c("white", "#86c0f7")) +
  scale_colour_manual(values = c("black", "black", "black")) +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("venn_set1.pdf", p_venn, width = 6, height = 6)
```

---

## 16. RDA Triplot

```r
# scaling=2 preserves correct arrow length interpretation
# Individual scores rescaled to fit SNP/arrow space
# Arrow labels added manually in Canva

ind_scores    <- scores(rda_partial, display = "sites",   scaling = 2, choices = 1:3)
snp_scores    <- scores(rda_partial, display = "species", scaling = 2, choices = 1:3)
env_scores_12 <- scores(rda_partial, display = "bp",      scaling = 2, choices = 1:2)
env_scores_13 <- scores(rda_partial, display = "bp",      scaling = 2, choices = c(1, 3))

ax1_var <- round(summary(eigenvals(rda_partial, model = "constrained"))["Proportion Explained", "RDA1"] * 100, 1)
ax2_var <- round(summary(eigenvals(rda_partial, model = "constrained"))["Proportion Explained", "RDA2"] * 100, 1)
ax3_var <- round(summary(eigenvals(rda_partial, model = "constrained"))["Proportion Explained", "RDA3"] * 100, 1)

pop_levels <- c("TZ", "MD", "SY", "AU", "RI", "MZ", "MA")
pop_cols   <- c("#E69F00",  # TZ - gold
                "#56B4E9",  # MD - light blue
                "#009E73",  # SY - green
                "#8B008B",  # AU - dark purple
                "#D55E00",  # RI - vermillion
                "#0072B2",  # MZ - dark blue
                "#FF69B4")  # MA - hot pink
names(pop_cols) <- pop_levels

pop_labels <- c("TZ" = "Tanzania",
                "MD" = "Madagascar",
                "SY" = "Seychelles",
                "AU" = "Australia",
                "RI" = "Reunion Island",
                "MZ" = "Mozambique",
                "MA" = "Malaysia")

pop_factor <- popmap$pop[match(coords$ID, popmap$ID)]

ind_scale_factor <- max(abs(snp_scores[, 1:2])) / max(abs(ind_scores[, 1:2])) * 0.9
ind_scores_plot  <- ind_scores * ind_scale_factor

arrow_scale <- 0.15

# Panel A: RDA1 vs RDA2
pdf("rda_triplot_set1_ax12.pdf", width = 7, height = 7)

  plot_range_12 <- max(abs(c(snp_scores[, 1:2],
                              env_scores_12 * arrow_scale))) * 1.3

  plot(snp_scores[, 1], snp_scores[, 2],
       pch      = 19,
       col      = adjustcolor("grey50", alpha.f = 0.25),
       cex      = 0.3,
       xlim     = c(-plot_range_12, plot_range_12),
       ylim     = c(-plot_range_12, plot_range_12),
       xlab     = paste0("RDA1 (", ax1_var, "%)"),
       ylab     = paste0("RDA2 (", ax2_var, "%)"),
       cex.lab  = 1.2,
       cex.axis = 1.1,
       las      = 1,
       bty      = "o")

  abline(h = 0, v = 0, lty = 2, col = "grey70", lwd = 0.8)

  for(pop in pop_levels){
    idx <- which(pop_factor == pop)
    if(length(idx) > 0)
      points(ind_scores_plot[idx, 1], ind_scores_plot[idx, 2],
             pch = 21, bg = pop_cols[pop],
             col = "grey30", cex = 1.5, lwd = 0.8)
  }

  arrows(0, 0,
         env_scores_12[, 1] * arrow_scale,
         env_scores_12[, 2] * arrow_scale,
         length = 0.08, col = "black", lwd = 1.2)

  legend("topright",
         legend    = pop_labels[pop_levels],
         pt.bg     = pop_cols,
         col       = "grey30",
         pch       = 21,
         pt.cex    = 1.5,
         pt.lwd    = 0.8,
         cex       = 0.85,
         bty       = "n",
         y.intersp = 0.8)

dev.off()

# Panel B: RDA1 vs RDA3
pdf("rda_triplot_set1_ax13.pdf", width = 7, height = 7)

  plot_range_13 <- max(abs(c(snp_scores[, c(1, 3)],
                              env_scores_13 * arrow_scale))) * 1.3

  plot(snp_scores[, 1], snp_scores[, 3],
       pch      = 19,
       col      = adjustcolor("grey50", alpha.f = 0.25),
       cex      = 0.3,
       xlim     = c(-plot_range_13, plot_range_13),
       ylim     = c(-plot_range_13, plot_range_13),
       xlab     = paste0("RDA1 (", ax1_var, "%)"),
       ylab     = paste0("RDA3 (", ax3_var, "%)"),
       cex.lab  = 1.2,
       cex.axis = 1.1,
       las      = 1,
       bty      = "o")

  abline(h = 0, v = 0, lty = 2, col = "grey70", lwd = 0.8)

  for(pop in pop_levels){
    idx <- which(pop_factor == pop)
    if(length(idx) > 0)
      points(ind_scores_plot[idx, 1], ind_scores_plot[idx, 3],
             pch = 21, bg = pop_cols[pop],
             col = "grey30", cex = 1.5, lwd = 0.8)
  }

  arrows(0, 0,
         env_scores_13[, 1] * arrow_scale,
         env_scores_13[, 2] * arrow_scale,
         length = 0.08, col = "black", lwd = 1.2)

  legend("topright",
         legend    = pop_labels[pop_levels],
         pt.bg     = pop_cols,
         col       = "grey30",
         pch       = 21,
         pt.cex    = 1.5,
         pt.lwd    = 0.8,
         cex       = 0.85,
         bty       = "n",
         y.intersp = 0.8)

dev.off()
```

---

## 17. Environmental Variable Distribution Maps

```r
library(terra)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sdmpredictors)
library(biooracler)

world_sp     <- as(ne_countries(scale = "medium", returnclass = "sf"), "Spatial")
study_extent <- extent(30, 145, -30, 15)

# Load MARSPEC rasters
options(sdmpredictors_datadir = ".")
marspec_rasters <- load_layers(c("MS_biogeo11_sss_range_5m", "MS_bathy_5m"))
sss_range_t     <- rast(crop(marspec_rasters[["MS_biogeo11_sss_range_5m"]], study_extent))
bathy_t         <- rast(crop(marspec_rasters[["MS_bathy_5m"]], study_extent))
bathy_display   <- clamp(bathy_t, lower = -200, upper = 0)

# Load Bio-ORACLE rasters
constraints <- list(latitude = c(-30, 15), longitude = c(30, 145))
phyc_mean_r <- mean(download_layers("phyc_baseline_2000_2020_depthsurf",
                                     variables = "phyc_mean", constraints = constraints))
ph_mean_r   <- mean(download_layers("ph_baseline_2000_2018_depthsurf",
                                     variables = "ph_mean",   constraints = constraints))
o2_mean_r   <- mean(download_layers("o2_baseline_2000_2018_depthsurf",
                                     variables = "o2_ltmin",  constraints = constraints))

# Population colours
pop_levels <- c("TZ", "MD", "SY", "AU", "RI", "MZ", "MA")
pop_cols   <- c("#E69F00", "#56B4E9", "#009E73", "#8B008B", "#D55E00", "#0072B2", "#FF69B4")
names(pop_cols) <- pop_levels
pop_factor <- popmap$pop[match(coords$ID, popmap$ID)]
pt_cols    <- pop_cols[pop_factor]

# Variable units: SSS range=PSU, Bathymetry=m, Phytoplankton=mmol/m3,
#                 pH=unitless, DO min=mmol/m3
# Legend and units added manually in Canva

make_env_map_base <- function(raster_layer, title){
  plot(raster_layer, main = title, axes = TRUE, box = FALSE, cex.main = 0.9)
  plot(world_sp, add = TRUE, border = "grey60", col = "grey85")
  points(coords$lon, coords$lat, pch = 21, bg = pt_cols,
         col = "black", cex = 1.2, lwd = 0.5)
}

for(map in list(
  list(r = sss_range_t,   title = "SSS range",    file = "map_sss_range"),
  list(r = bathy_display, title = "Bathymetry",    file = "map_bathymetry"),
  list(r = phyc_mean_r,   title = "Phytoplankton", file = "map_phytoplankton"),
  list(r = ph_mean_r,     title = "pH",            file = "map_ph"),
  list(r = o2_mean_r,     title = "DO min",        file = "map_do_min")
)){
  pdf(paste0(map$file, ".pdf"), width = 9, height = 5)
    make_env_map_base(map$r, map$title)
  dev.off()
}
cat("All maps exported\n")
```

---
