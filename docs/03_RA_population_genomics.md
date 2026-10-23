# 03. *R. australiae*: population genomics

Genetic diversity, differentiation, structure and isolation by distance and environment in the bottlenose wedgefish (*Rhynchobatus australiae*, RA), analysed separately for the four SNP datasets from [stage 2](02_RA_seascape_outliers.md): neutral, adaptive (any method), high confidence and environmental (RDA).

Everything is run in R from one script, [`scripts/03_RA_population_genomics.R`](../scripts/03_RA_population_genomics.R), which starts by loading the objects saved at the end of stage 2 (`RA_seascape_FINAL_checkpoint.RData`).

## Populations

MZ (Mozambique, n = 28), TZ (Tanzania, 8), MD (Madagascar, 2), RI (Réunion Island, 9), SY (Seychelles, 12), MA (Malaysia, 10), AU (Australia, 2). Madagascar and Australia are excluded from population-level diversity statistics and the site frequency spectra (SFS) because of their small sample sizes. SFS excluded in the manuscript.

## Sections

| Section | Analysis | Main settings | Main outputs |
|---|---|---|---|
| 1 | Diversity statistics: HO, HE, FIS | `hierfstat`; 95% confidence intervals for FIS from 500 bootstrap resamples of loci | `RA_set1_diversity_results.RData` |
| 2 | Allelic richness and private alleles | richness rarefied to the smallest sample (`hierfstat::allelic.richness`); private alleles counted as loci with `poppr::private_alleles` | `RA_set1_AR_PA_all_datasets.csv`, `RA_set1_AR_PA_results.RData` |
| 3 | Pairwise and global FST | Weir and Cockerham (`hierfstat`); 1000 bootstraps over loci; Benjamini-Hochberg correction | `RA_set1_FST_results.RData` |
| 4 | FST heatmaps | neutral above and adaptive below the diagonal (main text); high confidence and environmental (supplement) | `RA_set1_FST_heatmap*.pdf` |
| 5 | DAPC | PCs explaining 80% of variance; `find.clusters` for K = 1 to 10 (BIC); `xvalDapc` (100 replicates, 90% training set) and `optim.a.score` (100 simulations); with and without population priors | `RA_set1_DAPC_summary_table.csv`, one results file per dataset |
| 6 | Adaptive DAPC without Australia | Australia (n = 2) removed because its extreme scores compress the rest of the plot | `RA_set1_adaptive_dapc_noAU_results.RData` |
| 7 | sNMF | `LEA`; K = 1 to 10, 10 replicates, alpha = 100, K chosen by cross-entropy; barplots with `pophelper` and colours fixed to geographic groups | `RA_set1_snmf_summary.csv`, barplots |
| 8 | Site frequency spectra | folded (minor allele frequency), per population, four datasets (excluded from manuscript) | `RA_set1_SFS_plot.rds`, PDF and PNG |
| 9 | Isolation by distance and environment | see below | `RA_set1_IBD_IBE_results.RData`, `RA_set1_IBD_IBE_results.csv`, plots |
| 10 | Diversity statistics figure (neutral against adaptive) | reads `diversity.xlsx`, a table compiled from the section 1 and 2 results | `diversity_neutral_vs_adaptive.pdf/png` |
| 11 | Sampling map | `sf` and `rnaturalearth` | `RA_IndoPacific_zoomed_out.png/pdf` |

### Isolation by distance and environment (section 9)

Linearised pairwise FST (FST / (1 - FST)) was the response. Geographic distances are marine least-cost distances averaged to the population level, and environmental distances are Bray-Curtis, Manhattan or Canberra dissimilarities of the population mean values of the four environmental variables retained in the RDA (standardised), with the metric chosen by AIC (Canberra for three of the four datasets, so Canberra is used throughout). Three models (IBD only, IBE only, combined) were fitted with `gls()` from `nlme` by maximum likelihood to obtain AIC values for comparison. Significance of IBD and IBE was assessed with 10 000 permutations of the pairwise genetic values, with R² taken from the observed correlation. Pairwise comparisons share populations, so the p-values and AIC comparisons are indicative.

## Input

Objects from stage 2 (`RA_seascape_FINAL_checkpoint.RData`) and the stage 1 VCF (`RApop_set1.snps.vcf`). The diversity figure (section 10) also needs `diversity.xlsx`, a table compiled from the section 1 and 2 results that is not produced by code in this repository.
