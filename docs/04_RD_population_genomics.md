# 04. *R. djiddensis*: population genomics

Population structure, differentiation, diversity and isolation by distance in the whitespotted wedgefish (*Rhynchobatus djiddensis*, RD), using dataset A from [stage 1](01_assembly_and_filtering.md): 47 individuals and 27,366 SNPs. All steps are run in R. The marine distances used in section 6 come from a short separate script, [`scripts/04a_RD_spatial_data.R`](../scripts/04a_RD_spatial_data.R), which is run first; everything else is in [`scripts/04_RD_population_genomics.R`](../scripts/04_RD_population_genomics.R).

## Populations

Individuals were grouped as SA (South Africa: Sodwana Bay SB, Zinkwazi ZK and Durban DU; n = 29), MO (Mozambique; n = 8) and AR (Arabian region; n = 10, merged from sub-locations with one or two individuals each). The South African sites are analysed separately in the isolation by distance and divMigrate sections.

## Sections

| Section | Analysis | Main settings | Main outputs |
|---|---|---|---|
| 0 | Marine distances (`04a_RD_spatial_data.R`): marine least-cost distances between individuals with `marineDistances`, a PCoA and dbMEMs | resolution 0.5, 5 degree buffer; GSHHS coastline | `RD_set1_spatial_data.RData` (`marine_dist_RD`), `Contour _ Pairwise Marine Distances.txt` |
| 1 | Load VCF and popmap, build the genotype matrix | missing genotypes replaced by the most common genotype at each locus for DAPC and sNMF; FST uses the non-imputed genind | |
| 2 | Pairwise and global FST, heatmap | Weir and Cockerham (`hierfstat`); 1000 bootstraps; a pair is significant when the lower bound of the bootstrap 95% CI is above zero (p set to 0.001, otherwise 1), then Benjamini-Hochberg correction | `RD_set1_pairwiseFST.csv`, `RD_set1_FST_heatmap.pdf/png` |
| 3 | DAPC | PCs explaining 80% of variance; `find.clusters` for K = 1 to 10 (10 replicates); `xvalDapc`; `optim.a.score` (100 simulations); with and without population priors; final DAPC with 5 PCs (a-score) | `RD_set1_DAPC_summary.csv`, `RD_set1_dapc_noprior.pdf`, `RD_set1_dapc_5pc.pdf/png` |
| 4 | sNMF | `LEA`; K = 1 to 10, 10 replicates, alpha = 100, K chosen by cross-entropy; barplots with `pophelper` | `RD_set1_snmf_crossentropy.pdf` |
| 5 | Diversity: HO, HE, FIS, allelic richness, private alleles | `hierfstat`; 95% confidence intervals for FIS from 500 bootstraps of loci; richness rarefied to the smallest population; private alleles with `poppr` | `RD_set1_diversity.csv`, `RD_set1_diversity_results.RData` |
| 6 | Isolation by distance | five populations (SB, ZK, DU, MO, AR), ten pairs; see below | `RD_set1_5pop_pairwiseFST.csv`, `RD_set1_5pop_IBD_results.csv`, `RD_set1_5pop_IBD.pdf/png` |
| 7 | Unmerged FST and divMigrate (SB, ZK, DU, MO; Arabian region excluded) | pairwise FST and significance as in section 2; divMigrate with Jost's D, GST and Nm (see below) | `RD_set1_noRS_pairwiseFST.csv`, `RD_set1_noRS_divMigrate_results.RData`, `RD_set1_noRS_divMigrate_network_labelled.pdf/png` |
| 8 | Sampling map | `sf`, `rnaturalearth` | `RD_main_basemap.png/pdf`, `RD_inset_map.png/pdf` |

### Isolation by distance (section 6)

Linearised pairwise FST (FST / (1 - FST)) is the response and marine least-cost distance, averaged to the population level, is the predictor. The model is fitted with `gls()` from `nlme` by maximum likelihood, without a correlation structure, and significance comes from 10,000 permutations of the pairwise genetic values, with R² taken from the observed correlation. With five populations there are ten pairs that share populations and four of them involve the Arabian region, so the p-value is indicative.

### divMigrate (section 7)

Directional relative migration (Jost's D, GST and Nm; Sundqvist et al. 2016) among the South African and Mozambican sites (SB, ZK, DU, MO; the Arabian region is excluded) was estimated with `diveRsity` v1.9.90 (Keenan et al. 2013), as in Manuzzi et al. (2019). Bootstrapping fails in this version of the package, so the analysis is split in two:

1. **On the cluster:** the Arabian region individuals are removed, and a random subset of 5,000 SNPs is drawn.

```bash
vcftools --vcf RDpop_set1.snps.vcf --remove-indv <each Arabian region individual> \
  --recode --recode-INFO-all --out RD_set1_noRS

grep -v "^#" RD_set1_noRS.recode.vcf | shuf -n 5000 | sort -k1,1 -k2,2n > snps_5000.txt
awk '{print $1, $2}' snps_5000.txt > positions_5000.txt

vcftools --vcf RD_set1_noRS.recode.vcf --positions positions_5000.txt \
  --recode --recode-INFO-all --out RD_set1_noRS_5000
```

   `shuf` is not seeded, so the 5,000-SNP subset cannot be reproduced exactly.
2. **Convert to Genepop** on Galaxy Europe (usegalaxy.eu), for the full dataset and the subset. Galaxy writes all individuals into one population section, so both files are restructured in R by population.
3. **Relative migration** was calculated in R on the full dataset (27,366 SNPs) without bootstraps.
4. **Significance** of asymmetric migration was tested in the [divMigrate online application](https://popgen.shinyapps.io/divMigrate-online/) on the 5,000-SNP subset (500 bootstraps, all three statistics, alpha = 0.05).
5. **The network plot** was drawn in R from the 5,000-SNP subset. Its node labels show the first individual of each population and were corrected to population names manually after export.

References: Keenan et al. (2013) *Methods in Ecology and Evolution* 4:782-788; Sundqvist et al. (2016) *Ecology and Evolution* 6:3461-3475; Manuzzi et al. (2019) *Molecular Ecology* 28:2509-2526.

## Input

`RDpop_set1.snps.vcf` (stage 1), `RDpopmap.txt` (two columns: ID, population) and, for divMigrate, the two Genepop files from Galaxy Europe (`genepop.txt` for the full dataset and `genepop_5000.txt` for the subset). For the marine distances (`04a_RD_spatial_data.R`): `RD_locations.xlsx` (ID, location, latitude, longitude per individual), the GSHHS coastline shapefile `GSHHS_h_L1.shp` and the `marineDistances` script (Klein et al. 2024).
