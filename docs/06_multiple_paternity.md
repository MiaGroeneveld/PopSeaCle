# 06. *R. djiddensis*: multiple paternity and twinning

Analysis of a litter of eleven pups (SC1 to SC11) and their deceased mother (SC12) of the whitespotted wedgefish (*Rhynchobatus djiddensis*) from Mozambique: sibship reconstruction and number of sires (COLONY), pairwise kinship and relatedness, and identity-by-descent estimates for putative monozygotic twins.

## Dataset

Dataset D contains the 12 family individuals and a SWIO reference panel of 35 individuals (47 in total) and **1 378 SNPs with no missing data**. The family individuals were extracted from it for the sub-analyses (`fam_only.vcf`, the same 1 378 SNPs; see the cluster steps). The reference panel provides allele frequencies for the relatedness analyses.

### Cluster steps [HPC]

Dataset D comes from the same de novo assembly as dataset A ([stage 1](01_assembly_and_filtering.md)).

```bash
# combined popmap (family + reference panel, 48 individuals): loci present in all individuals
populations --in-path $ASSEMBLY --out-path $OUT --popmap combined_popmap.txt \
  --min-samples-overall 1.0 --max-obs-het 0.5 --min-mac 3 --write-single-snp --vcf

# genotype filters (2 523 SNPs before, 1 378 after)
vcftools --vcf $OUT/populations.snps.vcf --minGQ 20 --minDP 5 --max-missing 1 \
  --recode --recode-INFO-all --out MP_combined_GQ20_DP5_MISS1.0

# remove a PD1 (48 to 47 individuals)
vcftools --vcf MP_combined_GQ20_DP5_MISS1.0.recode.vcf --remove-indv PD1 \
  --recode --recode-INFO-all --out MP_combined_noPD1

# the 12 family individuals extracted from the filtered VCF (same 1 378 SNPs, no missing data);
# the local copy of this file is fam_only.vcf
vcftools --vcf MP_combined_GQ20_DP5_MISS1.0.recode.vcf --keep family_individuals.txt \
  --recode --recode-INFO-all --out MP_family_1378

# whitelist of the loci in the family VCF
grep -v '^#' MP_family_1378.recode.vcf | awk '{print $1}' | sort | uniq > family_1378_whitelist.txt

# STRUCTURE file for COLONY: the 12 family individuals, restricted to those loci
populations --in-path $ASSEMBLY --out-path $FAMILY_OUT --popmap family_individuals_popmap.txt \
  --whitelist family_1378_whitelist.txt --structure
```

## Scripts

| Script | What it does | Where |
|---|---|---|
| [`scripts/06a_MP_related_input.R`](../scripts/06a_MP_related_input.R) | converts dataset D to the `related` input format and extracts the reference panel allele frequencies | R |
| [`scripts/hpc/MP_related_simulation.R`](../scripts/hpc/MP_related_simulation.R) | relatedness simulation | HPC |
| [`scripts/hpc/MP_related_empirical.R`](../scripts/hpc/MP_related_empirical.R) | empirical relatedness (LynchRD, TrioML, DyadML) | HPC |
| [`scripts/06_multiple_paternity.R`](../scripts/06_multiple_paternity.R) | COLONY input files, heatmaps, relatedness plots, PISibs and IBD | R |

## Sections of `06_multiple_paternity.R`

| Section | What it does | Main settings |
|---|---|---|
| COLONY input | reads the STRUCTURE file (`skip = 2` for the two header lines of the Stacks output), and writes the offspring genotypes, the mother genotype, the marker file with the error rate and the known maternal sibship, with and without SC3 | genotype error rate 0.01, allelic dropout 0 |
| Kinship heatmaps | pairwise genetic distance (`adegenet::bitwise.dist`, pup-subset) and KING-robust kinship (`SNPRelate::snpgdsIBDKING` - which is calculated in R on all 47 individuals, so that the reference panel contributes to the allele frequencies, and then subset to the pups) | `MP_heatmaps_combined.pdf/png` |
| Relatedness | plots of the empirical LynchRD, TrioML and DyadML estimates from the cluster, with the simulated distributions | |
| PISibs | probability of identity for siblings (Waits et al. 2001), per locus and cumulative, for the family individuals | |
| IBD | `SNPRelate` method of moments (`snpgdsIBDMoM`) and maximum likelihood (`snpgdsIBDMLE`) for SC1 to SC12 | `SC1_to_SC12_IBD_matrix.csv`, `SC1_to_SC12_IBD_MLE_results.csv` |

### COLONY

COLONY (full-likelihood method) is run in its graphical interface with the input files written by the script. A second run excludes SC3.

### Relatedness simulation

The allele frequencies of the 35 reference individuals replace the frequencies in the combined input, and 400 loci are sampled at random (the memory limit of the Fortran backend of `related`). Then 150 dyads per class (parent-offspring, full siblings, half siblings, unrelated) are simulated and scored with seven estimators (TrioML, Wang, LynchLI, LynchRD, Ritland, QuellerGT and DyadML), with a genotype error rate of 0.01 and `trioml.num.reference = 100`. The empirical estimates are run on all 1 378 loci.

## Input

`MP_combined.vcf` (dataset D, from the cluster steps), `fam_only.vcf` (the family individuals extracted), `MP_family_1378.structure` (the Stacks STRUCTURE output renamed), and the `related` input files written by `06a_MP_related_input.R`.
