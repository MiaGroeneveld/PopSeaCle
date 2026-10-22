# 05. *R. djiddensis*: hybridisation

Detection of hybrids between the whitespotted wedgefish (*Rhynchobatus djiddensis*, RD) and two congeners, using two hybridisation datasets built from the same sequencing run:

- **Mozambique** (RD x *R. australiae*, RA): the putative hybrid MO3. A second individual, MZ16, was suspected as a hybrid from its spotting pattern.
- **Arabian region** (RD x *R. laevis*, RL): the putative hybrid RS2.

Each dataset has a version with and without the putative hybrid. The versions without it define the diagnostic loci from the parental species, and the versions with it are used to classify the individual. The analysis in R is in [`scripts/05_hybridisation.R`](../scripts/05_hybridisation.R).

## Datasets

| Dataset | Individuals | SNPs | File |
|---|---|---|---|
| Mozambique, without MO3 | 35 (RA 28, RD 7) | 36 612 | `hyb_noMO3.vcf` |
| Mozambique, with MO3 | 36 (RA 28, RD 8) | 36 317 | `hyb_withMO3.vcf` |
| Arabian region, without RS2 | 13 (RD 10, RL 3) | 40 089 | `rs_noRS2.vcf` |
| Arabian region, with RS2 | 14 (RD 11, RL 3) | 39 551 | `rs_withRS2.vcf` |

## Cluster steps

The steps follow stage 1 and are the same for both datasets ([HPC]):

1. **Assembly** (`denovo_map.pl`, m = 3, M = 3, n = 4, as in stage 1). The Mozambique datasets come from **one** joint assembly of the RA and RD individuals, and the individuals to exclude (MO5, and MO3 for the "without" version) are left out with a popmap in the next step. The two Arabian datasets were assembled **separately**, with and without RS2.
2. **Base filter** with `populations` and the popmap for each version: loci present in at least 80% of samples, maximum observed heterozygosity 0.5, minor allele count of at least 3, one SNP per locus.
3. **Genotype filters** with VCFtools: `--minGQ 20 --minDP 5 --max-missing 0.80`.
4. **Whitelist** of the loci that remain, then a final `populations` run on the assembly with the whitelist (`--write-single-snp --vcf --structure`).

```bash
populations --in-path $ASSEMBLY --out-path $OUT/base_R80_Het50_MAC3/ --popmap $POPMAP \
  --min-samples-overall 0.80 --max-obs-het 0.5 --min-mac 3 --write-single-snp --vcf

vcftools --vcf $OUT/base_R80_Het50_MAC3/populations.snps.vcf --minGQ 20 --minDP 5 --max-missing 0.80 \
  --recode --recode-INFO-all --out $OUT/hyb_GQ20_DP5_MISS0.80

grep -v '^#' $OUT/hyb_GQ20_DP5_MISS0.80.recode.vcf | cut -f1 | sort -n | uniq > whitelist.txt

populations --in-path $ASSEMBLY --out-path $FINAL --popmap $POPMAP \
  --whitelist whitelist.txt --write-single-snp --vcf --structure
```

## Analyses in R

| Step | What it does | Main settings |
|---|---|---|
| PCA | of each dataset including the putative hybrid (the script also draws the version without it, as a species-separation check) | `pcadapt`, K = 2 |
| Diagnostic loci | per-locus FST between the two parental species in the dataset without the putative hybrid; loci with FST = 1 are kept (FST > 0.9 if there are fewer than 200); the loci are then found in the dataset with the putative hybrid, and loci with any missing genotype are removed | `hierfstat::basic.stats` |
| sNMF | ancestry coefficients on the diagnostic loci | `LEA`; K = 1 to 5 (Mozambique) or 1 to 4 (Arabian region), 100 repetitions, seed 123; K = 2 |
| NewHybrids | random subsets of 200 diagnostic loci for three seeds (123, 456, 789), run with burn-in 10000 and 100000 sweeps | `parallelnewhybrid::parallelnh_WIN` |

NewHybrids runs outside R. It uses the Windows executable (the `newhybrids` program) through the `parallelnewhybrid` package, with one folder per run. NewHybrids is meant for about 200 loci, so random 200-locus subsets were used instead. For the Arabian region, NewHybrids did not converge because only three *R. laevis* reference individuals are available, so RS2 is classified from the PCA, the sNMF ancestry coefficients and its heterozygosity at the diagnostic loci. Heterozygosity is the number of heterozygous genotypes among the 200 loci of each subset (`het_MO3` and `het_RS2` in the script).

## Input

The four VCFs above (stage 1 and the cluster steps), the NewHybrids executable, and, for the NewHybrids runs, folders of the input files written by the script.
