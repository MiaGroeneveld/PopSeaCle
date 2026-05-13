# Population and Seascape Genomics of the whitespotted (*Rhynchobatus djiddensis*) and bottlenose (*R. australiae*) wedgefishes

## Overview

This file contains QC, de novo assembly, and filtering scripts used to process raw 3RADseq data generated using the NovaSeq S4 high-throughput sequencing platform with paired-end (PE) sequencing (2 × 150 bp). PBS and R scripts are referenced in [brackets]. The pipeline covers both species sequenced on the same run.

**Final datasets:**

| Dataset | Species | SNPs | n | Missingness | Notes |
|---|---|---|---|---|---|
| RD_set1 | *R. djiddensis* | 27,366 | 47 | 4.83% | MAC≥3, FS1 — **PRIMARY** |
| RD_set2 | *R. djiddensis* | 18,942 | 47 | 5.09% | MAC≥3, FS2 (DP8) — robustness |
| RA_set1 | *R. australiae* | 27,562 | 71 | 6.32% | MAC≥3, FS1 — **PRIMARY** |
| RA_set2 | *R. australiae* | 30,380 | 71 | 6.24% | MAC≥2, FS1 — robustness |

---

## Step 1: Preparing Raw Data

### 1.1 Download raw data via filesender onto the HPC

```bash
# Interactive job required (password protected)
qsub -I -l select=1:ncpus=2:mem=4gb -l walltime=20:00:00
cd $HOME/radseq
wget --user=maduna-radlibs11 --ask-password --accept "*.tar" \
  --recursive --no-directories --no-parent \
  https://hts-nonsecure.web.sigma2.no/260220_LH00770.A.Project_Maduna-RADlibs11-2026-02-06/

# Unzip and verify integrity
tar -xvf 260220_LH00770.A.Project_Maduna-RADlibs11-2026-02-06.tar
md5sum -c md5sum.txt
```

---

## *Rhynchobatus djiddensis* (RD)

### 1.2 Quality control of raw datasets [mergeA.sh] [fastqcA.sh]

```bash
# Merge reads — two dual-indexing combinations per plate
# R1
cat Rdji-3RAD-D2-Plate10-1_S3_L001_R1_001.fastq.gz \
    Rdji-3RAD-D2-Plate10-2_S4_L001_R1_001.fastq.gz \
  > Rdji-3RAD-D2-Plate10_L001_R1_001.fastq.gz

# R2
cat Rdji-3RAD-D2-Plate10-1_S3_L001_R2_001.fastq.gz \
    Rdji-3RAD-D2-Plate10-2_S4_L001_R2_001.fastq.gz \
  > Rdji-3RAD-D2-Plate10_L001_R2_001.fastq.gz

# FastQC
module load app/fastqc/0.11.5
fastqc Rdji-3RAD-D2-Plate10_L001_R1_001.fastq.gz
fastqc Rdji-3RAD-D2-Plate10_L001_R2_001.fastq.gz
```

---

## Step 2: Demultiplexing

```bash
# Create output directories:
# /new-home/mjgroeneveld/radseq/RD_demultiplex/
# /new-home/mjgroeneveld/radseq/RD_demultiplex_A/
# /new-home/mjgroeneveld/radseq/RD_demultiplex_T/
```

### 2.1 Adapter trimming with Cutadapt [cutadaptA.sh]

```bash
module load app/cutadapt/1.9.1

cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --pair-filter=any \
  -o RD.R1.trimmed.fastq.gz \
  -p RD.R2.trimmed.fastq.gz \
  Rdji-3RAD-D2-Plate10_L001_R1_001.fastq.gz \
  Rdji-3RAD-D2-Plate10_L001_R2_001.fastq.gz
```

### 2.2 Demultiplexing with process_radtags — three strategies

```bash
module load app/stacks/2.59

# (a) Untrimmed [process_radtagsA.sh]
/apps/stacks/2.59/bin/process_radtags \
  -1 Rdji-3RAD-D2-Plate10_L001_R1_001.fastq.gz \
  -2 Rdji-3RAD-D2-Plate10_L001_R2_001.fastq.gz \
  -i gzfastq \
  -b intags_d2_RD.txt --inline_inline \
  -o /new-home/mjgroeneveld/radseq/RD_demultiplex/reads/ \
  -c -q -r -t 140 -w 0.15 -s 10 \
  --renz_1 mspI --renz_2 bamHI

# (b) Trimmed (Cutadapt) [process_radtagsA_T.sh]
/apps/stacks/2.59/bin/process_radtags \
  -1 RD.R1.trimmed.fastq.gz \
  -2 RD.R2.trimmed.fastq.gz \
  -i gzfastq \
  -b intags_d2_RD.txt --inline_inline \
  -o /new-home/mjgroeneveld/radseq/RD_demultiplex_T/reads/ \
  -c -q -r -t 140 -w 0.15 -s 10 \
  --renz_1 mspI --renz_2 bamHI

# (c) Untrimmed with adapters specified — SELECTED [process_radtagsA2.sh]
/apps/stacks/2.59/bin/process_radtags \
  -1 Rdji-3RAD-D2-Plate10_L001_R1_001.fastq.gz \
  -2 Rdji-3RAD-D2-Plate10_L001_R2_001.fastq.gz \
  -i gzfastq \
  -b intags_d2_RD.txt --inline_inline \
  -o /new-home/mjgroeneveld/radseq/RD_demultiplex_A/reads/ \
  --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --adapter_mm 2 \
  -c -q -r -t 140 -w 0.15 -s 10 \
  --renz_1 mspI --renz_2 bamHI
```

### 2.3 Compare trimming strategies [extract_tsv_AT2.sh]

```bash
stacks-dist-extract process_radtagsA.log  per_barcode_raw_read_counts > persamplecounts_A.tsv
stacks-dist-extract process_radtagsA_T.log per_barcode_raw_read_counts > persamplecounts_AT.tsv
stacks-dist-extract process_radtagsA2.log per_barcode_raw_read_counts > persamplecounts_A2.tsv

# Compare in R [trimming_compare_three.txt]
# Decision: A2 selected — no cutadapt, adapters specified in process_radtags
```

---

## Step 3: Parameter Optimisation using RADstackshelpR

```bash
# Select 20 samples for optimisation (~equal representation per locality,
# close to average read depth)
mkdir rem_files && mv *.rem.*.fq.gz rem_files/

# Optimise m [radstackhelpr_optim_m.sh]
WORK=/new-home/mjgroeneveld/radseq/RD_assembly_A/optim/
READS_DIR=/new-home/mjgroeneveld/radseq/RD_demultiplex_A/reads/
POPMAP=$WORK/popmap_optim.txt

for m in {3..7}; do
  OUT=$WORK/denovo.m${m}
  mkdir -p $OUT
  /apps/stacks/2.59/bin/denovo_map.pl \
    --samples $READS_DIR --popmap $POPMAP \
    --out-path $OUT --paired \
    -m $m -T 12 \
    -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
done

# Optimise M [radstackshelpr_optim_bigM.sh]
for i in {1..8}; do
  OUT=$WORK/denovo.bigM${i}
  mkdir -p $OUT
  /apps/stacks/2.59/bin/denovo_map.pl \
    --samples $READS_DIR --popmap $POPMAP \
    --out-path $OUT --paired \
    -m 3 -M $i -T 12 \
    -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
done

# Optimise n [radstackshelpr_optim_n.sh]
for n in {1..8}; do
  OUT=$WORK/denovo.n${n}
  mkdir -p $OUT
  /apps/stacks/2.59/bin/denovo_map.pl \
    --samples $READS_DIR --popmap $POPMAP \
    --out-path $OUT --paired \
    -m 3 -M 3 -n $n -T 12 \
    -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
done

# Visualise in R [radstackshelpr_visualisation.txt]
# Optimised parameters: m=3, M=3, n=4
```

---

## Step 4: STACKS De Novo Assembly [denovo_assembly.sh]

```bash
module load app/stacks/2.59
WORK=/new-home/mjgroeneveld/radseq/RD_assembly_A/core/
READS_DIR=/new-home/mjgroeneveld/radseq/RD_demultiplex_A/reads/
POPMAP=$WORK/popmap_full.txt
OUT=$WORK/denovo_m3_M3_n4
mkdir -p $OUT

/apps/stacks/2.59/bin/denovo_map.pl \
  --samples $READS_DIR --popmap $POPMAP \
  --out-path $OUT --paired \
  -m 3 -M 3 -n 4 -T 12 \
  -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
```

### 4.1 Evaluate assembly effectiveness

```bash
# Coverage after ustacks
stacks-dist-extract denovo_map.log cov_per_sample | grep -v '^#'

# Coverage after gstacks
cat gstacks.log | grep -B 2 -A 3 '^Genotyped'
# RD: Genotyped 450,000+ loci
# RA: Genotyped 515,960 loci
#   effective per-sample coverage: mean=75.4x, stdev=36.3x, min=6.2x, max=140.8x
#   mean number of sites per locus: 267.8
#   consistent phasing: 70.2% of diploid loci

stacks-dist-extract gstacks.log.distribs effective_coverages_per_sample \
  | grep -v '^#' | cut -f 1-2,4-5,8
# Flag problematic samples in R [denovo_checks.txt]
```

---

## Step 5: Locus and Individual-Based Filtering

### 5.1 Base filtering with populations module [base_filter.sh]

```bash
# R80/Het50/MAC3 base filter (repeat for R60 and MAC2 variants)
populations \
  --in-path /new-home/mjgroeneveld/radseq/RD_assembly_A/core/denovo_m3_M3_n4/ \
  --out-path /new-home/mjgroeneveld/radseq/RD_assembly_A/core/filtering/base_R80_Het50_MAC3/ \
  --popmap popmap_onepop.txt \
  --min-samples-overall 0.80 \
  --max-obs-het 0.5 \
  --mac 3 \
  --vcf
```

### 5.2 VCFtools filter sets [filter_sets.sh]

```bash
# GQ = minimum genotype quality, DP = minimum depth, max-missing = site missingness
# Four filter sets tested:
run_filter_set FS1 20 5 0.8   # GQ20, DP5, MISS0.8
run_filter_set FS2 20 8 0.8   # GQ20, DP8, MISS0.8
run_filter_set FS3 20 8 0.9   # GQ20, DP8, MISS0.9
run_filter_set FS4 25 10 0.9  # GQ25, DP10, MISS0.9

# Evaluate in R [filter_sets_R.txt]
# Decision: FS1 selected for primary datasets (best balance of missingness and SNP retention)
# FS2 (DP8) used for RD_set2 robustness dataset
```

### 5.3 Remove replicates and low-quality individuals

```bash
# Remove replicates only [remove1.sh]
vcftools \
  --vcf FS1_GQ20_DP5_MISS0.8.recode.vcf \
  --remove remove1.txt \
  --recode --recode-INFO-all \
  --out FS1_R80_Het50_MAC3_GQ20_DP5_MISS0.8_remove1

# Remove replicates + flagged individuals (MZ13, TZ5t for RA) [remove2.sh]
vcftools \
  --vcf FS1_GQ20_DP5_MISS0.8.recode.vcf \
  --remove remove2.txt \
  --recode --recode-INFO-all \
  --out FS1_R80_Het50_MAC3_GQ20_DP5_MISS0.8_remove2
```

### 5.4 Final populations run with whitelist [final_clean1.sh] [final_clean2.sh]

```bash
# Write one SNP per locus, create VCF and additional formats
# in-path points to original assembly, not filtered VCF
populations \
  --in-path /new-home/mjgroeneveld/radseq/RD_assembly_A/core/denovo_m3_M3_n4 \
  --out-path /new-home/mjgroeneveld/radseq/RD_assembly_A/core/final/popgen1 \
  --popmap popmap_final1.txt \
  --whitelist RD1_whitelist.txt \
  --write-single-snp \
  --vcf --structure --radpainter
```

### 5.5 Final QC statistics

```bash
# Missingness per individual
module load app/vcftools
vcftools --vcf populations.snps.vcf --missing-indv
awk 'NR>1 {sum+=$5} END {print sum/(NR-1)}' out.imiss

# SNP count
grep -vc '^#' populations.snps.vcf > popgen_nSNPs.txt
```

---

## *Rhynchobatus australiae* (RA)

### 1.2 Quality control of raw datasets

```bash
# Merge reads — two dual-indexing combinations per plate
# R1
cat Raus-3RAD-D2-Plate11-1_S5_L001_R1_001.fastq.gz \
    Raus-3RAD-D2-Plate11-2_S6_L001_R1_001.fastq.gz \
  > Raus-3RAD-D2-Plate11_L001_R1_001.fastq.gz

# R2
cat Raus-3RAD-D2-Plate11-1_S5_L001_R2_001.fastq.gz \
    Raus-3RAD-D2-Plate11-2_S6_L001_R2_001.fastq.gz \
  > Raus-3RAD-D2-Plate11_L001_R2_001.fastq.gz

# FastQC
module load app/fastqc/0.11.5
fastqc Raus-3RAD-D2-Plate11_L001_R1_001.fastq.gz
fastqc Raus-3RAD-D2-Plate11_L001_R2_001.fastq.gz
```

### Steps 2–3: Demultiplexing and Parameter Optimisation

```bash
# Two strategies compared [process_radtagsBT.sh] [process_radtagsB2.sh]
# Decision: B2 selected — no cutadapt, adapters specified in process_radtags

# Optimised parameters: m=3, M=3, n=4 (same as RD)
```

### Step 4: Assembly [denovo_assembly.sh]

```bash
# Remove non-target species individuals before assembly
mv TZ1.1.fq.gz TZ4.1.fq.gz TZ8.1.fq.gz TZ9.1.fq.gz TZ11.1.fq.gz TZ12.1.fq.gz nonspecies/
mv TZ1.2.fq.gz TZ4.2.fq.gz TZ8.2.fq.gz TZ9.2.fq.gz TZ11.2.fq.gz TZ12.2.fq.gz nonspecies/

WORK=/new-home/mjgroeneveld/radseq/RA_assembly_B/core/
READS_DIR=/new-home/mjgroeneveld/radseq/RA_demultiplex_B/reads/
POPMAP=$WORK/popmap_full.txt
OUT=$WORK/denovo_m3_M3_n4
mkdir -p $OUT

/apps/stacks/2.59/bin/denovo_map.pl \
  --samples $READS_DIR --popmap $POPMAP \
  --out-path $OUT --paired \
  -m 3 -M 3 -n 4 -T 12 \
  -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
```

### Steps 5: Filtering

```bash
# Same filtering pipeline as RD, excep MAC≥2 to increase SNP retention
# FS1 selected for RA_set1 (MAC≥3) and RA_set2 (MAC≥2)

# Final populations run — one SNP per locus
populations \
  --in-path /new-home/mjgroeneveld/radseq/RA_assembly_B/core/denovo_m3_M3_n4 \
  --out-path /new-home/mjgroeneveld/radseq/RA_assembly_B/core/final/popgen1 \
  --popmap popmap_final1.txt \
  --whitelist RA1_whitelist.txt \
  --write-single-snp \
  --vcf --structure --radpainter
```

---
