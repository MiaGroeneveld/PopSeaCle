# PopSeaCle
**Pop**ulation and **sea**scape genetics **c**ol**le**ction - *Rhynchobatus* 3RADseq data analysis

Scripts and workflow notes for the RADseq analyses of the whitespotted wedgefish (*Rhynchobatus djiddensis*) and the bottlenose wedgefish (*R. australiae*) in the Western Indian Ocean and Indo-West Pacific.

> Companion code for [citation to be added on publication].

## Pipeline

| Stage | Contents | Where it runs | Document |
|---|---|---|---|
| 1 | Raw data QC, demultiplexing, de novo assembly (Stacks), filtering | HPC, with some checks in R | [`docs/01_assembly_and_filtering.md`](docs/01_assembly_and_filtering.md) |
| 2 | *R. australiae*: environmental data, variable selection, partial RDA, outlier detection, four-dataset preparation | R, with ordistep and partial RDA on HPC | [`docs/02_RA_seascape_outliers.md`](docs/02_RA_seascape_outliers.md), [`scripts/02_RA_seascape_outliers.R`](scripts/02_RA_seascape_outliers.R) |
| 3 | *R. australiae*: population genomics (diversity, FST, DAPC, sNMF, IBD/IBE) | R | [`docs/03_RA_population_genomics.md`](docs/03_RA_population_genomics.md), [`scripts/03_RA_population_genomics.R`](scripts/03_RA_population_genomics.R) |
| 4 | *R. djiddensis*: population genomics | R | to be added |
| 5 | *R. djiddensis*: hybridisation analyses | R | to be added |
| 6 | *R. djiddensis*: multiple paternity and twinning | R, with simulations on HPC | to be added |

Each document marks which steps ran on a cluster and which in R.

## Repository layout

```
.
├── README.md
├── docs/
│   ├── 01_assembly_and_filtering.md
│   ├── 02_RA_seascape_outliers.md
│   └── 03_RA_population_genomics.md
└── scripts/
    ├── 02_RA_seascape_outliers.R
    ├── 03_RA_population_genomics.R
    └── hpc/
        ├── RA_set1_RDA_ordistep11.R
        └── RA_set1_partialRDA.R
```

Data files (VCFs, popmaps, coordinates, environmental data) are not included.

## Software

| Software | Version |
|---|---|
| FastQC | 0.11.5 |
| Cutadapt | 1.9.1 |
| Stacks | 2.59 |
| VCFtools | 0.1.14 |
| RADstackshelpR | 0.1.0 |
| R | 4.4.1 locally (packages listed in each stage); 4.3.2 for the two cluster RDA jobs |

## Data availability

Raw sequence data: [accession to be added]. Sample metadata and locality information are not included in this repository.

## Contact

[Name, email or ORCID]

