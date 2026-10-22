# PopSeaCle
**Pop**ulation and **sea**scape genetics **c**ol**le**ction

Scripts and workflow notes for the 3RADseq analyses of the whitespotted wedgefish (*Rhynchobatus djiddensis*) and the bottlenose wedgefish (*R. australiae*) in the Western Indian Ocean and Indo-West Pacific.

> Companion code for [citations to be added once published].

## Pipeline

| Stage | Contents | Where it runs | Document |
|---|---|---|---|
| 1 | Raw data QC, demultiplexing, de novo assembly (Stacks), filtering | HPC, with some checks in R | [`docs/01_assembly_and_filtering.md`](docs/01_assembly_and_filtering.md) |
| 2 | *R. australiae*: environmental data, variable selection, partial RDA, outlier detection, four-dataset preparation | R, with ordistep and partial RDA on HPC | [`docs/02_RA_seascape_outliers.md`](docs/02_RA_seascape_outliers.md), [`scripts/02_RA_seascape_outliers.R`](scripts/02_RA_seascape_outliers.R) |
| 3 | *R. australiae*: population genomics (diversity, FST, DAPC, sNMF, IBD/IBE) | R | [`docs/03_RA_population_genomics.md`](docs/03_RA_population_genomics.md), [`scripts/03_RA_population_genomics.R`](scripts/03_RA_population_genomics.R) |
| 4 | *R. djiddensis*: population genomics (FST, DAPC, sNMF, diversity, IBD, divMigrate) | R | [`docs/04_RD_population_genomics.md`](docs/04_RD_population_genomics.md), [`scripts/04_RD_population_genomics.R`](scripts/04_RD_population_genomics.R), [`scripts/04a_RD_spatial_data.R`](scripts/04a_RD_spatial_data.R) |
| 5 | *R. djiddensis*: hybridisation (PCA, diagnostic loci, sNMF, NewHybrids) | HPC, then R | [`docs/05_hybridisation.md`](docs/05_hybridisation.md), [`scripts/05_hybridisation.R`](scripts/05_hybridisation.R) |
| 6 | *R. djiddensis*: multiple paternity and twinning (COLONY input, kinship, relatedness simulation, PISibs, IBD) | R, with simulations on HPC | [`docs/06_multiple_paternity.md`](docs/06_multiple_paternity.md), [`scripts/06_multiple_paternity.R`](scripts/06_multiple_paternity.R), [`scripts/06a_MP_related_input.R`](scripts/06a_MP_related_input.R) |

## Repository layout

```
.
├── README.md
├── docs/
│   ├── 01_assembly_and_filtering.md
│   ├── 02_RA_seascape_outliers.md
│   ├── 03_RA_population_genomics.md
│   ├── 04_RD_population_genomics.md
│   ├── 05_hybridisation.md
│   └── 06_multiple_paternity.md
└── scripts/
    ├── 02_RA_seascape_outliers.R
    ├── 03_RA_population_genomics.R
    ├── 04_RD_population_genomics.R
    ├── 04a_RD_spatial_data.R
    ├── 05_hybridisation.R
    ├── 06_multiple_paternity.R
    ├── 06a_MP_related_input.R
    └── hpc/
        ├── RA_set1_RDA_ordistep11.R
        ├── RA_set1_partialRDA.R
        ├── MP_related_simulation.R
        └── MP_related_empirical.R
```

Data files (VCFs, popmaps, coordinates, environmental data) are not included.
R version 4.4.1 used locally (packages are loaded at the top of each script); 4.3.2 used for the cluster R jobs (partial RDA, relatedness)

## Running the pipeline

- ***R. australiae*:** stage 1, then stage 2 (with the two cluster jobs in `scripts/hpc/`), then stage 3.
- ***R. djiddensis*:** stage 1, then `04a_RD_spatial_data.R` and stage 4. Stages 5 and 6 start from the cluster steps described in their documents.
- Each document lists the input files its scripts read under "Input". The R scripts are written to be run section by section, in order, from the top of the file.

## Data availability

Raw sequence data: [accession to be added]. Sample metadata and locality information are not included in this repository.

## Contact

ResearchGate: https://www.researchgate.net/profile/Mia-Groeneveld.

