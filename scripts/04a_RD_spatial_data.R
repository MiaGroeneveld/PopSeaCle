# ============================================================
# RD marine distances (input to the isolation by distance analysis in 04_RD_population_genomics.R)
# Species: Rhynchobatus djiddensis (RD set1: 47 individuals, 27 366 SNPs)
# Marine least-cost distances between individuals with marineDistances,
# averaged to the population level; same procedure as for R. australiae (stage 2, section 5)
# ============================================================

# ── Step 1: Setup and data loading ───────────────────────────────────────────

library(vcfR)
library(adegenet)
library(vegan)
library(readxl)

# ── 1. Load popmap ────────────────────────────────────────────────────────────
popmap <- read.table("RDpopmap.txt",
                     header = FALSE,
                     col.names = c("ID", "pop"))
cat("Popmap loaded:", nrow(popmap), "individuals\n")
cat("Populations:", levels(factor(popmap$pop)), "\n")
table(popmap$pop)

# ── 2. Load coordinates ───────────────────────────────────────────────────────
coords <- read_excel("RD_locations.xlsx")
colnames(coords) <- c("ID", "Location", "lat", "lon")

# Match to popmap order and assign population
coords <- coords[match(popmap$ID, coords$ID), ]
coords$pop <- popmap$pop[match(coords$ID, popmap$ID)]

cat("\nCoordinates loaded:", nrow(coords), "individuals\n")
cat("Any missing coords:\n")
print(coords[is.na(coords$lat) | is.na(coords$lon), "ID"])

# ── 3. Load VCF and build genind ──────────────────────────────────────────────
vcf       <- read.vcfR("RDpop_set1.snps.vcf", verbose = FALSE)
gind_full <- vcfR2genind(vcf)

# Fix locus names: "CHR:POS" (0-based) → "CHR_POS" (1-based)
loc_fixed <- gsub(":", "_", locNames(gind_full))
loc_split <- strsplit(loc_fixed, "_")
loc_fixed <- sapply(loc_split, function(x) paste0(x[1], "_", as.integer(x[2]) + 1))
locNames(gind_full) <- loc_fixed

# Assign populations
pop_factor <- popmap$pop[match(indNames(gind_full), popmap$ID)]
pop(gind_full) <- factor(pop_factor)

cat("\nGenind built:", nInd(gind_full), "individuals |",
    nLoc(gind_full), "loci\n")
cat("Populations:", levels(pop(gind_full)), "\n")
cat("Pop sizes:\n")
print(table(pop(gind_full)))

# ── 4. Load genotype matrix for RDA/LFMM ─────────────────────────────────────
# Read VCF dosage matrix (0/1/2)
geno_mat <- t(extract.gt(vcf, element = "GT",
                          as.numeric = FALSE))

# Convert GT to dosage
gt_to_dosage <- function(gt) {
  sapply(gt, function(x) {
    if (is.na(x) || x == "./.") return(NA)
    alleles <- as.integer(strsplit(x, "[/|]")[[1]])
    sum(alleles)
  })
}

G <- apply(geno_mat, c(1,2), gt_to_dosage)
rownames(G) <- indNames(gind_full)
colnames(G) <- locNames(gind_full)

cat("\nGenotype matrix:", nrow(G), "individuals x", ncol(G), "loci\n")
cat("Missing data rate:", round(mean(is.na(G)) * 100, 2), "%\n")

# ── 5. Impute missing data (mean imputation per locus) ────────────────────────
G_imp <- G
for (j in 1:ncol(G_imp)) {
  miss <- is.na(G_imp[, j])
  if (any(miss)) {
    G_imp[miss, j] <- mean(G_imp[, j], na.rm = TRUE)
  }
}
cat("Imputation complete — missing data remaining:",
    sum(is.na(G_imp)), "\n")

# ── 6. Verify individual order is consistent ──────────────────────────────────
stopifnot(all(rownames(G_imp) == popmap$ID))
stopifnot(all(coords$ID == popmap$ID))
cat("\nAll individual orders consistent ✓\n")

# ── 7. Save session objects ───────────────────────────────────────────────────
save(popmap, coords, gind_full, G, G_imp,
     file = "RD_set1_session_objects.RData")
cat("Session objects saved\n")


# ============================================================
# RD Seascape Genomics — Step 3: Marine distances + dbMEMs
# Using marineDistances
# Requires GSHHS_h_L1.shp in working directory
# ============================================================

library(adespatial)

# ── 1. Write coordinates file ─────────────────────────────────────────────────
write.table(
  data.frame(
    Name = coords$ID,
    Lon  = coords$lon,
    Lat  = coords$lat
  ),
  file      = "RD_individuals.txt",
  sep       = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

cat("Individual coordinates file written\n")
head(read.table("RD_individuals.txt"))

# ── 2. Run marineDistances ────────────────────────────────────────────────────
# Requires GSHHS_h_L1.shp + companions in working directory
# Download from: http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip

source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")

contour(
  global.polygon = "GSHHS_h_L1.shp",
  file           = "RD_individuals.txt",
  file.sep       = " ",
  file.dec       = ".",
  file.strucutre = 1,
  file.header    = FALSE,
  resolution     = 0.5,
  buffer         = c(5, 5, 5, 5),
  export.file    = TRUE
)

# ── 3. Load distance matrix ───────────────────────────────────────────────────
marine_dist_RD <- as.matrix(
  read.table("Contour _ Pairwise Marine Distances.txt",
             header = TRUE, row.names = 1,
             sep = " ", dec = ".")
)

cat("Marine distance matrix dimensions:", dim(marine_dist_RD), "\n")

# ── 4. Sanity checks — update pairs as appropriate for RD ────────────────────
cat("ZK1 to MO1 (expect ~800-1000 km):", marine_dist_RD["ZK1", "MO1"], "\n")
cat("SB1t to ZK1 (expect ~200-300 km):", marine_dist_RD["SB1t", "ZK1"], "\n")
cat("MO1 to RS1 (expect ~5000-7000 km):", marine_dist_RD["MO1", "RS1"], "\n")
cat("SB1t to RS1 (expect ~7000-9000 km):", marine_dist_RD["SB1t", "RS1"], "\n")

# ── 5. PCoA to Euclidean ──────────────────────────────────────────────────────
marine_pcoa_RD <- cmdscale(as.dist(marine_dist_RD),
                            k = nrow(marine_dist_RD) - 1,
                            eig = TRUE)

n_pos <- sum(marine_pcoa_RD$eig > 0)
cat("Positive eigenvalues:", n_pos, "\n")

marine_euclid_RD <- marine_pcoa_RD$points[, 1:n_pos]
rownames(marine_euclid_RD) <- coords$ID

marine_euclid_dist_RD <- dist(marine_euclid_RD)

# ── 6. dbMEMs ─────────────────────────────────────────────────────────────────
dbmems_RD <- dbmem(marine_euclid_dist_RD)
cat("dbMEMs computed:", ncol(dbmems_RD), "variables\n")
cat("MEM eigenvalues (first 5):\n")
print(attributes(dbmems_RD)$values[1:5])

# ── 7. Population-level distances for IBD/IBE ────────────────────────────────
pop_levels_RD <- unique(popmap$pop)

pop_coords_RD <- aggregate(cbind(lon, lat) ~ pop,
                            data = coords, FUN = mean)
rownames(pop_coords_RD) <- pop_coords_RD$pop
pop_coords_RD <- pop_coords_RD[match(pop_levels_RD,
                                      pop_coords_RD$pop), ]

# Extract population-level distances from individual matrix
# by averaging within population pairs
n_pops <- length(pop_levels_RD)
geo_dist_pop_RD <- matrix(NA, n_pops, n_pops,
                           dimnames = list(pop_levels_RD, pop_levels_RD))

for (i in 1:n_pops) {
  for (j in 1:n_pops) {
    idx_i <- which(popmap$pop == pop_levels_RD[i])
    idx_j <- which(popmap$pop == pop_levels_RD[j])
    geo_dist_pop_RD[i, j] <- mean(marine_dist_RD[idx_i, idx_j])
  }
}
diag(geo_dist_pop_RD) <- 0

cat("\nPopulation-level distances (km):\n")
print(round(geo_dist_pop_RD, 0))

# ── 8. Alignment checks ───────────────────────────────────────────────────────
cat("\n── Alignment checks ──\n")
cat("G_imp dimensions:", dim(G_imp), "\n")
cat("dbMEMs dimensions:", dim(dbmems_RD), "\n")
cat("Row order consistent:", all(rownames(G_imp) == coords$ID), "\n")

# ── 9. Save ───────────────────────────────────────────────────────────────────
save(marine_dist_RD, marine_euclid_dist_RD,
     dbmems_RD, geo_dist_pop_RD, pop_coords_RD,
     file = "RD_set1_spatial_data.RData")

cat("Spatial data saved: RD_set1_spatial_data.RData\n")
