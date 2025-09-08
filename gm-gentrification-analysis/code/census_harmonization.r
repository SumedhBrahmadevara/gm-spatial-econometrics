# ==============================================================================
# Census Data Harmonization: 1991 ED to 2011 MSOA
# ==============================================================================
# This script performs areal interpolation to harmonize 1991 Enumeration District
# data to 2011 Middle Super Output Area boundaries using area-weighted interpolation.
#
# Requirements:
# - 1991 ED shapefile
# - 1991 ED attribute data (CSV)
# - 2011 MSOA shapefile
# - Greater Manchester MSOA codes (CSV)
#
# Author: Sumedh Brahmadevara
# ==============================================================================

# Load required libraries
library(sf)
library(dplyr)
library(readr)

# ==============================================================================
# Configuration
# ==============================================================================

# Define file paths - update these paths to match your data location
DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- "path/to/your/output"

# Input file paths
ED_1991_SHAPEFILE <- file.path(DATA_DIR, "england_ed_1991.shp")
ED_1991_ATTRIBUTES <- file.path(DATA_DIR, "Scaled_1991_ED-Level_Dataset.csv")
MSOA_2011_SHAPEFILE <- file.path(DATA_DIR, "infuse_msoa_lyr_2011.shp")
GM_MSOA_CODES <- file.path(DATA_DIR, "gm_msoas.csv")

# Output file path
OUTPUT_FILE <- file.path(OUTPUT_DIR, "harmonized_1991_to_2011_census.csv")

# Target coordinate reference system (WGS84)
TARGET_CRS <- 4326

# ==============================================================================
# Load and prepare 1991 Enumeration District data
# ==============================================================================

message("Loading 1991 ED shapefile...")
ed_1991 <- st_read(ED_1991_SHAPEFILE)

message("Loading 1991 ED attribute data...")
ed_attributes <- read_csv(ED_1991_ATTRIBUTES, show_col_types = FALSE)

# Standardize join column name
colnames(ed_attributes)[1] <- "label"

# Join attributes to spatial data
ed_1991 <- left_join(ed_1991, ed_attributes, by = "label")

message(sprintf("Loaded %d enumeration districts", nrow(ed_1991)))

# ==============================================================================
# Load and prepare 2011 MSOA data
# ==============================================================================

message("Loading 2011 MSOA shapefile...")
msoa_2011 <- st_read(MSOA_2011_SHAPEFILE)

# Filter to Greater Manchester MSOAs only
message("Filtering to Greater Manchester MSOAs...")
gm_msoa_codes <- read_csv(GM_MSOA_CODES, show_col_types = FALSE)
msoa_2011 <- msoa_2011 %>% 
  filter(geo_code %in% gm_msoa_codes$msoa11cd)

message(sprintf("Filtered to %d Greater Manchester MSOAs", nrow(msoa_2011)))

# ==============================================================================
# Harmonize coordinate reference systems
# ==============================================================================

message("Harmonizing coordinate reference systems...")
ed_1991 <- st_transform(ed_1991, crs = TARGET_CRS)
msoa_2011 <- st_transform(msoa_2011, crs = TARGET_CRS)

# ==============================================================================
# Ensure valid geometries
# ==============================================================================

message("Validating geometries...")
ed_1991 <- st_make_valid(ed_1991)
msoa_2011 <- st_make_valid(msoa_2011)

# ==============================================================================
# Perform spatial intersection
# ==============================================================================

message("Computing spatial intersection...")
intersection <- st_intersection(ed_1991, msoa_2011)

# Ensure intersection geometries remain valid
intersection <- st_make_valid(intersection)

message(sprintf("Created %d intersection polygons", nrow(intersection)))

# ==============================================================================
# Calculate area weights for interpolation
# ==============================================================================

message("Computing area weights...")

# Calculate intersection areas
intersection$intersection_area <- st_area(intersection)

# Calculate original ED areas
intersection$ed_area <- st_area(ed_1991[match(intersection$label, ed_1991$label), ])

# Calculate area weight (proportion of ED area in each MSOA)
intersection$area_weight <- as.numeric(intersection$intersection_area / intersection$ed_area)

# ==============================================================================
# Define variables for interpolation
# ==============================================================================

# Population and demographic variables
demographic_vars <- c("pop", "age0_14", "age15_24", "age25_44", "age45_64", "age65plus")

# Education variables
education_vars <- c("total18", "educ_lev1", "educ_lev2", "educ_lev3", "educ_lev4")

# Occupation variables
occupation_vars <- c(
  "occ_total", "occ_managers_seniors", "occ_professional", "occ_associate",
  "occ_admin", "occ_skilled", "occ_personal", "occ_sales", "occ_machine", "occ_other"
)

# Employment sector variables
employment_vars <- c(
  "primary", "secondary", "tertiary", "quaternary", "total_emp",
  "primary_share", "secondary_share", "tertiary_share", "quaternary_share"
)

# Housing variables
housing_vars <- c(
  "dwelling_total", "hh_total", "hh_1person_65plus", "hh_1person_other",
  "hh_loneparent", "hh_married_nochildren", "hh_married_withchildren", "hh_other",
  "total_households"
)

# Housing tenure variables
tenure_vars <- c(
  "owned_outright", "owned_mortgage", "owned_total",
  "private_landlord", "private_other", "private_rented_total",
  "social_rented_other", "social_rented_council", "social_rented_total"
)

# Housing type variables
housing_type_vars <- c(
  "type_detach", "type_semi", "type_terra", "type_flat", 
  "type_conv", "type_other", "type_temp"
)

# Combine all variables
interpolation_vars <- c(
  demographic_vars, education_vars, occupation_vars, 
  employment_vars, housing_vars, tenure_vars, housing_type_vars
)

# ==============================================================================
# Validate variables and perform interpolation
# ==============================================================================

# Check which variables exist in the data
existing_vars <- interpolation_vars[interpolation_vars %in% colnames(intersection)]
missing_vars <- setdiff(interpolation_vars, existing_vars)

if (length(missing_vars) > 0) {
  warning("The following variables are missing and will be skipped: ",
          paste(missing_vars, collapse = ", "))
}

message(sprintf("Interpolating %d variables using area weights...", length(existing_vars)))

# Apply area weights to create weighted variables
for (var in existing_vars) {
  weighted_var_name <- paste0(var, "_weighted")
  intersection[[weighted_var_name]] <- intersection[[var]] * intersection$area_weight
}

# ==============================================================================
# Aggregate to MSOA level
# ==============================================================================

message("Aggregating to MSOA level...")

# Create weighted variable names
weighted_var_names <- paste0(existing_vars, "_weighted")

# Aggregate weighted values by MSOA
harmonized_data <- intersection %>%
  st_drop_geometry() %>%
  group_by(geo_code) %>%
  summarise(
    across(all_of(weighted_var_names), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Rename weighted variables back to original names
names(harmonized_data) <- gsub("_weighted$", "", names(harmonized_data))

message(sprintf("Harmonized data for %d MSOAs", nrow(harmonized_data)))

# ==============================================================================
# Export results
# ==============================================================================

message(sprintf("Exporting results to: %s", OUTPUT_FILE))

# Ensure output directory exists
dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

# Export to CSV
write_csv(harmonized_data, OUTPUT_FILE)

message("Data harmonization complete!")

# ==============================================================================
# Summary statistics
# ==============================================================================

message("\n--- Summary ---")
message(sprintf("Input EDs: %d", nrow(ed_1991)))
message(sprintf("Output MSOAs: %d", nrow(harmonized_data)))
message(sprintf("Variables interpolated: %d", length(existing_vars)))
message(sprintf("Total intersection polygons: %d", nrow(intersection)))

if (length(missing_vars) > 0) {
  message(sprintf("Variables skipped (missing): %d", length(missing_vars)))
}
