library(ebirdst)
library(terra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(tibble)
library(stringr)
library(glue)
library(fs)
library(fields)
library(smoothr)
library(readr)
library(nngeo)
library(scales)

# get correct projection (from eBird status)
crs <- get_species_path("example_data") %>% 
  dir_ls(glob = "*.tif", recurse = TRUE) %>% 
  head(1) %>% 
  rast() %>% 
  st_crs()

# get country borders
central_america <- ne_countries(scale = 10,
                                continent = "North America", 
                                returnclass = "sf") %>% 
  filter(subregion == "Central America") %>% 
  st_transform(crs = crs) %>% 
  select(country_code = iso_a2, country = name) %>%
  filter(country != "Mexico") %>% 
  filter(country != "Belize") %>% 
  st_geometry()
# project
central_america_proj <- st_transform(central_america, crs)
# vectorize
central_america_vect <- vect(central_america_proj)

# MIGRATORY FOREST BIRDS ----
# read in polygons showing areas of overlap for migratory landbirds ----
# primary suitability
overlap_mig_forest <- read_sf("data/final/overlap_nt_suitability_mig_forest_birds.gpkg") %>%
  vect()

# plot with country borders
png("plots/overlap_all_nt_increase_mig_forest-birds.png", width = 1500, height = 1500)
par(mar = c(0.25, 0.25, 4, 0.25))
plot(central_america)
plot(overlap_mig_forest, axes = FALSE, border = NA, col = alpha("blue", 0.5),
     add = TRUE)

# add title
title(main = "Overlap between migratory forest bird IBLs and \nincreased suitability for narco-trafficking",
      cex.main = 3, line = -2)

dev.off()


# save combined polygons for migratory forest birds 
writeVector(overlap_mig_forest, "data/final/combined_overlap_mig_forest_birds_nt.gpkg")

