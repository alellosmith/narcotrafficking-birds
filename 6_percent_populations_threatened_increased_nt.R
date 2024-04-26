
# This script estimates the % of each species' population that is at risk from increased
# narcotrafficking suitability

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
library(tidyr)
library(smoothr)
library(readr)
library(scales)
library(nngeo)
library(RColorBrewer)

# get country borders
crs <- get_species_path("example_data") %>% 
  dir_ls(glob = "*.tif", recurse = TRUE) %>% 
  head(1) %>% 
  rast() %>% 
  st_crs()

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
central_america_vect <- vect(central_america)


# set first part of path
raster_dir <- "data/raw/ebirdst_high-res/"

# read in list of species with eBird abundance data in Central America
# have to do it this complicated way because two species didn't have 
# high enough of a % population in the higher-resolution (smaller) cells
species_list <- dir_ls(raster_dir, 
                       glob = "*percent-population_median_hr_2021.tif") %>% 
  basename() %>% 
  str_extract("^[a-z0-9]+")


# select migratory species only
# get data table from eBird
species_info_all <- ebirdst_runs |> as.data.frame()
species_list_df <- as.data.frame(species_list) %>%
  rename(species_code = species_list)
# subset data table to just species with abundance maps in Central America
species_list_info <- species_list_df %>%
  inner_join(species_info_all, by = "species_code") 
# join to ACAD data to get breeding habitat
acad_data <- read_csv("data/raw/acad_global-scores_CA-2021.csv") %>%
  rename(common_name = "Common Name", acad_mig = "Mig Status",
         breeding_habitat = "Primary Breeding Habitat", nonbreeding_habitat = "Primary Nonbreeding Habitat") %>%
  select(common_name, acad_mig, nonbreeding_habitat)
# subset to migrants only
species_list_mig <- species_list_info %>%
  filter(resident == F) %>%
  select(species_code, common_name)
# filter to species with Forest as primary breeding habitat 
species_list_forest <- species_list_mig %>%
  left_join(acad_data, by = "common_name") %>%
  filter(grepl('Forest', nonbreeding_habitat))
# according to the ACAD this list of 261 species contains a few residents,
# a bunch of partial migrants, and "true" (NN) migrants
# subsetting JUST to M migrants
species_list_forest_mig <- species_list_forest %>%
  filter(acad_mig == "M") %>%
  left_join(species_info_all) %>%
  filter(nonbreeding_quality > 1) %>%
  select(species_code, common_name, nonbreeding_habitat)

# convert to character list
species_list_fd_mig <- paste(species_list_forest_mig$species_code)

# read in polygons with increased primary and secondary suitability
increased_nt_all <- read_sf("data/final/areas_primary_secondary_increased_nt.gpkg")
# make list to hold % of global pop values for all species
mig_spp_total_percents <- list()
# loops over all species
for (species in species_list_fd_mig) {
  message(species)
  # read in percent population raster
  prop_pop <- glue("{species}_percent-population_median_hr_2021.tif") %>% 
    path(raster_dir, .) %>% 
    rast() %>%
    # mask to area impacted by increased suitability for narco-trafficking
    mask(increased_nt_all)
  # for each week, sum the % population values across the area impacted by narco-trafficking
  # this gives the percent of the global population that is in that area in each week
  weekly_prop_pop <- global(prop_pop, "sum", na.rm = TRUE)
  # get the week with the highest % population value
  max_sum_pop <- max(weekly_prop_pop, na.rm = TRUE)
  # add single-species value to list
  mig_spp_total_percents[[species]] <- max_sum_pop
}

mig_forest_total_percents <- as_tibble(mig_spp_total_percents,
                                       .name_repair = "unique") %>%
  pivot_longer(everything(), names_to = "species_code") %>%
  mutate(value = round(value, digits = 2)) %>%
  rename(perc_pop = "value") 

mig_percents_table <- mig_forest_total_percents %>%
  mutate(perc_pop = perc_pop * 100) %>%
  inner_join(species_list_info) %>%
  select(common_name, scientific_name, perc_pop) %>%
  rename("English Name" = common_name, "Scientific Name" = scientific_name,
         "% Population" = perc_pop) %>%
  write_csv("data/final/perc_pop_migrants_threatened_incr_suitability.csv")

### END SCRIPT ###