
# This script identifies Important Bird Landscapes for resident forest birds
# in Central America

# install 2021 version of ebirdst package on which analysis was based
remotes::install_version("ebirdst", version = "2.2021.3")

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
library(RColorBrewer)
library(scales)


# get s&t raster template
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


# filter to resident forest birds
# get species data table from eBird
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
  select(common_name, acad_mig, breeding_habitat)
# subset to residents only
species_list_res <- species_list_info %>%
  filter(resident == T) %>%
  select(species_code, common_name)
# filter to species with Forest as primary breeding habitat
species_list_forest <- species_list_res %>%
  left_join(acad_data, by = "common_name") %>%
  filter(grepl('Forest', breeding_habitat))
species_list_res_forest <- species_list_forest %>%
  left_join(species_info_all) %>%
  filter(resident_quality > 1) %>%
  select(species_code, common_name, breeding_habitat)

# convert to character list
species_list_res <- paste(species_list_res_forest$species_code)

# make list to hold raster layers for all species
res_spp_rasters <- list()
# loops over all species
for (species in species_list_res) {
  message(species)
  # read in percent population raster
  prop_pop <- glue("{species}_percent-population_median_hr_2021.tif") %>% 
    path(raster_dir, .) %>% 
    rast() %>%
    # mask out Belize
    mask(central_america_vect)
  # get the max value for each cell across 52 weeks - will end up with single
  # layer instead of 52
  max_prop_pop <- app(prop_pop, "max", na.rm = TRUE)
  #  plot(max_prop_pop)
  # to threshold top 20% most important cells for each species, first get 
  # rid of zeros - want to prioritize WITHIN species' range
  max_prop_pop[max_prop_pop == 0] <- NA 
  # after removing zeros, get threshold value (20% of values above threshold)
  threshold <- global(max_prop_pop, quantile, probs = 0.8, na.rm = TRUE)
  # get numeric value
  threshold <- threshold[[1]][1]
  # turn all cells below threshold value in to zeros
  max_prop_pop[max_prop_pop < threshold] <- NA
  max_prop_pop_zeros_removed <- max_prop_pop
  #  plot(max_prop_pop_zeros_removed)
  # add single-species raster to list
  res_spp_rasters[[species]] <- max_prop_pop_zeros_removed
}

# create single SpatRaster with a layer for each species
res_spp_rasters_rast <- rast(res_spp_rasters)

# this gives for each cell an "importance" for birds by percent of the population across all species
prop_pop_res_spp <- app(res_spp_rasters_rast, "sum", na.rm = TRUE)
# read out
writeRaster(prop_pop_res_spp, "data/final/ibl_fd_residents_central-america.tif",
            overwrite = TRUE)

# Plot Important Bird Landscapes for residents

# first, drop the 0 values
prop_pop_res_spp_vector <- values(prop_pop_res_spp)
prop_pop_res_spp_vector <- prop_pop_res_spp_vector[prop_pop_res_spp_vector != 0]

breaks <- quantile(prop_pop_res_spp_vector, 
                   seq(0, 1, by = 0.2),
                   na.rm = TRUE)
pal <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#2b8cbe", "brown")


# plot -----
pdf("plots/ibl_prop_pop_forest_residents.pdf", height = 10, width = 10)
par(mar = c(1.25, 0.25, 0.25, 0.25))
plot(prop_pop_res_spp, 
     breaks = breaks, col = pal,
     legend = FALSE, axes = FALSE)
plot(central_america, col = NA, border = "black", lwd = 2, add = TRUE)

# add legend
par(new = TRUE, mar = c(0.5, 0, 0, 0))
title <- "Important areas for 129 Resident \nForest Birds in Central America"
image.plot(zlim = c(0, 1), legend.only = TRUE, 
           breaks = seq(0, 1, length.out = length(breaks)),
           col = pal,
           smallplot = c(0.20, 0.8, 0.05, 0.08),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), 
                            labels = c("low", "medium", "high"),
                            fg = "black", col.axis = "black",
                            cex.axis = 1.2, lwd.ticks = 0.5,
                            padj = 0.05),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1.5, line = 0.5,
                              padj = -0.1))
dev.off()

# ID important polygons for residents ----
# identify and map polygons showing areas of high importance for residents ----
threshold <- global(prop_pop_res_spp, quantile, probs = 0.7, na.rm = TRUE)
threshold <- threshold[[1]][1]
r <- prop_pop_res_spp > threshold
r <- subst(r, from = FALSE, to = NA)
pol <- as.polygons(r)
# includes many tiny polygons - drop small ones below certain area threshold
area_thresh <- units::set_units(50, km^2)
pol_dropped <- st_as_sf(pol) %>%
  drop_crumbs(threshold = area_thresh)
# fill in holes in polygons
pol_filled <- st_remove_holes(pol_dropped, as.numeric(units::set_units(area_thresh, "m2")))
# smooth edges (using Matt's package!)
pol_smoothed <- smooth(pol_filled, method = "ksmooth")
pol_smoothed_spatrast <- vect(pol_smoothed)

# plot these polygons with Central America country borders
# plot -----
pdf("plots/important_polygons_residents.pdf", width = 22, height = 22)
par(mar = c(0.25, 0.25, 5, 0.25))
# add polygons
pol_smoothed_spatrast_mask <- st_as_sf(pol_smoothed_spatrast) %>% 
  st_intersection(central_america_proj) %>%
  write_sf("data/final/important_polygons_residents.gpkg")
plot(central_america, col = NA, border = "black", lwd = 4)
plot(pol_smoothed_spatrast_mask, main = NULL, border = NA, col = alpha("blue", 0.5),
     add = TRUE)

# add title
title(main = "Important areas for 129 resident \nforest birds in Central America",
      cex.main = 3, line = -2)

dev.off()

### END SCRIPT ###

