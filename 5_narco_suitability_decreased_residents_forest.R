# This script maps the overlap between decreased narcotrafficking suitability and IBLs for resident birds

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

# read in important bird landscapes raster, match shapefile to
# this extent and resolution
prop_pop_forest_mig <- rast("data/final/ibl_forest-dependent_migrants_central-america.tif")

# read in shapefile with change in suitability for primary shipments
narco_suit_primary <- read_sf("data/raw/NarcoSuit/primv_chngcat.shp") %>%
  st_transform(crs) %>%
  select(suit_change = gridcode) %>%
  vect() 
# convert to raster
narco_suit_primary_rast <- terra::rasterize(narco_suit_primary, prop_pop_forest_mig,
                                            field = "suit_change") 

# read in shapefile with change in suitability for secondary movements
narco_suit_secondary <- read_sf("data/raw/NarcoSuit/secv_chngcat.shp") %>%
  st_transform(crs) %>%
  select(suit_change = gridcode) %>%
  vect()
# convert to raster
narco_suit_secondary_rast <- terra::rasterize(narco_suit_secondary, prop_pop_forest_mig,
                                              field = "suit_change")


# restrict both primary and secondary to areas where suitability DECREASED
# primary
suit_decreased_both <- narco_suit_primary_rast <= 0 & narco_suit_secondary_rast <= 0 &
  (narco_suit_primary_rast < 0 | narco_suit_secondary_rast < 0)
suit_decreased_both <- subst(suit_decreased_both, from = FALSE, to = NA)

# convert to polygons
pol_suit_decreased <- as.polygons(suit_decreased_both)

# turn back into multipolygon to intersect with important bird areas
multipol_suit_decreased <- st_as_sf(pol_suit_decreased) %>%
  st_intersection(central_america_proj)

# OVERLAY DECREASED SUIT AND RESIDENT FOREST IBLs
# read in migratory tipping-point important polygons
polygons_resident <- read_sf("data/final/important_polygons_residents.gpkg")

# intersect - note that sf is intelligent with attribute data!
overlap_decrease <- st_intersection(polygons_resident, multipol_suit_decreased) 
overlap_decrease %>%
  st_geometry() %>%
  write_sf("data/final/overlap_suitability_decrease_resident_forest_birds.gpkg")

# now, plot both areas of nt suitability increase and bird importance, and intersection of the layers
# plot -----
png("plots/4_overlap_ibls_resident_forest_birds_nt_decrease.png", width = 1500, height = 1500)
par(mar = c(2, 0.25, 4, 0.25))
# plot country borders
plot(central_america, lwd = 4)
# plot important bird areas
plot(st_geometry(polygons_resident), col = alpha("blue", alpha = 0.5), 
     border = NA, add = TRUE)
# plot areas where suitability for primary movements increase
plot(st_geometry(multipol_suit_decreased), col = alpha("darkorange", alpha = 0.5), border = NA, main = NULL,  add = TRUE)
# plot intersection of bird ibls and primary mvmt increase
plot(st_geometry(overlap_decrease), add = TRUE, col = 'darkred', border = NA)

# # add title
# title(main = "Footprint of decreased suitability for narco-trafficking",
#       cex.main = 4, line = -1)
# 
# legend(x = "bottom", legend=c("IBLs for resident forest birds", "Areas of decreased suitability for narco-trafficking",
#                               "Overlap between IBLs and decreased suitability"), 
#        fill = c(alpha("blue", alpha = 0.5), alpha("darkorange", alpha = 0.5), "darkred"),
#        cex = 3)

dev.off()


# GET AREA OF BOTH IMPORTANT POLYGONS FOR BIRDS and DECREASED SUIT - then % overlap
# add in area of intersections in km2
area_intersect <- overlap_decrease %>% 
  mutate(area = st_area(.) %>% units::set_units("km2") %>% as.numeric()) %>%
  as_tibble() %>% 
  summarize(area_overlap = sum(area))

# get total areas of important polygons for resident birds
area_fd_res_ibl <- polygons_resident %>% 
  mutate(area = st_area(.) %>% units::set_units("km2") %>% as.numeric()) %>% 
  as_tibble() %>% 
  summarize(area_total = sum(area))

# identify which bird group and shipment level these stats are for, then save to data/intermediate
area_stats <- area_intersect %>%
  cbind(area_fd_res_ibl) %>%
  mutate(percent_overlap = area_overlap/area_total) %>%
  mutate(id = "res_forest_birds_decrease") %>%
  write_csv("data/intermediate/forest-resident_stats_nt_decrease.csv")

### END SCRIPT ###
