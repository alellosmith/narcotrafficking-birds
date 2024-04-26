
# This script maps the overlap between increased narcotrafficking suitability and IBLs for resident birds

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

# restrict both primary and secondary to areas where suitability INCREASED
# primary
narco_suit_prim_increased <- narco_suit_primary_rast > 0
narco_suit_prim_increased <- subst(narco_suit_prim_increased, from = FALSE, to = NA) 
# convert to polygons
pol_prim_increase <- as.polygons(narco_suit_prim_increased)
pol_sf_prim <- st_as_sf(pol_prim_increase) %>%
  # convert to SpatVect
  vect()
# secondary
narco_suit_sec_increased <- narco_suit_secondary_rast > 0
narco_suit_sec_increased <- subst(narco_suit_sec_increased, from = FALSE, to = NA)
# convert to polygons
pol_sec_increase <- as.polygons(narco_suit_sec_increased)
pol_sf_sec <- st_as_sf(pol_sec_increase) %>%
  # convert to SpatVect
  vect()

# combine two SpatVectors of polygons (for primary and secondary) into one
narco_suit_combined <- combineGeoms(pol_sf_prim, 
                                    pol_sf_sec,
                                    distance = FALSE)

# turn back into multipolygon to intersect with important bird areas
narco_suit_combined <- st_as_sf(narco_suit_combined) %>%
  st_intersection(central_america_proj)


# OVERLAY IMPORTANT BIRD LANDSCAPES AND INCREASED SUITABILITY FOR NT ----

# read in forest resident important polygons
polygons_resident <- read_sf("data/final/important_polygons_residents.gpkg")

# intersect - note that sf is intelligent with attribute data!
overlap_increase <- st_intersection(polygons_resident, narco_suit_combined) 

overlap_increase %>%
  st_geometry() %>%
  write_sf("data/final/overlap_suitability_forest-residents.gpkg")

# now, plot both areas of nt suitability increase and bird importance, and intersection of the layers
# plot -----
png("plots/3_overlap_ibls_forest-residents_nt_increase.png", width = 1500, height = 1500)
par(mar = c(2, 0.25, 4, 0.25))
# plot country borders
plot(central_america, lwd = 4)
# plot important bird areas
plot(st_geometry(polygons_resident), col = alpha("blue", alpha = 0.5), 
     border = NA, add = TRUE)
# plot areas where suitability for primary movements increase
plot(st_geometry(narco_suit_combined), col = alpha("darkorange", alpha = 0.5), border = NA, main = NULL,  add = TRUE)
# plot intersection of bird ibls and primary mvmt increase
plot(st_geometry(overlap_increase), add = TRUE, col = 'darkred', border = NA)

# # add title
# title(main = "Footprint of increased suitability for narco-trafficking",
#       cex.main = 4, line = -1)
# 
# legend(x = "bottom", legend=c("IBLs for resident forest species", "Areas of increased suitability for narco-trafficking",
#                               "Overlap between IBLs and increased suitability"), 
#        fill = c(alpha("blue", alpha = 0.5), alpha("darkorange", alpha = 0.5), "darkred"),
#        cex = 3)

dev.off()


# GET AREA OF BOTH IMPORTANT POLYGONS FOR BIRDS and INCREASED SUIT - then % overlap
# add in area of intersections in km2
area_intersect <- overlap_increase %>% 
  mutate(area = st_area(.) %>% units::set_units("km2") %>% as.numeric()) %>%
  as_tibble() %>% 
  summarize(area_overlap = sum(area))

# get total areas of important polygons for migratory birds
area_res_ibl <- polygons_resident %>% 
  mutate(area = st_area(.) %>% units::set_units("km2") %>% as.numeric()) %>% 
  as_tibble() %>% 
  summarize(area_total = sum(area))

# identify which bird group and shipment level these stats are for, then save to data/intermediate
area_stats <- area_intersect %>%
  cbind(area_res_ibl) %>%
  mutate(percent_overlap = area_overlap/area_total) %>%
  mutate(id = "forest_resident") %>%
  write_csv("data/intermediate/forest-resident_stats.csv")

### END SCRIPT ###