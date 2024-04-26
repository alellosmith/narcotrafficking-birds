# This script makes the map of change in suitability for narcotrafficking

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

# get lakes to mask out
ne_lakes <- ne_download(scale = 10, category = "physical",
                        type = "lakes", returnclass = "sf")
ne_lakes_proj <- st_transform(ne_lakes, crs) %>%
  # crop to Central America
  st_crop(central_america) %>%
  st_intersection(central_america) %>%
  st_geometry()

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
                                              field = "suit_change") %>%
  as.factor()
cls <- data.frame(id = -1:1, change=c("decrease", "none", "increase"))
levels(narco_suit_primary_rast) <- cls

# read in shapefile with change in suitability for secondary movements
narco_suit_secondary <- read_sf("data/raw/NarcoSuit/secv_chngcat.shp") %>%
  st_transform(crs) %>%
  select(suit_change = gridcode) %>%
  vect()
# convert to raster
narco_suit_secondary_rast <- terra::rasterize(narco_suit_secondary, prop_pop_forest_mig,
                                            field = "suit_change") %>%
  as.factor()
cls <- data.frame(id = -1:1, change=c("decrease", "none", "increase"))
levels(narco_suit_secondary_rast) <- cls

plot(narco_suit_secondary_rast)
plot(narco_suit_primary_rast)

# concatenate layers together to get all possible combos of increase, no change, or decrease for each
# kind of narco suitability (primary and secondary)
narco_suit_all_change <- concats(narco_suit_primary_rast, narco_suit_secondary_rast)
# plot combined raster
plot(narco_suit_all_change)
levels(narco_suit_all_change)[[1]]
# want to use same color for multiple categories; create dataframe assigning proper colors
# to each change class
coltb = data.frame(value = c(0, 2:6, 8), col = c("gold", "skyblue3", "gold",
                                  "lightgray", "skyblue3", "lightskyblue1",
                                  "royalblue4"))
coltab(narco_suit_all_change) <- coltb

# plot narco suit map with all categories
pdf("plots/changes_narco_suit_5cats.pdf", width = 22, height = 22)
par(mar = c(0.25, 0.25, 6, 0.25))
plot(narco_suit_all_change, legend = FALSE, axes = FALSE)
plot(central_america, lwd = 3, add = TRUE)
plot(ne_lakes_proj, add = TRUE, col = "white")

# add title
title(main = "Predicted change in suitability for primary and secondary \n cocaine shipments after peak interdiction pressure",
      cex.main = 3, line = -4)

# add legend
legend("bottomleft",
       legend = c("primary & secondary increased", "secondary increased",
                  "primary increased", "no change", "suitability decreased"),
       col = c("royalblue4", "skyblue3", "lightskyblue1", "lightgray", "gold"),
       pch = c(19),
       bty = "n",
       pt.cex = 4,
       cex = 3,
       text.col = "black",
       horiz = F ,
       inset = c(0.01, 0.01))

# save plot
dev.off()

### END SCRIPT ###