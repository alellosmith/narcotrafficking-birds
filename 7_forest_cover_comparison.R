library(ebirdst)
library(terra)
library(sf)
library(rnaturalearth)
library(tibble)
library(dplyr)
library(stringr)
library(glue)
library(fs)
library(exactextractr)
library(fields)
library(smoothr)
library(readr)
library(scales)  
library(nngeo)
library(ggplot2)
library(tidyr)
library(landscapemetrics)

# read in netCDF file (global map of landcover classes 2018) using terra
lc_2018_global <- terra::rast("data/raw/C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.nc")
lc_2018_global

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
# get crs of central america vector
crs_ca <- crs(central_america_vect)
# project central am vect to lc_2018_global, bc other way around takes too long (global too big)
central_america_vect_ll <- project(central_america_vect, lc_2018_global)

# crop giant global landcover map to just central america
lc_2018_ca <- subset(lc_2018_global, "lccs_class") %>% 
  crop(central_america_vect_ll) %>% 
  project(crs_ca, method = "near")

# crop to CA borders
lc_2018_ca_masked <- mask(lc_2018_ca, central_america_vect)
 # writeRaster("data/intermediate/landcover_ca_2018.tif")

# reclassify to forest vs nonforest
# read in csv to manually reclassify
rcl_mat <- read_csv("data/raw/reclass_forest_als.csv")

lc_2018_ca_reclassed <- classify(lc_2018_ca_masked, rcl = rcl_mat) 
lc_2018_ca_forest <- ifel(lc_2018_ca_reclassed == 0, NA, lc_2018_ca_reclassed) %>%
  writeRaster("data/final/forest_cover_ca.tif", overwrite = TRUE)

plot(lc_2018_ca_forest, col = "black", legend = FALSE, axes = FALSE)



# get polygons for areas where ALL narcotrafficking has become (1) more suitable and (2) less suitable

# less suitable

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

# more suitable

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

# convert to SpatRaster to write to file
narco_suit_combined_rast <- rasterize(narco_suit_combined, suit_decreased_both) 


# turn back into multipolygon to intersect with important bird areas
narco_suit_combined <- st_as_sf(narco_suit_combined) %>%
  st_intersection(central_america_proj)


# mask forest cover layer to areas where NT decreased
lc_2018_nt_decreased <- mask(lc_2018_ca_reclassed, multipol_suit_decreased)
# plot
plot(lc_2018_nt_decreased, col = c("blue", "orange"))

# mask forest cover layer to areas where NT increased
lc_2018_nt_increased <- mask(lc_2018_ca_reclassed, narco_suit_combined)
# plot
plot(lc_2018_nt_increased, col = c("blue", "orange"))

# GET CORE AREA METRICS FOR areas where NT DECREASED ----

# METRICS ----
# number of forest patches in areas of decreased NT suit
patch_metrics_decreased <- dplyr::bind_rows(
  lsm_p_area(lc_2018_nt_decreased),
  lsm_p_cai(lc_2018_nt_decreased),
  lsm_p_para(lc_2018_nt_decreased)
) 

patch_metrics_decreased_forest <- patch_metrics_decreased %>%
  filter(class == 1)

patch_area_forest_decreased <- patch_metrics_decreased_forest %>%
  filter(metric == "area") %>%
  mutate(landscape = "Decreased") 

forest_patch_area_top20_decreased <- patch_area_forest_decreased %>%
  arrange(desc(value)) %>% 
  slice(1:20)

mean(forest_patch_area_top20_decreased$value)
sd(forest_patch_area_top20_decreased$value)

# number of forest patches in areas of increased NT suit
patch_metrics_increased <- dplyr::bind_rows(
  lsm_p_area(lc_2018_nt_increased),
  lsm_p_cai(lc_2018_nt_increased),
  lsm_p_para(lc_2018_nt_increased)
) 

patch_metrics_increased_forest <- patch_metrics_increased %>%
  filter(class == 1)

patch_area_forest_increased <- patch_metrics_increased_forest %>%
  filter(metric == "area") %>%
  mutate(landscape = "Increased") 

forest_patch_area_top20_increased <- patch_area_forest_increased %>%
  arrange(desc(value)) %>% 
  slice(1:20) 

mean(forest_patch_area_top20_increased$value)
sd(forest_patch_area_top20_increased$value)

patch_area_forest_both <- patch_area_forest_increased %>%
  rbind(patch_area_forest_decreased)


# plot density of area
ggplot(patch_area_forest_both, aes(x = log(value), fill = landscape)) + 
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("gold", "skyblue3")) +
  labs(title = "Distribution of forest patch size\nin landscapes of increased and decreased\nsuitability for narco-trafficking",
       x = "log(Area(ha))", y = "Density",
       fill = "Change in\nnarco-trafficking\nsuitability") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

top20_forest_both <- forest_patch_area_top20_increased %>%
  rbind(forest_patch_area_top20_decreased)

top20_forest_both <- top20_forest_both %>%
  group_by(landscape) %>%
  mutate(mean = mean(value),
         sd = sd(value))
  

# make function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

g <- ggplot(top20_forest_both, aes(x = landscape, y = (value/1000000), fill = landscape), color = "black") + 
  geom_violin() +
  scale_fill_manual(values = c("gold", "skyblue3")) +
  stat_summary(fun.data = data_summary,
               colour = "black") +
  labs(title = "Areas of 20 largest forest patches", x = "Change in narco-trafficking suitability", y = "Area (1,000,000 ha)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = "none")
ggsave(g, filename = "plots/forest_patch_size_comparison.pdf",
       width = 10, height = 9) 
      
  

landscape_metrics_decreased <- dplyr::bind_rows(
  lsm_l_core_mn(lc_2018_nt_decreased),
  lsm_l_core_sd(lc_2018_nt_decreased)
)

landscape_metrics_increased <- dplyr::bind_rows(
  lsm_l_core_mn(lc_2018_nt_increased),
  lsm_l_core_sd(lc_2018_nt_increased)
)



# mean core area and SD core area
mean_core_decreased <- lsm_l_core_mn(lc_2018_nt_decreased,
              directions = 8,
              consider_boundary = FALSE,
              edge_depth = 2
)

sd_core_decreased <- lsm_l_core_sd(lc_2018_nt_decreased,
              directions = 8,
              consider_boundary = FALSE,
              edge_depth = 2
)




# calculate % forest cover within areas where NT INcreased
# first, mask forest cover raster to areas where NT increased
lc_2018_nt_increased <- mask(lc_2018_ca_reclassed, narco_suit_combined)
# plot
plot(lc_2018_nt_increased)

forest_increase_perc <- global(lc_2018_nt_increased, fun = "mean", na.rm = TRUE)
# 63.2% of area of increased suitability for NT forested in 2019

# calculate % forest cover within areas where NT DEcreased
forest_decrease_perc <- global(lc_2018_nt_decreased, fun = "mean", na.rm = TRUE)
# 57.7% of area of decreased suitability for NT

# area of increased suit 5.5% more forested than area of decreased suit
central_america_crop <- st_bbox(
  c(xmin = -9928529.53568959, ymin = 803249.84054262, 
    xmax = -8497933.31263572, ymax = 1981053.2246607)) %>% 
  st_crop(central_america, .)
# plot forest cover in areas of decreased NT
png("forest_cover_narco-decrease.png", width = 1500, height = 1500)
plot(central_america_crop, 
     border = NA, col = "grey90",
     axes = FALSE, bty = "n", reset = FALSE)
plot(lc_2018_nt_decreased, 
     breaks = c(-0.5, 0.5, 1.5), col = c("orange", "forestgreen"),
     axes = FALSE, legend = FALSE, add = TRUE)
plot(central_america, col = NA, border = "black", lwd = 4, add = TRUE)
dev.off()

# plot forest cover in areas of increased NT
png("forest_cover_narco-increase.png", width = 1500, height = 1500)
plot(central_america_crop, 
     border = NA, col = "grey90",
     axes = FALSE, bty = "n", reset = FALSE)
plot(lc_2018_nt_increased, 
     breaks = c(-0.5, 0.5, 1.5), col = c("orange", "forestgreen"),
     axes = FALSE, legend = FALSE, add = TRUE)
plot(central_america_crop, col = NA, border = "black", lwd = 4, add = TRUE)
dev.off()


  




