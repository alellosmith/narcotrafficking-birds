# This script generates the migratory bird connections map (Figure 2 in paper)
# need to get API from eBird "download data" page

library(tidyverse)
library(glue)
library(fs)
library(viridis)
library(terra)
library(sf)
library(rnaturalearth)
library(ebirdst)
library(ebirdstwf)
library(foreach)
library(doParallel)
registerDoParallel(cores = 15)

group <- "mig-forest-birds"
group_name <- "Migratory Forest Birds"
region_boundary <- glue("overlap_nt_suitability_{group}.gpkg") %>% 
  path("data", .) %>% 
  read_sf() %>% 
  st_combine()
region_boundary_vect <- vect(region_boundary)


# focal region ----

hawaii <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  filter(postal == "HI") %>% 
  st_buffer(10000)
wh <- st_bbox(c(xmin = -170, xmax = 0, ymin = 0, ymax = 90), crs = 4326) %>% 
  st_as_sfc() %>% 
  smoothr::densify(n = 1000)
north_america <- ne_countries(scale = 50, returnclass = "sf") %>% 
  filter(iso_a2 %in% c("US", "CA", "MX")) %>% 
  st_difference(hawaii) %>% 
  st_intersection(wh)
north_america_vect <- north_america %>% 
  st_transform(st_crs(region_boundary)) %>% 
  vect()

# species list ----

# migrants passing review
migrants <- ebirdst_runs %>% 
  filter(!resident, breeding_quality > 0, nonbreeding_quality > 0) %>% 
  select(species_code, common_name)

# species list
species_list <- glue("species_list_{group}.csv") %>% 
  path("data", .) %>% 
  read_csv() %>% 
  transmute(species_code = value) %>% 
  left_join(migrants,by = "species_code")
if (any(is.na(species_list$common_name))) {
  dropped <- filter(species_list, is.na(common_name)) %>% 
    pull(species_code)
  message("Missing species:\n", paste(dropped, collapse = "\n"))
}
species_list <- filter(species_list, !is.na(common_name))

# download associated raster layer
dl_tifs <- character()
for (s in species_list$species_code) {
  path <- ebirdst_download(s, pattern = "percent-population_seasonal_mean_mr")
  path <- path(path, "seasonal", 
               glue("{s}_percent-population_seasonal_mean_mr_2021.tif"))
  stopifnot(file_exists(path))
  dl_tifs <- c(dl_tifs, path)
}
species_list$tif <- dl_tifs


# weights ----

species_list$prop_pop <- foreach (f_tif = species_list$tif, 
                                  .combine = c) %dopar% {
  r <- rast(f_tif, lyrs = "nonbreeding")
  w <- extract(r, region_boundary_vect, 
               fun = "sum", weights = TRUE, 
               na.rm = TRUE)
  w[["nonbreeding"]]
}
species_list <- filter(species_list, prop_pop > 0.005)


# connectivity raster ----

# breeding season prop population layers
prop_pop_species <- map(species_list$tif, rast, lyrs = "breeding") %>% 
  rast() %>% 
  setNames(species_list$species_code)

# weighted sum of prop population
s_num <- sum(prop_pop_species * species_list$prop_pop, na.rm = TRUE)
# number of species with values for each cell
s_cnt <- sum(!is.na(prop_pop_species), na.rm = TRUE)
# weighted prop pop
s_wpp <- s_num / sum(species_list$prop_pop, na.rm = TRUE)
# mask out any cells with models only for some species
connectivity <- mask(s_wpp, s_cnt >= nlyr(prop_pop_species), 
                     maskvalue = c(0, NA))
# mask out regions outside of north america
filename <- glue("connectivity_{str_to_lower(group)}_mr_2021.tif") %>% 
  path("data", "connectivity", .)
connectivity <- mask(connectivity, north_america_vect) %>% 
  trim(filename = filename, overwrite = TRUE,
       gdal = c("COMPRESS=DEFLATE",
                "TILED=YES",
                "COPY_SRC_OVERVIEWS=YES"))


# connectivity map ----

# projection and extent
# lat/lng bounding box
combined_regions <- st_bbox(connectivity) %>% 
  st_as_sfc() %>% 
  c(region_boundary)
bb_crs <- combined_regions %>% 
  st_bbox()
bb_ll <- combined_regions %>% 
  st_transform(crs = 4326) %>% 
  st_bbox() %>% 
  round(3)
# custom crs
crs <- define_projection(bb_ll, force_global = TRUE)
# bounding box in custom crs
ext_order <- c("xmin", "xmax", "ymin", "ymax")
ext_crs <- ext(bb_crs[ext_order])

# basemap data
basemap <- projected_basemap(crs)

# define raster template and map extent
rt <- generate_raster_template(crs = crs$wkt, 
                               ext = ext_crs, 
                               res = res(connectivity))
proj_template <- rt$template
map_bbox <- rt$map_bbox

# project data
connectivity_proj <- project_raster(connectivity, proj_template, method = "near")
region_proj <- st_transform(region_boundary, crs = crs) %>% 
  st_geometry()

# bins
v <- values(connectivity_proj, na.rm = TRUE, mat = FALSE) %>% 
  as.numeric() %>% 
  keep(~ . > 0)
v_cutoff <- quantile(v, probs = 0.5)
v <- v[v >= v_cutoff]
q <- c(1 - 1 / 1.5^(0:9), 1)
bins <- quantile(v, probs = q, na.rm = TRUE)
bins <- c(0, bins)

# palette
pal <- c(map_theme$col_zero, rev(plasma(length(bins) - 2, end = 0.9)))

# map
glue("connectivity_{str_to_lower(group)}_mr_2021.png") %>% 
  path("maps", .) %>% 
  png(width = 2400, height = 2400)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

create_basemap(basemap, bbox = st_bbox(map_bbox), zeros = connectivity_proj)
plot(connectivity_proj, col = pal, breaks = bins,
     maxcell = ncell(connectivity_proj), 
     legend = FALSE, axes = FALSE, bty = "n",
     reset = FALSE, add = TRUE)
add_boundaries(basemap)
plot(region_proj, border = NA, col = "black", lwd = map_theme$lwd_countries, 
     add = TRUE)

# annotations
add_map_legend(palette = pal,
               title = "Strength of connection",
               labels = c("Low", "Moderate", "High"),
               labels_at = c(0, 0.1, 1),
               x = 0.90, y = 0.33, width = 0.015, height = 0.33,
               less_than = FALSE)
add_rectangle(labels = "No\nPrediction", cex = 2.5, pos = 4, x = 0.90, y = 0.28)

dev.off()
