library(ebirdst)
library(terra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(glue)
library(fs)

key <- ebirdst:::get_ebirdst_access_key()
api_url <- "https://st-download.ebird.org/v1/fetch"
dl_dir <- path("data", "raw", "ebirdst_high-res")
dir_create(dl_dir)

# get s&t projection
# path <- ebirdst_download("example_data")
crs <- get_species_path("example_data") %>% 
  dir_ls(glob = "*.tif", recurse = TRUE) %>% 
  head(1) %>% 
  rast() %>% 
  st_crs()

# define focal region
central_america <- ne_countries(scale = 10,
                                continent = "North America", 
                                returnclass = "sf") %>% 
  filter(subregion == "Central America") %>% 
  st_transform(crs = crs) %>% 
  select(country_code = iso_a2, country = name) %>%
  filter(country != "Mexico")
central_america_vect <- vect(central_america)

# read in list of species with eBird abundance data in Central America
# THIS WILL SPEED UP THE PROCESS - got this list from downloading the 
# low resolution maps, and can just download the high-res maps for this
# subset of ~400 species
species_list <- readLines("data/raw/species-list.txt")

# getting error for large files because download time >60 secs; increase
# timeout 
options(timeout=500)

# download all low resolution weekly abundance raster stacks
for (species in sort(species_list)) {
  message("Downloading ", species)
  
  # download percent of population
  filename <- glue("{species}_percent-population_median_hr_2021.tif")
  src_url <- glue("{api_url}?objKey=2021/{species}/weekly/{filename}&key={key}")
  temp_pop <- tempfile(fileext = ".tif")
  dst_pop <- path(dl_dir, filename)
  if (file_exists(dst_pop)) {
    next
  }
  download.file(src_url, temp_pop, quiet = TRUE)
  
  # is the species found in this region
  r_pop <- rast(temp_pop) %>% 
    crop(central_america_vect) %>% 
    mask(central_america_vect)
  # for each week, calculate total % of pop in region
  weekly_pop <- global(r_pop, sum, na.rm = TRUE)
  if (max(weekly_pop$sum, na.rm = TRUE) < 0.0001) {
    file_delete(temp_pop)
    next
  }
  
  # save to data/ directory
  writeRaster(r_pop, dst_pop, 
              gdal = c("COMPRESS=DEFLATE"),
              overwrite = TRUE)
}