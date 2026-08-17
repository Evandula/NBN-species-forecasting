library(dplyr)
library(readr)
library(sf)
library(terra)
library(tidyr)
library(rnaturalearth)
library(lubridate)
library(spThin)

setwd("~/Desktop/SPHERE/NBN-LPI-forecasting/")

keep_cols <- c("Occurrence.ID", "Scientific.name", "Event.Date", "Day", "Month",
               "Year", "Latitude..WGS84.", "Longitude..WGS84.")

extract = terra::extract

#------------------------------------------------------------------------------#
# Load data
#------------------------------------------------------------------------------#

uk_grid_5km <- st_read("data/uk_grid_5km.gpkg")
uk_regions <- st_read("data/uk_regions.gpkg")
uk_1km_master <- rast(crs = "EPSG:27700", ext = ext(0, 700000, 0, 1250000), res = 1000)


# Elevation
elevation <- rast("data/env_predictors/gb_dem/uk_srtm/w001001.adf")
  crs(elevation) <- "EPSG:27700"

elev_1km_mean <- resample(elevation, uk_1km_master, method = "average")
elev_1km_min   <- resample(elevation, uk_1km_master, method = "min")
elev_1km_max   <- resample(elevation, uk_1km_master, method = "max")
elev_1km_range <- elev_1km_max - elev_1km_min

elevation_1km_uk_list <- c(mean_elevation = elev_1km_mean, elevation_range = elev_1km_range)
elevation_1km_uk <- rast(elevation_1km_uk_list)


# Soil properties
soil_file_list <- list.files(path = "data/env_predictors/soil/")
soil_file_list <- paste0("data/env_predictors/soil/", soil_file_list)
soil_raster_stack <- rast(soil_file_list)

uk_bbox_wgs84 <- ext(-10, 3, 49, 61)
uk_bbox_homolosine <- project(uk_bbox_wgs84, from = "EPSG:4326", to = crs(soil_raster_stack))

soil_uk_cropped <- crop(soil_raster_stack, uk_bbox_homolosine)
soil_1km_uk <- project(soil_uk_cropped, uk_1km_master, method = "bilinear") # now 1km, EPSG:27700


#################################################################################
#LUH2 data
load("~/Desktop/luh2_data_hamonized.RData")
anthrop_layers <- rast(anthrop_layers)
# anthrop_layers <- aggregate(anthrop_layers, fact = 2, fun = mean)

nature_layers <- rast(nature_layers)
# nature_layers <- aggregate(nature_layers, fact = 2, fun = mean)

layer_years <- as.integer(gsub("X", "", names(anthrop_layers)))


# Bioclim data (COPERNICUS)
# bioclim_files = list.files("data/env_predictors/bioclim/", pattern = ".nc")
  bioclim_files = list.files("/Users/alexrabeau/Desktop/SPHERE/ADAS-pest-forecasting/data/bioclim_data/2c867a29093da8fad8e3b6fcf88d895a/", pattern = ".nc")

bioclim_vars = gsub("-es_rcp45_r2i1p1_1950-2100_v1.0.nc", "", bioclim_files)

# bioclim_dir = "data/env_predictors/bioclim"
  bioclim_dir = "/Users/alexrabeau/Desktop/SPHERE/ADAS-pest-forecasting/data/bioclim_data/2c867a29093da8fad8e3b6fcf88d895a"
#################################################################################
 
  
#------------------------------------------------------------------------------#
# Load functions
#------------------------------------------------------------------------------#

clean_species <- function(df) {
  df <- df %>%
    filter(Coordinate.Uncertainty..m. <= 5000) %>%
    filter(Year >= 1970) %>%
    filter(Identification.verification.status == "Accepted" |
             Identification.verification.status == "Accepted - considered correct" |
             Identification.verification.status == "Accepted - correct") %>%
    distinct(Occurrence.ID, .keep_all = TRUE) %>%
    select(all_of(keep_cols))
  # df <- thin(loc.data = df, lat.col = "Latitude..WGS84.", long.col = "Longitude..WGS84.",
  #            spec.col = "Scientific.name", thin.par = 0.5, reps = 1,
  #            locs.thinned.list.return = TRUE, write.files = FALSE)
  
  return(df)
  
}


extract_bioclim <- function(df, bioclim_dir, bioclim_files, bioclim_vars,
                            lon_col = "Longitude..WGS84.",
                            lat_col = "Latitude..WGS84.") {
  
  pts <- vect(df, geom = c(lon_col, lat_col), crs = "EPSG:4326")
  
  for (i in seq_along(bioclim_files)) {
    var_name <- bioclim_vars[i]
    r <- rast(file.path(bioclim_dir, bioclim_files[i]))
    layer_years <- year(time(r))
    
    vals <- raster::extract(r, pts, ID = FALSE)  
    
    col_idx <- match(df[["Year"]], layer_years)
    df[[var_name]] <- vals[cbind(seq_len(nrow(vals)), col_idx)]
  }
  
  df
  
}


extract_lu <- function(df, r, layer_years, var_name,
                       lon_col = "Longitude..WGS84.",
                       lat_col = "Latitude..WGS84.") {
  
  df <- df %>% filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
  pts <- vect(df, geom = c(lon_col, lat_col), crs = "EPSG:4326")
  
  vals <- extract(r, pts, ID = FALSE)
  col_idx <- match(df$Year, layer_years) #check this
  df[[var_name]] <- vals[cbind(seq_len(nrow(vals)), col_idx)]
  
  df
}


extract_static <- function(df, lon_col = "Longitude..WGS84.", lat_col = "Latitude..WGS84.") {
  
  pts <- vect(df, geom = c(lon_col, lat_col), crs = "EPSG:4326")
  elevation_df <- extract(elevation_1km_uk, pts, ID = FALSE)
  df$mean_elevation <- elevation_df$mean_elevation
  df$elevation_range <- elevation_df$elevation_range
  df_soil_extract <- extract(soil_1km_uk, pts, ID = FALSE)
  df <- cbind(df, df_soil_extract)
  
  df
  
}



extract_grid <- function(df){
  
  occurrence_sf <- df %>% 
    st_as_sf(coords = c("Longitude..WGS84.", "Latitude..WGS84."), crs = 4326, remove = FALSE) %>% 
    st_transform(st_crs(uk_grid_5km)) %>%
    st_join(
      uk_grid_5km %>%
        select(grid_id, Region),
      join = st_within
    ) %>%
    st_drop_geometry()
  
  return(occurrence_sf)
  
}




#------------------------------------------------------------------------------#
# Data Processing
#------------------------------------------------------------------------------#

#Euplagia_quadripunctaria
euplagia_quadripunctaria <- read.csv("data/NBN-atlas/species/Euplagia_quadripunctaria/records-2026-08-04.csv")
euplagia_quadripunctaria_v1 <- clean_species(euplagia_quadripunctaria)
euplagia_quadripunctaria_v2 <- extract_bioclim(euplagia_quadripunctaria_v1, bioclim_dir, bioclim_files, bioclim_vars)
euplagia_quadripunctaria_v3 <- extract_lu(euplagia_quadripunctaria_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
euplagia_quadripunctaria_v4 <- extract_lu(euplagia_quadripunctaria_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
euplagia_quadripunctaria_v5 <- extract_static(euplagia_quadripunctaria_v4)
euplagia_quadripunctaria_v6 <- extract_grid(euplagia_quadripunctaria_v5)
euplagia_quadripunctaria_v7 <- euplagia_quadripunctaria_v6 %>% drop_na()
write.csv(euplagia_quadripunctaria_v7, 
          "data/NBN-atlas/species/Euplagia_quadripunctaria/euplagia_quadripunctaria_processed.csv", row.names = FALSE)

#Fratercula_arctica
fratercula_arctica <- read.csv("data/NBN-atlas/species/Fratercula_arctica/records-2026-08-04.csv")
fratercula_arctica_v1 <- clean_species(fratercula_arctica)
fratercula_arctica_v2 <- extract_bioclim(fratercula_arctica_v1, bioclim_dir, bioclim_files, bioclim_vars)
fratercula_arctica_v3 <- extract_lu(fratercula_arctica_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
fratercula_arctica_v4 <- extract_lu(fratercula_arctica_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
fratercula_arctica_v5 <- extract_static(fratercula_arctica_v4)
fratercula_arctica_v6 <- extract_grid(fratercula_arctica_v5)
fratercula_arctica_v7 <- fratercula_arctica_v6 %>% drop_na()
write.csv(fratercula_arctica_v7, 
          "data/NBN-atlas/species/Fratercula_arctica/fratercula_arctica_processed.csv", row.names = FALSE)

#Haliaeetus_albicilla
haliaeetus_albicilla <- read.csv("data/NBN-atlas/species/Haliaeetus_albicilla/records-2026-08-04.csv")
haliaeetus_albicilla_v1 <- clean_species(haliaeetus_albicilla)
haliaeetus_albicilla_v2 <- extract_bioclim(haliaeetus_albicilla_v1, bioclim_dir, bioclim_files, bioclim_vars)
haliaeetus_albicilla_v3 <- extract_lu(haliaeetus_albicilla_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
haliaeetus_albicilla_v4 <- extract_lu(haliaeetus_albicilla_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
haliaeetus_albicilla_v5 <- extract_static(haliaeetus_albicilla_v4)
haliaeetus_albicilla_v6 <- extract_grid(haliaeetus_albicilla_v5)
haliaeetus_albicilla_v7 <- haliaeetus_albicilla_v6 %>% drop_na()
write.csv(haliaeetus_albicilla_v7, 
          "data/NBN-atlas/species/Haliaeetus_albicilla/haliaeetus_albicilla_processed.csv", row.names = FALSE)

#Lutra_lutra
lutra_lutra <- read.csv("data/NBN-atlas/species/Lutra_lutra/records-2026-08-04.csv")
lutra_lutra_v1 <- clean_species(lutra_lutra)
lutra_lutra_v2 <- extract_bioclim(lutra_lutra_v1, bioclim_dir, bioclim_files, bioclim_vars)
lutra_lutra_v3 <- extract_lu(lutra_lutra_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
lutra_lutra_v4 <- extract_lu(lutra_lutra_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
lutra_lutra_v5 <- extract_static(lutra_lutra_v4)
lutra_lutra_v6 <- extract_grid(lutra_lutra_v5)
lutra_lutra_v7 <- lutra_lutra_v6 %>% drop_na()
write.csv(lutra_lutra_v7, 
          "data/NBN-atlas/species/Lutra_lutra/lutra_lutra_processed.csv", row.names = FALSE)

#Pandion_haliaetus
pandion_haliaetus <- read.csv("data/NBN-atlas/species/Pandion_haliaetus/records-2026-08-04.csv")
pandion_haliaetus_v1 <- clean_species(pandion_haliaetus)
pandion_haliaetus_v2 <- extract_bioclim(pandion_haliaetus_v1, bioclim_dir, bioclim_files, bioclim_vars)
pandion_haliaetus_v3 <- extract_lu(pandion_haliaetus_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
pandion_haliaetus_v4 <- extract_lu(pandion_haliaetus_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
pandion_haliaetus_v5 <- extract_static(pandion_haliaetus_v4)
pandion_haliaetus_v6 <- extract_grid(pandion_haliaetus_v5)
pandion_haliaetus_v7 <- pandion_haliaetus_v6 %>% drop_na()
write.csv(pandion_haliaetus_v7, 
          "data/NBN-atlas/species/Pandion_haliaetus/pandion_haliaetus_processed.csv", row.names = FALSE)

#Pipistrellus_pipistrellus
pipistrellus_pipistrellus <- read.csv("data/NBN-atlas/species/Pipistrellus_pipistrellus/records-2026-08-04.csv")
pipistrellus_pipistrellus_v1 <- clean_species(pipistrellus_pipistrellus)
pipistrellus_pipistrellus_v2 <- extract_bioclim(pipistrellus_pipistrellus_v1, bioclim_dir, bioclim_files, bioclim_vars)
pipistrellus_pipistrellus_v3 <- extract_lu(pipistrellus_pipistrellus_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
pipistrellus_pipistrellus_v4 <- extract_lu(pipistrellus_pipistrellus_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
pipistrellus_pipistrellus_v5 <- extract_static(pipistrellus_pipistrellus_v4)
pipistrellus_pipistrellus_v6 <- extract_grid(pipistrellus_pipistrellus_v5)
pipistrellus_pipistrellus_v7 <- pipistrellus_pipistrellus_v6 %>% drop_na()
write.csv(pipistrellus_pipistrellus_v7, 
          "data/NBN-atlas/species/Pipistrellus_pipistrellus/pipistrellus_pipistrellus_processed.csv", row.names = FALSE)

#Sciurus_vulgaris
sciurus_vulgaris <- read.csv("data/NBN-atlas/species/Sciurus_vulgaris/records-2026-08-04.csv")
sciurus_vulgaris_v1 <- clean_species(sciurus_vulgaris)
sciurus_vulgaris_v2 <- extract_bioclim(sciurus_vulgaris_v1, bioclim_dir, bioclim_files, bioclim_vars)
sciurus_vulgaris_v3 <- extract_lu(sciurus_vulgaris_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
sciurus_vulgaris_v4 <- extract_lu(sciurus_vulgaris_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
sciurus_vulgaris_v5 <- extract_static(sciurus_vulgaris_v4)
sciurus_vulgaris_v6 <- extract_grid(sciurus_vulgaris_v5)
sciurus_vulgaris_v7 <- sciurus_vulgaris_v6 %>% drop_na()
write.csv(sciurus_vulgaris_v7, 
          "data/NBN-atlas/species/Sciurus_vulgaris/sciurus_vulgaris_processed.csv", row.names = FALSE)


#Triturus_cristatus
triturus_cristatus <- read.csv("data/NBN-atlas/species/Triturus_cristatus/records-2026-08-04.csv")
triturus_cristatus_v1 <- clean_species(triturus_cristatus)
triturus_cristatus_v2 <- extract_bioclim(triturus_cristatus_v1, bioclim_dir, bioclim_files, bioclim_vars)
triturus_cristatus_v3 <- extract_lu(triturus_cristatus_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
triturus_cristatus_v4 <- extract_lu(triturus_cristatus_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
triturus_cristatus_v5 <- extract_static(triturus_cristatus_v4)
triturus_cristatus_v6 <- extract_grid(triturus_cristatus_v5)
triturus_cristatus_v7 <- triturus_cristatus_v6 %>% drop_na()
write.csv(triturus_cristatus_v7, 
          "data/NBN-atlas/species/Triturus_cristatus/triturus_cristatus_processed.csv", row.names = FALSE)

#Vipera_berus
vipera_berus <- read.csv("data/NBN-atlas/species/Vipera_berus/records-2026-08-04.csv")
vipera_berus_v1 <- clean_species(vipera_berus)
vipera_berus_v2 <- extract_bioclim(vipera_berus_v1, bioclim_dir, bioclim_files, bioclim_vars)
vipera_berus_v3 <- extract_lu(vipera_berus_v2, r = anthrop_layers, layer_years = layer_years, var_name = "anthrop_lu")
vipera_berus_v4 <- extract_lu(vipera_berus_v3, r = nature_layers, layer_years = layer_years, var_name = "nature_lu")
vipera_berus_v5 <- extract_static(vipera_berus_v4)
vipera_berus_v6 <- extract_grid(vipera_berus_v5)
vipera_berus_v7 <- vipera_berus_v6 %>% drop_na()
write.csv(vipera_berus_v7, 
          "data/NBN-atlas/species/Vipera_berus/vipera_berus_processed.csv", row.names = FALSE)
