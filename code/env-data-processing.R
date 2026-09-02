library(dplyr)
library(readr)
library(sf)
library(terra)
library(tidyr)
library(lubridate)


setwd("~/Desktop/SPHERE/NBN-LPI-forecasting/")

set.seed(123)

keep_cols <- c("Occurrence.ID", "Scientific.name", "Event.Date", "Day", "Month",
               "Year", "Latitude..WGS84.", "Longitude..WGS84.")

extract = terra::extract
project = terra::project
crs = terra::crs

output_dir <- "data/env_predictors/gridded/"


# ============================================================
# Load data
# ============================================================

uk_grid_5km <- st_read("data/uk_grid_5km.gpkg")
uk_5km_master <- rast(crs = "EPSG:27700", ext = ext(0, 700000, -100000, 1250000), res = 5000)

# Elevation
elevation <- rast("data/env_predictors/gb_dem/uk_srtm/w001001.adf")
crs(elevation) <- "EPSG:27700"

elev_5km_mean <- resample(elevation, uk_5km_master, method = "average")
elev_5km_min   <- resample(elevation, uk_5km_master, method = "min")
elev_5km_max   <- resample(elevation, uk_5km_master, method = "max")
elev_5km_range <- elev_5km_max - elev_5km_min

elevation_5km_uk_list <- c(mean_elevation = elev_5km_mean, elevation_range = elev_5km_range)
elevation_5km_uk <- rast(elevation_5km_uk_list)


# Soil properties
soil_file_list <- list.files(path = "data/env_predictors/soil/")
soil_file_list <- paste0("data/env_predictors/soil/", soil_file_list)
soil_raster_stack <- rast(soil_file_list)

uk_bbox_wgs84 <- ext(-10, 3, 49, 61)
uk_bbox_homolosine <- project(uk_bbox_wgs84, from = "EPSG:4326", to = crs(soil_raster_stack))

soil_uk_cropped <- crop(soil_raster_stack, uk_bbox_homolosine)
soil_5km_uk <- project(soil_uk_cropped, uk_5km_master, method = "average") # now 5km, EPSG:27700


#LUH2 data
load("~/Desktop/luh2_data_hamonized.RData")

anthrop_layers <- rast(anthrop_layers)
anthrop_uk <- crop(anthrop_layers, uk_bbox_wgs84)
anthrop_uk <- project(anthrop_uk, uk_5km_master, method = "average") # now 5km, EPSG:27700

nature_layers <- rast(nature_layers)
nature_uk <- crop(nature_layers, uk_bbox_wgs84)
nature_uk <- project(nature_uk, uk_5km_master, method = "average") # now 5km, EPSG:27700

layer_years_lu <- as.integer(gsub("X", "", names(anthrop_layers)))


# Bioclim data
bioclim_dir = "/Volumes/ajr221/ephemeral/UK_bioclim_1km/dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.2.ceda/1km/"
bioclim_dir_2 = "/ann/v20260512/"
  
bioclim_vars = c("groundfrost", "hurs", "pv", "rainfall", "sun", "tasmax", "tasmin")

layer_years_bioclim = c(1970:2025)
  # layer_years_bioclim = names(bioclim_var1_5km) %>% as.integer()

# var 1: groundfrost
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[1], bioclim_dir_2), pattern = ".nc")
bioclim_var1_stack = rast(paste0(bioclim_dir, bioclim_vars[1], bioclim_dir_2, bioclim_files))
bioclim_var1_5km <- project(bioclim_var1_stack, uk_5km_master, method = "average")
  names(bioclim_var1_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var1_stack))

# var 2: hurs
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[2], bioclim_dir_2), pattern = ".nc")
bioclim_var2_stack = rast(paste0(bioclim_dir, bioclim_vars[2], bioclim_dir_2, bioclim_files))
bioclim_var2_5km <- project(bioclim_var2_stack, uk_5km_master, method = "average")
names(bioclim_var2_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var2_stack))
  
# var 3: pv
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[3], bioclim_dir_2), pattern = ".nc")
bioclim_var3_stack = rast(paste0(bioclim_dir, bioclim_vars[3], bioclim_dir_2, bioclim_files))
bioclim_var3_5km <- project(bioclim_var3_stack, uk_5km_master, method = "average")
names(bioclim_var3_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var3_stack))

# var 4: rainfall
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[4], bioclim_dir_2), pattern = ".nc")
bioclim_var4_stack = rast(paste0(bioclim_dir, bioclim_vars[4], bioclim_dir_2, bioclim_files))
bioclim_var4_5km <- project(bioclim_var4_stack, uk_5km_master, method = "average")
names(bioclim_var4_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var4_stack))

# var 5: sun
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[5], bioclim_dir_2), pattern = ".nc")
bioclim_var5_stack = rast(paste0(bioclim_dir, bioclim_vars[5], bioclim_dir_2, bioclim_files))
bioclim_var5_5km <- project(bioclim_var5_stack, uk_5km_master, method = "average")
names(bioclim_var5_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var5_stack))

# var 6: tasmax
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[6], bioclim_dir_2), pattern = ".nc")
bioclim_var6_stack = rast(paste0(bioclim_dir, bioclim_vars[6], bioclim_dir_2, bioclim_files))
bioclim_var6_5km <- project(bioclim_var6_stack, uk_5km_master, method = "average")
names(bioclim_var6_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var6_stack))

# var 7: tasmin
bioclim_files = list.files(paste0(bioclim_dir, bioclim_vars[7], bioclim_dir_2), pattern = ".nc")
bioclim_var7_stack = rast(paste0(bioclim_dir, bioclim_vars[7], bioclim_dir_2, bioclim_files))
bioclim_var7_5km <- project(bioclim_var7_stack, uk_5km_master, method = "average")
names(bioclim_var7_5km) <- sub(".*uk_1km_ann_([0-9]{4}).*", "\\1", sources(bioclim_var7_stack))

bioclim_vars_list <- list(groundfrost = bioclim_var1_5km, hurs = bioclim_var2_5km, pv = bioclim_var3_5km,
     rainfall = bioclim_var4_5km, sun = bioclim_var5_5km, tasmax = bioclim_var6_5km, tasmin = bioclim_var7_5km)


# ============================================================
# Load functions
# ============================================================

clean_species <- function(df) {
  df <- df %>%
    filter(Coordinate.Uncertainty..m. <= 5000) %>%
    filter(Year >= 1970) %>%
    filter(Identification.verification.status == "Accepted" |
             Identification.verification.status == "Accepted - considered correct" |
             Identification.verification.status == "Accepted - correct") %>%
    distinct(Occurrence.ID, .keep_all = TRUE) %>%
    select(all_of(keep_cols))

  return(df)
  
}


extract_bioclim <- function(df, lon_col = "Longitude..WGS84.", lat_col = "Latitude..WGS84.") {
  
  pts <- vect(df, geom = c(lon_col, lat_col))
  crs(pts) <- "EPSG:4326"
  pts <- project(pts, "EPSG:27700")
  
  for (i in seq_along(bioclim_vars)) {
    
    var_name <- bioclim_vars[i]
    
    r <- bioclim_vars_list[[var_name]]
    
    vals <- extract(r, pts, ID = FALSE)  
    
    env_year <- df[["Year"]] - 1
    
    col_idx <- match(env_year, layer_years_bioclim)
    
    df[[var_name]] <- vals[cbind(seq_len(nrow(vals)), col_idx)]
    
    df[["env_year"]] <- env_year
    
  }
  
  df
  
}


extract_lu <- function(df, r, var_name, lon_col = "Longitude..WGS84.", lat_col = "Latitude..WGS84.") {
  
  pts <- vect(df, geom = c(lon_col, lat_col))
  crs(pts) <- "EPSG:4326"
  pts <- project(pts, "EPSG:27700")
  
  vals <- extract(r, pts, ID = FALSE)
  
  env_year <- df[["Year"]] - 1
  
  col_idx <- match(env_year, layer_years_lu)
  
  df[[var_name]] <- vals[cbind(seq_len(nrow(vals)), col_idx)]
  
  df[["env_year"]] <- env_year
  
  
  df
}


extract_static <- function(df, lon_col = "Longitude..WGS84.", lat_col = "Latitude..WGS84.") {
  
  pts <- vect(df, geom = c(lon_col, lat_col))
  crs(pts) <- "EPSG:4326"
  pts <- project(pts, "EPSG:27700")
  
  elevation_df <- extract(elevation_5km_uk, pts, ID = FALSE)
  df$mean_elevation <- elevation_df$mean_elevation
  df$elevation_range <- elevation_df$elevation_range
  df_soil_extract <- extract(soil_5km_uk, pts, ID = FALSE)
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
    st_drop_geometry() %>%
    filter(!is.na(grid_id)) %>%
    group_by(grid_id) %>%
    mutate(.random = runif(n())) %>%
    arrange(.random) %>%
    slice_head(n = 100) %>%
    ungroup() %>%
    select(-.random)
    
  return(occurrence_sf)
  
}


# ============================================================
# 1. Occurrence Data Processing
# ============================================================

#Euplagia_quadripunctaria
euplagia_quadripunctaria <- read.csv("data/NBN-atlas/species/Euplagia_quadripunctaria/records-2026-08-04.csv")
euplagia_quadripunctaria_v1 <- clean_species(euplagia_quadripunctaria)
euplagia_quadripunctaria_v2 <- extract_bioclim(euplagia_quadripunctaria_v1)
euplagia_quadripunctaria_v3 <- extract_lu(euplagia_quadripunctaria_v2, r = anthrop_uk, var_name = "anthrop_lu")
euplagia_quadripunctaria_v4 <- extract_lu(euplagia_quadripunctaria_v3, r = nature_uk, var_name = "nature_lu")
euplagia_quadripunctaria_v5 <- extract_static(euplagia_quadripunctaria_v4)
euplagia_quadripunctaria_v6 <- extract_grid(euplagia_quadripunctaria_v5)
euplagia_quadripunctaria_v7 <- euplagia_quadripunctaria_v6 %>% drop_na()
write.csv(euplagia_quadripunctaria_v7, 
          "data/NBN-atlas/species/Euplagia_quadripunctaria/euplagia_quadripunctaria_processed.csv", row.names = FALSE)

#Fratercula_arctica
fratercula_arctica <- read.csv("data/NBN-atlas/species/Fratercula_arctica/records-2026-08-04.csv")
fratercula_arctica_v1 <- clean_species(fratercula_arctica)
fratercula_arctica_v2 <- extract_bioclim(fratercula_arctica_v1)
fratercula_arctica_v3 <- extract_lu(fratercula_arctica_v2, r = anthrop_uk, var_name = "anthrop_lu")
fratercula_arctica_v4 <- extract_lu(fratercula_arctica_v3, r = nature_uk, var_name = "nature_lu")
fratercula_arctica_v5 <- extract_static(fratercula_arctica_v4)
fratercula_arctica_v6 <- extract_grid(fratercula_arctica_v5)
fratercula_arctica_v7 <- fratercula_arctica_v6 %>% drop_na()
write.csv(fratercula_arctica_v7, 
          "data/NBN-atlas/species/Fratercula_arctica/fratercula_arctica_processed.csv", row.names = FALSE)

#Haliaeetus_albicilla
haliaeetus_albicilla <- read.csv("data/NBN-atlas/species/Haliaeetus_albicilla/records-2026-08-04.csv")
haliaeetus_albicilla_v1 <- clean_species(haliaeetus_albicilla)
haliaeetus_albicilla_v2 <- extract_bioclim(haliaeetus_albicilla_v1)
haliaeetus_albicilla_v3 <- extract_lu(haliaeetus_albicilla_v2, r = anthrop_uk, var_name = "anthrop_lu")
haliaeetus_albicilla_v4 <- extract_lu(haliaeetus_albicilla_v3, r = nature_uk, var_name = "nature_lu")
haliaeetus_albicilla_v5 <- extract_static(haliaeetus_albicilla_v4)
haliaeetus_albicilla_v6 <- extract_grid(haliaeetus_albicilla_v5)
haliaeetus_albicilla_v7 <- haliaeetus_albicilla_v6 %>% drop_na()
write.csv(haliaeetus_albicilla_v7, 
          "data/NBN-atlas/species/Haliaeetus_albicilla/haliaeetus_albicilla_processed.csv", row.names = FALSE)

#Lutra_lutra
lutra_lutra <- read.csv("data/NBN-atlas/species/Lutra_lutra/records-2026-08-04.csv")
lutra_lutra_v1 <- clean_species(lutra_lutra)
lutra_lutra_v2 <- extract_bioclim(lutra_lutra_v1)
lutra_lutra_v3 <- extract_lu(lutra_lutra_v2, r = anthrop_uk, var_name = "anthrop_lu")
lutra_lutra_v4 <- extract_lu(lutra_lutra_v3, r = nature_uk, var_name = "nature_lu")
lutra_lutra_v5 <- extract_static(lutra_lutra_v4)
lutra_lutra_v6 <- extract_grid(lutra_lutra_v5)
lutra_lutra_v7 <- lutra_lutra_v6 %>% drop_na()
write.csv(lutra_lutra_v7, 
          "data/NBN-atlas/species/Lutra_lutra/lutra_lutra_processed.csv", row.names = FALSE)

#Pandion_haliaetus
pandion_haliaetus <- read.csv("data/NBN-atlas/species/Pandion_haliaetus/records-2026-08-04.csv")
pandion_haliaetus_v1 <- clean_species(pandion_haliaetus)
pandion_haliaetus_v2 <- extract_bioclim(pandion_haliaetus_v1)
pandion_haliaetus_v3 <- extract_lu(pandion_haliaetus_v2, r = anthrop_uk,  var_name = "anthrop_lu")
pandion_haliaetus_v4 <- extract_lu(pandion_haliaetus_v3, r = nature_uk, var_name = "nature_lu")
pandion_haliaetus_v5 <- extract_static(pandion_haliaetus_v4)
pandion_haliaetus_v6 <- extract_grid(pandion_haliaetus_v5)
pandion_haliaetus_v7 <- pandion_haliaetus_v6 %>% drop_na()
write.csv(pandion_haliaetus_v7, 
          "data/NBN-atlas/species/Pandion_haliaetus/pandion_haliaetus_processed.csv", row.names = FALSE)

#Pipistrellus_pipistrellus
pipistrellus_pipistrellus <- read.csv("data/NBN-atlas/species/Pipistrellus_pipistrellus/records-2026-08-04.csv")
pipistrellus_pipistrellus_v1 <- clean_species(pipistrellus_pipistrellus)
pipistrellus_pipistrellus_v2 <- extract_bioclim(pipistrellus_pipistrellus_v1)
pipistrellus_pipistrellus_v3 <- extract_lu(pipistrellus_pipistrellus_v2, r = anthrop_uk, var_name = "anthrop_lu")
pipistrellus_pipistrellus_v4 <- extract_lu(pipistrellus_pipistrellus_v3, r = nature_uk, var_name = "nature_lu")
pipistrellus_pipistrellus_v5 <- extract_static(pipistrellus_pipistrellus_v4)
pipistrellus_pipistrellus_v6 <- extract_grid(pipistrellus_pipistrellus_v5)
pipistrellus_pipistrellus_v7 <- pipistrellus_pipistrellus_v6 %>% drop_na()
write.csv(pipistrellus_pipistrellus_v7, 
          "data/NBN-atlas/species/Pipistrellus_pipistrellus/pipistrellus_pipistrellus_processed.csv", row.names = FALSE)

#Sciurus_vulgaris
sciurus_vulgaris <- read.csv("data/NBN-atlas/species/Sciurus_vulgaris/records-2026-08-04.csv")
sciurus_vulgaris_v1 <- clean_species(sciurus_vulgaris)
sciurus_vulgaris_v2 <- extract_bioclim(sciurus_vulgaris_v1)
sciurus_vulgaris_v3 <- extract_lu(sciurus_vulgaris_v2, r = anthrop_uk, var_name = "anthrop_lu")
sciurus_vulgaris_v4 <- extract_lu(sciurus_vulgaris_v3, r = nature_uk, var_name = "nature_lu")
sciurus_vulgaris_v5 <- extract_static(sciurus_vulgaris_v4)
sciurus_vulgaris_v6 <- extract_grid(sciurus_vulgaris_v5)
sciurus_vulgaris_v7 <- sciurus_vulgaris_v6 %>% drop_na()
write.csv(sciurus_vulgaris_v7, 
          "data/NBN-atlas/species/Sciurus_vulgaris/sciurus_vulgaris_processed.csv", row.names = FALSE)


#Triturus_cristatus
triturus_cristatus <- read.csv("data/NBN-atlas/species/Triturus_cristatus/records-2026-08-04.csv")
triturus_cristatus_v1 <- clean_species(triturus_cristatus)
triturus_cristatus_v2 <- extract_bioclim(triturus_cristatus_v1)
triturus_cristatus_v3 <- extract_lu(triturus_cristatus_v2, r = anthrop_uk, var_name = "anthrop_lu")
triturus_cristatus_v4 <- extract_lu(triturus_cristatus_v3, r = nature_uk, var_name = "nature_lu")
triturus_cristatus_v5 <- extract_static(triturus_cristatus_v4)
triturus_cristatus_v6 <- extract_grid(triturus_cristatus_v5)
triturus_cristatus_v7 <- triturus_cristatus_v6 %>% drop_na()
write.csv(triturus_cristatus_v7, 
          "data/NBN-atlas/species/Triturus_cristatus/triturus_cristatus_processed.csv", row.names = FALSE)

#Vipera_berus
vipera_berus <- read.csv("data/NBN-atlas/species/Vipera_berus/records-2026-08-04.csv")
vipera_berus_v1 <- clean_species(vipera_berus)
vipera_berus_v2 <- extract_bioclim(vipera_berus_v1)
vipera_berus_v3 <- extract_lu(vipera_berus_v2, r = anthrop_uk, var_name = "anthrop_lu")
vipera_berus_v4 <- extract_lu(vipera_berus_v3, r = nature_uk, var_name = "nature_lu")
vipera_berus_v5 <- extract_static(vipera_berus_v4)
vipera_berus_v6 <- extract_grid(vipera_berus_v5)
vipera_berus_v7 <- vipera_berus_v6 %>% drop_na()
write.csv(vipera_berus_v7, 
          "data/NBN-atlas/species/Vipera_berus/vipera_berus_processed.csv", row.names = FALSE)




# ============================================================
# 2. Gridded Environmental Predictors
# ============================================================

grid_centres <- st_centroid(uk_grid_5km)
grid_centres <- st_transform(grid_centres, 4326)

coords <- st_coordinates(grid_centres)

grid_centres$longitude <- coords[, 1]
grid_centres$latitude  <- coords[, 2]


# ============================================================
# Load functions
# ============================================================

extract_raster <- function(points, raster_file, variable_name) {
  
  r <- raster_file
  
  # Convert points to terra
  pts <- vect(points)
  
  # Reproject points to raster CRS
  pts <- project(pts, crs(r))
  
  # Extract raster values
  values <- extract(r, pts)
  
  # First column is ID generated by terra
  values <- values[, -1, drop = FALSE]
  
  # If raster has one layer
  if (ncol(values) == 1) {
    
    points[[variable_name]] <- values[[1]]
    
  } else {
    
    # Multiple layers
    names(values) <- paste0(
      variable_name,
      "_",
      seq_len(ncol(values))
    )
    
    points <- bind_cols(
      points,
      values
    )
  }
  
  return(points)
}


# ============================================================
# Static predictors
# ============================================================

static_rasters <- c(
  setNames(as.list(elevation_5km_uk), names(elevation_5km_uk)),
  setNames(as.list(soil_5km_uk), names(soil_5km_uk))
)


static_grid <- extract_raster(grid_centres, static_rasters[[1]], names(static_rasters)[1])

for (variable in names(static_rasters)[2:length(static_rasters)]) {
  
  cat("Extracting static:", variable, "\n")
  
  static_grid <- extract_raster(static_grid, static_rasters[[variable]], variable)
  
}

# ============================================================
# Dynamic predictors
# ============================================================
# year = 1970

cols_to_edit = c("clay_0.5cm_mean_1000", "clay_15.30cm_mean_1000", "clay_5.15cm_mean_1000", "nitrogen_0.5cm_mean_1000", "nitrogen_15.30cm_mean_1000", 
                 "nitrogen_5.15cm_mean_1000", "phh2o_0.5cm_mean_1000", "phh2o_15.30cm_mean_1000", "phh2o_5.15cm_mean_1000", "sand_0.5cm_mean_1000", 
                 "sand_15.30cm_mean_1000", "sand_5.15cm_mean_1000", "soc_0.5cm_mean_1000", "soc_15.30cm_mean_1000", "soc_5.15cm_mean_1000")

for (year in years) {
  
  cat("\nProcessing", year, "\n")
  
  static_grid_year = static_grid
  
  # ----------------------------------------------------------
  # i. Extract bioclim to grid centroids
  # ----------------------------------------------------------
  bioclim_vars_list <- list(groundfrost = bioclim_var1_5km, hurs = bioclim_var2_5km, pv = bioclim_var3_5km,
                            rainfall = bioclim_var4_5km, sun = bioclim_var5_5km, tasmax = bioclim_var6_5km, tasmin = bioclim_var7_5km)
  
  for(j in 1:length(bioclim_vars_list)){
    
    bioclim <- bioclim_vars_list[[j]]
    bioclim_v2 <- bioclim[[as.character(year)]]
    
    pts <- vect(grid_centres)
    pts <- project(pts, crs(bioclim_v2))
    
    bioclim_values <- extract(bioclim_v2, pts)
    bioclim_values <- bioclim_values[, -1, drop = FALSE]
    
    static_grid_year[[names(bioclim_vars_list)[j]]] <- bioclim_values[[1]]
    
  }
  
  cat("  bioclim processed\n")
  
  # ----------------------------------------------------------
  # ii. Extract LUH2 to grid centroids
  # ----------------------------------------------------------
  
  anthrop_luh2 <- anthrop_uk[[grep(year, names(anthrop_uk))]]
  anthrop_luh2_values <- extract(anthrop_luh2, grid_centres)
  anthrop_luh2_values <- anthrop_luh2_values[, -1, drop = FALSE]
  static_grid_year[["anthrop_luh2"]] <- anthrop_luh2_values[[1]]
  
  nature_luh2 <- nature_uk[[grep(year, names(nature_uk))]]
  nature_luh2_values <- extract(nature_luh2, grid_centres)
  nature_luh2_values <- nature_luh2_values[, -1, drop = FALSE]
  static_grid_year[["nature_luh2"]] <- nature_luh2_values[[1]]
  
  cat("  LUH2 processed\n")
  
  # ----------------------------------------------------------
  # iii. Export
  # ----------------------------------------------------------
  static_grid_year_export[['env_year']] <- year - 1
  colnames(static_grid_year_export)[colnames(static_grid_year_export) %in% cols_to_edit] = gsub('\\.', "-", colnames(static_grid_year_export)[colnames(static_grid_year_export) %in% cols_to_edit])
  
  static_grid_year_export <- static_grid_year %>% st_drop_geometry()
  static_grid_year_export <- static_grid_year_export %>% dplyr::rename('Longitude..WGS84.' = 'longitude', 
                                                                       'Latitude..WGS84.' = 'latitude', 
                                                                       'anthrop_lu' = 'anthrop_luh2', 
                                                                       'nature_lu' = 'nature_luh2')
  
  # colnames(static_grid_year_export) = gsub("\\.", "-", colnames(static_grid_year_export))
  
  write.csv(static_grid_year_export, paste0(output_dir, "gridden_env_", year, ".csv"), row.names = FALSE)
  
}

