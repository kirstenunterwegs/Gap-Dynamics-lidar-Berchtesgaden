########################################################
#
# Script to generate Canopy Height Model from Lidar data
# for the Berchtesgaden National Park (year 2021)
#
# author: Dirk Pflugmacher (Humboldt Universität zu Berlin)
# last revised: 20.10.2023
#
#######################################################

# --- libaries ---

library(lidR)
library(sf)
library(dplyr)
library(terra)
require(future)

#--- paths --------------

root <- "z:/Fonda/workspace/berchtesgaden/lidar/2021"
root <- "f:/Spring/Berchtesgaden/lidar/2021"
p_las <- file.path(root, "las")
p_norm <- file.path(root, "las_norm")
p_dtm <- file.path(root, "dtm")
p_chm <- file.path(root, "chm")
p_int <- file.path(root, "intensity")
dtm_file <- file.path(root, "berchtesgaden_2021_dtm.tif")
chm_file <- file.path(root, "berchtesgaden_2021_chm.tif")
chm1_file <- file.path(root, "berchtesgaden_2021_chm1.tif")
density_file <- file.path(root, "berchtesgaden_2021_lidar_density.tif")
int_file <- file.path(root, "berchtesgaden_2021_lidar_intensity.tif")

#--- setup future --------------

plan(multisession, workers = 15L)

#--- density ----------------------
ctg <- readLAScatalog(p_las)
lidR:::catalog_laxindex(ctg)

opt_output_files(ctg) <- paste0(file.path(root, "density"), "/{*}_density")
dens <- rasterize_density(ctg, 1)

terra::writeRaster(dens, density_file)


#--- terrain-----------------------

ctg <- readLAScatalog(p_las)
opt_output_files(ctg) <- paste0(p_dtm, "/{*}_dtm")
opt_filter(ctg) = "-keep_class 2 3 -drop_scan_angle_below 75 -drop_scan_angle_above 105"

dtm <- rasterize_terrain(ctg, res = 0.5, algorithm = tin())

terra::writeRaster(dtm, dtm_file)


#--- normalize height --------------

ctg <- readLAScatalog(p_las)
dtm <- raster::raster(dtm_file)

opt_output_files(ctg) <-  paste0(p_norm, "/{*}_norm")
opt_filter(ctg) = "-keep_class 2 3 -drop_scan_angle_below 75 -drop_scan_angle_above 105"

ctg_norm <- normalize_height(ctg, algorithm=dtm)


# suffix <- strsplit(basename(opt_output_files(ctg)), split="_")[[1]][2]
# tile_processed <- sub(paste0("_", suffix, ".las"), "", list.files(dirname(opt_output_files(ctg)), ".las$"))
# in_tiles <- sub(".las", "", basename(ctg$filename))
# ctg_new <- ctg[! in_tiles %in% tile_processed,]
# length((ctg_new$filename))
# ctg_norm <- normalize_height(ctg_new, dtm)


#--- canopy height model 1 m -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(file.path(root, "chm1"), "/{*}_chm1")
opt_filter(ctg_norm) = "-drop_z_above 50"

chm1 <- rasterize_canopy(ctg_norm, res = 1, algorithm = pitfree())
terra::writeRaster(chm1, chm1_file)


#--- canopy height model -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(p_chm, "/{*}_chm")
opt_filter(ctg_norm) = "-drop_z_above 50"

# opt_restart(ctg_norm) <- 73
# chm <- rasterize_canopy(ctg_norm, res = 0.5, algorithm = p2r(subcircle = 0.15, na.fill = tin()))

chm <- rasterize_canopy(ctg_norm, res = 0.5, algorithm = pitfree())
terra::writeRaster(chm, chm_file)


#--- intensity raster -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(p_int, "/{*}_int")
opt_filter(ctg_norm) = "-drop_z_above 50"

intensity <- pixel_metrics(ctg_norm, func = ~mean(Intensity), res=0.5, filter = ~ReturnNumber == 1L)
intensity <- pixel_metrics(ctg_norm, func = ~mean(Intensity), res=1)
terra::writeRaster(intensity, int_file)

#--- intensity raster -----------

ctg <- readLAScatalog(p_las)
opt_filter(ctg) = "-keep_class 2 3 -drop_scan_angle_below 75 -drop_scan_angle_above 105"

# sensor <- track_sensor(ctg, Roussel2020())
# write_sf(sensor, file.path(root, "sensor_track.gpkg"))
sensor <- read_sf(file.path(root, "sensor_track.gpkg"))

opt_output_files(ctg) <-  paste0(file.path(root, "rangecor"), "/{*}_rangecor")
ctg@data <- ctg@data[ctg@data$Number.of.point.records > 10e+06, ]
ctg_rangec <- normalize_intensity(ctg, range_correction(sensor, Rs = 2000))


fns = list.files(p_las, "las$")
for (fn in fns) {
  message(fn)
  ofn <- file.path(root, "rangecor", fn)
  if (!file.exists(ofn)) {
    las <- readLAS(file.path(p_las, fn), filter="-keep_class 2 3 -drop_scan_angle_below 75 -drop_scan_angle_above 105")
    las_rc <- try(normalize_intensity(las, range_correction(sensor, Rs = 2000)))
    if (class(las_rc)=="LAS") {
      writeLAS(las_rc, ofn, index=T)
    }
  }
}


ctg <- readLAScatalog(file.path(root, "rangecor"))
opt_output_files(ctg) <-  paste0(file.path(root, "intensity_rangecor"), "/{*}_int")
opt_filter(ctg) = "-keep_class 2 3 -drop_scan_angle_below 75 -drop_scan_angle_above 105"

intensity <- pixel_metrics(ctg, func = ~mean(Intensity), res=0.5, filter = ~ReturnNumber == 1L)
# intensity <- pixel_metrics(ctg, func = ~mean(Intensity), res=1)
terra::writeRaster(intensity, file.path(root, "berchtesgaden_2021_lidar_intensity_rangecor.tif"))



