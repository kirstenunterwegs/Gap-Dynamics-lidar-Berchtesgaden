########################################################
#
# Script to generate Canopy Height Model from Lidar data
# for the Berchtesgaden National Park (year 2009)
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


# --- paths ---

root <- "i:/Fonda/workspace/berchtesgaden/lidar/2009"
root <- "f:/Spring/Berchtesgaden/lidar/2009"
p_las <- file.path(root, "las")
p_norm <- file.path(root, "las_norm")
p_dtm <- file.path(root, "dtm")
p_chm <- file.path(root, "chm")
dtm_file <- file.path(root, "berchtesgaden_2009_dtm.tif")
chm_file <- file.path(root, "berchtesgaden_2009_chm.tif")
density_file <- file.path(root, "berchtesgaden_2009_lidar_density.tif")


#--- setup future -----------------

plan(multisession, workers = 20L)


#--- density ----------------------

ctg <- readLAScatalog(p_las)
dens <- rasterize_density(ctg, 10)

terra::writeRaster(dens, density_file)


#--- terrain ----------------------

ctg <- readLAScatalog(p_las)
opt_output_files(ctg) <- paste0(p_dtm, "/{*}_dtm")

dtm <- rasterize_terrain(ctg, res = 0.5, algorithm = tin())

terra::writeRaster(dtm, dtm_file)


#--- normalize height --------------

ctg <- readLAScatalog(p_las)
dtm <- raster::raster(dtm_file)

opt_output_files(ctg) <-  paste0(p_norm, "/{*}_norm")
ctg_norm <- normalize_height(ctg, dtm)


#--- canopy height model 1m -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(file.path(root, "chm1"), "/{*}_chm")
opt_filter(ctg_norm) = "-drop_z_above 50"

# chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = p2r(subcircle = 0.1, na.fill = knnidw()))
chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = pitfree(subcircle=0.2))

terra::writeRaster(chm, file.path(root, "berchtesgaden_2009_chm_1m.tif"))


#--- canopy height model -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(p_chm, "/{*}_chm")
opt_filter(ctg_norm) = "-drop_z_above 50"

# chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = p2r(subcircle = 0.1, na.fill = knnidw()))
chm <- rasterize_canopy(ctg_norm, res = 0.5, algorithm = pitfree(subcircle=0.2))

terra::writeRaster(chm, chm_file)


#--- debugging -----------

chunk <- readRDS("C:/Users/pflugmad/AppData/Local/Temp/6/RtmpYlrjGP/chunk35.rds")
las <- readLAS(chunk)
las_check(las)

las$Z[las$Z < 0] <- 0
las_check(las)
las2 <- filter_duplicates(las)
chm <- rasterize_canopy(las, res = 1, algorithm = p2r(subcircle = 0.1, na.fill = knnidw()))

