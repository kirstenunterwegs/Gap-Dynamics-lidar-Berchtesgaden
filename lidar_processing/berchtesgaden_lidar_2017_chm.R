########################################################
#
# Script to generate Canopy Height Model from Lidar data
# for the Berchtesgaden National Park (year 2017)
#
# author: Dirk Pflugmacher (Humboldt Universität zu Berlin)
# last revised: 20.10.2023
#
#######################################################

# --- libaries ----

library(lidR)
library(sf)
library(dplyr)
library(terra)
require(future)

#--- paths --------------

root <- "i:/Fonda/workspace/berchtesgaden/lidar/2017"
root <- "f:/Spring/Berchtesgaden/lidar/2017"
p_las <- file.path(root, "las")
p_norm <- file.path(root, "las_norm")
p_dtm <- file.path(root, "dtm")
p_chm <- file.path(root, "chm")
dtm_file <- file.path(root, "berchtesgaden_2017_dtm.tif")
chm_file <- file.path(root, "berchtesgaden_2017_chm.tif")
density_file <- file.path(root, "berchtesgaden_2017_lidar_density.tif")


#--- setup future --------------

plan(multisession, workers = 20L)

#--- density ----------------------
ctg <- readLAScatalog(p_las)
dens <- rasterize_density(ctg, 10)

terra::writeRaster(dens, density_file)


#--- terrain-----------------------

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
opt_output_files(ctg_norm) <-  paste0(file.path(root, "chm1"), "/{*}_chm1")

chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = pitfree(subcircle=0.2))
terra::writeRaster(chm, file.path(root, "berchtesgaden_2017_chm1.tif"))


#--- canopy height model -----------

ctg_norm <- readLAScatalog(p_norm)
opt_output_files(ctg_norm) <-  paste0(p_chm, "/{*}_chm")

#opt_restart(ctg_norm) <- 73

# chm <- rasterize_canopy(ctg_norm, res = 0.5, algorithm = p2r(subcircle = 0.15, na.fill = tin()))
chm <- rasterize_canopy(ctg_norm, res = 0.5, algorithm = pitfree(subcircle=0.2))
terra::writeRaster(chm, chm_file)


#--- debugging -----------

missing <- c("790000_5282000", "789000_5282000", "790000_5281000", "790000_5283000")

chunk <- readRDS("C:/Users/pflugmad/AppData/Local/Temp/6/Rtmp2zQFUb/chunk3.rds")

las_norm <- readLAS(file.path(p_norm, "789000_5282000_norm.las"))# , filter="-first_only")
las_check(las_norm)
las_norm2 <- filter_duplicates(las_norm)
las_norm2 <- las_norm[!duplicated(las_norm$gpstime),]
chm1 <- rasterize_canopy(las_norm2, res = 1, algorithm = p2r(subcircle = 0.2, na.fill = tin()))

las_check(las_norm2)
table(las_norm2$ReturnNumber)

las3 <- readLAS(file.path(p_las, "785000_5274000.las"))
las_check(las3)



las <- readLAS(file.path(p_las, "789000_5282000.las"))
las_check(las)

track_sensor(las)

las2 <- LAS(las, check=T)
las_check(las2)

las0 <- readLAS("i:/Fonda/workspace/berchtesgaden/lidar/2017/las_epsg5678/4564_5277_all.las")
las_check(las0)
dups <- which(duplicated(las0$gpstime))
tmp <- las0[las0$gpstime==las0$gpstime[10],]
tmp@data
las_check(tmp)

table(las0$UserData)
table(las0$Classification)

table(las0$ReturnNumber)

las0 <- readLAS("i:/Fonda/workspace/berchtesgaden/lidar/2017/las_epsg5678/4563_5277_all.las")
las_check(las0)

chunk <- readRDS("C:/Users/pflugmad/AppData/Local/Temp/6/RtmpUb9ouy/chunk15.rds")
las <- readLAS(chunk)
las_check(las)
hist(las$Z)

las$Z[las$Z < 0] <- 0
las <- filter_duplicates(las)
chm0 <- rasterize_canopy(las, res = 1, algorithm = p2r(subcircle = 0.2))#, na.fill = tin()))


