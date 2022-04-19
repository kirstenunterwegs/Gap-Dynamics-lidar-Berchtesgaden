################################################################
# creating metrics of lidar to identify forest structure
############################################################

library(lidR)
library(sf)
library(dplyr)
library(terra)
require(future)


#--- load LAS data and derive metrics for every year --------------
#--- Grid aggregated - 200m red --------------

#2009
root <- "i:/Fonda/workspace/berchtesgaden/lidar/2009"
p_norm9 <- file.path(root, "las_norm")
metrics_file9 <- file.path(root, "berchtesgaden_2009_metrics_30.tif")

ctg9 <- readLAScatalog(p_norm9, filter = "-drop_z_below 0 -drop_z_above 55")

metrics9 <- grid_metrics(ctg9, .stdmetrics_z, res = 30)
terra::writeRaster(metrics9, metrics_file9)

#metrics <- raster::as.data.frame(metrics, xy=TRUE, na.rm = TRUE)
#write.csv(metrics,"i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgarden_2009_metrics.csv")


#2017
root <- "i:/Fonda/workspace/berchtesgaden/lidar/2017"
p_norm17 <- file.path(root, "las_norm")
metrics_file17 <- file.path(root, "berchtesgaden_2017_metrics_30.tif")

ctg17 <- readLAScatalog(p_norm17, filter = "-drop_z_below 0 -drop_z_above 55")

metrics17 <- grid_metrics(ctg17, .stdmetrics_z, res = 30)
terra::writeRaster(metrics17, metrics_file17)

#2021
root <- "i:/Fonda/workspace/berchtesgaden/lidar/2021"
p_norm21 <- file.path(root, "las_norm")
metrics_file21 <- file.path(root, "berchtesgaden_2021_metrics_30.tif")

ctg21 <- readLAScatalog(p_norm21, filter = "-drop_z_below 0 -drop_z_above 55")

metrics21 <- grid_metrics(ctg21, .stdmetrics_z, res = 30)
terra::writeRaster(metrics21, metrics_file21)

#--- option to calculate metrics on plot level with LAS catalouge
