# === Evaluating lidar metrics to learn about forest structure ===

library(lidR)
library(sf)
library(dplyr)
library(terra)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)


#--- load lidar metrics --------------

metrics_name.list <-  c("zmax","zmean", "zsd","zskew", "zkurt" , "zentropy" ,  "pzabovezmean", "pzabove2","zq5" , "zq10",
                        "zq15", "zq20","zq25", "zq30", "zq35", "zq40", "zq45","zq50" ,"zq55","zq60" ,"zq65","zq70",
                        "zq75","zq80", "zq85","zq90", "zq95", "zpcum1", "zpcum2", "zpcum3", "zpcum4", "zpcum5", "zpcum6",
                        "zpcum7", "zpcum8","zpcum9") #layer names taken from original output

#2009
metrics2009 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_metrics_crop.tif") 
names(metrics2009) <- metrics_name.list
plot(metrics2009[[c("zmean","zsd", "pzabovezmean", "pzabove2","zq25", "zq50","zq75")]])

#2017
metrics2017 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_metrics_crop.tif") 
names(metrics2017) <- metrics_name.list
plot(metrics2017[[c("zmean","zsd", "pzabovezmean", "pzabove2","zq25", "zq50","zq75")]])

#2017
metrics2021 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_metrics_30.tif") 
names(metrics2021) <- metrics_name.list
plot(metrics2021[[c("zmean","zsd", "pzabovezmean", "pzabove2","zq25", "zq50","zq75")]])

#DEM 
DEM2017 <- raster::raster("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_dtm.tif")

#--- adjust 2009 and 2017 stacks to 2021 extent ----------

#extent <- extent(metrics2021)
#metrics2009_crop <- crop(metrics2009, extent)
#metrics2017_crop <- crop(metrics2017, extent)

#metrics2009_crop[is.na(metrics2021)] <- NA 
#metrics2017_crop[is.na(metrics2021)] <- NA

#par(mfrow = c(1,3))

#plot(metrics2021[["zmean"]], main='2021 zmean')
#plot(metrics2017_crop[["zmean"]], main='2017 zmean')
#plot(metrics2009_crop[["zmean"]], main='2009 zmean')

#writeRaster(metrics2009_crop,  "i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_metrics_crop.tif")
#writeRaster(metrics2017_crop,  "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_metrics_crop.tif")


#--- reclassify DEM ---
#assign the DEM elevation values to class intervals

m = c(500, 900, 1, 900, 1300, 2, 1300, 1800, 3, 1800, 3000, 4)# Generate a reclass matrix (begin, end, value ...)
rclmat = matrix(m, ncol=3, byrow=TRUE)
dem <- rast(DEM2017)
reclassified_dem = terra::classify(dem, rclmat, include.lowest=TRUE) # Reclassify the raster layer
terra::writeRaster(reclassified_dem, "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_classified_dtm.tif")

#--- reclassify Metrics for 2017---
# not using 2009, as only first and last return, which makes metrics less comparable with 2017 and 2021

#zmean
zmean2017 <- raster::raster(metrics2017_crop, layer="zmean") #extract metrics layer
m_mean = c(0, 8, 1, 8, 16, 2, 16, 23, 3)# Generate a reclass matrix (begin, end, value ...)
clmat_mean = matrix(m_mean, ncol=3, byrow=TRUE)
zmean2017.rast <- rast(zmean2017)
reclassified_mean = terra::classify(zmean2017.rast, clmat_mean, include.lowest=TRUE) 
terra::writeRaster(reclassified_mean, "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_mean_classified.tif")

#zsd
zsd2017 <- raster::raster(metrics2017_crop, layer="zsd") #extract metrics layer
m_zsd = c(0, 4, 1, 4, 8, 2, 8, 13, 3)# Generate a reclass matrix (begin, end, value ...)
clmat_zsd = matrix(m_zsd, ncol=3, byrow=TRUE)
zsd2017.rast <- rast(zsd2017)
reclassified_zsd = terra::classify(zsd2017.rast, clmat_zsd, include.lowest=TRUE) 
terra::writeRaster(reclassified_zsd, "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_zsd_classified.tif")

#pzabovezmean
pzabovezmean2017 <- raster::raster(metrics2017_crop, layer="pzabovezmean") #extract metrics layer
m_pzabovezmean = c(0, 40, 1, 40, 60, 2, 60, 70, 3)# Generate a reclass matrix (begin, end, value ...)
clmat_pzabovezmean = matrix(m_pzabovezmean, ncol=3, byrow=TRUE)
pzabovezmean2017.rast <- rast(pzabovezmean2017)
reclassified_pzabovezmean = terra::classify(pzabovezmean2017.rast, clmat_pzabovezmean, include.lowest=TRUE) 
terra::writeRaster(reclassified_pzabovezmean, "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_pzabovezmean_classified.tif")

#--- plotting different metrics ---

par(mfrow = c(3,4))

plot(metrics2021[[c("zmean")]], main='2021 zmean')
plot(metrics2021[[c("zsd")]], main='2021 zsd')
plot(metrics2021[[c("pzabovezmean")]], main='2021 pzabovezmean')
plot(metrics2021[[c("pzabove2")]], main='2021 pzabove2')
plot(metrics2017[["zmean"]], main='2017 zmean')
plot(metrics2017[[c("zsd")]], main='2017 zsd')
plot(metrics2017[[c("pzabovezmean")]], main='2017 pzabovezmean')
plot(metrics2017[[c("pzabove2")]], main='2017 pzabove2')
plot(metrics2009[["zmean"]], main='2009 zmean')
plot(metrics2009[[c("zsd")]], main='2009 zsd')
plot(metrics2009[[c("pzabovezmean")]], main='2009 pzabovezmean')
plot(metrics2009[[c("pzabove2")]], main='2009 pzabove2')

par(mfrow=c(3,2))
plot(metrics2021[[c("zmean")]], main='2021 zmean')
plot(metrics2021[[c("pzabovezmean")]], main='2021 pzabovezmean')
plot(metrics2017[["zmean"]], main='2017 zmean')
plot(metrics2017[[c("pzabovezmean")]], main='2017 pzabovezmean')
plot(metrics2009[["zmean"]], main='2009 zmean')
plot(metrics2009[[c("pzabovezmean")]], main='2009 pzabovezmean')

par(mfrow=c(1,1))
plot(metrics2021[[c("zmean")]], main='2021 zmean')
plot(metrics2017[["zmean"]], main='2017 zmean')
plot(metrics2009[["zmean"]], main='2009 zmean')

#------ aggregate CHM to 30m res ------
#chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_pf-sc02.tif")
#chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_pf-sc02.tif")
#chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_fix.tif")

#chm9_30 <- aggregate(chm9, fact=60, fun="mean")
#chm17_30 <- aggregate(chm17, fact=60, fun="mean")
#chm21_30 <- aggregate(chm21, fact=60, fun="mean")

#writeRaster(chm9_30,"i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_pf-sc02_30.tif" )
#writeRaster(chm17_30, "i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_pf-sc02_30.tif")
#writeRaster(chm21_30, "i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_fix_30.tif")

chm9_30 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_pf-sc02_30.tif")
chm17_30<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_pf-sc02_30.tif")
chm21_30 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_fix_30.tif")



focus_sites <-as(focus_sites, "Spatial")
levelplot(metrics2021[["zmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 21 aggregated 30m")) + layer(sp.polygons(focus_sites, col='green'))

