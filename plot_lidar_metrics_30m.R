############################################################
# plotting derived lidar metrics
###########################################################

library(dplyr)
library(terra)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)


#--- load lidar metrics --------------

metrics_name.list <-  c("zmax","zmean", "zsd","zskew", "zkurt" , "zentropy" ,  "pzabovezmean", "pzabove2","zq5" , "zq10",
                        "zq15", "zq20","zq25", "zq30", "zq35", "zq40", "zq45","zq50" ,"zq55","zq60" ,"zq65","zq70",
                        "zq75","zq80", "zq85","zq90", "zq95", "zpcum1", "zpcum2", "zpcum3", "zpcum4", "zpcum5", "zpcum6",
                        "zpcum7", "zpcum8","zpcum9") #layer names taken from original output

#2009
metrics2009 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_metrics_crop.tif") 
names(metrics2009) <- metrics_name.list

#2017
metrics2017 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_metrics_crop.tif") 
names(metrics2017) <- metrics_name.list

#2017
metrics2021 <- raster::stack("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_metrics_30.tif") 
names(metrics2021) <- metrics_name.list

#DEM 
DEM2017 <- raster::raster("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_dtm.tif")

#focus sites
focus_sites <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/focus_sites_buffer.shp")
focus_sites <-as(focus_sites, "Spatial")

# ----- plot metrics

zmean9 <- levelplot(metrics2009[["zmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 9 zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))
zmean17 <- levelplot(metrics2017[["zmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 17 zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))
zmean21 <- levelplot(metrics2021[["zmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 21 zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))

pzabovezmean9 <- levelplot(metrics2009[["pzabovezmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 9 percentage z above zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))
pzabovezmean17 <-levelplot(metrics2017[["pzabovezmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 17 percentage z above zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))
pzabovezmean21 <-levelplot(metrics2021[["pzabovezmean"]], margin = FALSE, par.setting =magmaTheme(), main= list("CHM 21 percentage z above zmean 30m")) + layer(sp.polygons(focus_sites, col='green'))
