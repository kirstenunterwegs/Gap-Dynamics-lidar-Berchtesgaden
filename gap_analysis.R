######################################
# Analysing Gaps
#####################################

library(sf)
library(plyr)
library(dplyr)
library(dplyr)
library(terra)
library(ForestGapR)
library(ggplot2)
library(ForestTools)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)
library(landscapemetrics)
require(scales)
library(tidyverse)

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

# --- load CHM and Gap layers ----

chm_stack <- rast("chm_berchtesgaden_stack_1m_maskedartifacts.tif")
chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]

#chm17 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_1m.tif")

gaps2009 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/gaps/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/gaps/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/gaps/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop gaps to CHM 
gaps2009 <- crop(gaps2009, chm21, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, chm21, snap="near",mask=TRUE) 
#gaps2021 <-crop(gaps2021, chm21, snap="near",mask=TRUE) 


# --- load NP information 
wd <- "i:/Fonda/workspace/berchtesgaden/focus_sites/NP_data"
setwd(wd)

foresttype <- vect("forest_type_map.shp")
management <- vect("npb_zonierung_22_epsg25832.shp")
aspect<-  rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/aspect_2021_classified_1m.tif")
elevation <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_classified_200steps_dtm_1m.tif")

#resample elevation and aspect
#aspect_cre <- resample(aspect, chm21, method="bilinear")
#elevation_cre <- resample(elevation, chm21, method="bilinear")
#terra::writeRaster(aspect_cre, "i:/Fonda/workspace/berchtesgaden/lidar/2021/aspect_2021_classified_1m.tif",overwrite=TRUE)
#terra::writeRaster(elevation_cre, "i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_classified_200steps_dtm_1m.tif",overwrite=TRUE)

# -- rasterize NP information

#crate empty raster for rasterization
r_1 <- rast()
ext(r_1) <- ext(chm9)
terra::res(r_1) <- terra::res(chm9)  
terra::crs(r_1) <- terra::crs(chm9)

# rasterize forest type Information
ftype <- terra::rasterize(foresttype, r_1, field="type")

# rasterize Management Information
mtype <- terra::rasterize(management, r_1, field="zone_id")

# --- define functions -----

#gap_layer <- gap_stack_fs1$gaps_17_fs1
#chm_layer <- chm17_fs1

Gap_Stats <- function (gap_layer, chm_layer, year) 
 {  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

gap_polygon <- as.polygons(gap_layer)
gap_list$perimeter <- perim(gap_polygon)#add perimeter  
#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter  
print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "year")
return(gap_list)
}


# ------- management type

#function to assign management type
getManagementType <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, management) %>% #count pixels per mtype per gap
    summarize(count = n())
  #identify dominating management type in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one management type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several mtypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no management type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other mtype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(management = as.factor(recode(management,
                                              `2`="buffer zone",
                                              `4`="core zone")))
  return(xx)
}



Gap_Stats_management <- function (gap_layer, chm_layer, management, year) 
{  t <- Sys.time()
  gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, management), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "management")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

gap_polygon <- as.polygons(gap_layer)
gap_list$perimeter <- perim(gap_polygon)#add perimeter  
#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter  
print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$management <- as.data.frame(getManagementType(gap_chm))[,2]
print("get management type: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "management", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


# ---- forest type

#function to assign forest type
getForestType <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, ftype) %>% #count pixels per ftype per gap
    summarize(count = n())
  #identify dominating forest type in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  return(xx)
}


Gap_Stats_ftype <- function (gap_layer, chm_layer, foresttype, year) 
{  t <- Sys.time()
  gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, foresttype), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "ftype")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

gap_polygon <- as.polygons(gap_layer)
gap_list$perimeter <- perim(gap_polygon)#add perimeter  
#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter  
print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$ftype <- as.data.frame(getForestType(gap_chm))[,2]
print("get forest type: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "forest_type", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#---elevation 

#function to assign elevation class
getElevation <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, elevation) %>% #count pixels per elevation class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(elevation = as.factor(recode(elevation,
                                                   `1`="600-800",
                                                  `2`="800-1000",
                                                  `3`="1000-1200",
                                                  `4`="1200-1400",
                                                  `5`="1400-1600",
                                                  `6`="1600-1800",
                                                  `7`="1800-2000",
                                                   `8`="2000-2800")))
  return(xx)
}


Gap_Stats_elevation <- function (gap_layer, chm_layer, dtm_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, dtm_class), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "elevation")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

gap_polygon <- as.polygons(gap_layer)
gap_list$perimeter <- perim(gap_polygon)#add perimeter  
#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter  
print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$elevation <- as.data.frame(getElevation(gap_chm))[,2]
print("get elevation: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "elevation", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


# ---- aspect 

#function to assign elevation class
getAspect <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, aspect) %>% #count pixels per aspect class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one aspect type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several aspect classes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no aspect info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other aspect info assign that one to ID
    }
  }
  xx<- xx %>% mutate(aspect = as.factor(recode(aspect,
                                                  `1`="North",
                                                  `2`="East",
                                                  `3`="South",
                                                  `4`="West")))
  return(xx)
}


Gap_Stats_aspect <- function (gap_layer, chm_layer, ascpect_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, ascpect_class), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "aspect")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

gap_polygon <- as.polygons(gap_layer)
gap_list$perimeter <- perim(gap_polygon)#add perimeter  
#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter  
print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$aspect <- as.data.frame(getAspect(gap_chm))[,2]
print("get aspect: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "aspect", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}

#------- calculate gap stats per forest type
stats_ftype_2009<- Gap_Stats_ftype(gaps2009, chm9, ftype, "2009")
stats_ftype_2017 <- Gap_Stats_ftype(gaps2017, chm17, ftype, "2017")
stats_ftype_2021 <- Gap_Stats_ftype(gaps2021, chm21, ftype, "2021")
#saveRDS(stats_ftype_2021, "i:/Fonda/workspace/berchtesgaden/gaps/stats_ftype_2021.rds")

stats_all_ftype <- rbind(stats_ftype_2009, stats_ftype_2017, stats_ftype_2021)
stats_all_ftype$pa_ratio <- stats_all_ftype$perimeter/stats_all_ftype$gap_area #calculate perimeter/area ratio
stats_all_ftype$shape.index <- stats_all_ftype$perimeter/ (2*(pi*stats_all_ftype$gap_area)^0.5)
stats_all_ftype$gap_area_ha <- stats_all_ftype$gap_area/10000 #convert area to ha

saveRDS(stats_all_ftype, "i:/Fonda/workspace/berchtesgaden/gaps/stats_all_ftype.rds")
# change year label
#stats_all_ftype$year <- factor(stats_all_ftype$year, levels=c("chm9", "chm17", "chm21"), labels=c("2009", "2017", "2021"))

#------- calculate gap stats per management type
stats_mtype_2009<- Gap_Stats_management(gaps2009, chm9, mtype, "2009")
stats_mtype_2017 <- Gap_Stats_management(gaps2017, chm17, mtype, "2017")
stats_mtype_2021 <- Gap_Stats_management(gaps2021, chm21, mtype, "2021")


stats_all_mtype <- rbind(stats_mtype_2009, stats_mtype_2017, stats_mtype_2021)
stats_all_mtype$pa_ratio <- stats_all_mtype$perimeter/stats_all_mtype$gap_area #calculate perimeter/area ratio
stats_all_mtype$shape.index <- stats_all_mtype$perimeter/ (2*(pi*stats_all_mtype$gap_area)^0.5)
stats_all_mtype$gap_area_ha <- stats_all_mtype$gap_area/10000 #convert area to ha

saveRDS(stats_all_mtype, "i:/Fonda/workspace/berchtesgaden/gaps/stats_all_mtype.rds")

# change year label
#stats_all_mtype$year <- factor(stats_all_mtype$year, levels=c("chm9", "chm17", "chm21"), labels=c("2009", "2017", "2021"))

#stats_all_ftype[stats_all_ftype == ""] <- NA #replace empty fields for forest type (any info) with NA


#------- calculate gap stats per elevation
stats_elevation_2009<- Gap_Stats_elevation(gaps2009, chm9, elevation, "2009")
stats_elevatione_2017 <- Gap_Stats_elevation(gaps2017, chm17, elevation, "2017")
stats_elevation_2021 <- Gap_Stats_elevation(gaps2021, chm21, elevation, "2021")

stats_all_elevation <- rbind(stats_elevation_2009, stats_elevatione_2017, stats_elevation_2021)
stats_all_elevation$pa_ratio <- stats_all_elevation$perimeter/stats_all_elevation$gap_area #calculate perimeter/area ratio
stats_all_elevation$shape.index <- stats_all_elevation$perimeter/ (2*(pi*stats_all_elevation$gap_area)^0.5)
stats_all_elevation$gap_area_ha <- stats_all_elevation$gap_area/10000 #convert area to ha

# stats_all_elevation <- stats_all_elevation %>% mutate(elevation = recode(elevation,
#                                            `1800-1000`="1800-2000"))
stats_all_elevation$elevation <- factor(stats_all_elevation$elevation , levels=c("600-800", "800-1000", "1000-1200",
                                                                                 "1200-1400", "1400-1600", "1600-1800",
                                                                                 "1800-2000", "2000-2800"))

saveRDS(stats_all_elevation, "i:/Fonda/workspace/berchtesgaden/gaps/stats_all_elevation.rds")

#------- calculate gap stats per aspect
stats_aspect_2009<- Gap_Stats_aspect(gaps2009, chm9, aspect, "2009")
stats_aspect_2017 <- Gap_Stats_aspect(gaps2017, chm17, aspect, "2017")
stats_aspect_2021 <- Gap_Stats_aspect(gaps2021, chm21, aspect, "2021")

stats_all_aspect <- rbind(stats_aspect_2009, stats_aspect_2017, stats_aspect_2021)
stats_all_aspect$pa_ratio <- stats_all_aspect$perimeter/stats_all_aspect$gap_area #calculate perimeter/area ratio
stats_all_aspect$shape.index <- stats_all_aspect$perimeter/ (2*(pi*stats_all_aspect$gap_area)^0.5)
stats_all_aspect$gap_area_ha <- stats_all_aspect$gap_area/10000 #convert area to ha

saveRDS(stats_all_aspect, "i:/Fonda/workspace/berchtesgaden/gaps/stats_all_aspect.rds")

#-------calculate gap stats across whole Nationalpark

stats_2009<- Gap_Stats(gaps2009, chm9, "2009")
stats_2017 <- Gap_Stats(gaps2017, chm17,"2017")
stats_2021 <- Gap_Stats(gaps2021, chm21,"2021")

stats_all <- rbind(stats_2009, stats_2017, stats_2021)
stats_all$pa_ratio <- stats_all$perimeter/stats_all$gap_area #calculate perimeter/area ratio
stats_all$shape.index <- stats_all$perimeter/ (2*(pi*stats_all$gap_area)^0.5)
stats_all$gap_area_ha <- stats_all$gap_area/10000 #convert area to ha

saveRDS(stats_all, "i:/Fonda/workspace/berchtesgaden/gaps/stats_all_berchtesgaden.rds")

###########################################################################################################################################
###########################################################################################################################################

# load stats

stats_all <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_berchtesgaden.rds")
stats_all_aspect <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_aspect.rds")
stats_all_elevation <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_elevation.rds")
stats_all_ftype <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_ftype.rds")
stats_all_mtype <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_mtype.rds")

# --- plot gap size vs. total disturbance area ---(cumulative distribution function)

#subset elevation by excluding 2000-2800 elevation level, as there is only one stable, wrongly classified gap
stats_all_elevation <- stats_all_elevation[stats_all_elevation$elevation != "2000-2800", ]  

#--ftype
#relabel no information forest type
stats_all_ftype$forest_type [stats_all_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
stats_all_ftype$forest_type <- as.character(stats_all_ftype$forest_type)
stats_all_ftype[is.na(stats_all_ftype)] <- "no information" # exchange NAs <- "no information"
stats_all_ftype$forest_type <- as.factor(stats_all_ftype$forest_type)

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 19,angle = 45, hjust=1),
  axis.text.y = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=20),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"))
 # legend.position="bottom")

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_analysis/"
setwd(wd)

tiff("ecdf_aspect_gap.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_aspect, aes(x=gap_area_ha, colour = aspect)) +
  stat_ecdf(geom = "step", pad = FALSE, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_minimal() +
  scale_color_brewer(palette="Dark2") + My_Theme + facet_wrap(~year) +
  labs(x= "gap area [ha]", y= "cumulative % of disturbed area")
dev.off()

tiff("ecdf_ftype_gap.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_ftype, aes(x=gap_area_ha, colour = forest_type)) +
  stat_ecdf(geom = "step", pad = FALSE, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_minimal() +
  scale_color_brewer(palette="Dark2") + My_Theme + facet_wrap(~year) +
  labs(x= "gap area [ha]", y= "cumulative % of disturbed area", color= "forest type")
dev.off()

tiff("ecdf_elevation_gap.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_elevation, aes(x=gap_area_ha, colour = elevation)) +
  stat_ecdf(geom = "step", pad = FALSE, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_minimal() +
  scale_color_brewer(palette="GnBu") + My_Theme + facet_wrap(~year) +
  labs(x= "gap area [ha]", y= "cumulative % of disturbed area", color = "elevation [m]")
dev.off()

tiff("ecdf_mtype_gap.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_mtype, aes(x=gap_area_ha, colour = management)) +
  stat_ecdf(geom = "step", pad = FALSE, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_minimal() +
  scale_color_brewer(palette="Dark2") + My_Theme + facet_wrap(~year) +
  labs(x= "gap area [ha]", y= "cumulative % of disturbed area")
dev.off()

# --- meta stats ---

#reorder elevation labels
# stats_all_elevation$elevation <- factor(stats_all_elevation$elevation , levels=c("600-800", "800-1000", "1000-1200",
#                                                                                  "1200-1400", "1400-1600", "1600-1800",
#                                                                                  "1800-2000", "2000-2800"))
#stats overview

summary_stats_all <- stats_all %>% group_by(year)%>% 
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = sum(gap_area),
            max_gap_size = max(gap_area),
            mean_gap_area = mean(gap_area),
            max_perimeter = max(perimeter),
            mean_perimeter = mean(perimeter),
            min_perimeter = min(perimeter),
            max_pa = max(pa_ratio),
            mean_pa = mean(pa_ratio),
            min_pa = min(pa_ratio))

summary_stats_all_t <- t(summary_stats_all)

summary_stats_mtype <- stats_all_mtype %>% group_by(year, management)%>% 
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = sum(gap_area),
            max_gap_size = max(gap_area),
            mean_gap_area = mean(gap_area),
            max_perimeter = max(perimeter),
            mean_perimeter = mean(perimeter),
            min_perimeter = min(perimeter),
            max_pa = max(pa_ratio),
            mean_pa = mean(pa_ratio),
            min_pa = min(pa_ratio))
summary_stats_mtype_t <- t(summary_stats_mtype)

# summary_stats_mtype <- stats_all_mtype %>% group_by(year, management)%>% 
#   summarise(n_gaps = n_distinct(gap_id),
#             gap_area_total = round(sum(gap_area),2),
#             max_gap_size = round(max(gap_area),2),
#             mean_gap_area = round(mean(gap_area),2),
#             max_perimeter = round(max(perimeter),2),
#             mean_perimeter = round(mean(perimeter),2),
#             min_perimeter = round(min(perimeter),2),
#             max_pa = round(max(pa_ratio),2),
#             mean_pa = round(mean(pa_ratio),2),
#             min_pa = round(min(pa_ratio),2))
# summary_stats_mtype$gap_area_total <- round(summary_stats_mtype$gap_area_total/10000,2)
# summary_stats_mtype$mean_gap_area <- round(summary_stats_mtype$mean_gap_area/10000,2)
# summary_stats_mtype$max_gap_size <- round(summary_stats_mtype$max_gap_size/10000,2)
# summary_stats_mtype.sub <- summary_stats_mtype[,-c(7:9)]
# xtable(summary_stats_mtype.sub)



summary_stats_ftype <- stats_all_ftype %>% group_by(year, forest_type)%>% 
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = sum(gap_area),
            max_gap_size = max(gap_area),
            mean_gap_area = mean(gap_area),
            max_perimeter = max(perimeter),
            mean_perimeter = mean(perimeter),
            min_perimeter = min(perimeter),
            max_pa = max(pa_ratio),
            mean_pa = mean(pa_ratio),
            min_pa = min(pa_ratio))

# summary_stats_ftype <- stats_all_ftype %>% group_by(year, forest_type)%>% 
#   summarise(n_gaps = n_distinct(gap_id),
#             gap_area_total = round(sum(gap_area),2),
#             max_gap_size = round(max(gap_area),2),
#             mean_gap_area = round(mean(gap_area),2),
#             max_perimeter = round(max(perimeter),2),
#             mean_perimeter = round(mean(perimeter),2),
#             min_perimeter = round(min(perimeter),2),
#             max_pa = round(max(pa_ratio),2),
#             mean_pa = round(mean(pa_ratio),2),
#             min_pa = round(min(pa_ratio),2))
# summary_stats_ftype$gap_area_total <- round(summary_stats_ftype$gap_area_total/10000,2)
# summary_stats_ftype$mean_gap_area <- round(summary_stats_ftype$mean_gap_area/10000,2)
# summary_stats_ftype$max_gap_size <- round(summary_stats_ftype$max_gap_size/10000,2)
# summary_stats_ftype.sub <- summary_stats_ftype[,-c(7:9)]
# xtable(summary_stats_ftype.sub)


summary_stats_elevation <- stats_all_elevation %>% group_by(year, elevation)%>%
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = sum(gap_area),
            max_gap_size = max(gap_area),
            mean_gap_area = mean(gap_area),
            max_perimeter = max(perimeter),
            mean_perimeter = mean(perimeter),
            min_perimeter = min(perimeter),
            max_pa = max(pa_ratio),
            mean_pa = mean(pa_ratio),
            min_pa = min(pa_ratio))

# summary_stats_elevation <- stats_all_elevation %>% group_by(year, elevation)%>%
#   summarise(n_gaps = n_distinct(gap_id),
#             gap_area_total = round(sum(gap_area),2),
#             max_gap_size = round(max(gap_area),2),
#             mean_gap_area = round(mean(gap_area),2),
#             max_perimeter = round(max(perimeter),2),
#             mean_perimeter = round(mean(perimeter),2),
#             min_perimeter = round(min(perimeter),2),
#             max_pa = round(max(pa_ratio),2),
#             mean_pa = round(mean(pa_ratio),2),
#             min_pa = round(min(pa_ratio),2))
# summary_stats_elevation$gap_area_total <- round(summary_stats_elevation$gap_area_total/10000,2)
# summary_stats_elevation$mean_gap_area <- round(summary_stats_elevation$mean_gap_area/10000,2)
# summary_stats_elevation$max_gap_size <- round(summary_stats_elevation$max_gap_size/10000,2)
# summary_stats_elevation.sub <- summary_stats_elevation[,-c(7:9)]
# xtable(summary_stats_elevation.sub)


summary_stats_aspect <- stats_all_aspect %>% group_by(year, aspect)%>% 
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = sum(gap_area),
            max_gap_size = max(gap_area),
            mean_gap_area = mean(gap_area),
            max_perimeter = max(perimeter),
            mean_perimeter = mean(perimeter),
            min_perimeter = min(perimeter),
            max_pa = max(pa_ratio),
            mean_pa = mean(pa_ratio),
            min_pa = min(pa_ratio))

# summary_stats_aspect <- stats_all_aspect %>% group_by(year, aspect)%>% 
#   summarise(n_gaps = n_distinct(gap_id),
#             gap_area_total = round(sum(gap_area),2),
#             max_gap_size = round(max(gap_area),2),
#             mean_gap_area = round(mean(gap_area),2),
#             max_perimeter = round(max(perimeter),2),
#             mean_perimeter = round(mean(perimeter),2),
#             min_perimeter = round(min(perimeter),2),
#             max_pa = round(max(pa_ratio),2),
#             mean_pa = round(mean(pa_ratio),2),
#             min_pa = round(min(pa_ratio),2))
# summary_stats_aspect$gap_area_total <- round(summary_stats_aspect$gap_area_total/10000,2)
# summary_stats_aspect$mean_gap_area <- round(summary_stats_aspect$mean_gap_area/10000,2)
# summary_stats_aspect$max_gap_size <- round(summary_stats_aspect$max_gap_size/10000,2)
# summary_stats_aspect.sub <- summary_stats_aspect[,-c(7:9)]
# xtable(summary_stats_aspect.sub)

#------ plotting summary stats n vs. area----

# cumulative distribution functions (ECDFs)
ggplot(stats_all, aes(x=gap_area_ha, colour= year)) + stat_ecdf() +theme_classic()

ggplot(stats_all_aspect, aes(x=gap_area_ha, colour= aspect)) + stat_ecdf() + facet_grid(~year)

ggplot(stats_all_elevation, aes(x=gap_area_ha, colour= elevation)) + stat_ecdf() + facet_grid(~year)

ggplot(stats_all_mtype, aes(x=gap_area_ha, colour= management)) + stat_ecdf() + facet_grid(~year)

ggplot(stats_all_ftype, aes(x=gap_area_ha, colour= forest_type)) + stat_ecdf() + facet_grid(~year)

#prepare dfs for plotting

###############################!!!!!!!!!!!!!!!!!!!!! merge with area share information
area_shares <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/area_shares_subcategory_NP.rds")

#--aspect
summary_aspect_sub <- summary_stats_aspect[,c("year", "aspect", "n_gaps", "gap_area_total")]
summary_aspect_sub$gap_area_total <- summary_aspect_sub$gap_area_total/10000 #convert to ha
summary_aspect_sub_long <- gather(summary_aspect_sub,category, value, -c("year", "aspect"))

summary_aspect_sub_long<- summary_aspect_sub_long %>%mutate(
  value = round(value,0),
  category =  recode(category,
                    `gap_area_total`="total gap area [ha]",
                    `n_gaps`="number of gaps"),
  plot_value = ifelse(category=="number of gaps", value*(-1),
                      value*1)) 

#merge with area share information
summary_aspect_sub_long <- merge(summary_aspect_sub_long, area_shares[,c("individual", "area_share_subcategory")],
                                 by.x = "aspect", by.y = "individual", all.x = TRUE)


#scale gap area and number by area share of subcategory
summary_aspect_sub_long$total.area <- 8139 #add total research area in ha
#get area per subcaterory and get factor to scale gap size and area towards 1000 ha
summary_aspect_sub_long$area.scaling.factor <- (summary_aspect_sub_long$total.area*summary_aspect_sub_long$area_share_subcategory)/1000 

summary_aspect_sub_long$value_ascaled <- round(summary_aspect_sub_long$value / summary_aspect_sub_long$area.scaling.factor,0)
summary_aspect_sub_long$plot_value_ascaled <- summary_aspect_sub_long$plot_value / summary_aspect_sub_long$area.scaling.factor

#--elevation
summary_elevation_sub <- summary_stats_elevation[,c("year", "elevation", "n_gaps", "gap_area_total")]
summary_elevation_sub$gap_area_total <- summary_elevation_sub$gap_area_total/10000
summary_elevation_sub_long <- gather(summary_elevation_sub,category, value, -c("year", "elevation"))

summary_elevation_sub_long<- summary_elevation_sub_long %>%mutate(
  value = round(value,0),
  category =  recode(category,
                     `gap_area_total`="total gap area [ha]",
                     `n_gaps`="number of gaps"),
  plot_value = ifelse(category=="number of gaps", value*(-1),
                      value*1)) 

# there is one gap in 2000-2800 class, as closed forest mask identifies it as forest, but actually there are only
# rocks and gap is in same position

#merge with area share information
summary_elevation_sub_long <- merge(summary_elevation_sub_long, area_shares[,c("individual", "area_share_subcategory")],
                                 by.x = "elevation", by.y = "individual", all.x = TRUE)

#scale gap area and number by area share of subcategory
summary_elevation_sub_long$total.area <- 8139 #add total research area in ha
#get area per subcaterory and get factor to scale gap size and area towards 1000 ha
summary_elevation_sub_long$area.scaling.factor <- (summary_elevation_sub_long$total.area*summary_elevation_sub_long$area_share_subcategory)/1000

summary_elevation_sub_long$value_ascaled <- round(summary_elevation_sub_long$value / summary_elevation_sub_long$area.scaling.factor,0)
summary_elevation_sub_long$plot_value_ascaled <- summary_elevation_sub_long$plot_value / summary_elevation_sub_long$area.scaling.factor

#subset elevation by excluding 2000-2800 elevation level, as there is only one stable, wrongly classified gap
summary_elevation_sub_long <- summary_elevation_sub_long[summary_elevation_sub_long$elevation != "2000-2800", ]  

#--ftype
#relabel no information forest type
summary_stats_ftype$forest_type [summary_stats_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
summary_stats_ftype$forest_type <- as.character(summary_stats_ftype$forest_type)
summary_stats_ftype[is.na(summary_stats_ftype)] <- "no information" # exchange NAs <- "no information"
summary_stats_ftype$forest_type <- as.factor(summary_stats_ftype$forest_type)
#merge no information labels
summary_stats_ftype <- ddply(summary_stats_ftype,c("forest_type","year"),numcolwise(sum))


summary_ftype_sub <- summary_stats_ftype[,c("year", "forest_type", "n_gaps", "gap_area_total")]
summary_ftype_sub$gap_area_total <- summary_ftype_sub$gap_area_total/10000

#modify NA values and empty category outside of R
summary_ftype_sub <- as.data.frame(summary_ftype_sub)
summary_ftype_sub$gap_area_total <- round(summary_ftype_sub$gap_area_total,0)

# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! not needed when code above runs!
# #modify forest type labels in excel
# write.csv(summary_ftype_sub,"i:/Fonda/workspace/berchtesgaden/gaps/summary_ftype_sub.csv" )
# # modify table in excel by combing empty forest_type with NA
# summary_ftype_sub <- read.csv("i:/Fonda/workspace/berchtesgaden/gaps/summary_ftype_sub.csv", sep=";")
# summary_ftype_sub <- summary_ftype_sub[,-1]

summary_ftype_sub_long <- gather(summary_ftype_sub,category, value, -c("year", "forest_type"))

summary_ftype_sub_long<- summary_ftype_sub_long %>%mutate(
  value = round(value,0),
  category =  recode(category,
                     `gap_area_total`="total gap area [ha]",
                     `n_gaps`="number of gaps"),
  plot_value = ifelse(category=="number of gaps", value*(-1),
                      value*1)) 

#merge with area share information
summary_ftype_sub_long <- merge(summary_ftype_sub_long, area_shares[,c("individual", "area_share_subcategory")],
                                    by.x = "forest_type", by.y = "individual", all.x = TRUE)

#scale gap area and number by area share of subcategory
summary_ftype_sub_long$total.area <- 8139 #add total research area in ha
#get area per subcaterory and get factor to scale gap size and area towards 1000 ha
summary_ftype_sub_long$area.scaling.factor <- (summary_ftype_sub_long$total.area*summary_ftype_sub_long$area_share_subcategory)/1000

summary_ftype_sub_long$value_ascaled <- round(summary_ftype_sub_long$value / summary_ftype_sub_long$area.scaling.factor,0)
summary_ftype_sub_long$plot_value_ascaled <- summary_ftype_sub_long$plot_value / summary_ftype_sub_long$area.scaling.factor



#--mtype
summary_mtype_sub <- summary_stats_mtype[,c("year", "management", "n_gaps", "gap_area_total")]
summary_mtype_sub$gap_area_total <- summary_mtype_sub$gap_area_total/10000
summary_mtype_sub_long <- gather(summary_mtype_sub,category, value, -c("year", "management"))

summary_mtype_sub_long<- summary_mtype_sub_long %>%mutate(
  value = round(value,0),
  category =  recode(category,
                     `gap_area_total`="total gap area [ha]",
                     `n_gaps`="number of gaps"),
  plot_value = ifelse(category=="number of gaps", value*(-1),
                      value*1)) 

#merge with area share information
summary_mtype_sub_long <- merge(summary_mtype_sub_long, area_shares[,c("individual", "area_share_subcategory")],
                                by.x = "management", by.y = "individual", all.x = TRUE)

#scale gap area and number by area share of subcategory
summary_mtype_sub_long$total.area <- 8139 #add total research area in ha
#get area per subcaterory and get factor to scale gap size and area towards 1000 ha
summary_mtype_sub_long$area.scaling.factor <- (summary_mtype_sub_long$total.area*summary_mtype_sub_long$area_share_subcategory)/1000

summary_mtype_sub_long$value_ascaled <- round(summary_mtype_sub_long$value / summary_mtype_sub_long$area.scaling.factor,0)
summary_mtype_sub_long$plot_value_ascaled <- summary_mtype_sub_long$plot_value / summary_mtype_sub_long$area.scaling.factor



#-----

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_analysis/"
setwd(wd)

Theme_area_n = theme(axis.text.x=element_blank(),  #remove axis labels
                    axis.ticks.x=element_blank(),
                    title = element_text(size = 18),
                    axis.title.x = element_blank(),
                    axis.text.y = element_text(size = 24),
                    axis.title.y = element_blank(),
                    legend.key.height = unit(1, 'cm'),
                    legend.title = element_text(size=17),
                    legend.text = element_text(size=24),
                    strip.text.x = element_text(size = 20),
                    legend.position="bottom")

# https://ggplot2.tidyverse.org/reference/geom_text.html
# zum spielen soft gecaded - Lars 
b <- "outward" # label position (Achtung durch Rotation nicht intuitiv)
c <- 90 # label rotation 

#aspect
tiff("n_vs_area_aspect.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_aspect_sub_long, aes(x = aspect,y = plot_value, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value), size = 4, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps    vs.    total gap area [ha]") 
dev.off()

tiff("n_vs_area_aspect_ascaled.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_aspect_sub_long, aes(x = aspect,y = plot_value_ascaled, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value_ascaled), size = 8, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps  vs.  total gap area [ha]", fill="") 
dev.off()

#elevation
tiff("n_vs_area_elevation.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_elevation_sub_long, aes(x = elevation,y = plot_value, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value), size = 8, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps    vs.    total gap area [ha]") 
dev.off()

tiff("n_vs_area_elevation_ascaled.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_elevation_sub_long, aes(x = elevation,y = plot_value_ascaled, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value_ascaled), size = 7, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps  vs.  total gap area [ha]", fill="") 
dev.off()

# forest type
tiff("n_vs_area_ftype.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_ftype_sub_long, aes(x = forest_type,y = plot_value, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value), size = 4, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps    vs.    total gap area [ha]") 
dev.off()

tiff("n_vs_area_ftype_ascaled.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_ftype_sub_long, aes(x = forest_type,y = plot_value_ascaled, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value_ascaled), size = 8, hjust = b, angle = c) +
  facet_wrap(~year)+
  labs( y = "number of gaps  vs.  total gap area [ha] ", fill="") 
dev.off()


#management
tiff("n_vs_area_mtype.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_mtype_sub_long, aes(x = management, y = plot_value, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value), size = 6, hjust = b, angle = c) +
  facet_grid(~year)+
  labs( y = "number of gaps    vs.    total gap area [ha]") 
dev.off()

tiff("n_vs_area_mtype_ascaled.tiff", units="in", width=12, height=8, res=300)  
ggplot(summary_mtype_sub_long, aes(x = management, y = plot_value_ascaled, fill=category)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_brewer(type = "seq",palette = 7) +
  theme_classic()+
  Theme_area_n +
  geom_text(aes(label = value_ascaled), size = 8, hjust = b, angle = c) +
  facet_grid(~year)+
  labs( y = "number of gaps  vs.  total gap area [ha]", fill="") 
dev.off()


#------ plotting summary stats single ----
#change ftype labels
summary_stats_ftype$forest_type [summary_stats_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
summary_stats_ftype$forest_type <- as.character(summary_stats_ftype$forest_type)
summary_stats_ftype[is.na(summary_stats_ftype)] <- "no information" # exchange NAs <- "no information"
summary_stats_ftype$forest_type <- as.factor(summary_stats_ftype$forest_type)



wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_analysis/"
setwd(wd)

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 17),
  axis.title.y = element_text(size = 22),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=17),
  legend.text = element_text(size=16))


# n -forest type
tiff("n_gaps_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_ftype, aes(x=forest_type, y=n_gaps, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme + labs( x="forest type", y="number of gaps")
dev.off()

# total gap forest type
tiff("totalarea_gaps_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_ftype, aes(x=forest_type, y=gap_area_total, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="forest type", y="gap area [m2]")
dev.off()

# n -aspect
tiff("n_gaps_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_aspect, aes(x=aspect, y=n_gaps, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="aspect", y="number of gaps")
dev.off()

# total gap aspect
tiff("totalarea_gaps_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_aspect, aes(x=aspect, y=gap_area_total, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="aspect", y="gap area [m2]")
dev.off()

#n -elevation
tiff("n_gaps_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_elevation, aes(x=elevation, y=n_gaps, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="elevation [m]", y="number of gaps")
dev.off()
# total gap area-elevation
tiff("totalarea_gaps_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_elevation, aes(x=elevation, y=gap_area_total, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2") + My_Theme+ labs( x="elevation [m]", y="gap area [m2]")
dev.off()


# n -management type
tiff("n_gaps_mtype.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_mtype, aes(x=management, y=n_gaps, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="management type", y="number of gaps")
dev.off()
# total gap managmenet type
tiff("totalarea_gaps_mtype.tiff", units="in", width=12, height=8, res=300)
ggplot(summary_stats_mtype, aes(x=management, y=gap_area_total, fill=year))+ 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme + My_Theme+ labs( x="management type", y="gap area [m2]")
dev.off()



# --- gaps numbers per sizes ---

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=20),
  strip.text.x = element_text(size = 20))

require(scales)

# as histogram
tiff("gaps_size_histograml.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all, aes(x=gap_area)) + geom_histogram(color = "black", fill = "white",bins = 50)+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + facet_grid(~year) + My_Theme +
  labs(x=expression ("gap area in log"(m^2)))
dev.off()



# subdivide area into bins 
stats_all_bins<-stats_all %>% 
  mutate(gap_area_bins = cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000)))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                            `(399,500]`="400-500",
                            `(500,1e+03]`="500-1000",
                            `(1e+03,5e+03]`="1000-5000",
                            `(5e+03,1e+04]`="5000-10,000",
                            `(1e+04,5e+04]`="10,000-50,000",
                            `(5e+04,4.5e+05]`="50,000-450,000")))
head(stats_all_bins,10)
tiff("gaps_size_bins_all.png", units="in", width=12, height=8, res=300)
ggplot(stats_all_bins, aes(y=gap.size)) + geom_bar() + theme_classic()+ 
  labs( y="gap size class [m2]")+facet_grid(~year) +My_Theme +
  geom_text(stat='count', aes(label=..count..), size=8, hjust = "inward")
dev.off()


stats_all_aspect_bins<-stats_all_aspect %>% 
    mutate(gap_area_bins = cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000)))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))
tiff("gaps_size_bins_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_aspect_bins, aes(y=gap.size, fill=aspect)) + geom_bar(color="black") + theme_classic()+ 
  scale_fill_brewer(palette = "Dark2") +facet_grid(~year)+ 
  labs( y="gap size class [m2]")+My_Theme
dev.off()
  
stats_all_elevation_bins<-stats_all_elevation %>% 
  mutate(gap_area_bins = cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000)))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))
tiff("gaps_size_bins_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_elevation_bins, aes(y=gap.size, fill=elevation)) + geom_bar(color="black") + theme_classic()+ 
  scale_fill_brewer(palette = "Dark2")+facet_grid(~year)+ 
  labs( y="gap size class [m2]") +My_Theme
dev.off()

stats_all_mtype_bins<-stats_all_mtype %>% 
  mutate(gap_area_bins = cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000)))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))
tiff("gaps_size_bins_mtype.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_mtype_bins, aes(y=gap.size, fill=management)) + geom_bar(color="black") + theme_classic()+ 
  scale_fill_brewer(palette = "Dark2")+facet_grid(~year)+ 
  labs( y="gap size class [m2]")+My_Theme
dev.off()

stats_all_ftype_bins<-stats_all_ftype %>% 
  mutate(gap_area_bins = cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000)))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))
tiff("gaps_size_bins_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_ftype_bins, aes(y=gap.size, fill=forest_type)) + geom_bar( color="black") + theme_classic()+ 
  scale_fill_brewer(palette = "Dark2")+facet_grid(~year)+ 
  labs( y="gap size class [m2]")+My_Theme
dev.off()
  


# --- plot gap characteristics according to forest type
#change ftype labels
stats_all_ftype$forest_type [stats_all_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
stats_all_ftype$forest_type <- as.character(stats_all_ftype$forest_type)
stats_all_ftype[is.na(stats_all_ftype)] <- "no information" # exchange NAs <- "no information"
stats_all_ftype$forest_type <- as.factor(stats_all_ftype$forest_type)

# #gap area below 0.5 ha
# stats_all_ftype %>% filter(gap_area_ha <= 0.5) %>% ggplot(aes(gap_area_ha, fill = year)) +  
#   #geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
#     geom_density(col="#FFBB00",size=0.2, alpha=0.5) +
#   labs(x="Area of gap size in ha", y="Density") + facet_wrap(~forest_type) + theme_minimal()
# #gap area above btw 0.5 and 2 ha
# stats_all_ftype %>% filter(gap_area_ha <= 2 & gap_area_ha >0.5) %>% ggplot(aes(gap_area_ha, fill = year)) +  
#   #geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
#   geom_density(col="#FFBB00",size=0.2, alpha=0.5) +
#   labs(x="Area of gap size in ha", y="Density") + facet_wrap(~forest_type) + theme_minimal()
# #gap area above 2 ha
# stats_all_ftype %>% filter(gap_area_ha > 2) %>% ggplot(aes(gap_area_ha, fill = year)) +  
#   #geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
#   geom_density(col="#FFBB00",size=0.2, alpha=0.5) +
#   labs(x="Area of gap size in ha", y="Density") + facet_wrap(~forest_type) + theme_minimal()



wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_analysis/"
setwd(wd)

# plot perimeter/area ratio
tiff("pa_ratio_management.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_mtype, aes(x=management, y=pa_ratio, fill=year)) + geom_boxplot() + 
  theme_minimal() + labs( x="management type", y="perimeter/area ratio") + 
  My_Theme  + scale_fill_brewer(palette = "Dark2")
dev.off()


tiff("pa_ratio_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_ftype, aes(x=forest_type, y=pa_ratio, fill=year)) + geom_boxplot() + 
  theme_minimal() + labs( x="forest type", y="perimeter/area ratio") + 
  My_Theme + scale_fill_brewer(palette = "Dark2")
dev.off()

tiff("pa_ratio_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_elevation, aes(x=elevation, y=pa_ratio, fill=year)) + geom_boxplot() + 
  theme_minimal() + labs( x="elevation class [m]", y="perimeter/area ratio") + 
  My_Theme + scale_fill_brewer(palette = "Dark2")
dev.off()

tiff("pa_ratio_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_all_aspect, aes(x=aspect, y=pa_ratio, fill=year)) + geom_boxplot() + 
  theme_minimal() + labs( x="aspect", y="perimeter/area ratio") + 
  My_Theme + scale_fill_brewer(palette = "Dark2")
dev.off()


### ---- check correlation between elevation and forest type per gap ----

stats_all_elevation <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_elevation.rds")
stats_all_ftype <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_ftype.rds")

#relabel no information forest type
stats_all_ftype$forest_type [stats_all_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
stats_all_ftype$forest_type <- as.character(stats_all_ftype$forest_type)
stats_all_ftype[is.na(stats_all_ftype)] <- "no information" # exchange NAs <- "no information"
stats_all_ftype$forest_type <- as.factor(stats_all_ftype$forest_type)

# merge elevation and forest type information
stats_elev_ftype <- merge(stats_all_elevation, stats_all_ftype[,c("gap_id","forest_type")], by= "gap_id", all.x=TRUE)
# delete information for gaps above 2000 m elevation 
stats_elev_ftype <- stats_elev_ftype[! stats_elev_ftype$elevation== "2000-2800",]  

#merge no information labels
#stats_elev_ftype <- ddply(stats_elev_ftype,c("forest_type","year"),numcolwise(sum))

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 22),
  axis.text.y = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=20),
  strip.text.x = element_text(size = 20))

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_analysis/"
setwd(wd)

tiff("gaps_elevation_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(stats_elev_ftype, aes(y=elevation,  fill=forest_type)) + geom_bar( width=.5, position = "dodge") +
  facet_grid(~year) +theme_minimal() + My_Theme + scale_fill_brewer(palette = "Paired") + labs(fill='forest type') 
dev.off()
