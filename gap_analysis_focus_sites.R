######################################
# Analysing Gaps
#####################################

library(sf)
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

#wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
wd <- "~/Documents/global_change_geography/masterthesis/data/processed_data_from_server/data"
setwd(wd)

# --- load CHM and Gap layers ----

chm_names <- list("chm9_fs1", "chm17_fs1", "chm21_fs1", "chm9_fs2", "chm17_fs2", "chm21_fs2",
                  "chm9_fs3", "chm17_fs3", "chm21_fs3","chm9_fs4", "chm17_fs4", "chm21_fs4")

# load CHM croped to focus sites + 500m Buffer
chm_fs1_crop <- rast("chm_focus_site1_large_stack.tif")
chm9_fs1 <- chm_fs1_crop[[1]]
chm17_fs1<- chm_fs1_crop[[2]]
chm21_fs1<- chm_fs1_crop[[3]]
chm_fs2_crop <- rast("chm_focus_site2_large_stack.tif")
chm9_fs2 <- chm_fs2_crop[[1]]
chm17_fs2<- chm_fs2_crop[[2]]
chm21_fs2<- chm_fs2_crop[[3]]
chm_fs3_crop <- rast("chm_focus_site3_large_stack.tif")
chm9_fs3 <- chm_fs3_crop[[1]]
chm17_fs3<- chm_fs3_crop[[2]]
chm21_fs3<- chm_fs3_crop[[3]]
chm_fs4_crop <- rast("chm_focus_site4_large_stack.tif")
chm9_fs4 <- chm_fs4_crop[[1]]
chm17_fs4<- chm_fs4_crop[[2]]
chm21_fs4<- chm_fs4_crop[[3]]

chm_list <- list(chm9_fs1, chm17_fs1, chm21_fs1,
                 chm9_fs2, chm17_fs2, chm21_fs2,
                 chm9_fs3, chm17_fs3, chm21_fs3,
                 chm9_fs4, chm17_fs4, chm21_fs4)
names(chm_list) <- chm_names


# load gap layer stacks

# min 100
gap_stack_fs1 <- rast("gap_layers_fs1_100.tif")
gap_stack_fs2 <- rast("gap_layers_fs2_100.tif")
gap_stack_fs3 <- rast("gap_layers_fs3_100.tif")
gap_stack_fs4 <- rast("gap_layers_fs4_100.tif")


gap_list_100 <- list(gap_stack_fs1$gaps_9_fs1,gap_stack_fs1$gaps_17_fs1,gap_stack_fs1$gaps_21_fs1,
                     gap_stack_fs2$gaps_9_fs2,gap_stack_fs2$gaps_17_fs2,gap_stack_fs2$gaps_21_fs2,
                     gap_stack_fs3$gaps_9_fs3,gap_stack_fs3$gaps_17_fs3,gap_stack_fs3$gaps_21_fs3,
                     gap_stack_fs4$gaps_9_fs4,gap_stack_fs4$gaps_17_fs4,gap_stack_fs4$gaps_21_fs4 )
names(gap_list_100) <- chm_names



# min 250
gap_stack_fs1 <- rast("gap_layers_fs1_250.tif")
gap_stack_fs2 <- rast("gap_layers_fs2_250.tif")
gap_stack_fs3 <- rast("gap_layers_fs3_250.tif")
gap_stack_fs4 <- rast("gap_layers_fs4_250.tif")


gap_list_250 <- list(gap_stack_fs1$gaps_9_fs1,gap_stack_fs1$gaps_17_fs1,gap_stack_fs1$gaps_21_fs1,
                     gap_stack_fs2$gaps_9_fs2,gap_stack_fs2$gaps_17_fs2,gap_stack_fs2$gaps_21_fs2,
                     gap_stack_fs3$gaps_9_fs3,gap_stack_fs3$gaps_17_fs3,gap_stack_fs3$gaps_21_fs3,
                     gap_stack_fs4$gaps_9_fs4,gap_stack_fs4$gaps_17_fs4,gap_stack_fs4$gaps_21_fs4 )
names(gap_list_250) <- chm_names



# min 400
gap_stack_fs1 <- rast("gap_layers_fs1_400.tif")
gap_stack_fs2 <- rast("gap_layers_fs2_400.tif")
gap_stack_fs3 <- rast("gap_layers_fs3_400.tif")
gap_stack_fs4 <- rast("gap_layers_fs4_400.tif")


gap_list_400 <- list(gap_stack_fs1$gaps_9_fs1,gap_stack_fs1$gaps_17_fs1,gap_stack_fs1$gaps_21_fs1,
                     gap_stack_fs2$gaps_9_fs2,gap_stack_fs2$gaps_17_fs2,gap_stack_fs2$gaps_21_fs2,
                     gap_stack_fs3$gaps_9_fs3,gap_stack_fs3$gaps_17_fs3,gap_stack_fs3$gaps_21_fs3,
                     gap_stack_fs4$gaps_9_fs4,gap_stack_fs4$gaps_17_fs4,gap_stack_fs4$gaps_21_fs4 )
names(gap_list_400) <- chm_names



# load gap polygons

polygons_100 <- list()
polygons_100 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_100_", n ,"/","gaps_polygons_100_", n ,".shp", sep=""))
})
names(polygons_100) <- chm_names

polygons_250 <- list()
polygons_250 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_250_", n ,"/","gaps_polygons_250_", n ,".shp", sep=""))
})
names(polygons_250) <- chm_names

polygons_400 <- list()
polygons_400 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_400_", n ,"/","gaps_polygons_400_", n ,".shp", sep=""))
})
names(polygons_400) <- chm_names

# --- load NP information 
wd <- "~/Documents/global_change_geography/masterthesis/data/berchtesgaden_data"
setwd(wd)

foresttype <- vect("forest_type_map.shp")
management <- vect("npb_zonierung_22_epsg25832.shp")

# --- define functions -----

# rasterize NP information

r_1 <- rast()
ext(r_1) <- ext(chm9_fs1)
terra::res(r_1) <- terra::res(chm9_fs1)  
terra::crs(r_1) <- terra::crs(chm9_fs1)
r_2 <- rast()
ext(r_2) <- ext(chm9_fs2)
terra::res(r_2) <- terra::res(chm9_fs2)  
terra::crs(r_2) <- terra::crs(chm9_fs2)
r_3 <- rast()
ext(r_3) <- ext(chm9_fs3)
terra::res(r_3) <- terra::res(chm9_fs3)  
terra::crs(r_3) <- terra::crs(chm9_fs3)
r_4 <- rast()
ext(r_4) <- ext(chm9_fs4)
terra::res(r_4) <- terra::res(chm9_fs4)  
terra::crs(r_4) <- terra::crs(chm9_fs4)

ftype_1 <- terra::rasterize(foresttype, r_1, field="type")
ftype_2 <- terra::rasterize(foresttype, r_2, field="type")
ftype_3 <- terra::rasterize(foresttype, r_3, field="type")
ftype_4 <- terra::rasterize(foresttype, r_4, field="type")

# rasterize Management Information
mtype_1 <- terra::rasterize(management, r_1, field="zone_id")
mtype_2 <- terra::rasterize(management, r_2, field="zone_id")
mtype_3 <- terra::rasterize(management, r_3, field="zone_id")
mtype_4 <- terra::rasterize(management, r_4, field="zone_id")

# function to calulate gap statistics

# Gap_Stats <- function (gap_layer, chm_layer, gap_polygon) 
# { gap_list <- data.frame(terra::freq(gap_layer))
# gap_list$count <- gap_list$count * raster::res(chm_layer)[1]^2
# gap_list <- gap_list[!is.na(gap_list[, 1]), ]
# # extract raster values per gap
# gap_chm <- terra::extract(chm_layer, gap_polygon) #versucen hier einfach gap_layer und chm_layer zu stacken und in df überführen
# names(gap_chm)[2] <- "chm_values"
# gap_chm<- gap_chm %>%
#   mutate(chm_values = ifelse(chm_values > 5, NA, chm_values)) #filter out values above 5
# # calculate gap statistics
# gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
#                                     summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]
# 
# gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
#                                     summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]
# 
# gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
#                                      summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]
# 
# gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
#                                    summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]
# 
# gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
# gap_list$perimeter <- perim(gap_polygon) #hierfür bräuchte man ggf. dann nicht das exakte Polygon, sondern könnte das von as-polygons nehmen
# 
# gap_list$year <- sub("_.*", "", names(chm_layer))
# gap_list$site <- sub(".*_", "", names(chm_layer))
# 
# gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]
# # add perimeter to calculate perimeter/area ratio
# colnames(gap_list) <- c("gap_id", "gap_area", 
#                         "chm_max", "chm_min", "chm_mean", "chm_sd"
#                         ,"chm_range", "perimeter", "year", "site")
# return(gap_list)
# }

# function to calulate gap statistics without polygon

gap_layer <- gap_stack_fs1$gaps_17_fs1
chm_layer <- chm17_fs1

Gap_Stats <- function (gap_layer, chm_layer) 
{ gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * raster::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer))
names(gap_chm) <- c("ID", "chm_values")
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

gap_list <- cbind(gap_list, lsm_p_perim(gap_layer)[,6] ) #add perimeter

gap_list$year <- sub("_.*", "", names(chm_layer))
gap_list$site <- sub(".*_", "", names(chm_layer))

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "year", "site")
return(gap_list)
}


# --------------------------------- function to calculate gap statistics without polygon + NP information

gap_layer <- gap_stack_fs1$gaps_17_fs1
chm_layer <- chm17_fs1
foresttype <- ftype_1

#do the same with management types!!!!!!

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


Gap_Stats <- function (gap_layer, chm_layer, foresttype) 
{ gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, foresttype), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "ftype")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
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

#gap_list <- cbind(gap_list, lsm_p_perim(raster::raster(gap_layer))[,6] ) #add perimeter   !!!reactivate function on server

gap_list$ftype <- as.data.frame(getForestType(gap_chm))[,2]

gap_list$year <- sub("_.*", "", names(chm_layer))
gap_list$site <- sub(".*_", "", names(chm_layer))

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "forest_type", "year", "site") #add , "perimeter" after range!
return(gap_list)
}



# ---  calculate gap statistics ---

#stats_list <- mapply(function(x,y,z){
#  Gap_Stats_apply(gap_layer=x, chm_layer=y, gap_polygon=z)}, x= gap_list , y= chm_list, z= polygons_erosion) 

# min 100
stats9_1_100 <- Gap_Stats(gap_list_100$chm9_fs1, chm_list$chm9_fs1, polygons_100$chm9_fs1)
stats_17_1_100 <- Gap_Stats(gap_list_100$chm17_fs1, chm_list$chm17_fs1, polygons_100$chm17_fs1)
stats_21_1_100 <- Gap_Stats(gap_list_100$chm21_fs1, chm_list$chm21_fs1, polygons_100$chm21_fs1)

stats9_2_100 <- Gap_Stats(gap_list_100$chm9_fs2, chm_list$chm9_fs2, polygons_100$chm9_fs2)
stats_17_2_100 <- Gap_Stats(gap_list_100$chm17_fs2, chm_list$chm17_fs2, polygons_100$chm17_fs2)
stats_21_2_100 <- Gap_Stats(gap_list_100$chm21_fs2, chm_list$chm21_fs2, polygons_100$chm21_fs2)

stats9_3_100 <- Gap_Stats(gap_list_100$chm9_fs3, chm_list$chm9_fs3, polygons_100$chm9_fs3)
stats_17_3_100 <- Gap_Stats(gap_list_100$chm17_fs3, chm_list$chm17_fs3, polygons_100$chm17_fs3)
stats_21_3_100 <- Gap_Stats(gap_list_100$chm21_fs3, chm_list$chm21_fs3, polygons_100$chm21_fs3)

stats9_4_100 <- Gap_Stats(gap_list_100$chm9_fs4, chm_list$chm9_fs4, polygons_100$chm9_fs4)
stats_17_4_100 <- Gap_Stats(gap_list_100$chm17_fs4, chm_list$chm17_fs4, polygons_100$chm17_fs4)
stats_21_4_100 <- Gap_Stats(gap_list_100$chm21_fs4, chm_list$chm21_fs4, polygons_100$chm21_fs4)

stats_all_100 <- rbind(stats9_1_100, stats9_2_100, stats9_3_100, stats9_4_100, 
                       stats_17_1_100, stats_17_2_100, stats_17_3_100, stats_17_4_100, 
                       stats_21_1_100, stats_21_2_100, stats_21_3_100, stats_21_4_100)
stats_all_100$pa_ratio <- stats_all_100$perimeter/stats_all_100$gap_area

# min 250
stats9_1_250 <- Gap_Stats(gap_list_250$chm9_fs1, chm_list$chm9_fs1, polygons_250$chm9_fs1)
stats_17_1_250 <- Gap_Stats(gap_list_250$chm17_fs1, chm_list$chm17_fs1, polygons_250$chm17_fs1)
stats_21_1_250 <- Gap_Stats(gap_list_250$chm21_fs1, chm_list$chm21_fs1, polygons_250$chm21_fs1)

stats9_2_250 <- Gap_Stats(gap_list_250$chm9_fs2, chm_list$chm9_fs2, polygons_250$chm9_fs2)
stats_17_2_250 <- Gap_Stats(gap_list_250$chm17_fs2, chm_list$chm17_fs2, polygons_250$chm17_fs2)
stats_21_2_250 <- Gap_Stats(gap_list_250$chm21_fs2, chm_list$chm21_fs2, polygons_250$chm21_fs2)

stats9_3_250 <- Gap_Stats(gap_list_250$chm9_fs3, chm_list$chm9_fs3, polygons_250$chm9_fs3)
stats_17_3_250 <- Gap_Stats(gap_list_250$chm17_fs3, chm_list$chm17_fs3, polygons_250$chm17_fs3)
stats_21_3_250 <- Gap_Stats(gap_list_250$chm21_fs3, chm_list$chm21_fs3, polygons_250$chm21_fs3)

stats9_4_250 <- Gap_Stats(gap_list_250$chm9_fs4, chm_list$chm9_fs4, polygons_250$chm9_fs4)
stats_17_4_250 <- Gap_Stats(gap_list_250$chm17_fs4, chm_list$chm17_fs4, polygons_250$chm17_fs4)
stats_21_4_250 <- Gap_Stats(gap_list_250$chm21_fs4, chm_list$chm21_fs4, polygons_250$chm21_fs4)

stats_all_250 <- rbind(stats9_1_250, stats9_2_250, stats9_3_250, stats9_4_250, 
                       stats_17_1_250, stats_17_2_250, stats_17_3_250, stats_17_4_250, 
                       stats_21_1_250, stats_21_2_250, stats_21_3_250, stats_21_4_250)
stats_all_250$pa_ratio <- stats_all_250$perimeter/stats_all_250$gap_area

# min 400
stats9_1_400 <- Gap_Stats(gap_list_400$chm9_fs1, chm_list$chm9_fs1, polygons_400$chm9_fs1)
stats_17_1_400 <- Gap_Stats(gap_list_400$chm17_fs1, chm_list$chm17_fs1, polygons_400$chm17_fs1)
stats_21_1_400 <- Gap_Stats(gap_list_400$chm21_fs1, chm_list$chm21_fs1, polygons_400$chm21_fs1)

stats9_2_400 <- Gap_Stats(gap_list_400$chm9_fs2, chm_list$chm9_fs2, polygons_400$chm9_fs2)
stats_17_2_400 <- Gap_Stats(gap_list_400$chm17_fs2, chm_list$chm17_fs2, polygons_400$chm17_fs2)
stats_21_2_400 <- Gap_Stats(gap_list_400$chm21_fs2, chm_list$chm21_fs2, polygons_400$chm21_fs2)

stats9_3_400 <- Gap_Stats(gap_list_400$chm9_fs3, chm_list$chm9_fs3, polygons_400$chm9_fs3)
stats_17_3_400 <- Gap_Stats(gap_list_400$chm17_fs3, chm_list$chm17_fs3, polygons_400$chm17_fs3)
stats_21_3_400 <- Gap_Stats(gap_list_400$chm21_fs3, chm_list$chm21_fs3, polygons_400$chm21_fs3)

stats9_4_400 <- Gap_Stats(gap_list_400$chm9_fs4, chm_list$chm9_fs4, polygons_400$chm9_fs4)
stats_17_4_400 <- Gap_Stats(gap_list_400$chm17_fs4, chm_list$chm17_fs4, polygons_400$chm17_fs4)
stats_21_4_400 <- Gap_Stats(gap_list_400$chm21_fs4, chm_list$chm21_fs4, polygons_400$chm21_fs4)

stats_all_400 <- rbind(stats9_1_400, stats9_2_400, stats9_3_400, stats9_4_400, 
                       stats_17_1_400, stats_17_2_400, stats_17_3_400, stats_17_4_400, 
                       stats_21_1_400, stats_21_2_400, stats_21_3_400, stats_21_4_400)
stats_all_400$pa_ratio <- stats_all_400$perimeter/stats_all_400$gap_area

# min 400 without polygon
stats9_1_400 <- Gap_Stats(gap_list_400$chm9_fs1, chm_list$chm9_fs1)
stats_17_1_400 <- Gap_Stats(gap_list_400$chm17_fs1, chm_list$chm17_fs1)
stats_21_1_400 <- Gap_Stats(gap_list_400$chm21_fs1, chm_list$chm21_fs1)

stats9_2_400 <- Gap_Stats(gap_list_400$chm9_fs2, chm_list$chm9_fs2)
stats_17_2_400 <- Gap_Stats(gap_list_400$chm17_fs2, chm_list$chm17_fs2)
stats_21_2_400 <- Gap_Stats(gap_list_400$chm21_fs2, chm_list$chm21_fs2)

stats9_3_400 <- Gap_Stats(gap_list_400$chm9_fs3, chm_list$chm9_fs3)
stats_17_3_400 <- Gap_Stats(gap_list_400$chm17_fs3, chm_list$chm17_fs3)
stats_21_3_400 <- Gap_Stats(gap_list_400$chm21_fs3, chm_list$chm21_fs3)

stats9_4_400 <- Gap_Stats(gap_list_400$chm9_fs4, chm_list$chm9_fs4)
stats_17_4_400 <- Gap_Stats(gap_list_400$chm17_fs4, chm_list$chm17_fs4)
stats_21_4_400 <- Gap_Stats(gap_list_400$chm21_fs4, chm_list$chm21_fs4)

stats_all_400 <- rbind(stats9_1_400, stats9_2_400, stats9_3_400, stats9_4_400, 
                       stats_17_1_400, stats_17_2_400, stats_17_3_400, stats_17_4_400, 
                       stats_21_1_400, stats_21_2_400, stats_21_3_400, stats_21_4_400)
stats_all_400$pa_ratio <- stats_all_400$perimeter/stats_all_400$gap_area

# min 400 without polygon with forest type info
stats9_1_400 <- Gap_Stats(gap_list_400$chm9_fs1, chm_list$chm9_fs1, ftype_1)
stats_17_1_400 <- Gap_Stats(gap_list_400$chm17_fs1, chm_list$chm17_fs1, ftype_1)
stats_21_1_400 <- Gap_Stats(gap_list_400$chm21_fs1, chm_list$chm21_fs1, ftype_1)

stats9_2_400 <- Gap_Stats(gap_list_400$chm9_fs2, chm_list$chm9_fs2, ftype_2)
stats_17_2_400 <- Gap_Stats(gap_list_400$chm17_fs2, chm_list$chm17_fs2, ftype_2)
stats_21_2_400 <- Gap_Stats(gap_list_400$chm21_fs2, chm_list$chm21_fs2, ftype_2)

stats9_3_400 <- Gap_Stats(gap_list_400$chm9_fs3, chm_list$chm9_fs3, ftype_3)
stats_17_3_400 <- Gap_Stats(gap_list_400$chm17_fs3, chm_list$chm17_fs3, ftype_3)
stats_21_3_400 <- Gap_Stats(gap_list_400$chm21_fs3, chm_list$chm21_fs3, ftype_3)

stats9_4_400 <- Gap_Stats(gap_list_400$chm9_fs4, chm_list$chm9_fs4, ftype_4)
stats_17_4_400 <- Gap_Stats(gap_list_400$chm17_fs4, chm_list$chm17_fs4, ftype_4)
stats_21_4_400 <- Gap_Stats(gap_list_400$chm21_fs4, chm_list$chm21_fs4, ftype_4)

stats_all_400 <- rbind(stats9_1_400, stats9_2_400, stats9_3_400, stats9_4_400, 
                       stats_17_1_400, stats_17_2_400, stats_17_3_400, stats_17_4_400, 
                       stats_21_1_400, stats_21_2_400, stats_21_3_400, stats_21_4_400)
stats_all_400$pa_ratio <- stats_all_400$perimeter/stats_all_400$gap_area


# change year label
stats_all_100$year <- factor(stats_all_100$year, levels=c("chm9", "chm17", "chm21"), labels=c("2009", "2017", "2021"))
stats_all_250$year <- factor(stats_all_250$year, levels=c("chm9", "chm17", "chm21"), labels=c("2009", "2017", "2021"))
stats_all_400$year <- factor(stats_all_400$year, levels=c("chm9", "chm17", "chm21"), labels=c("2009", "2017", "2021"))

# --- meta stats ---
stats_100_summary <- stats_all_100 %>% count(site, year)
stats_250_summary <- stats_all_250 %>% count(site, year)
stats_400_summary <- stats_all_400 %>% count(site, year)

n_100 <- ggplot(stats_100_summary, aes(x=year, y=n, color=site, group=site)) + geom_point() + geom_line() + theme_minimal() +
  labs(title = "number of gaps per site and year with min 100")
n_250 <- ggplot(stats_250_summary, aes(x=year, y=n, color=site, group=site)) + geom_point() + geom_line() + theme_minimal() +
  labs(title = "number of gaps per site and year with min 250")
n_400 <- ggplot(stats_400_summary, aes(x=year, y=n, color=site, group=site)) + geom_point() + geom_line() + theme_minimal() +
  labs(title = "number of gaps per site and year with min 400")

# --- plotting stats ---

#min 100
area_100 <- ggplot(stats_all_100, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 100")

max_100 <- ggplot(stats_all_100, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 100")

range_100 <- ggplot(stats_all_100, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 100")

mean_100 <- ggplot(stats_all_100, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 100")

pa_100 <- ggplot(stats_all_100, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 100")

#min 250
area_250 <- ggplot(stats_all_250, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 250")

max_250 <- ggplot(stats_all_250, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 250")

range_250 <- ggplot(stats_all_250, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 250")

mean_250 <- ggplot(stats_all_250, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 250")

pa_250 <- ggplot(stats_all_250, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 250")

#min 400
area_400 <- ggplot(stats_all_400, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 400")

max_400 <- ggplot(stats_all_400, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 400")

range_400 <- ggplot(stats_all_400, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 400")

mean_400 <- ggplot(stats_all_400, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 400")

pa_400 <- ggplot(stats_all_400, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with min 400")

#plotting gap size in histogram and density

stats_all_400$gap_area_ha <- stats_all_400$gap_area/10000

stats_all_400 %>% ggplot(aes(gap_area, fill=site)) +
  geom_histogram(alpha=0.25)+
  labs(x="Area of gap size in m2")


stats_all_400 %>% ggplot(aes(gap_area_ha, fill = year)) +  
  #  geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
  geom_density(col="#FFBB00",size=0.5, alpha =0.5) +
  labs(x="Area of gap size in ha", y="Density") + facet_wrap(~site) + theme_minimal()

stats_all_400 %>% filter(gap_area_ha < 1) %>% ggplot(aes(gap_area_ha, fill = year)) +  
  #  geom_histogram(aes(y=..density..),alpha=0.5, bins =100) + #FFBB00
  geom_density(col="#FFBB00",size=1, alpha=0.5) +
  labs(x="Area of gap size in ha", y="Density") + facet_wrap(~site) + theme_minimal()

stats_all_400 %>% filter(gap_area_ha < 1) %>% ggplot(aes(gap_area_ha, fill = year)) +  
  geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
  #  geom_density(col="#FFBB00",size=1, alpha=0.5) +
  labs(x="Area of gap size in ha", y="Density") + facet_wrap(~site) + theme_minimal()

stats_all_400 %>% filter(gap_area_ha < 2) %>% 
  ggplot(aes(x=site, y=gap_area_ha, color=year)) + geom_boxplot() + theme_minimal() 

# plot gap area according to forest type
setwd("~/Documents/global_change_geography/masterthesis/data")
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=16),
  legend.text = element_text(size=14))

tiff("gap_area_ftype.tiff", units="in", width=12, height=8, res=500)
ggplot(stats_all_400, aes(x=forest_type, y=gap_area_ha, color=year)) + geom_boxplot() + theme_minimal() + labs(title="Gap area per forest type", x="forest type", y="gap area in ha") + My_Theme  + scale_color_brewer(palette = "Dark2")
dev.off()


stats_all_400 %>% ggplot(aes(gap_area_ha, fill = year)) +  
  #  geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
  geom_density(col="#FFBB00",size=0.5, alpha =0.5) +
  labs(x="Area of gap size in ha", y="Density") + facet_wrap(~forest_type) + theme_minimal()

stats_all_400 %>% filter(gap_area_ha < 2) %>% ggplot(aes(gap_area_ha, fill = year)) +  
  #  geom_histogram(aes(y=..density..),alpha=0.5, bins =100) + #FFBB00
  geom_density(col="#FFBB00",size=1, alpha=0.5) +
  labs(x="Area of gap size in ha", y="Density") + facet_wrap(~forest_type) + theme_minimal()

stats_all_400 %>% filter(gap_area_ha < 2) %>% ggplot(aes(gap_area_ha, fill = year)) +  
  geom_histogram(aes(y=..density..),alpha=0.5, bins =50) + #FFBB00
  #  geom_density(col="#FFBB00",size=1, alpha=0.5) +
  labs(x="Area of gap size in ha", y="Count") + facet_wrap(~forest_type) + theme_minimal()

##############################################################
# calculate Gap-size-frequency-distribution
##############################################################


Gap_Stats2 <- function (gap_layer, chm_layer) 
{
  GiniCoeff <- function(x, finite.sample = TRUE, na.rm = TRUE) {
    if (!na.rm && any(is.na(x))) {
      return(NA_real_)
    }
    x <- as.numeric(stats::na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- 2 * sum(x * 1L:n)/sum(x) - (n + 1L)
    if (finite.sample) {
      GC <- G/(n - 1L)
    }
    else {
      GC <- G/n
    }
    return(GC)
  }
  gap_list <- data.frame(raster::freq(gap_layer))
  gap_list$count <- gap_list$count * raster::res(chm_layer)[1]^2
  gap_list <- gap_list[!is.na(gap_list[, 1]), ]
  gap_list$chm_max <- tapply(chm_layer[], gap_layer[], max)
  gap_list$chm_min <- tapply(chm_layer[], gap_layer[], min)
  gap_list$chm_mean <- round(tapply(chm_layer[], gap_layer[], 
                                    mean), 2)
  gap_list$chm_sd <- round(tapply(chm_layer[], gap_layer[], 
                                  stats::sd), 2)
  gap_list$chm_gini <- round(tapply(chm_layer[], gap_layer[], 
                                    GiniCoeff), 2)
  gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 
                              2)
  gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]
  colnames(gap_list) <- c("gap_id", "gap_area", 
                          "chm_max", "chm_min", "chm_mean", "chm_sd", 
                          "chm_gini", "chm_range")
  return(gap_list)
}


gap_stack_fs1$gaps_9_fs1[is.na(gap_stack_fs1$gaps_9_fs1)] <- 0 #change NAs to 0
gap_stack_fs1$gaps_17_fs1[is.na(gap_stack_fs1$gaps_17_fs1)] <- 0
gap_stack_fs1$gaps_21_fs1[is.na(gap_stack_fs1$gaps_21_fs1)] <- 0
#stats_fs1_9 <- GapStats(gap_stack_fs1$gaps_9_fs1, chm9_fs1)
stats_fs1_9 <- Gap_Stats2(gap_stack_fs1$gaps_9_fs1, chm9_fs1)
stats_fs1_17<- Gap_Stats2(gap_stack_fs1$gaps_17_fs1, chm17_fs1)
stats_fs1_21<- Gap_Stats2(gap_stack_fs1$gaps_21_fs1, chm21_fs1)
Stats_fs1_9<- stats_fs1_9[!(stats_fs1_9$gap_id == 0),] # exclude non gap area
Stats_fs1_17<- stats_fs1_17[!(stats_fs1_17$gap_id == 0),] # exclude non gap area
Stats_fs1_21<- stats_fs1_21[!(stats_fs1_21$gap_id == 0),] # exclude non gap area


# Gap-size Frequency Distributions
GapSizeFDist(
  gaps_stats = gap_list, method = "Hanel_2017", col = "forestgreen", pch = 16, cex = 1,
  axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
)
axis(1)
axis(2)
grid(4, 4)

# Gap-size Frequency Distributions
GapSizeFDist(
  gaps_stats = Stats_fs1_17, method = "Hanel_2017", col = "forestgreen", pch = 16, cex = 1,
  axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
)
axis(1)
axis(2)
grid(4, 4)

# Gap-size Frequency Distributions
GapSizeFDist(
  gaps_stats = Stats_fs1_21, method = "Hanel_2017", col = "forestgreen", pch = 16, cex = 1,
  axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
)
axis(1)
axis(2)
grid(4, 4)




######################################################################################
# Subset gaps per elevation level, management, other regional criteria
######################################################################################


