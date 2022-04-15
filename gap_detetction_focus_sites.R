######################################################
# gap detection, filtering and export
#####################################################

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
library(fieldRS)


#---- load CHMs ------


# load CHM croped to focus sites + 500m Buffer
chm_fs1_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site1_large_stack.tif")
chm9_fs1_l <- chm_fs1_crop[[1]]
chm17_fs1_l<- chm_fs1_crop[[2]]
chm21_fs1_l<- chm_fs1_crop[[3]]
chm_fs2_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site2_large_stack.tif")
chm9_fs2_l <- chm_fs2_crop[[1]]
chm17_fs2_l<- chm_fs2_crop[[2]]
chm21_fs2_l<- chm_fs2_crop[[3]]
chm_fs3_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site3_large_stack.tif")
chm9_fs3_l <- chm_fs3_crop[[1]]
chm17_fs3_l<- chm_fs3_crop[[2]]
chm21_fs3_l<- chm_fs3_crop[[3]]
chm_fs4_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site4_large_stack.tif")
chm9_fs4_l <- chm_fs4_crop[[1]]
chm17_fs4_l<- chm_fs4_crop[[2]]
chm21_fs4_l<- chm_fs4_crop[[3]]

# load CHM croped to focus sites + Buffer for average Height detetction for filtering
chm_fs1_crop_buffer <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site1_lbuffer_stack.tif")
chm9_fs1_lbuffer <- chm_fs1_crop_buffer[[1]]
chm17_fs1_lbuffer<- chm_fs1_crop_buffer[[2]]
chm21_fs1_lbuffer<- chm_fs1_crop_buffer[[3]]
chm_fs2_crop_buffer <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site2_lbuffer_stack.tif")
chm9_fs2_lbuffer <- chm_fs2_crop_buffer[[1]]
chm17_fs2_lbuffer<- chm_fs2_crop_buffer[[2]]
chm21_fs2_lbuffer<- chm_fs2_crop_buffer[[3]]
chm_fs3_crop_buffer <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site3_lbuffer_stack.tif")
chm9_fs3_lbuffer <- chm_fs3_crop_buffer[[1]]
chm17_fs3_lbuffer<- chm_fs3_crop_buffer[[2]]
chm21_fs3_lbuffer<- chm_fs3_crop_buffer[[3]]
chm_fs4_crop_buffer <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_focus_site4_lbuffer_stack.tif")
chm9_fs4_lbuffer <- chm_fs4_crop_buffer[[1]]
chm17_fs4_lbuffer<- chm_fs4_crop_buffer[[2]]
chm21_fs4_lbuffer<- chm_fs4_crop_buffer[[3]]

# ---- Define functions

# --- filter gaps via bounding box

crs <- "epsg:25832" #define crs outside of function
filter_bbox20x20 <- function(gap_layer){
  
  gaps_poly_spat <- as.polygons(gap_layer)
  bboxfilter <- c()
  
  for(i in 1:length(gaps_poly_spat)) {
    single_spat <- gaps_poly_spat[gaps_poly_spat$gaps == i] #extract single polygon
    s_spat_ext <- ext(single_spat)  #get extent of single polygon
    # create border points for bounding box
    ext_points<- vect(cbind(c(s_spat_ext[1],s_spat_ext[1], s_spat_ext[2]), c( s_spat_ext[3],s_spat_ext[4],s_spat_ext[3])),crs=crs)
    dist_ext <- distance(ext_points) #calculate distance between points in m (1 is y dimension, 1 is x dimension)
    #check bounding box distance for 20m threshold and if below add ID to bboxfilter  
    ifelse(dist_ext[1] <20 | dist_ext[2] <20, bboxfilter <- append(bboxfilter, i), next)
  }
  
  bboxfilter <- as.data.frame(bboxfilter)
  if (dim(bboxfilter)[1] == 0) {gaps_filtered <- gap_layer    #assign original gap layer when no gaps filtered out
  } 
  else{   #else filter out gaps
    bboxfilter$replace <- NA
    gaps_filtered <- subst(gap_layer, from=bboxfilter$bboxfilter, to=bboxfilter$replace)
    return(gaps_filtered)
  }
}


#filter function canopy height around gap

#filter_gaps_20m_75quant <- function(gap_layer, chm) { #need bigger extent of chm to cover the whole area of the buffer!
#  if(is.na(raster::minValue(gap_layer)) == TRUE & is.na(raster::maxValue(gap_layer)) ==TRUE) {  # check if there are gaps at all, if no keep NA raster
#    gaps_filtered <- gap_layer
#  }
#  else{
#  gaps__polygon <- GapSPDF(gap_layer)
#  gaps_buffer <- buffer(vect(gaps__polygon), width = 20) 
#  buffer_area <- mask(chm, gaps_buffer)
#  buffer_area <- mask(buffer_area, vect(gaps__polygon), inverse=TRUE)   #extract CHM values for only buffer region without gap

#  mean_canopy_around_gap <-  terra::extract(buffer_area, gaps_buffer)
#  mean_canopy_around_gap <- mean_canopy_around_gap %>%                   #extract 75th quantile of CHM height in buffer area
#                                         group_by(ID) %>% 
#                                         summarize(quant75 = quantile(mean_canopy_around_gap[,2], probs =0.75, na.rm=TRUE))
#  mean_canopy_filter <- mean_canopy_around_gap[!(mean_canopy_around_gap[,2] > 10),]                         #check if buffers mean < 10m
#  if (dim(mean_canopy_filter)[1] == 0) {gaps_filtered <- gap_layer                                         #assign original gap layer when no gaps filtered out
#  } 
#  else{                                                                                                   #else filter out gaps
#    mean_canopy_filter$replace <- NA
#    gaps_filtered <- raster::subs(gap_layer, mean_canopy_filter, by=1, which="replace",subsWithNA=FALSE )
#    return(gaps_filtered)
#  }
# }
#}

filter_gaps_20m_75quant <- function(gap_layer, chm) { #need bigger extent of chm to cover the whole area of the buffer!
  if(is.na(minmax(gap_layer))[1,] == TRUE & is.na(minmax(gap_layer))[2,] ==TRUE) {  # check if there are gaps at all, if no keep NA raster
    gaps_filtered <- gap_layer
  }
  else{
    gaps__polygon <- as.polygons(gap_layer)
    gaps_buffer <- buffer(gaps__polygon, width = 20) 
    buffer_area <- mask(chm, gaps_buffer)
    buffer_area <- mask(buffer_area, gaps__polygon, inverse=TRUE)   #extract CHM values for only buffer region without gap area
    
    chm_buffer <-  terra::extract(buffer_area, gaps_buffer)
    chm_buffer <- chm_buffer %>%                   #extract 75th quantile of CHM height in buffer area
      group_by(ID) %>% 
      summarize(quant75 = quantile(chm_buffer[,2], probs =0.75, na.rm=TRUE))
    canopy_filter <- chm_buffer[!(chm_buffer[,2] > 10),]                     #check if buffers 75th percentile < 10m, if yes write in filter df
    if (dim(canopy_filter)[1] == 0) {gaps_filtered <- gap_layer              #assign original gap layer when no gaps filtered out
    } 
    else{   #else filter out gaps
      canopy_filter$replace <- NA
      gaps_filtered <- subst(gap_layer, from=canopy_filter$value, to=canopy_filter$replace) 
      return(gaps_filtered)
    }
  }
}


# function to retrieve gap changes

GapChangeDecTerra <- function (gap_layer1, gap_layer2)  # adapted from ForestRGap GapChangeDec function
{
  gap_layer1[!is.na(gap_layer1)] <- 1
  gap_layer2[!is.na(gap_layer2)] <- 2
  gap_layer1[is.na(gap_layer1)] <- 0
  gap_layer2[is.na(gap_layer2)] <- 0
  gap_diff <- gap_layer2 - gap_layer1
  gap_diff[gap_diff != 2] <- NA
  return(gap_diff)
}

getGapChanges <- function(gap_1, gap2){ #gap 1 = 2009 gaps, gap2 = 2021 gaps
  change1 <- GapChangeDecTerra(gap_layer1 = gap_1, gap_layer2 = gap2)
  change2 <- GapChangeDecTerra(gap_layer1 = gap2, gap_layer2 = gap_1)
  gap_changes <- raster::merge(change1, change2)
  return(gap_changes)   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!muss gucken das output SpatRaster ist oder nächste Func raster nimmmt
}


#get Gap changes function (detecting gains and losses!)
#getGapChanges <- function(gap_1, gap2){ #gap 1 = 2009 gaps, gap2 = 2021 gaps
#change1 <- GapChangeDec(gap_layer1 = gap_1, gap_layer2 = gap2)
#change2 <- GapChangeDec(gap_layer1 = gap2, gap_layer2 = gap_1)
#gap_changes <- raster::merge(change1, change2)
#return(gap_changes)
#}

# function to get gap change direction
m = c(-200, -1, 1, -1, 1, 2, 1, 200, 3)# 1= <-1 (loss); 2= -1/1 (steady); 3= >1 (gain)
rclmat = matrix(m, ncol=3, byrow=TRUE)

GetChangeDir <- function(chm1, chm2, gapchange, rclmat) { #chm1 ist chm21, chm2 is chm9
  diff <- chm1 - chm2  # create simple difference
  diff_class = terra::classify(diff, rclmat, include.lowest=TRUE) # Reclassify the raster layer
  # gap_change_dir <- raster::mask(raster::raster(diff_class), gapchange)
  gap_change_dir <- mask(diff_class, gapchange)
  return(gap_change_dir)
}

# function to get polygons from gaps
getGapPolygons <- function(gap_layer) {
  if(is.na(raster::minValue(gap_layer)) == TRUE & is.na(raster::maxValue(gap_layer)) ==TRUE) {  # check if there are gaps at all, if no keep NA raster
    print("no gaps")
  }
  else{GapSPDF(gap_layer)}
}

#---- GAP DETECTION -------------------------------

threshold <- 5
size <- c(100,10^10) #m2 min and max gap size  

# prepare lists

chm_l_list <- list(chm9_fs1_l, chm17_fs1_l, chm21_fs1_l, chm9_fs2_l, chm17_fs2_l, chm21_fs2_l,
                   chm9_fs3_l, chm17_fs3_l, chm21_fs3_l,chm9_fs4_l, chm17_fs4_l, chm21_fs4_l)
chm_l_list2 <- list(chm9_fs1_lbuffer, chm17_fs1_lbuffer, chm21_fs1_lbuffer, chm9_fs2_lbuffer, chm17_fs2_lbuffer, chm21_fs2_lbuffer,
                    chm9_fs3_lbuffer, chm17_fs3_lbuffer, chm21_fs3_lbuffer,chm9_fs4_lbuffer, chm17_fs4_lbuffer, chm21_fs4_lbuffer)
chm_names <- list("chm9_fs1", "chm17_fs1", "chm21_fs1", "chm9_fs2", "chm17_fs2", "chm21_fs2",
                  "chm9_fs3", "chm17_fs3", "chm21_fs3","chm9_fs4", "chm17_fs4", "chm21_fs4")
names(chm_l_list) <- chm_names

# small for trial 
#chm_l_list <- list(chm9_fs3_l, chm17_fs3_l, chm21_fs3_l)
#chm_l_list2 <- list(chm9_fs3_lbuffer, chm17_fs3_lbuffer, chm21_fs3_lbuffer)
#chm_names <- list("chm9_fs3", "chm17_fs3", "chm21_fs3")
#names(chm_l_list) <- chm_names


#--- detect gaps ----
gaps_nofil <- lapply(chm_l_list, function(n){
  rast(getForestGaps(chm_layer = raster::raster(n), threshold = threshold, size = size))
})

# --- loop bounding box filter over gaps

gaps_bbox_fil <- lapply(gaps_nofil, function(n){
  filter_bbox20x20(gap_layer = n)
})

#--- looping canopy height and bounding box filter over gaps ---

gaps_filtered <- mapply(function(G,C){
  filter_gaps_20m_75quant(gap_layer = G, chm = C)}, G=gaps_bbox_fil, C=chm_l_list2) #important to chose chm with buffer for canopy estimation around gaps

# --- identify Gap changes and change direction per site ---

gap_change_fs1_917 <- getGapChanges(gaps_filtered$chm9_fs1, gaps_filtered$chm17_fs1)
gap_change_dir1_917 <- GetChangeDir(chm17_fs1_l, chm9_fs1_l, gap_change_fs1_917, rclmat) 
gap_change_fs1_1721 <- getGapChanges(gaps_filtered$chm21_fs1, gaps_filtered$chm17_fs1)
gap_change_dir1_1721 <- GetChangeDir(chm21_fs1_l, chm17_fs1_l, gap_change_fs1_1721, rclmat)

gap_change_fs2_917 <- getGapChanges(gaps_filtered$chm9_fs2, gaps_filtered$chm17_fs2)
gap_change_dir2_917 <- GetChangeDir(chm17_fs2_l, chm9_fs2_l, gap_change_fs2_917, rclmat)
gap_change_fs2_1721 <- getGapChanges(gaps_filtered$chm21_fs2, gaps_filtered$chm17_fs2)
gap_change_dir2_1721 <- GetChangeDir(chm21_fs2_l, chm17_fs2_l, gap_change_fs2_1721, rclmat)

gap_change_fs3_917 <- getGapChanges(gaps_filtered$chm9_fs3, gaps_filtered$chm17_fs3)
gap_change_dir3_917 <- GetChangeDir(chm17_fs3_l, chm9_fs3_l, gap_change_fs3_917, rclmat)
gap_change_fs3_1721 <- getGapChanges(gaps_filtered$chm21_fs3, gaps_filtered$chm17_fs3)
gap_change_dir3_1721 <- GetChangeDir(chm21_fs3_l, chm17_fs3_l, gap_change_fs3_1721, rclmat)

gap_change_fs4_917 <- getGapChanges(gaps_filtered$chm9_fs4, gaps_filtered$chm17_fs4)
gap_change_dir4_917 <- GetChangeDir(chm17_fs4_l, chm9_fs4_l, gap_change_fs4_917, rclmat)
gap_change_fs4_1721 <- getGapChanges(gaps_filtered$chm21_fs4, gaps_filtered$chm17_fs4)
gap_change_dir4_1721 <- GetChangeDir(chm21_fs4_l, chm17_fs4_l, gap_change_fs4_1721, rclmat)


# --- export gap layers ---

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

#focus site 1
#gap_layers_fs1 <- c(rast(gaps_filtered$chm9_fs1), rast(gaps_filtered$chm17_fs1), rast(gaps_filtered$chm21_fs1), 
#                    rast(gap_change_fs1_917 ), rast(gap_change_dir1_917),rast(gap_change_fs1_1721 ), rast(gap_change_dir1_1721))
gap_layers_fs1 <- c(gaps_filtered$chm9_fs1, gaps_filtered$chm17_fs1,gaps_filtered$chm21_fs1, 
                    gap_change_fs1_917,gap_change_dir1_917,gap_change_fs1_1721,gap_change_dir1_1721)
names(gap_layers_fs1) <- c("gaps_9_fs1", "gaps_17_fs1", "gaps_21_fs1", "gap_changes_fs1_917", "gap_change_dir1_917", "gap_changes_fs1_1721", "gap_change_dir1_1721")
terra::writeRaster(gap_layers_fs1, "gap_layers_fs1_100.tif",overwrite=TRUE)

#focus site 2
#gap_layers_fs2 <- c(rast(gaps_filtered$chm9_fs2), rast(gaps_filtered$chm17_fs2), rast(gaps_filtered$chm21_fs2), 
#                    rast(gap_change_fs2_917 ), rast(gap_change_dir2_917),rast(gap_change_fs2_1721 ), rast(gap_change_dir2_1721))
gap_layers_fs2 <- c(gaps_filtered$chm9_fs2, gaps_filtered$chm17_fs2, gaps_filtered$chm21_fs2, 
                    gap_change_fs2_917, gap_change_dir2_917,gap_change_fs2_1721 ,gap_change_dir2_1721)
names(gap_layers_fs2) <- c("gaps_9_fs2", "gaps_17_fs2", "gaps_21_fs2", "gap_changes_fs2_917", "gap_change_dir2_917", "gap_changes_fs2_1721", "gap_change_dir2_1721")
terra::writeRaster(gap_layers_fs2, "gap_layers_fs2_100.tif",overwrite=TRUE)

#focus site 3
#gap_layers_fs3 <- c(rast(gaps_filtered$chm9_fs3), rast(gaps_filtered$chm17_fs3), rast(gaps_filtered$chm21_fs3), 
#                    rast(gap_change_fs3_917 ), rast(gap_change_dir3_917),rast(gap_change_fs3_1721 ), rast(gap_change_dir3_1721))
gap_layers_fs3 <- c(gaps_filtered$chm9_fs3, gaps_filtered$chm17_fs3,gaps_filtered$chm21_fs3, 
                    gap_change_fs3_917 , gap_change_dir3_917,gap_change_fs3_1721,gap_change_dir3_1721)
names(gap_layers_fs3) <- c("gaps_9_fs3", "gaps_17_fs3", "gaps_21_fs3", "gap_changes_fs3_917", "gap_change_dir3_917", "gap_changes_fs3_1721", "gap_change_dir3_1721")
terra::writeRaster(gap_layers_fs3, "gap_layers_fs3_100.tif",overwrite=TRUE)

#focus site 4
#gap_layers_fs4 <- c(rast(gaps_filtered$chm9_fs4), rast(gaps_filtered$chm17_fs4), rast(gaps_filtered$chm21_fs4), 
#                    rast(gap_change_fs4_917 ), rast(gap_change_dir4_917),rast(gap_change_fs4_1721 ), rast(gap_change_dir4_1721))
gap_layers_fs4 <- c(gaps_filtered$chm9_fs4, gaps_filtered$chm17_fs4, gaps_filtered$chm21_fs4, 
                    gap_change_fs4_917 , gap_change_dir4_917,gap_change_fs4_1721, gap_change_dir4_1721)
names(gap_layers_fs4) <- c("gaps_9_fs4", "gaps_17_fs4", "gaps_21_fs4", "gap_changes_fs4_917", "gap_change_dir4_917", "gap_changes_fs4_1721", "gap_change_dir4_1721")
terra::writeRaster(gap_layers_fs4, "gap_layers_fs4_100.tif",overwrite=TRUE)


# --- convert gaps to polygons ----

gaps_filtered_polygons <- lapply(gaps_filtered, function(n){
  as.polygons(n)
})
polygon_names <- as.list(c(names(gaps_filtered_polygons)))  #create list of names for export

###### for raster input---
#gaps_filtered_polygons <- lapply(gaps_filtered, function(n){
#  getGapPolygons(n)
#})
#polygon_names <- as.list(c(names(gaps_filtered_polygons)))
# delete layers without polygons/gaps
#list.condition <- sapply(gaps_filtered_polygons, function(x) class(x) == "SpatialPolygonsDataFrame")
#gap_polygons <- gaps_filtered_polygons[list.condition]

#gap_polygons_vect <- lapply(gap_polygons, function(n) vect(n)) # transform to SpatVector for export
#polygon_names <- as.list(c(names(gap_polygons_vect)))  #create list of names for export


#gaps_filtered_polygons_stack <- list(vect(gaps_filtered_polygons$chm9_fs1), vect(gaps_filtered_polygons$chm17_fs1), vect(gaps_filtered_polygons$chm21_fs1),
#                                  vect(gaps_filtered_polygons$chm9_fs2), vect(gaps_filtered_polygons$chm17_fs2), vect(gaps_filtered_polygons$chm21_fs2),
#                                  vect(gaps_filtered_polygons$chm9_fs3), vect(gaps_filtered_polygons$chm17_fs3), vect(gaps_filtered_polygons$chm21_fs3),
#                                  vect(gaps_filtered_polygons$chm9_fs4), vect(gaps_filtered_polygons$chm17_fs4), vect(gaps_filtered_polygons$chm21_fs4))
#names(gaps_filtered_polygons_stack) <- chm_names

# write polygons to file (can I write the collection to file?)
#mapply(function(x,y){
#  writeVector(x, paste("gaps_polygons", y , "max1ha.shp", sep="_"), overwrite=TRUE)}, x = gaps_filtered_polygons_stack, y= chm_names)

# write polygons to file 
mapply(function(x,y){
  writeVector(x, paste("gaps_polygons_100", y , sep="_"), overwrite=TRUE)}, x = gaps_filtered_polygons, y= polygon_names)


#----------------- Code which has been used in trials -----------------#

#--- loop erosion filter over gap layer ----

# filter function applying erosion filter
GapErosionFilter <- function(gaps) {
  gaps_erode <- pixFilter(gaps, 2, method = "erode")
  gaps_erode_patch <- patches(rast(gaps_erode), allowGaps=FALSE)
  filter_size <- as.data.frame(freq(gaps_erode_patch))%>% select(count, value) %>% group_by(value) %>% summarize(patch_m2 = count/4)
  filter_size <- filter_size[!(filter_size[,2] > 400),]                         #check if patch size > 400 m2                                                                                               #else filter out gaps
  filter_size$replace <- NA
  gaps_erode_filtered <- subst(gaps_erode_patch, from=filter_size$value, to=filter_size$replace) #replace gap ids of gaps < 400m2 with NA
  
  gaps_dilate <- pixFilter(gaps_erode_filtered, 2, method = "dilate") # add back 2 pixels to gap which I deduced with erosion filter
  gaps_dilate_patch <- patches(gaps_dilate, allowGaps=FALSE)
  return(gaps_dilate_patch)
}

gaps_erosion_fil <- lapply(gaps_nofil, function(n){
  GapErosionFilter(n)
})

# add layer name
for(i in chm_names){ 
  gaps_erosion_fil[[i]]@data@names <- paste("gaps", i, sep="_")
}

# add layer name
#for(i in chm_names){ 
#  gaps_nofil[[i]]@data@names <- paste("gaps", i, sep="_")
#}

# BBox filter Erarbeitung
gap_layer <- gaps_nofil$chm9_fs3 
chm <- chm9_fs3_l #need bigger extent of chm to cover the whole area of the buffer!
crs <- "epsg:25832"
bboxfilter <- c()
gaps_poly_spat <- as.polygons(gap_layer)


for(i in 1:length(gaps_poly_spat)) {
  single_spat <- gaps_poly_spat[gaps_poly_spat$gaps == i] #extract single polygon
  s_spat_ext <- ext(single_spat)  #get extent of single polygon
  # create border points for bounding box
  ext_points<- vect(cbind(c(s_spat_ext[1],s_spat_ext[1], s_spat_ext[2]), c( s_spat_ext[3],s_spat_ext[4],s_spat_ext[3])),crs=crs)
  # distance(ext_points) 
  dist_ext <- distance(ext_points) #calculate distance between points in m (1 is y dimension, 1 is x dimension)
  #dist_ext[1] # length y dimension in m
  #dist_ext[2] # length x dimension in m
  # those must be min 20, here 5 for filter trial 
  ifelse(dist_ext[1] <20 | dist_ext[2] <20, bboxfilter <- append(bboxfilter, i), next)
}

bboxfilter <- as.data.frame(bboxfilter)
bboxfilter$replace <- NA
gaps_filtered <- subst(gap_layer, from=bboxfilter$bboxfilter, to=bboxfilter$replace)
