######################################################
# gap detection, filtering and export
#####################################################

library(sf)
library(dplyr)
library(terra)
library(ForestGapR)
library(ForestTools)
library(sp)
library(fieldRS)
library(maptools)


#---- load CHMs ------

# load CHM croped to focus sites + 500m Buffer
chm_valley <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m_valley.tif")
chm9 <- chm_valley[[1]]
chm17<- chm_valley[[2]]
chm21<- chm_valley[[3]]


# ---- Define functions ------

#gap detetction for terra

getForestGaps_terra <- function (chm_layer, threshold = 10, size = c(1, 10^4))  # adapted from ForestRGap getForestGaps function
{
  chm_layer[chm_layer > threshold] <- NA
  chm_layer[chm_layer <= threshold] <- 1
  gaps <- terra::patches(chm_layer, directions = 8, allowGaps=FALSE)
  rcl <- terra::freq(gaps)
  rcl <- rcl[ ,-1]
  rcl[, 2] <- rcl[, 2] * raster::res(chm_layer)[1]^2
  rcl <- cbind(rcl[, 1], rcl)
  z <- raster::reclassify(gaps, rcl = rcl, right = NA)
  z <- classify(gaps, rcl = rcl, right = NA)
  z[is.na(gaps)] <- NA
  gaps[z > size[2]] <- NA
  gaps[z < size[1]] <- NA
  gaps <- terra::patches(gaps, directions = 8, allowGaps=FALSE)
  names(gaps) <- "gaps"
  return(gaps)
}

# filter gaps via bounding box

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

# filter gaps via min height of surrounding dominant trees

filter_gaps_20m_75quant <- function(gap_layer, chm) { #need bigger extent of chm to cover the whole area of the buffer!
  
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

getGapChanges <- function(gap_1, gap2){ #gap 1 = t1 gaps, gap2 = t2 gaps
  change1 <- GapChangeDecTerra(gap_layer1 = gap_1, gap_layer2 = gap2)
  change2 <- GapChangeDecTerra(gap_layer1 = gap2, gap_layer2 = gap_1) #change detetction twice to get gains and losses
  gap_changes <- terra::merge(change1, change2)
  return(gap_changes)   
}


# function to get gap change direction
m = c(-200, -1, 1, -1, 1, 2, 1, 200, 3)# 1= <-1 (loss); 2= -1/1 (steady); 3= >1 (gain)
rclmat = matrix(m, ncol=3, byrow=TRUE)

GetChangeDir <- function(chm1, chm2, gapchange, rclmat) { #chm1 at t1, chm2 at t2
  diff <- chm1 - chm2  # create simple difference
  diff_class = terra::classify(diff, rclmat, include.lowest=TRUE) # Reclassify the raster layer
  gap_change_dir <- mask(diff_class, gapchange)
  return(gap_change_dir)
}


#---- GAP DETECTION -------------------------------

threshold <- 5
size <- c(400,10^100) #m2 min and max gap size  

# prepare lists

chm_l_list <- list(chm9, chm17, chm21)
chm_names <- list("chm9", "chm17", "chm21")
names(chm_l_list) <- chm_names


#--- detect gaps ----
gaps_nofil <- lapply(chm_l_list, function(n){
  getForestGaps_terra(chm_layer = n, threshold = threshold, size = size)
})

# --- loop bounding box filter over gaps

gaps_bbox_fil <- lapply(gaps_nofil, function(n){
  filter_bbox20x20(gap_layer = n)
})

#--- looping canopy height filter over gaps ---

gaps_filtered <- mapply(function(G,C){
  filter_gaps_20m_75quant(gap_layer = G, chm = C)}, G=gaps_bbox_fil, C=chm_l_list) #important to chose chm with buffer for canopy estimation around gaps

# --- identify Gap changes and change direction per site ---

gap_change_917 <- getGapChanges(gaps_filtered$chm9, gaps_filtered$chm17)
gap_change_dir_917 <- GetChangeDir(chm17, chm9, gap_change_917, rclmat) 
gap_change_1721 <- getGapChanges(gaps_filtered$chm21, gaps_filtered$chm17)
gap_change_dir_1721 <- GetChangeDir(chm21, chm17, gap_change_1721, rclmat)



# --- export gap layers ---

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

#stack layer
gap_layers<- c(gaps_filtered$chm9, gaps_filtered$chm17,gaps_filtered$chm21, 
               gap_change_917,gap_change_dir_917,gap_change_1721,gap_change_dir_1721)
names(gap_layers_fs1) <- c("gaps_9", "gaps_17", "gaps_21", "gap_changes_917", "gap_change_dir_917", "gap_changes_1721", "gap_change_dir_1721")
terra::writeRaster(gap_layers, "gap_layers_valley.tif",overwrite=TRUE)



