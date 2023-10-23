##########################################################
#
# crop Berchtesgaden CHMs to closed Forest & mask artifacts 
#
#########################################################

# --- libraries ---

library(terra)

# --- load CHM and artifacts mask ---

chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_1m.tif")
chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_1m.tif")
chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_1m.tif")
waldmaske <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/Waldmaske.shp")
artifacts <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/artifacts_epsg25832.shp") # manually delineated 

#np_border <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/berchtesgaden.gpkg")


# --- crop to closed forest areas and mask artifacs --- 

############################################################# check if I really cropped to closed forests and artifacts

# closed_forest <- subset(waldmaske, waldmaske$ORIG_FID == 1) # subset to closed forest areas (=1)
# 
# # mask out artifacts from closed forest mask: get mask for CHM
# 
# crs(artifacts) <- crs(closed_forest)
# f_artifacts_mask <- terra::erase(closed_forest, artifacts)
# 
# 
# #--- create one mask per extent
# 
# f_artifacts_mask_2009<- as(f_artifacts_mask, "Spatial") #convert to SpatialPolygonsDataFrame to access bbox
# f_artifacts_mask_2017<- as(f_artifacts_mask, "Spatial")
# f_artifacts_mask_2021<- as(f_artifacts_mask, "Spatial")
# 
# f_artifacts_mask_2009@bbox <- as.matrix(extent(raster::raster(chm9))) #set bbox to extent of respective CHM
# f_artifacts_mask_2017@bbox <- as.matrix(extent(raster::raster(chm17)))
# f_artifacts_mask_2021@bbox <- as.matrix(extent(raster::raster(chm21)))
# 
# # export artifacts masks
# 
# wd <- "i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/"
# setwd(wd)
# terra::writeVector(vect(f_artifacts_mask_2009), "f_artifacts_mask_2009.shp",overwrite=TRUE)
# terra::writeVector(vect(f_artifacts_mask_2017), "f_artifacts_mask_2017.shp",overwrite=TRUE)
# terra::writeVector(vect(f_artifacts_mask_2021), "f_artifacts_mask_2021.shp",overwrite=TRUE)
# 
# # rasterize mask
# 
# mask9 <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2009.shp")
# mask17<- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2017.shp")
# mask21<- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2021.shp")
# 
# mask9$ID<- 1 #change ID to 1 for rasterization (1=valid, 0=invalid)
# mask17$ID<- 1
# mask21$ID<- 1
# 
# # crate empty raster for rasterization
# r_9 <- rast()
# ext(r_9) <- ext(chm9)
# terra::res(r_9) <- terra::res(chm9)  
# terra::crs(r_9) <- terra::crs(chm9)
# 
# r_17 <- rast()
# ext(r_17) <- ext(chm17)
# terra::res(r_17) <- terra::res(chm17)  
# terra::crs(r_17) <- terra::crs(chm17)
# 
# r_21 <- rast()
# ext(r_21) <- ext(chm21)
# terra::res(r_21) <- terra::res(chm21)  
# terra::crs(r_21) <- terra::crs(chm21)
# 
# 
# # rasterize 
# 
# f_artifacts_mask_2009 <- terra::rasterize(mask9, r_9, field="ID",background=0)
# f_artifacts_mask_2017 <- terra::rasterize(mask17, r_17, field="ID",background=0)
# f_artifacts_mask_2021 <- terra::rasterize(mask21, r_21, field="ID",background=0)
# 
# #export
# 
# wd <- "i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/"
# setwd(wd)
# terra::writeRaster(f_artifacts_mask_2009, "f_artifacts_mask_2009.tif",overwrite=TRUE)
# terra::writeRaster(f_artifacts_mask_2017, "f_artifacts_mask_2017.tif",overwrite=TRUE)
# terra::writeRaster(f_artifacts_mask_2021, "f_artifacts_mask_2021.tif",overwrite=TRUE)


### ---- create CHM without artifacts ----

artifacts <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/artifacts_outclosedforest_epsg25832.shp")
artifacts$replace <- 0

artifacts_2009 <- terra::rasterize(artifacts, r_9, field="replace",background=NA)
artifacts_2017 <- terra::rasterize(artifacts, r_17, field="replace",background=NA)
artifacts_2021 <- terra::rasterize(artifacts, r_21, field="replace",background=NA)

chm9_masked <- terra::mask(chm9, artifacts,inverse=TRUE, updatevalue=0)
chm17_masked <- terra::mask(chm17, artifacts,inverse=TRUE,updatevalue=0) 
chm21_masked <- terra::mask(chm21, artifacts,inverse=TRUE, updatevalue=0)

terra::writeRaster(chm21_masked, "i:/Fonda/workspace/berchtesgaden/lidar/2021/chm21_artifacts_masked.tif",overwrite=TRUE)
terra::writeRaster(chm17_masked, "i:/Fonda/workspace/berchtesgaden/lidar/2017/chm17_artifacts_masked.tif",overwrite=TRUE)
terra::writeRaster(chm9_masked, "i:/Fonda/workspace/berchtesgaden/lidar/2009/chm9_artifacts_masked.tif",overwrite=TRUE)


chm21 <- mask(chm21, closed_forest)
#crop and adjust extent
chm17 <- crop(chm17, chm21, snap="near",mask=TRUE) 
ext(chm17) <- ext(chm21)
chm9 <- crop(chm9, chm21, snap="near",mask=TRUE) 
ext(chm9) <- ext(chm21)
#mask 2017 and 2009 CHM
chm17 <- mask(chm17, chm21) # use masked CHM21, to ensure same extent
chm9 <- mask(chm9, chm21)


wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

chm_berchtesgaden <- c(chm9, chm17, chm21)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_1m.tif",overwrite=TRUE)


# ---- mask out artifacts ----

############################################# check if I have really done that
# 
# chm_stack <- rast("chm_berchtesgaden_stack_1m.tif")
# chm9 <- chm_stack[[1]]
# chm17<- chm_stack[[2]]
# chm21<- chm_stack[[3]]
# 
# artifacts <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/artifacts_epsg25832.shp")
# 
# chm9m <- mask(chm9, artifacts, inverse=TRUE)
# chm17m <-  mask(chm17, artifacts, inverse=TRUE)
# chm21m <-  mask(chm21, artifacts, inverse=TRUE)
# 
# chm_berchtesgaden <- c(chm9m, chm17m, chm21m)
# names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
# terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_1m_maskedartifacts.tif",overwrite=TRUE)




# -- crop to Klausbach valley for tree segmentation and growth height extraction ---

valley <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/valley_utm32N.shp") # manually delineated layer
chm_stack <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m.tif")

chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]


chm9_valley <- crop(chm9, valley, snap="near",mask=TRUE) 
chm17_valley <- crop(chm17, valley, snap="near",mask=TRUE) 
chm21_valley <- crop(chm21, valley, snap="near",mask=TRUE) 

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

chm_berchtesgaden <- c(chm9_valley, chm17_valley, chm21_valley)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_1m_valley_crop.tif",overwrite=TRUE)

