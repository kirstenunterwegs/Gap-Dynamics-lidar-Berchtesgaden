##########################################################
#
# crop Berchtesgaden CHMs to closed Forest & mask artifacts 
#
#########################################################

# --- libraries ---

library(terra)

# --- load CHM, closed forest and artifacts mask ---

chm9 <- rast("data/processed/CHM_data/lidar/2009/berchtesgaden_2009_chm1.tif")
chm17<- rast("data/processed/CHM_data/lidar/2017/berchtesgaden_2017_chm1.tif")
chm21<- rast("data/processed/CHM_data/lidar/2021/berchtesgaden_2021_chm1.tif")
closed_forest <- vect("data/raw/closed_forest.gpkg")
artifacts <- vect("data/processed/artifacts_mask/artifacts_mask.gpkg") # manually masked out artifacts from closed forest layer in QGis


# --- crop to closed forest areas --- 


chm21 <- mask(chm21, closed_forest)

#crop and adjust extent

chm17 <- crop(chm17, chm21, snap="near",mask=TRUE) 
ext(chm17) <- ext(chm21)
chm9 <- crop(chm9, chm21, snap="near",mask=TRUE) 
ext(chm9) <- ext(chm21)

#mask 2017 and 2009 CHM

chm17 <- mask(chm17, chm21) # use masked CHM21, to ensure same extent
chm9 <- mask(chm9, chm21)


chm_berchtesgaden <- c(chm9, chm17, chm21)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "data/processed/CHM_data/chm_berchtesgaden_stack_1m.tif")


# ---- mask out artifacts ----

chm_stack <- rast("data/processed/CHM_data/chm_berchtesgaden_stack_1m.tif")
chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]


chm9m <- mask(chm9, artifacts)
chm17m <-  mask(chm17, artifacts)
chm21m <-  mask(chm21, artifacts)

chm_berchtesgaden <- c(chm9m, chm17m, chm21m)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "data/processed/CHM_data/chm_berchtesgaden_stack_1m_maskedartifacts.tif")

# store individual layer 

chm9 <- chm_berchtesgaden[[1]]
chm17<- chm_berchtesgaden[[2]]
chm21<- chm_berchtesgaden[[3]]

terra::writeRaster(chm9, "data/processed/CHM_data/chm9_artifacts_masked.tif")
terra::writeRaster(chm17, "data/processed/CHM_data/chm17_artifacts_masked.tif")
terra::writeRaster(chm21, "data/processed/CHM_data/chm21_artifacts_masked.tif")

# -- crop to Klausbach valley for tree segmentation and growth height extraction ---

valley <- vect("data/raw/klausbach_valley.gpkg") # manually created layer covering the Klausbach (most western) valley of study area
chm_stack <- rast("data/processed/CHM_data/chm_berchtesgaden_stack_1m_maskedartifacts.tif")

chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]


chm9_valley <- crop(chm9, valley, snap="near",mask=TRUE) 
chm17_valley <- crop(chm17, valley, snap="near",mask=TRUE) 
chm21_valley <- crop(chm21, valley, snap="near",mask=TRUE) 


chm_stack_valley <- c(chm9_valley, chm17_valley, chm21_valley)
names(chm_stack_valley) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_stack_valley, "data/processed/CHM_data/chm_berchtesgaden_stack_1m_valley_crop.tif")

