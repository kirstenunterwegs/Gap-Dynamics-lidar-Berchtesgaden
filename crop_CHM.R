##########################################################
# crop Berchtesgaden CHMs to closed Forest and focus sites
#########################################################

library(terra)

chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_1m.tif")
chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_1m.tif")
chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_1m.tif")
waldmaske <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/Waldmaske.shp")
#np_border <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/berchtesgaden.gpkg")
focus_sites_large <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/focus_sites_buffer.shp")
focus_sites <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/walddynamik_intensiv_fp.shp")

# --- crop to closed forest areas --- 

closed_forest <- subset(waldmaske, waldmaske$ORIG_FID == 1)
artifacts <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/artifacts_epsg25832.shp")

crs(artifacts) <- crs(closed_forest)
f_artifacts_mask <- terra::erase(closed_forest, artifacts)

#f_artifacts_mask <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_smallext_EPSG25832.shp")

#--- create one mask per extent
f_artifacts_mask_2009<- as(f_artifacts_mask, "Spatial") #convert to SpatialPolygonsDataFrame to access bbox
f_artifacts_mask_2017<- as(f_artifacts_mask, "Spatial")
f_artifacts_mask_2021<- as(f_artifacts_mask, "Spatial")

f_artifacts_mask_2009@bbox <- as.matrix(extent(raster::raster(chm9))) #set bbox to extent of respective CHM
f_artifacts_mask_2017@bbox <- as.matrix(extent(raster::raster(chm17)))
f_artifacts_mask_2021@bbox <- as.matrix(extent(raster::raster(chm21)))

#export
wd <- "i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/"
setwd(wd)
terra::writeVector(vect(f_artifacts_mask_2009), "f_artifacts_mask_2009.shp",overwrite=TRUE)
terra::writeVector(vect(f_artifacts_mask_2017), "f_artifacts_mask_2017.shp",overwrite=TRUE)
terra::writeVector(vect(f_artifacts_mask_2021), "f_artifacts_mask_2021.shp",overwrite=TRUE)

#rasterize mask
mask9 <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2009.shp")
mask17<- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2017.shp")
mask21<- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/f_artifacts_mask_2021.shp")

mask9$ID<- 1 #change ID to 1 for rasterization (1=valid, 0=invalid)
mask17$ID<- 1
mask21$ID<- 1

#crate empty raster for rasterization
r_9 <- rast()
ext(r_9) <- ext(chm9)
terra::res(r_9) <- terra::res(chm9)  
terra::crs(r_9) <- terra::crs(chm9)

r_17 <- rast()
ext(r_17) <- ext(chm17)
terra::res(r_17) <- terra::res(chm17)  
terra::crs(r_17) <- terra::crs(chm17)

r_21 <- rast()
ext(r_21) <- ext(chm21)
terra::res(r_21) <- terra::res(chm21)  
terra::crs(r_21) <- terra::crs(chm21)


# rasterize 
f_artifacts_mask_2009 <- terra::rasterize(mask9, r_9, field="ID",background=0)
f_artifacts_mask_2017 <- terra::rasterize(mask17, r_17, field="ID",background=0)
f_artifacts_mask_2021 <- terra::rasterize(mask21, r_21, field="ID",background=0)

#export
wd <- "i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/"
setwd(wd)
terra::writeRaster(f_artifacts_mask_2009, "f_artifacts_mask_2009.tif",overwrite=TRUE)
terra::writeRaster(f_artifacts_mask_2017, "f_artifacts_mask_2017.tif",overwrite=TRUE)
terra::writeRaster(f_artifacts_mask_2021, "f_artifacts_mask_2021.tif",overwrite=TRUE)

### ---- create CHM without artifacts for canopy height filter ----
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

#chm_fs1_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack.tif")
#chm9 <- chm_fs1_crop[[1]]
#chm17<- chm_fs1_crop[[2]]
#chm21<- chm_fs1_crop[[3]]

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

chm_stack <- rast("chm_berchtesgaden_stack_1m.tif")
chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]

artifacts <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/artifacts_epsg25832.shp")

chm9m <- mask(chm9, artifacts, inverse=TRUE)
chm17m <-  mask(chm17, artifacts, inverse=TRUE)
chm21m <-  mask(chm21, artifacts, inverse=TRUE)

chm_berchtesgaden <- c(chm9m, chm17m, chm21m)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_1m_maskedartifacts.tif",overwrite=TRUE)

#----- create inner buffer 

chm <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m_maskedartifacts.tif")
chm9 <- chm[[1]]

chm_poly <- as.polygons(chm9)
chm_poly_simply <- simplifyGeom(chm_poly, tolerance = 0.8)
terra::writeVector(chm_poly_simply,"i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_poly_geomsimply.shp" )
chm_inner_buffer <- buffer(chm_poly_simply, width = -50)
terra::writeVector(chm_inner_buffer,"i:/Fonda/workspace/berchtesgaden/gaps/clsoed_forest_inner_buffer.shp" )
#---- crop CHMs to Berchtesgaden extent ----- 

#np_border <- project(np_border, chm9)
#np_border_mngt <- subset(np_border, np_border$Id == 6) 
#np_border_noin <- subset(np_border, np_border$Id == 5) 
#np_border_all <- subset(np_border, np_border$Id == 5 | np_border$Id == 6) 

#chm21 <- crop(chm21, np_border_all, snap="near",mask=TRUE)
#chm9 <- crop(chm9, chm21)
#chm17 <- crop(chm17, chm21)

# --- crop CHMS to focus sites ------

#extract single focus sites
fs1 <- subset(focus_sites, focus_sites$Id == 1)
fs2 <- subset(focus_sites, focus_sites$Id == 2)
fs3 <- subset(focus_sites, focus_sites$Id == 3)
fs4 <- subset(focus_sites, focus_sites$Id == 4)

fs1_l <- subset(focus_sites_large, focus_sites$Id == 1)
fs2_l <- subset(focus_sites_large, focus_sites$Id == 2)
fs3_l <- subset(focus_sites_large, focus_sites$Id == 3)
fs4_l <- subset(focus_sites_large, focus_sites$Id == 4)

#crop to focus sites
chm9_fs1 <- crop(chm9, fs1, snap="near",mask=TRUE) # 2009
chm9_fs2 <- crop(chm9, fs2, snap="near",mask=TRUE)
chm9_fs3 <- crop(chm9, fs3, snap="near",mask=TRUE)
chm9_fs4 <- crop(chm9, fs4, snap="near",mask=TRUE)
chm17_fs1 <- crop(chm17, fs1, snap="near",mask=TRUE) # 2017
chm17_fs2 <- crop(chm17, fs2, snap="near",mask=TRUE)
chm17_fs3 <- crop(chm17, fs3, snap="near",mask=TRUE)
chm17_fs4 <- crop(chm17, fs4, snap="near",mask=TRUE)
chm21_fs1 <- crop(chm21, fs1, snap="near",mask=TRUE) # 2021
chm21_fs2 <- crop(chm21, fs2, snap="near",mask=TRUE)
chm21_fs3 <- crop(chm21, fs3, snap="near",mask=TRUE)
chm21_fs4 <- crop(chm21, fs4, snap="near",mask=TRUE)

#--- export 
chm_focus_site1 <- c(chm9_fs1, chm17_fs1, chm21_fs1)
names(chm_focus_site1) <- c("chm9_fs1","chm17_fs1", "chm21_fs1")
terra::writeRaster(chm_focus_site1, "chm_focus_site1_stack_1m.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2, chm17_fs2, chm21_fs2)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_stack_1m.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3, chm17_fs3, chm21_fs3)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_stack_1m.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4, chm17_fs4, chm21_fs4)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_stack_1m.tif",overwrite=TRUE)

#crop to focus sites large (with 500 m buffer)
chm9_fs1_l <- crop(chm9, fs1_l, snap="near",mask=TRUE) # 2009
chm9_fs2_l <- crop(chm9, fs2_l, snap="near",mask=TRUE)
chm9_fs3_l <- crop(chm9, fs3_l, snap="near",mask=TRUE)
chm9_fs4_l <- crop(chm9, fs4_l, snap="near",mask=TRUE)
chm17_fs1_l <- crop(chm17, fs1_l, snap="near",mask=TRUE) # 2017
chm17_fs2_l <- crop(chm17, fs2_l, snap="near",mask=TRUE)
chm17_fs3_l <- crop(chm17, fs3_l, snap="near",mask=TRUE)
chm17_fs4_l <- crop(chm17, fs4_l, snap="near",mask=TRUE)
chm21_fs1_l <- crop(chm21, fs1_l, snap="near",mask=TRUE) # 2021
chm21_fs2_l <- crop(chm21, fs2_l, snap="near",mask=TRUE)
chm21_fs3_l <- crop(chm21, fs3_l, snap="near",mask=TRUE)
chm21_fs4_l <- crop(chm21, fs4_l, snap="near",mask=TRUE)

#--- export 
chm_focus_site1 <- c(chm9_fs1_l, chm17_fs1_l, chm21_fs1_l)
names(chm_focus_site1) <- c("chm9_fs1","chm17_fs1", "chm21_fs1")
terra::writeRaster(chm_focus_site1, "chm_focus_site1_large_stack_1m.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2_l, chm17_fs2_l, chm21_fs2_l)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_large_stack_1m.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3_l, chm17_fs3_l, chm21_fs3_l)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_large_stack_1m.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4_l, chm17_fs4_l, chm21_fs4_l)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_large_stack_1m.tif",overwrite=TRUE)

#--- crop to large focus site with Buffer ----
# using CHM not clipped to Buffer, in order to retrieve the canopy height surrounding the gaps

chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_1m.tif")
chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_1m.tif")
chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_1m.tif")
 
fs1_buffer<- buffer(fs1_l, width = 30) #30 m Buffer, so 20m Buffer around gaps is covered with CHM
fs2_buffer<- buffer(fs2_l, width = 30) 
fs3_buffer<- buffer(fs3_l, width = 30) 
fs4_buffer<- buffer(fs4_l, width = 30) 

chm9_fs1_lbuffer <- crop(chm9, fs1_buffer, snap="near",mask=TRUE) # 2009
chm9_fs2_lbuffer <- crop(chm9, fs2_buffer, snap="near",mask=TRUE)
chm9_fs3_lbuffer <- crop(chm9, fs3_buffer, snap="near",mask=TRUE)
chm9_fs4_lbuffer <- crop(chm9, fs4_buffer, snap="near",mask=TRUE)
chm17_fs1_lbuffer <- crop(chm17, fs1_buffer, snap="near",mask=TRUE) # 2017
chm17_fs2_lbuffer <- crop(chm17, fs2_buffer, snap="near",mask=TRUE)
chm17_fs3_lbuffer <- crop(chm17, fs3_buffer, snap="near",mask=TRUE)
chm17_fs4_lbuffer <- crop(chm17, fs4_buffer, snap="near",mask=TRUE)
chm21_fs1_lbuffer <- crop(chm21, fs1_buffer, snap="near",mask=TRUE) # 2021
chm21_fs2_lbuffer <- crop(chm21, fs2_buffer, snap="near",mask=TRUE)
chm21_fs3_lbuffer <- crop(chm21, fs3_buffer, snap="near",mask=TRUE)
chm21_fs4_lbuffer <- crop(chm21, fs4_buffer, snap="near",mask=TRUE)

#--- export 
chm_focus_site1 <- c(chm9_fs1_lbuffer, chm17_fs1_lbuffer, chm21_fs1_lbuffer)
names(chm_focus_site1) <- c("chm9_fs1","chm17_fs1", "chm21_fs1")
terra::writeRaster(chm_focus_site1, "chm_focus_site1_lbuffer_stack_1m.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2_lbuffer, chm17_fs2_lbuffer, chm21_fs2_lbuffer)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_lbuffer_stack_1m.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3_lbuffer, chm17_fs3_lbuffer, chm21_fs3_lbuffer)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_lbuffer_stack_1m.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4_lbuffer, chm17_fs4_lbuffer, chm21_fs4_lbuffer)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_lbuffer_stack_1m.tif",overwrite=TRUE)

#---- crop to valley

valley <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/valley_utm32N.shp")
chm_fs1_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m.tif")
chm9 <- chm_fs1_crop[[1]]
chm17<- chm_fs1_crop[[2]]
chm21<- chm_fs1_crop[[3]]

# chm9_valley <- mask(chm9, valley) 
# chm17_valley<- mask(chm17, valley) 
# chm21_valley<- mask(chm21, valley) 

chm9_valley <- crop(chm9, valley, snap="near",mask=TRUE) 
chm17_valley <- crop(chm17, valley, snap="near",mask=TRUE) 
chm21_valley <- crop(chm21, valley, snap="near",mask=TRUE) 

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

chm_berchtesgaden <- c(chm9_valley, chm17_valley, chm21_valley)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_1m_valley_crop.tif",overwrite=TRUE)


#---- crop chm to tile + buffer -- trial running time gap detetction

tile <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/chm/789000_5278000_norm_chm.tif")
tile_ext <- ext(tile)

tile_ext <- as.polygons(tile_ext)
crs(tile_ext) <- crs(tile)
tile_ext_buffer<- buffer(tile_ext, width = 30) 

chm_valley <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m_valley.tif")
chm9 <- chm_valley[[1]]

tile_buffer <- crop(chm9, tile_ext_buffer, snap="near",mask=TRUE) # 2009

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)
terra::writeRaster(tile_buffer, "tile_buffer.tif",overwrite=TRUE)


#----- 16.01.23 ---- crop artifacts masked CHM to subarea for sensitivity analysis

chm9 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm9_artifacts_masked.tif")
chm17<- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm21_artifacts_masked.tif")

#subarea <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/sensitivity_analysis_subarea.gpkg")
subarea <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/variability_analysis_subarea.shp")
closed.forest <- vect("F:/Projects/CanopyDynamicsBDG/data/data/Waldmaske/closed_forest.shp")
closed.forest <- project(closed.forest, subarea)

chm9.sub <- crop(mask(chm9, subarea), subarea)
chm9.sub.f <- crop(chm9.sub, closed.forest, mask =T)

chm17.sub <- crop(mask(chm17, subarea), subarea)
chm17.sub <- crop(chm17.sub, closed.forest, mask =T)

chm21.sub <- crop(mask(chm21, subarea), subarea)
chm21.sub <- crop(chm21.sub, closed.forest, mask =T)

writeRaster(chm9.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm9_sub_sensitivity.tif")
writeRaster(chm17.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm17_sub_sensitivity.tif")
writeRaster(chm21.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm21_sub_sensitivity.tif")