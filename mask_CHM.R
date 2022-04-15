##########################################################
# crop Berchtesgaden CHMs to closed Forest and focus sites
#########################################################

library(terra)

chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_pf-sc02.tif")
chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_pf-sc02.tif")
chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_fix.tif")
waldmaske <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/Waldmaske/Waldmaske.shp")
#np_border <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/berchtesgaden.gpkg")
focus_sites_large <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/focus_sites_buffer.shp")
focus_sites <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/walddynamik_intensiv_fp.shp")

# --- crop to closed forest areas --- 

closed_forest <- subset(waldmaske, waldmaske$ORIG_FID == 1)

chm_fs1_crop <- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack.tif")
chm9 <- chm_fs1_crop[[1]]
chm17<- chm_fs1_crop[[2]]
chm21<- chm_fs1_crop[[3]]

chm21 <- mask(chm21, closed_forest)
chm17 <- mask(chm17, chm21) # use masked CHM21, to ensure same extent
chm9 <- mask(chm9, chm21)


wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

chm_berchtesgaden <- c(chm9, chm17, chm21)
names(chm_berchtesgaden) <- c("chm9", "chm17", "chm21")
terra::writeRaster(chm_berchtesgaden, "chm_berchtesgaden_stack_0.5.tif",overwrite=TRUE)

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
terra::writeRaster(chm_focus_site1, "chm_focus_site1_stack.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2, chm17_fs2, chm21_fs2)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_stack.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3, chm17_fs3, chm21_fs3)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_stack.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4, chm17_fs4, chm21_fs4)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_stack.tif",overwrite=TRUE)

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
terra::writeRaster(chm_focus_site1, "chm_focus_site1_large_stack.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2_l, chm17_fs2_l, chm21_fs2_l)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_large_stack.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3_l, chm17_fs3_l, chm21_fs3_l)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_large_stack.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4_l, chm17_fs4_l, chm21_fs4_l)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_large_stack.tif",overwrite=TRUE)

# crop to large focus site with Buffer

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
terra::writeRaster(chm_focus_site1, "chm_focus_site1_lbuffer_stack.tif",overwrite=TRUE)

chm_focus_site2 <- c(chm9_fs2_lbuffer, chm17_fs2_lbuffer, chm21_fs2_lbuffer)
names(chm_focus_site2) <- c( "chm9_fs2", "chm17_fs2","chm21_fs2")
terra::writeRaster(chm_focus_site2, "chm_focus_site2_lbuffer_stack.tif",overwrite=TRUE)

chm_focus_site3 <- c( chm9_fs3_lbuffer, chm17_fs3_lbuffer, chm21_fs3_lbuffer)
names(chm_focus_site3) <- c("chm9_fs3",  "chm17_fs3", "chm21_fs3")
terra::writeRaster(chm_focus_site3, "chm_focus_site3_lbuffer_stack.tif",overwrite=TRUE)

chm_focus_site4 <- c(chm9_fs4_lbuffer, chm17_fs4_lbuffer, chm21_fs4_lbuffer)
names(chm_focus_site4) <- c( "chm9_fs4", "chm17_fs4", "chm21_fs4")
terra::writeRaster(chm_focus_site4, "chm_focus_site4_lbuffer_stack.tif",overwrite=TRUE)
