#==== check focus sites ====
# collection of code snippets for gap analysis
#--------------------------

library(lidR)
library(sf)
library(dplyr)
library(terra)
require(future)
library(ForestGapR)
library(ggplot2)
library(ForestTools)

#---- load chms and focus sites ------

focus_sites_large <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/focus_sites_buffer.shp")
focus_sites <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/walddynamik_intensiv_fp.shp")
chm9 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/berchtesgaden_2009_chm_pf-sc02.tif")
chm17<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_pf-sc02.tif")
chm21<- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/berchtesgaden_2021_chm_fix.tif")


#---- clip CHMS to focus sites ------

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

#crop to focus sites large
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


# plot chms on focus sites through years
par(mfrow=c(4,3))
plot(chm9_fs1, main ="CHM 9 focus site 1") #site1
plot(chm17_fs1, main ="CHM 17 focus site 1")
plot(chm21_fs1, main ="CHM 21 focus site 1")
plot(chm9_fs2, main ="CHM 9 focus site 2") #site2
plot(chm17_fs2, main ="CHM 17 focus site 2")
plot(chm21_fs2, main ="CHM 21 focus site 2")
plot(chm9_fs3, main ="CHM 9 focus site 3") #site3
plot(chm17_fs3, main ="CHM 17 focus site 3")
plot(chm21_fs3, main ="CHM 21 focus site 3")
plot(chm9_fs4, main ="CHM 9 focus site 4") #site4
plot(chm17_fs4, main ="CHM 17 focus site 4")
plot(chm21_fs4, main ="CHM 21 focus site 4")

#---- Forest Gap detection ----


#---- Using ForestGap package -------------------------------
#convert chm from Spat into raster

threshold <- 5
size <- c(25,10^4) #m2
gaps_9_fs1 <- getForestGaps(chm_layer = raster::raster(chm9_fs1), threshold = threshold, size = size)
gaps_17_fs1 <- getForestGaps(chm_layer = raster::raster(chm17_fs1), threshold = threshold, size = size)
gaps_21_fs1 <- getForestGaps(chm_layer = raster::raster(chm21_fs1), threshold = threshold, size = size)

gap_changes_fs1_1 <- GapChangeDec(gap_layer1 = gaps_9_fs1, gap_layer2 = gaps_17_fs1)
gap_changes_fs1_2 <- GapChangeDec(gap_layer1 = gaps_17_fs1, gap_layer2 = gaps_21_fs1)

par(mfrow=c(3,1))
plot(chm9_fs1, main="Forest Canopy Gap 2009")
grid()
plot(gaps_9_fs1, col="red", add=TRUE, legend=FALSE)
plot(chm17_fs1, main="Forest Canopy Gap 2017")
grid()
plot(gaps_17_fs1, col="red", add=TRUE, legend=FALSE)
plot(chm21_fs1, main="Forest Canopy Gap 2021")
grid()
plot(gaps_21_fs1, col="red", add=TRUE, legend=FALSE)

par(mfrow=c(2,1))
plot(gap_changes_fs1_1, col="red",  main="Forest Canopy Gap Change 2009-2017")
plot(gap_changes_fs1_2, col="red",  main="Forest Canopy Gap Change 2017-2021")
par(mfrow=c(1,1))

gaps_9_fs1_polygon <- GapSPDF(gaps_9_fs1)
gaps_9_fs1_buffer <- buffer(vect(gaps_9_fs1_polygon), width = 20) #do buffers have the same id as pixel values?
# yes buffer ID and pixel values of gaps are the same

#writeVector(gaps_9_fs1_buffer, "i:/Fonda/workspace/berchtesgaden/focus_sites/buffer_trial.shp")
#writeRaster(gaps_9_fs1, "i:/Fonda/workspace/berchtesgaden/focus_sites/gaps_trial.tif")

buffer_area <- mask(chm9_fs1_l, gaps_9_fs1_buffer) #need bigger extent of chm to cover the whole area of the buffer!
buffer_area <- mask(buffer_area, vect(gaps_9_fs1_polygon), inverse=TRUE)
plot(buffer_area)

mean_canopy_gaps9_fs1 <-  terra::extract(buffer_area, gaps_9_fs1_buffer, list=FALSE, fun=mean, na.rm=TRUE)
head(mean_canopy_gaps9_fs1)

plot(gaps_9_fs1)
lines(gaps_9_fs1_buffer)

# when mean of surrounding canopy is <10, than I can use the ID to delete gap from gap layer


# try with focus site 2
threshold <- 5
size <- c(25,10000) #m2

gaps_9_fs2 <- getForestGaps(chm_layer = raster::raster(chm9_fs2), threshold = threshold, size = size)
gaps_17_fs2 <- getForestGaps(chm_layer = raster::raster(chm17_fs2), threshold = threshold, size = size)
gaps_21_fs2 <- getForestGaps(chm_layer = raster::raster(chm21_fs2), threshold = threshold, size = size)

gap_changes_fs2_1 <- GapChangeDec(gap_layer1 = gaps_9_fs2, gap_layer2 = gaps_17_fs2)
gap_changes_fs2_2 <- GapChangeDec(gap_layer1 = gaps_17_fs2, gap_layer2 = gaps_21_fs2)


par(mfrow=c(3,1))
plot(chm9_fs2, main="Forest Canopy Gap 2009")
grid()
plot(gaps_9_fs2, col="red", add=TRUE, legend=FALSE)
plot(chm17_fs2, main="Forest Canopy Gap 2017")
grid()
plot(gaps_17_fs2, col="red", add=TRUE, legend=FALSE)
plot(chm21_fs2, main="Forest Canopy Gap 2021")
grid()
plot(gaps_21_fs2, col="red", add=TRUE, legend=FALSE)

par(mfrow=c(2,1))
plot(gap_changes_fs2_1, col="red",  main="Forest Canopy Gap Change 2009-2017")
plot(gap_changes_fs2_2, col="red",  main="Forest Canopy Gap Change 2017-2021")

par(mfrow=c(1,1))

# check surrounding canopy height of gaps
gaps_17_fs2_polygon <- GapSPDF(gaps_17_fs2)
gaps_17_fs2_buffer <- buffer(vect(gaps_17_fs2_polygon), width = 20)
buffer_area <- mask(chm17_fs2_l, gaps_17_fs2_buffer) #need bigger extent of chm to cover the whole area of the buffer!
buffer_area <- mask(buffer_area, vect(gaps_17_fs2_polygon), inverse=TRUE)
plot(buffer_area)

mean_canopy_gaps17_fs2 <-  terra::extract(buffer_area, gaps_17_fs2_buffer, list=FALSE, fun=mean, na.rm=TRUE)
head(mean_canopy_gaps17_fs2)

# exclude gaps with mean z <10
#when z in mean_canopy_gaps < 10 than id_value in Raster = NA
# create loop which checks if there are surrounding canopys < 10 , if no stop, if yes continue #####!!!!!!!!!!!!!!!!!!!!!!!!

mean_canopy_filter <- mean_canopy_gaps17_fs2[!(mean_canopy_gaps17_fs2$Z > 10),]
mean_canopy_filter$replace <- NA
gaps_17_fs2_filtered <- raster::subs(gaps_17_fs2, mean_canopy_filter, by=1, which="replace" )

par(mfrow=c(2,1))
plot(gaps_17_fs2, col="red",  main="Forest Canopy Gaps 2017 non-filtered")
plot(gaps_17_fs2_filtered, col="red",  main="Forest Canopy Gaps 2017 filtered")
par(mfrow=c(1,1))

#--- Gap detection own approach -----

chm9_fs3_below5 <- chm9_fs3
chm9_fs3_below5[chm9_fs3_below5 > 5] <- 0
patches_below5_fs3 <- patches(chm9_fs3_below5, directions=8, zeroAsNA=TRUE)

plot(chm9_fs3)
plot(patches_below5_fs3, col="red", add=TRUE, legend=FALSE)


chm9_fs2_below5 <- chm9_fs2
chm9_fs2_below5[chm9_fs2_below5 > 5] <- 0
patches_below5_fs2 <- patches(chm9_fs2_below5, directions=8, zeroAsNA=TRUE, allowGaps=TRUE)

patches_below5_fs2_sub <- patches_below5_fs2
patches_below5_fs2_sub[patches_below5_fs2_sub != 1] <- 0 
plot(patches_below5_fs2_sub) # there is one huge gap -> is this still a gap?

plot(chm9_fs2)
plot(patches_below5_fs2, col="red", add=TRUE, legend=FALSE)

#------------------------------------------------------------------------------
#--- following Leitold et al 2021 -----

par(mfrow=c(1,1))

#simple difference

diff_917_fs2 <- chm9_fs2 - chm17_fs2
diff_1721_fs2 <- chm17_fs2 - chm21_fs2
#diff_921_fs2 <- chm9_fs2 - chm21_fs2 #falsch da reverse!

diff_917_fs2 <- chm17_fs2 - chm9_fs2
diff_1721_fs2 <- chm21_fs2 - chm17_fs2
diff_921_fs2 <- chm21_fs2 - chm9_fs2

# mask differences with gap layer!!!!

# classify forest canopy change direction
m = c(-50, -1, 1, -1, 1, 2, 1, 50, 3)# 1= <-1 (loss); 2= -1/1 (steady); 3= >1 (gain)
rclmat = matrix(m, ncol=3, byrow=TRUE)
diff_917_fs2_class = terra::classify(diff_917_fs2, rclmat, include.lowest=TRUE) # Reclassify the raster layer
diff_1721_fs2_class = terra::classify(diff_1721_fs2, rclmat, include.lowest=TRUE)
diff_921_fs2_class = terra::classify(diff_921_fs2, rclmat, include.lowest=TRUE) 

par(mfrow=c(3,2))
hist(diff_917_fs2, main="Canopy Change 2009-2017", breaks = seq(-50,50,1),)
plot(diff_917_fs2_class, main="Canopy Change 2009-2017")
hist(diff_1721_fs2, main="Canopy Change 2017-2021", breaks = seq(-55,50,1),)
plot(diff_1721_fs2_class,  main="Canopy Change 2017-2021",)
hist(diff_921_fs2, main="Canopy Change 2009-2021", breaks = seq(-50,50,1),)
plot(diff_921_fs2_class,main="Canopy Change 2009-2021" )

par(mfrow=c(1,1))

#--- identify change trajectory------------------
# mask and boolean?


#--- identify closure mechanism 

#identify treetops

winFunction <- function(x){x * 0.06 + 0.5} #function for variable window radius
minHgt <- 5 #minimum tree height (treetops below not detected)

trees_9_fs2 <- vwf(raster::raster(chm9_fs2), winFunction, minHgt)
trees_21_fs2 <- vwf(raster::raster(chm21_fs2), winFunction, minHgt)

#delineate crown objects
crowns_9_fs2 <- mcws(trees_9_fs2, raster::raster(chm9_fs2), minHgt, format ="raster",verbose = TRUE)
crowns_21_fs2 <- mcws(trees_21_fs2, raster::raster(chm21_fs2), minHgt, format ="raster",verbose = TRUE)

plot(chm9_fs2)
plot(crowns_9_fs2, col="purple", add= TRUE, legend=FALSE)
plot(chm21_fs2)
plot(crowns_21_fs2, col="purple", add= TRUE, legend=FALSE)

#terra::writeRaster(crowns_9_fs2, "i:/Fonda/workspace/berchtesgaden/lidar/2009/crowns_2009_focussite2.tif")
#terra::writeRaster(crowns_21_fs2, "i:/Fonda/workspace/berchtesgaden/lidar/2021/crowns_2021_focussite2.tif")

boundaries_crown_21_fs2 <- boundaries(crowns_21_fs2, inner = FALSE, classes = TRUE, directions = 8) 
boundaries_crown_9_fs2 <- boundaries(crowns_9_fs2, inner = FALSE, classes = TRUE, directions = 8) 
# using outer boundary!!!

#extract crown growth pixels
crown_growth_921 <- diff_921_fs2_class
crown_growth_921[crown_growth_921 < 3] <- NA

# mask to growth within, at edge and out of crowns
growth_in_crowns <- mask(crown_growth_921,rast(boundaries_crown_9_fs2))
growth_out_crowns <- mask(crown_growth_921, rast(boundaries_crown_21_fs2), inverse=TRUE)
growth_edge_crowns <- mask(crown_growth_921, rast(boundaries_crown_21_fs2))
growth_edge_crowns <- mask(growth_edge_crowns,rast(boundaries_crown_9_fs2), inverse=TRUE )

#plotting
plot(chm9_fs2,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
#plot(chm21_fs2,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
plot(boundaries_crown_9_fs2, col="purple", add= TRUE, legend=FALSE)
plot(boundaries_crown_21_fs2, add=TRUE, legend=FALSE, alpha = 0.5)

#plot growth in, at edge and out of crowns
plot(chm9_fs2,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
plot(crowns_9_fs2, add=TRUE, col="purple", legend=FALSE)
plot(boundaries_crown_21_fs2, add=TRUE, legend=FALSE, alpha = 0.5)
#plot(growth_out_crowns, add=TRUE, col="green", alpha=0.5, legend=FALSE)
plot(growth_edge_crowns, add=TRUE, col="red", alpha=0.5, legend=FALSE)
plot(growth_in_crowns, add=TRUE, col="white", alpha=0.8, legend=FALSE)

#plot only growth categories
plot(growth_out_crowns, col="grey", legend=FALSE)
plot(growth_edge_crowns, add=TRUE, col="red", alpha=0.8, legend=FALSE)
plot(growth_in_crowns, add=TRUE, col="orange", alpha=0.8, legend=FALSE)

# ------  repeat with growth numbers themselves! -------------------
#extract crown growth pixels
crown_growth_921 <- diff_921_fs2
crown_growth_921[crown_growth_921 < 1] <- NA

# mask to growth within, at edge and out of crowns
growth_in_crowns <- mask(crown_growth_921,rast(boundaries_crown_9_fs2))
growth_out_crowns <- mask(crown_growth_921, rast(boundaries_crown_21_fs2), inverse=TRUE)
growth_edge_crowns <- mask(crown_growth_921, rast(boundaries_crown_21_fs2))
growth_edge_crowns <- mask(growth_edge_crowns,rast(boundaries_crown_9_fs2), inverse=TRUE )

#plot growth in, at edge and out of crowns in full numbers
plot(chm9_fs2,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), legend=FALSE)
plot(crowns_9_fs2, add=TRUE, col="purple", legend=FALSE)
plot(boundaries_crown_21_fs2, add=TRUE, legend=FALSE, alpha = 0.5)
#plot(growth_out_crowns, add=TRUE, col="green", alpha=0.5, legend=FALSE)
plot(growth_edge_crowns, add=TRUE, alpha=0.5, legend=FALSE)
plot(growth_in_crowns, add=TRUE, alpha=0.8)

# --- get Histogram of Growth-Pixel per Height class

growth_pixel_height_dist <- c(crown_growth_921, growth_in_crowns, growth_edge_crowns, growth_out_crowns)
names(growth_pixel_height_dist) <- c("growth all", "crown center", "crown edge", "non-crown")
growth_pixel_height_dist <- as.data.frame(growth_pixel_height_dist)

crown_center <- as.data.frame(growth_in_crowns)
names(crown_center) <- "height_gain"
crown_center$type <- "crown center"
hist(crown_center$height_gain)

crown_edge <- as.data.frame(growth_edge_crowns)
names(crown_edge) <- "height_gain"
crown_edge$type <- "crown edge"
hist(crown_edge$height_gain)

non_crown <- as.data.frame(growth_out_crowns)
names(non_crown) <- "height_gain"
non_crown$type <- "non-crown"
hist(non_crown$height_gain)

growth_pixel_height_dist <- rbind(crown_center, crown_edge, non_crown)

ggplot(growth_pixel_height_dist, aes(x=height_gain, fill=type)) + geom_bar(color = "white" )+
  scale_x_binned(n.breaks = 28)+
  #  scale_color_manual(values=c( "#E69F00", "#EE2525","#999999"))+
  scale_fill_manual(values=c("#E69F00","#EE2525","#999999"))+
  theme_minimal()


#Overlapping 2018 and 2020 canopy tree objects were used to estimate height gains within and
# adjacent to the 2018 crown extents (Figure 2). Gains within the 2018 extent of crown objects were
# attributed to vertical growth, while height gains around crown edges were attributed to lateral expansion.

#height gain threshold of 4 m was used to differentiate canopy closure via vertical growth (height gain £
#        4 m) versus horizontal infilling by adjacent trees (height gain > 4 m),

#--- identify trajectories ----------------------

#--- group change pixels into vertical canopy height bins (here all)

change_vertical_bins <- c(chm9_fs2, diff_921_fs2)
names(change_vertical_bins) <- c("Canopy_height", "height_change")
change_vertical_bins <- as.data.frame(change_vertical_bins)

change_vertical_bins<-change_vertical_bins %>% 
  mutate(Canopy_height_bins = cut(Canopy_height, breaks = c(0,5,10,15,20,25,30, 35,40,45,50)))
head(change_vertical_bins,10)

ggplot(change_vertical_bins, aes(x=height_change, y=Canopy_height_bins)) + geom_boxplot() + theme_minimal()

#-- classify canopy cover per height bin ----------------------LATER - use dataframes and ggplot

library(RColorBrewer)

hist(chm9_fs2)
hist(chm17_fs2)
hist(chm21_fs2)

chm9_fs2[chm9_fs2 <5] <- NA # exclude all pixels below 5m (=threshold of canopy)
chm17_fs2[chm17_fs2 <5] <- NA 
chm21_fs2[chm21_fs2 <5] <- NA

fs2_9hist <- hist(chm9_fs2, breaks = seq(5,50,1), plot = FALSE)
fs2_17hist <- hist(chm17_fs2, breaks = seq(5,50,1), plot = FALSE)
fs2_21hist <- hist(chm21_fs2, breaks = seq(5,50,1), plot = FALSE)

blue <- rgb(0, 0, 1, alpha=0.5)
red <- rgb(1, 0, 0, alpha=0.5)
green <- rgb(0, 1, 0, alpha=0.5)

barplot(fs2_9hist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", col=blue)
barplot(fs2_17hist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", add=TRUE, col=red)
barplot(fs2_21hist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", add=TRUE, col=green)

