###############################################################################
#
# Extract growth rate based on CHMs (based on lidar data) within tree crowns
# in the Klausbach valley of the Berchtesgaden National Park
#
##############################################################################

library(dplyr)
library(tidyr)
library(terra)
library(lidR)
library(ggplot2)


# original artifacts masked CHM was cropped to a subarea - in this case the Klausbach valley of the reserach region
chm_valley<- rast("data/processed/closure/height_growth_in_crowns/chm_berchtesgaden_stack_1m_valley_crop.tif")
chm9 <- chm_valley[[1]]
chm21 <- chm_valley[[3]]
chm17 <- chm_valley[[2]]


# calculate vegetation height gain

crown_growth_1721<- chm21-chm17
crown_growth_1721[crown_growth_1721 <= 0] <- NA

crown_growth_917 <- chm17 - chm9
crown_growth_917[crown_growth_917 <= 0] <- NA


# smoothing of the CHM, which is than more robust for later segmentation  

kernel <- matrix(1,3,3) 
chm_smooth9 <- terra::focal(chm9, w = kernel, fun = median, na.rm = TRUE)
chm_smooth17 <- terra::focal(chm17, w = kernel, fun = median, na.rm = TRUE)
chm_smooth21 <- terra::focal(chm21, w = kernel, fun = median, na.rm = TRUE)

# single tree detection:

ttops9<- locate_trees(chm_smooth9, lmf(5)) # local maxima filter with a window size of 5 
ttops17<- locate_trees(chm_smooth17, lmf(5)) 
ttops21<- locate_trees(chm_smooth21, lmf(5))

# delineate tree crowns 

crowns9 <- dalponte2016(chm_smooth9, ttops9, th_tree = 10, max_cr=10)() # dalponte algorithm: https://rdrr.io/cran/lidR/man/its_dalponte2016.html
crowns17 <- dalponte2016(chm_smooth17, ttops17, th_tree = 10, max_cr=10)() 
crowns21 <- dalponte2016(chm_smooth21, ttops21, th_tree = 10, max_cr=10)() 
#plot(crowns9, col = pastel.colors(5000))


# terra::writeRaster(crowns9, "data/processed/closure/height_growth_in_crowns/crowns9_klausbachtal.tif")
# terra::writeRaster(crowns17, "data/processed/closure/height_growth_in_crowns/crowns17_klausbachtal.tif")
# terra::writeRaster(crowns21, "data/processed/closure/height_growth_in_crowns/crowns21_klausbachtal.tif")


# --- mask to growth within and outside of tree crowns

# 2009 - 2017

growth_in_crowns917 <- mask(crown_growth_917,crowns9)
growth_out_crowns917 <- mask(crown_growth_917, crowns9, inverse=TRUE)

# 2017 - 2021

growth_in_crowns1721 <- mask(crown_growth_1721,crowns17)
growth_out_crowns1721 <- mask(crown_growth_1721, crowns17, inverse=TRUE)


#save growth layer 

#2009-2017
terra::writeRaster(growth_in_crowns917, "data/processed/closure/height_growth_in_crowns/growth_in_crowns_917_klausbachtal.tif")
terra::writeRaster(growth_out_crowns917, "data/processed/closure/height_growth_in_crowns/growth_out_crowns_917_klausbachtal.tif")

#2017-2021
terra::writeRaster(growth_in_crowns1721, "data/processed/closure/height_growth_in_crowns/growth_in_crowns_1721_klausbachtal.tif")
terra::writeRaster(growth_out_crowns1721, "data/processed/closure/height_growth_in_crowns/growth_out_crowns_1721_klausbachtal.tif")


# --- data wrangling for plotting height growth within and outside of tree crowns

# 2009 - 2017

crown_center917 <- as.data.frame(growth_in_crowns917)
names(crown_center917) <- "height_gain"
crown_center917$type <- "crown center"
hist(crown_center917$height_gain)

quantile(crown_center917$height_gain, 0.80)/8 # 0.480875 
quantile(crown_center917$height_gain, 0.90)/8 # 0.6229999

non_crown917 <- as.data.frame(growth_out_crowns917)
names(non_crown917) <- "height_gain"
non_crown917$type <- "non-crown"
hist(non_crown917$height_gain)

growth_pixel_height_dist917 <- rbind(crown_center917, non_crown917)

# 2017 - 2021
crown_center1721 <- as.data.frame(growth_in_crowns1721)
names(crown_center1721) <- "height_gain"
crown_center1721$type <- "crown center"
hist(crown_center1721$height_gain)

quantile(crown_center1721$height_gain, 0.80)/4 # 0.5240002 
quantile(crown_center1721$height_gain, 0.90)/4 # 0.7205 


non_crown1721 <- as.data.frame(growth_out_crowns1721)
names(non_crown1721) <- "height_gain"
non_crown1721$type <- "non-crown"
hist(non_crown1721$height_gain)

growth_pixel_height_dist1721 <- rbind(crown_center1721, non_crown1721) # combine in - and out-crown growth data


saveRDS(growth_pixel_height_dist1721, file = "data/processed/closure/height_growth_in_crowns/growth_pixel_height_dist_1721.rds")
saveRDS(growth_pixel_height_dist917, file = "data/processed/closure/height_growth_in_crowns/growth_pixel_height_dist_917.rds")

# --- get Histogram of Growth-Pixel per Height class

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=17),
  legend.text = element_text(size=16),
  strip.text.x = element_text(size = 16),
  legend.position="bottom")

tiff("data/results/height_gain_within_out_crown_1721_Klausbachtal.tiff", units="in", width=12, height=8, res=300)   
ggplot(growth_pixel_height_dist1721, aes(x=height_gain, fill=type)) + geom_bar(color = "white" )+
  scale_x_binned(n.breaks = 28)+
  scale_fill_manual(values=c("#E69F00","#EE2525","#999999"))+
  theme_minimal() + labs(x = "height gain in [m]", y = "n gain pixels")+ My_Theme
dev.off()

tiff("data/results/height_gain_within_out_crown_917_Klausbachtal.tiff", units="in", width=12, height=8, res=300)   
ggplot(growth_pixel_height_dist917, aes(x=height_gain, fill=type)) + geom_bar(color = "white" )+
  scale_x_binned(n.breaks = 28)+
  scale_fill_manual(values=c("#E69F00","#EE2525","#999999"))+
  theme_minimal()+ labs(x = "height gain in [m]", y = "n gain pixels")+ My_Theme
dev.off()
