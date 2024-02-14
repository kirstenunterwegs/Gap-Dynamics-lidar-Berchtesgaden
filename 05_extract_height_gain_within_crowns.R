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

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)


chm_valley<- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m_valley_crop.tif")
chm9 <- chm_valley[[1]]
chm21 <- chm_valley[[3]]
chm17 <- chm_valley[[2]]

#extract crown growth pixels
crown_growth_1721<- chm21-chm17
crown_growth_1721[crown_growth_1721 <= 0] <- NA

crown_growth_917 <- chm17 - chm9
crown_growth_917[crown_growth_917 <= 0] <- NA

# smoothing of the CHM, which wis than more robust for later segmentation  

kernel <- matrix(1,3,3) ###
chm_smooth9 <- terra::focal(chm9, w = kernel, fun = median, na.rm = TRUE)
chm_smooth17 <- terra::focal(chm17, w = kernel, fun = median, na.rm = TRUE)
chm_smooth21 <- terra::focal(chm21, w = kernel, fun = median, na.rm = TRUE)

# single tree detetction:

ttops9<- locate_trees(chm_smooth9, lmf(5)) # local maxima filter with a window size of 5 
ttops17<- locate_trees(chm_smooth17, lmf(5)) 
ttops21<- locate_trees(chm_smooth21, lmf(5))


crowns9 <- dalponte2016(chm_smooth9, ttops9, th_tree = 10, max_cr=10)() # dalponte algorithm: https://rdrr.io/cran/lidR/man/its_dalponte2016.html
crowns17 <- dalponte2016(chm_smooth17, ttops17, th_tree = 10, max_cr=10)() 
crowns21 <- dalponte2016(chm_smooth21, ttops21, th_tree = 10, max_cr=10)() 
#plot(crowns9, col = pastel.colors(5000))

# terra::writeRaster(crowns9, "i:/Fonda/workspace/berchtesgaden/gaps/crowns9_klausbachtal.tif")
# terra::writeRaster(crowns17, "i:/Fonda/workspace/berchtesgaden/gaps/crowns17_klausbachtal.tif")
# terra::writeRaster(crowns21, "i:/Fonda/workspace/berchtesgaden/gaps/crowns21_klausbachtal.tif")

# - load crown layer

crowns9<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns9_klausbachtal.tif")
crowns17<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns17_klausbachtal.tif")
crowns21<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns21_klausbachtal.tif")

# boundaries_crown9 <- boundaries(crowns9, inner = TRUE, classes = TRUE, directions = 8) 
# boundaries_crown17 <- boundaries(crowns17, inner = TRUE, classes = TRUE, directions = 8) 
# boundaries_crown21 <- boundaries(crowns21, inner = TRUE, classes = TRUE, directions = 8) 

# --- mask to growth within, at edge and out of crowns

# 2009 - 2017

growth_in_crowns917 <- mask(crown_growth_917,crowns9)
growth_out_crowns917 <- mask(crown_growth_917, crowns9, inverse=TRUE)
# growth_edge_crowns917 <- mask(crown_growth_917, boundaries_crown17)
# growth_edge_crowns917 <- mask(growth_edge_crowns917,boundaries_crown9 )

# 2017 - 2021

growth_in_crowns1721 <- mask(crown_growth_1721,crowns17)
growth_out_crowns1721 <- mask(crown_growth_1721, crowns17, inverse=TRUE)
# growth_edge_crowns1721 <- mask(crown_growth_1721, rast(boundaries_crown21))
# growth_edge_crowns1721 <- mask(growth_edge_crowns1721,rast(boundaries_crown17), inverse=TRUE )


#save growth layer 
#2009-2017
terra::writeRaster(growth_in_crowns917, "i:/Fonda/workspace/berchtesgaden/gaps/growth_in_crowns_917_klausbachtal.tif",overwrite=TRUE)
terra::writeRaster(growth_out_crowns917, "i:/Fonda/workspace/berchtesgaden/gaps/growth_out_crowns_917_klausbachtal.tif",overwrite=TRUE)
#terra::writeRaster(growth_edge_crowns917, "i:/Fonda/workspace/berchtesgaden/gaps/growth_edge_crowns_917_klausbachtal.tif",overwrite=TRUE)
#2017-2021
terra::writeRaster(growth_in_crowns1721, "i:/Fonda/workspace/berchtesgaden/gaps/growth_in_crowns_1721_klausbachtal.tif",overwrite=TRUE)
terra::writeRaster(growth_out_crowns1721, "i:/Fonda/workspace/berchtesgaden/gaps/growth_out_crowns_1721_klausbachtal.tif",overwrite=TRUE)
#terra::writeRaster(growth_edge_crowns1721, "i:/Fonda/workspace/berchtesgaden/gaps/growth_edge_crowns_1721_klausbachtal.tif",overwrite=TRUE)


# --- get Histogram of Growth-Pixel per Height class
# 
# growth_pixel_height_dist1721 <- c(crown_growth_1721, growth_in_crowns1721, growth_out_crowns1721, growth_edge_crowns1721)
# names(growth_pixel_height_dist1721) <- c("growth all", "crown center", "non-crown", "crown edge")
# growth_pixel_height_dist1721 <- as.data.frame(growth_pixel_height_dist1721, na.rm=FALSE)

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

# crown_edge <- as.data.frame(growth_edge_crowns1721)
# names(crown_edge) <- "height_gain"
# crown_edge$type <- "crown edge"
# hist(crown_edge$height_gain)

non_crown1721 <- as.data.frame(growth_out_crowns1721)
names(non_crown1721) <- "height_gain"
non_crown1721$type <- "non-crown"
hist(non_crown1721$height_gain)

growth_pixel_height_dist1721 <- rbind(crown_center1721, non_crown1721)

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)
saveRDS(growth_pixel_height_dist1721, file = "growth_pixel_height_dist_1721.rds")
saveRDS(growth_pixel_height_dist917, file = "growth_pixel_height_dist_917.rds")


wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes"
setwd(wd)

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

tiff("height_gain_within_out_crown_1721_Klausbachtal.tiff", units="in", width=12, height=8, res=300)   
ggplot(growth_pixel_height_dist1721, aes(x=height_gain, fill=type)) + geom_bar(color = "white" )+
  scale_x_binned(n.breaks = 28)+
  #  scale_color_manual(values=c( "#E69F00", "#EE2525","#999999"))+
  scale_fill_manual(values=c("#E69F00","#EE2525","#999999"))+
  theme_minimal() + labs(x = "height gain in [m]", y = "n gain pixels")+ My_Theme
dev.off()

tiff("height_gain_within_out_crown_917_Klausbachtal.tiff", units="in", width=12, height=8, res=300)   
ggplot(growth_pixel_height_dist917, aes(x=height_gain, fill=type)) + geom_bar(color = "white" )+
  scale_x_binned(n.breaks = 28)+
  #  scale_color_manual(values=c( "#E69F00", "#EE2525","#999999"))+
  scale_fill_manual(values=c("#E69F00","#EE2525","#999999"))+
  theme_minimal()+ labs(x = "height gain in [m]", y = "n gain pixels")+ My_Theme
dev.off()

### Nach Bayrische Staatsforste Kiefernrichtlinie Zuw?che bei Fickte bis 58 cm/yr, Kiefer 49 cm/yr, Buche 42 cm/yr
### Nach Coates schwankend je nach gap opening und Baumart f?r seedlings in medium gaps bis full open: zw. 12 cm - 55 cm