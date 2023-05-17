###############################################
# Identify closure mechanism
##############################################
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(lidR)
library(ForestGapR)
library(ggplot2)
library(ForestTools)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)
library(ggExtra)


# --- load layers ----

gaps2009 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/gaps/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/gaps/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2021/gaps/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE) 


wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

chm_stack <- rast("chm_berchtesgaden_stack_1m_maskedartifacts.tif")
chm9 <- chm_stack[[1]]
chm17<- chm_stack[[2]]
chm21<- chm_stack[[3]]

exp_clo_917 <- rast("expansion_closure_917_cn2cr2_mmu800.tif") 
exp_clo_1721<- rast("expansion_closure_1721_cn2cr2_mmu800.tif") 

# subset to closure areas
clo_1721 <- subst(exp_clo_1721, 2, NA)
clo_917 <- subst(exp_clo_917, 2, NA)

# --- explore chm difference ---

diff917 <- chm17 - chm9
diff1721 <- chm21 - chm17

# hist(diff917)
# hist(diff1721)

diff917_df<- as.data.frame(diff917)
diff1721_df<- as.data.frame(diff1721)

diff917_df$timestep <- "9-17"
diff1721_df$timestep <- "17-21"

names(diff1721_df) <- c("diff", "timestep")
names(diff917_df) <- c("diff", "timestep")

diff_all <- rbind(diff917_df, diff1721_df)

hist_diff_100 <- ggplot(diff_all, aes(diff, fill=timestep)) + geom_histogram(aes(bins= 100)) + theme_classic() + My_Theme

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes"
setwd(wd)

tiff("hist_chm_diffs_all.tiff", units="in", width=12, height=8, res=300)   
ggplot(diff_all, aes(diff, fill=timestep)) + geom_histogram(bins= 120 ) + 
  theme_classic() + My_Theme + geom_vline(xintercept = 0) +labs(x= "CHM difference", y= "No of pixels")
dev.off()


#convert into df and analyze difference for closure areas
growth_917.df <- subset(diff917_df, chm17 > 0)
grwoth_1721.df<- subset(diff1721_df, chm21 > 0)

# ggplot(growth_917, aes(chm17)) + geom_histogram()
# ggplot(growth_917, aes(chm17)) + geom_boxplot() 
# 
# ggplot(growth_1721, aes(chm17)) + geom_histogram()
# ggplot(growth_1721, aes(chm17)) + geom_boxplot() 


################################ growth extraction Klausbachtal ######################################################

chm_valley<- rast("i:/Fonda/workspace/berchtesgaden/gaps/chm_berchtesgaden_stack_1m_valley_crop.tif")
chm9 <- chm_valley[[1]]
chm21 <- chm_valley[[3]]
chm17 <- chm_valley[[2]]

#extract crown growth pixels
crown_growth_1721<- chm21-chm17
crown_growth_1721[crown_growth_1721 <= 0] <- NA

crown_growth_917 <- chm17 - chm9
crown_growth_917[crown_growth_917 <= 0] <- NA

#    Ich würde das CHM smoothen, dann ist es robuster was die spätere Segmentierung anbelangt: 
  
kernel <- matrix(1,3,3) ###
chm_smooth9 <- terra::focal(chm9, w = kernel, fun = median, na.rm = TRUE)
chm_smooth17 <- terra::focal(chm17, w = kernel, fun = median, na.rm = TRUE)
chm_smooth21 <- terra::focal(chm21, w = kernel, fun = median, na.rm = TRUE)
  
# Einzelbaumdetektion:
  
ttops9<- locate_trees(chm_smooth9, lmf(5)) # local maxima filter with a window size of 5 (play around with windows size)
ttops17<- locate_trees(chm_smooth17, lmf(5)) 
ttops21<- locate_trees(chm_smooth21, lmf(5))
  
#  Einzelbaumdetektion plotten

# col <- height.colors(50)
# plot(chm17, col = col)
# plot(sf::st_geometry(ttops), add = T, pch =3)

#  Einzelbaumsegmentierung
#ggf. downsampling im source code ändern!!!

crowns9 <- dalponte2016(chm_smooth9, ttops9, th_tree = 10, max_cr=10)() # dalponte algorithm: https://rdrr.io/cran/lidR/man/its_dalponte2016.html
crowns17 <- dalponte2016(chm_smooth17, ttops17, th_tree = 10, max_cr=10)() 
crowns21 <- dalponte2016(chm_smooth21, ttops21, th_tree = 10, max_cr=10)() 
#plot(crowns9, col = pastel.colors(5000))

# terra::writeRaster(crowns9, "i:/Fonda/workspace/berchtesgaden/gaps/crowns9_klausbachtal.tif")
# terra::writeRaster(crowns17, "i:/Fonda/workspace/berchtesgaden/gaps/crowns17_klausbachtal.tif")
# terra::writeRaster(crowns21, "i:/Fonda/workspace/berchtesgaden/gaps/crowns21_klausbachtal.tif")

#load crown layer

crowns9<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns9_klausbachtal.tif")
crowns17<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns17_klausbachtal.tif")
crowns21<- rast("i:/Fonda/workspace/berchtesgaden/gaps/crowns21_klausbachtal.tif")

# boundaries_crown9 <- boundaries(crowns9, inner = TRUE, classes = TRUE, directions = 8) 
# boundaries_crown17 <- boundaries(crowns17, inner = TRUE, classes = TRUE, directions = 8) 
# boundaries_crown21 <- boundaries(crowns21, inner = TRUE, classes = TRUE, directions = 8) 

# mask to growth within, at edge and out of crowns
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
### Nach Bayrische Staatsforste Kiefernrichtlinie Zuwäche bei Fickte bis 58 cm/yr, Kiefer 49 cm/yr, Buche 42 cm/yr
### Nach Coates schwankend je nach gap opening und Baumart für seedlings in medium gaps bis full open: zw. 12 cm - 55 cm

####################################################### classify vertical and horizontal closure ##################################

#extract vegetation difference per closure area
clo_growth_917 <- mask(diff917, clo_917, inverse=FALSE, maskvalues=NA, updatevalue=NA)
clo_growth_1721 <- mask(diff1721, clo_1721, inverse=FALSE, maskvalues=NA, updatevalue=NA)

 terra::writeRaster(clo_growth_917, "closure_area_growth_917.tif", overwrite=TRUE)
 terra::writeRaster(clo_growth_1721, "closure_area_growth_1721.tif", overwrite=TRUE)

clo_growth_917 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/closure_area_growth_917.tif")
clo_growth_1721 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/closure_area_growth_1721.tif")

# --- reclassify into vertical and horizontal closure---

# need to differ between both timesteps due to differnet growing periods

gap_closure_mechanism917 <- function(diff_closure_layer){       # 0.5 m * 8 = 4m height gain
  closure_mechanism <- rast() #create empty raster to classify
  ext(closure_mechanism) <- ext(diff_closure_layer)
  res(closure_mechanism) <- res(diff_closure_layer)
  crs(closure_mechanism) <- crs(diff_closure_layer)
  # classify change group
  closure_mechanism[diff_closure_layer > 4  ] <- 1 #horizontal closure (crown plasticity)
  closure_mechanism[diff_closure_layer <= 4 & diff_closure_layer >0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism1721 <- function(diff_closure_layer){      # 0.5 m * 4 = 2m height gain
  closure_mechanism <- rast() #create empty raster to classify
  ext(closure_mechanism) <- ext(diff_closure_layer)
  res(closure_mechanism) <- res(diff_closure_layer)
  crs(closure_mechanism) <- crs(diff_closure_layer)
  # classify change group
  closure_mechanism[diff_closure_layer > 2 ] <- 1 #horizontal closure (crown plasticity)
  closure_mechanism[diff_closure_layer <= 2 & diff_closure_layer >0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism917 <- gap_closure_mechanism917(clo_growth_917)
gap_closure_mechanism1721 <- gap_closure_mechanism1721(clo_growth_1721)

 wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
 setwd(wd)
terra::writeRaster(gap_closure_mechanism917, "gap_closure_mechanism917.tif", overwrite=TRUE)
terra::writeRaster(gap_closure_mechanism1721, "gap_closure_mechanism1721.tif", overwrite=TRUE)


###----------- analyze closure per gap (size), elevation, aspect, management and forest type ----------- ###

# #crop to subarea for trial with df aggregation
# focus_sites_large <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/focus_sites_buffer.shp")
# fs3_l <- subset(focus_sites_large, focus_sites_large$Id == 3)
# 
# gap_closure_mechanism917_fs3 <- crop(gap_closure_mechanism917, fs3_l, snap="near", mask=TRUE)
# gaps2009_fs3 <- crop(gaps2009, fs3_l, snap="near", mask=TRUE)

# #### trial focus site ####
# 
# gap_closure_mechanism_stack <- c(gap_closure_mechanism917_fs3, gaps2009_fs3)
# gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
# #exclude pixels without gap (and hence closure):
# gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
# names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")
# 
# # aggregate closure and gap information
# 
# gap_clo_per_id <-  gap_closure_mechanism_stack.df %>% group_by(gap_id) %>%
#   count(closure_mechanism) %>% 
#   mutate(gap_area = sum(n))
# 
# gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) # calculate share of closure mechanism
# 
# gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
# gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism)
# 
# ggplot(gap_clo_per_id_nona, aes(x=gap_area, y=closure_share, col=closure_mechanism)) + geom_point() # plot gap size vs. closure share
# 
# #aggregate closure share
# gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
#   mutate(closure_share_sum = sum(closure_share))
# 
# ggplot(gap_clo_per_id_nona, aes(x=gap_area, y=closure_share_sum)) + geom_point() # plot gap size vs. closure share
# 
# #join with aspect and gap complexity information 
# stats_aspect <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_aspect.rds")
# stats_aspect_2009 <- subset(stats_aspect, year == 2009)
# #merge p:a ratio and aspect to closure information
# gap_clo_per_id_nona_aspect <- merge(x = gap_clo_per_id_nona, y = stats_aspect_2009[ , c("gap_id","aspect", "pa_ratio")], by = "gap_id", all.x=TRUE)
# 
# ggplot(gap_clo_per_id_nona_aspect, aes(x=gap_area, y=closure_share_sum, colour=aspect)) + geom_point()

### --- for whole NP --- ###

####################################################################################################################
# prepare closure mechanism dfs per timestep
###################################################################################################################

# 2009 - 2017
###################################################################################################################

#merge closure mechanism with gaps
gap_closure_mechanism917 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism917.tif")
gaps2009 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2009/gaps/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

gap_closure_mechanism_stack <- c(gap_closure_mechanism917, gaps2009)
gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

saveRDS(gap_closure_mechanism_stack.df,"i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism_pergap_917.rds" )
gap_closure_mechanism_stack.df <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism_pergap_917.rds")

# aggregate closure and gap information - prepare df for plotting
gap_clo_per_id <-  gap_closure_mechanism_stack.df %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 7 gaps do not experience any closure from 2009-2017

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
        closure_mechanism = as.factor(recode(closure_mechanism,
                                             `0`="no closure", 
                                             `1`="lateral closure",
                                             `2`="vertical closure")))

##### if I want to include the no closure pixels, I have to disable following line: !!!!!!!!!!!!!!
#exclude no closure shares
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))

gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

#gap_clo_per_id_nona$gap_area <- as.numeric(gap_clo_per_id_nona$gap_area)

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_917<-gap_clo_per_id_nona %>% 
  mutate(gap_area_bins = (cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000))))%>% 
           mutate(gap.size = as.factor(recode(gap_area_bins,
                                              `(399,500]`="400-500",
                                              `(500,1e+03]`="500-1000",
                                              `(1e+03,5e+03]`="1000-5000",
                                              `(5e+03,1e+04]`="5000-10,000",
                                              `(1e+04,5e+04]`="10,000-50,000",
                                              `(5e+04,4.5e+05]`="50,000-450,000")))


# 2017- 2021
###################################################################################################################

#merge closure mechanism with gaps
gap_closure_mechanism1721 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism1721.tif")
gaps2017 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/gaps/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721, gaps2017)
gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")

saveRDS(gap_closure_mechanism_stack.df_1721,"i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism_pergap_1721.rds" )
gap_closure_mechanism_stack.df_1721 <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism_pergap_1721.rds" )

# aggregate closure and gap information - prepare df for plotting
gap_clo_per_id <-  gap_closure_mechanism_stack.df_1721 %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify number of gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 162 gaps do not experience any closure from 2009-2017

gaps_no_closing <- subset(gap_clo_per_id, contraction %in% 1)
range(gaps_no_closing$gap_area) # range 404-9446
#join with forest information
stats_ftype <- subset(readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_ftype.rds"), year %in% 2017)
gaps_no_closing.forest <- merge(x = gaps_no_closing, y = stats_ftype[ , c("gap_id","forest_type")], by = c("gap_id"), all.x=TRUE)
gaps_no_closing.forest %>% group_by(forest_type) %>% summarise(n = n()) # ~ 65% in coniferous or larch-dominated gaps
#---

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
         closure_mechanism = as.factor(recode(closure_mechanism,
                                              `0`="no closure", 
                                              `1`="lateral closure",
                                              `2`="vertical closure")))





#exclude no closure shares
##### if I want to include the no closure pixels, I have to disable following line: !!!!!!!!!!!!!!
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))

gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))



# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_1721<-gap_clo_per_id_nona %>% 
  mutate(gap_area_bins = (cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))


### -------------------------------merge 9-17 and 17-21 dfs for comparison ----------------------- ###
gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"

#define quantiles of interest
q = c(.25, .5, .75)

#calculate quantiles by grouping variable
gap_clo_per_id_nona_1721 %>%
  group_by(closure_mechanism) %>%
  summarize(quant25 = quantile(closure_share, probs = q[1]), 
            quant50 = quantile(closure_share, probs = q[2]),
            quant75 = quantile(closure_share, probs = q[3]))


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
#rearrange timestep labels
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))


# prepare plotting
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 28),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=20),
  strip.text.x = element_text(size = 16),
  legend.position="bottom")

require(scales)

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes/"
setwd(wd)


# ----- closure across NP ------------

clo_NP <- gap_clo_NP_91721 %>% group_by(timestep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

library(viridis)
#plot closure according to gap size bins
tiff("gap_closure-no_clo_all_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=closure_share, fill=timestep)) + geom_boxplot() + 
  theme_minimal()+ coord_flip() +
  theme(legend.position="bottom") +My_Theme  +    scale_fill_viridis(discrete = TRUE) 
dev.off()

tiff("gap_closure-no_clo_mechanism_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + geom_boxplot() +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism") +
  facet_grid(~timestep)
dev.off()


# plot gap size vs. closure share per mechanism

# 9-17
pm917 <- ggplot(subset(gap_clo_NP_91721, timestep %in% c("9-17")), aes(x=gap_area, y=closure_share, col=closure_mechanism)) + geom_point()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic()+  scale_color_brewer(palette="Dark2") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
tiff("gap_closure-no_clo_mechanism_917.tiff", units="in", width=12, height=8, res=300)
ggMarginal(pm917, type="boxplot", groupColour = TRUE, groupFill = TRUE)
dev.off()

# 17-21
pm1721 <- ggplot(subset(gap_clo_NP_91721, timestep %in% c("17-21")), aes(x=gap_area, y=closure_share, col=closure_mechanism)) + geom_point()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic()+  scale_color_brewer(palette="Dark2") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
tiff("gap_closure-no_clo_mechanism1721.tiff", units="in", width=12, height=8, res=300)
ggMarginal(pm1721, type="boxplot", groupColour = TRUE, groupFill = TRUE)
dev.off()

# complete closure share vs. gap size

#9-17
p917<- ggplot(subset(gap_clo_NP_91721, timestep %in% c("9-17")), aes(x=gap_area, y=closure_share_sum)) + geom_point() + 
                theme_classic() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + My_Theme+
  labs(x=expression ("gap size in"(m^2)), y= "% of gap area closing") 


tiff("gap_closure-no_clo_all_917.tiff", units="in", width=12, height=8, res=300)
ggMarginal(p917, type="histogram")
#ggMarginal(p, type="boxplot")
dev.off()

#17-21
p1721<- ggplot(subset(gap_clo_NP_91721, timestep %in% c("17-21")), aes(x=gap_area, y=closure_share_sum)) + geom_point() + theme_classic() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + My_Theme+
  labs(x=expression ("gap size in"(m^2)), y= "% of gap area closing") #+ facet_grid(~timestep)
tiff("gap_closure-no_clo_all_1721.tiff", units="in", width=12, height=8, res=300)
ggMarginal(p1721, type="histogram")
dev.off()



# --- closure by aspect---

#join with aspect information
stats_aspect <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_aspect.rds")
# recode years to timesteps for merges
stats_aspect <- stats_aspect %>% mutate( timestep = as.factor(recode(year,
                                                     `2009`="9-17", 
                                                     `2017`="17-21")))


#merge aspect to closure information
gap_clo_NP_91721_aspect <- merge(x = gap_clo_NP_91721, y = stats_aspect[ , c("gap_id","aspect", "pa_ratio", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)


# closure share per aspect 
tiff("gap_closure_aspect-no_clo_all.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721_aspect, aes(x=gap.size , y=closure_share_sum, fill=timestep)) + geom_boxplot() + facet_wrap(~aspect) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing") +    scale_fill_viridis(discrete = TRUE) 
dev.off()

#calculate average closure share per aspect
clo_per_aspect <- gap_clo_NP_91721_aspect %>% group_by(aspect, timestep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

# closure mechanism according to aspect 
tiff("gap_closure_mechanism_aspect-no_clo.917.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_aspect, timestep %in% c("9-17")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + geom_boxplot() + facet_wrap(~aspect) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

tiff("gap_closure_Mechanism_aspect-no_clo-1721.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_aspect, timestep %in% c("17-21")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + geom_boxplot() + facet_wrap(~aspect) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()



# --- closure by forest type ---

#join with forest information
stats_ftype <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_ftype.rds")
# recode years to timesteps for merges
stats_ftype <- stats_ftype %>% mutate( timestep = as.factor(recode(year,
                                                                     `2009`="9-17", 
                                                                     `2017`="17-21")))

#merge aspect to closure information
gap_clo_NP_91721_ftype <- merge(x = gap_clo_NP_91721, y = stats_ftype[ , c("gap_id","forest_type", "pa_ratio", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)


# add empty forest type to NA category and label no information
gap_clo_NP_91721_ftype$forest_type [gap_clo_NP_91721_ftype$forest_type == "" ] <- "no information" # trigger NA for blank fields
gap_clo_NP_91721_ftype$forest_type <- as.character(gap_clo_NP_91721_ftype$forest_type)
gap_clo_NP_91721_ftype$forest_type[is.na(gap_clo_NP_91721_ftype$forest_type)] <- "no information" # exchange NAs <- "no information"
gap_clo_NP_91721_ftype$forest_type <- as.factor(gap_clo_NP_91721_ftype$forest_type)

ggplot(gap_clo_NP_91721_ftype, aes(x=gap_area, y=closure_share, colour=closure_mechanism)) + geom_point() + facet_wrap(~forest_type) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_minimal() +
  theme(legend.position="bottom") 

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes/"
setwd(wd)
library(viridis)

# closure share per forest type 
tiff("gap_closure_ftype-no_clo_all.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721_ftype, aes(x=gap.size , y=closure_share_sum, fill=timestep)) + geom_boxplot() + facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing") +   scale_fill_viridis(discrete = TRUE) 
dev.off()

#calculate average closure share per forest type
clo_per_ftype <- gap_clo_NP_91721_ftype %>% group_by(forest_type, timestep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

clo_per_ftype.sizebins <- gap_clo_NP_91721_ftype %>% group_by(forest_type, timestep, closure_mechanism, gap.size) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

# prepare plotting
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 24),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 18),
  axis.title.y = element_text(size = 24),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=24),
  strip.text.x = element_text(size = 24),
  legend.position="bottom")

# closure mechanism according to aspect 
tiff("gap_closure_mechanism_ftype-no_clo.917.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_ftype, timestep %in% c("9-17")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x=expression ("gap area bins in"(m^2)), y= "% of gap area closing", colour= "closure mechanism")
dev.off()

+
  labs(x=expression ("gap area in log"(m^2)))

tiff("gap_closure_mechanism_ftype-no_clo-1721.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_ftype, timestep %in% c("17-21")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()


# --- closure by management type ---

#join with management information
stats_mtype <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_mtype.rds")
# recode years to timesteps for merges
stats_mtype <- stats_mtype %>% mutate( timestep = as.factor(recode(year,
                                                                     `2009`="9-17", 
                                                                     `2017`="17-21")))

#merge management to closure information
gap_clo_NP_91721_mtype <- merge(x = gap_clo_NP_91721, y = stats_mtype[ , c("gap_id","management", "pa_ratio", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)


# closure share per management 
tiff("gap_closure_management-no_clo_all.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721_mtype, aes(x=gap.size , y=closure_share_sum, fill=timestep)) + geom_boxplot() + facet_wrap(~management) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing") +    scale_fill_viridis(discrete = TRUE) 
dev.off()

#calculate average closure share per management
clo_per_management <- gap_clo_NP_91721_mtype %>% group_by(management, timestep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))


# closure mechanism according to management 
tiff("gap_closure_mechanism_management-no_clo.917.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_mtype, timestep %in% c("9-17")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~management) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

tiff("gap_closure_Mechanism_management-no_clo-1721.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_mtype, timestep %in% c("17-21")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~management) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

# --- closure by elevation ---

#join with management information
stats_elevation <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/stats_all_elevation.rds")
# recode years to timesteps for merges
stats_elevation <- stats_elevation %>% mutate( timestep = as.factor(recode(year,
                                                                   `2009`="9-17", 
                                                                   `2017`="17-21")))

#merge elevation to closure information
gap_clo_NP_91721_elevation <- merge(x = gap_clo_NP_91721, y = stats_elevation[ , c("gap_id","elevation", "pa_ratio", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)


# closure share per elevation 
tiff("gap_closure_elevation-no_clo_all.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721_elevation, aes(x=gap.size , y=closure_share_sum, fill=timestep)) + geom_boxplot() + facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing") +    scale_fill_viridis(discrete = TRUE) 
dev.off()

#calculate average closure share per elevation
clo_per_elevation <- gap_clo_NP_91721_elevation %>% group_by(elevation, timestep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

# closure mechanism according to elevation 
tiff("gap_closure_mechanism_elevation-no_clo.917.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_elevation, timestep %in% c("9-17")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

tiff("gap_closure_Mechanism_elevation-no_clo-1721.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo_NP_91721_elevation, timestep %in% c("17-21")), aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

#--------------------------------

######################################### different plotting options
#-------
#elevation
pm <- ggplot(gap_clo_NP_91721_elevation, aes(x=gap_area , y=closure_share_sum, col = elevation)) + geom_point() +
  theme_minimal()+ coord_flip()  + My_Theme + #theme(legend.key.size = unit(5,"point")) +
  labs(x = "gap size [m2]", y= "% of gap area closing") +  scale_color_brewer(palette="Dark2") + geom_jitter()
ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#aspect
pm <- ggplot(gap_clo_per_id_nona_aspect, aes(x=gap_area , y=closure_share_sum, col = aspect)) + geom_point() +
  theme_minimal()+ coord_flip()  + My_Theme + #theme(legend.key.size = unit(5,"point")) +
  labs(x = "gap size [m2]", y= "% of gap area closing") +  scale_color_brewer(palette="Dark2") + geom_jitter()
ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#management
pm <- ggplot(gap_clo_per_id_nona_mtype, aes(x=gap_area , y=closure_share_sum, col = management)) + geom_point() +
  theme_minimal()+ coord_flip()  + My_Theme + #theme(legend.key.size = unit(5,"point")) +
  labs(x = "gap size [m2]", y= "% of gap area closing") +  scale_color_brewer(palette="Dark2") + geom_jitter()
ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#forest type
pm <- ggplot(gap_clo_per_id_nona_ftype, aes(x=gap_area , y=closure_share_sum, col = forest_type)) + geom_point() +
  theme_minimal()+ coord_flip()  + My_Theme + #theme(legend.key.size = unit(5,"point")) +
  labs(x = "gap size [m2]", y= "% of gap area closing") +  scale_color_brewer(palette="Dark2") + geom_jitter()
ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)
#-------

######################################## analyze gap closure by Waldentwicklungsplan ########################################

#load layer WEP

wep <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/NP_data/NP_data_Bernd/Waldentwicklungsplan_epsg25832.shp")
gap_closure_mechanism1721 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism1721.tif")
gap_closure_mechanism917 <- rast("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_mechanism917.tif")

#2017-2021

#crate empty raster for rasterization
r_1 <- rast()
ext(r_1) <- ext(gap_closure_mechanism1721)
terra::res(r_1) <- terra::res(gap_closure_mechanism1721)  
terra::crs(r_1) <- terra::crs(gap_closure_mechanism1721)
# rasterize wep
wep.r <- terra::rasterize(wep, r_1, field="wep_massna")

wep_clo_mechanism <- c(wep.r, gap_closure_mechanism1721, clo_growth_1721, gaps2017)
wep_clo_mechanism.df <- as.data.frame(wep_clo_mechanism, na.rm=FALSE)
names(wep_clo_mechanism.df) <- c("wep", "closure_mechanism", "growth", "gap_id")

#calculate area share of each wep

#Wep area shares cropped to research area in QGis
# wep_area_share <- vect("i:/Fonda/workspace/berchtesgaden/focus_sites/NP_data/NP_data_Bernd/wep_massnahme_reseracharea_crop.shp")
# #rasterize croped wep
# r_1 <- rast()
# ext(r_1) <- ext(gap_closure_mechanism1721)
# terra::res(r_1) <- terra::res(gap_closure_mechanism1721)  
# terra::crs(r_1) <- terra::crs(gap_closure_mechanism1721)
# # rasterize wep
# wep.r <- terra::rasterize(wep_area_share, r_1, field="wep_massna")
# wep.df <- as.data.frame(wep.r)
# wep.df <- wep.df %>% group_by(wep_massna) %>% count()
# wep.df <- wep.df %>% mutate(n_ha = n/10000)
# wep.df$total_area <- sum(wep.df$n_ha)
# wep.df$wep_share <- round(wep.df$n_ha/wep.df$total_area,2)
# wep.df <- wep.df %>% mutate(wep_massna = as.factor(recode(wep_massna,
#                                  `1`="no regeneration measures", 
#                                  `2`="regeneration safeguarding",
#                                  `3`="planting after natural disturbance",
#                                  `4`="regeneration preparation", 
#                                  `5`="medium-term reg. preparation",
#                                  `6`="core zone")))
# wep.df<- wep.df[-c(1),]
# saveRDS(wep.df,"i:/Fonda/workspace/berchtesgaden/gaps/wep_area_share.rds" )
wep_area_share <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/wep_area_share.rds")


#exclude pixels without gap (and hence closure):
wep_clo_mechanism.df <- wep_clo_mechanism.df %>% drop_na(gap_id)


saveRDS(wep_clo_mechanism.df,"i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_wep_1721.rds" )
wep_clo_mechanism.df <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_wep_1721.rds")

#function to assign wep type, when gaps are across several area types
getWEPType <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap_id, wep) %>% #count pixels per mtype per gap
    summarize(count = n())
  #identify dominating management type in gap area
  xx <- data_frame()
  for (i in unique(x$gap_id)) {
    a <- x[x$gap_id == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one management type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several mtypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no management type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other mtype info assign that one to ID
    }
  }
  #xx<- xx %>% mutate(wep = as.factor(recode(wep,
   #                                                `2`="buffer zone",
  #                                                 `4`="core zone")))
  return(xx)
}

wep_clo_mechanism.df.clear <- getWEPType(wep_clo_mechanism.df)
names(wep_clo_mechanism.df.clear) <- c("gap_id", "wep.clear", "count")

#replace wep id with unique wep id from loop

target_df <- wep_clo_mechanism.df %>% 
  left_join(wep_clo_mechanism.df.clear,  
            by = "gap_id") %>% 
  mutate(wep.clear = ifelse(is.na(wep.clear), wep, wep.clear)) %>% 
  select(wep = wep.clear, gap_id, closure_mechanism)

# aggregate closure and gap information - prepare df for plotting
gap_clo_per_id <-  target_df %>% group_by(gap_id,wep) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

# 
# gap_clo_per_id <-  wep_clo_mechanism.df %>% group_by(gap_id,wep) %>%
#   count(closure_mechanism) %>% 
#   mutate(gap_area = sum(n))

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
         closure_mechanism = as.factor(recode(closure_mechanism,
                                              `0`="no closure", 
                                              `1`="horizontal closure",
                                              `2`="vertical closure")),
         wep = as.factor(recode(wep,
                                              `1`="no regeneration measures", 
                                              `2`="regeneration safeguarding",
                                              `3`="planting after natural disturbance",
                                              `4`="regeneration preparation", 
                                              `5`="medium-term reg. preparation",
                                              `6`="core zone")))

#delete gaps without wep information
gap_clo_per_id_nona <- gap_clo_per_id_nona %>% drop_na(wep)

##### if I want to include the no closure pixels, I have to disable following line: !!!!!!!!!!!!!!
#exclude no closure shares
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "horizontal closure"))

gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

#aggregate closure share > df without dropping no closure areas
# gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
#   mutate(closure_share_sum = sum(closure_share))

#gap_clo_per_id_nona$gap_area <- as.numeric(gap_clo_per_id_nona$gap_area)

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona<-gap_clo_per_id_nona %>% 
  mutate(gap_area_bins = (cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))


# prepare plotting
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 18),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=16),
  strip.text.x = element_text(size = 16),
  legend.position="bottom")

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes/"
setwd(wd)


#calculate average closure share per wep
clo_per_wep <- gap_clo_per_id_nona %>% group_by(wep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

# closure share per management type 
tiff("gap_closure_wep-no_clo_all.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_per_id_nona, aes(x=gap.size , y=closure_share_sum)) + geom_boxplot() + facet_wrap(~wep) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing")
dev.off()

#calculate average closure share per management type
clo_per_wep<- gap_clo_per_id_nona %>% group_by(wep) %>%
  summarise(avg_closure_share = mean(closure_share_sum))

tiff("gap_closure_wep-no_clo.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_per_id_nona, aes(x=gap.size , y=closure_share, fill=closure_mechanism)) + geom_boxplot() + facet_wrap(~wep) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2") + My_Theme +labs(x = "gap size [m2]", y= "% of gap area closing", colour= "closure mechanism")
dev.off()

pm <- ggplot(gap_clo_per_id_nona, aes(x=gap_area , y=closure_share_sum, col = wep)) + geom_point() +
  theme_minimal()+ coord_flip()  + My_Theme + #theme(legend.key.size = unit(5,"point")) +
  labs(x = "gap size [m2]", y= "% of gap area closing") +  scale_color_brewer(palette="Dark2") + geom_jitter()

tiff("gap_closure_wep-no_clo_marginal.tiff", units="in", width=12, height=8, res=300)
ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)
dev.off()

# ------ 2009-2017 ------ 

#crate empty raster for rasterization
r_1 <- rast()
ext(r_1) <- ext(gap_closure_mechanism917)
terra::res(r_1) <- terra::res(gap_closure_mechanism917)  
terra::crs(r_1) <- terra::crs(gap_closure_mechanism917)
# rasterize wep
wep.r <- terra::rasterize(wep, r_1, field="wep_massna")

wep_clo_mechanism <- c(wep.r, gap_closure_mechanism917, gaps2009)
wep_clo_mechanism.df <- as.data.frame(wep_clo_mechanism, na.rm=FALSE)
names(wep_clo_mechanism.df) <- c("wep", "closure_mechanism", "gap_id")

#exclude pixels without gap (and hence closure):
wep_clo_mechanism.df <- wep_clo_mechanism.df %>% drop_na(gap_id)


saveRDS(wep_clo_mechanism.df,"i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_wep_917.rds" )
wep_clo_mechanism.df <- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/gap_closure_wep_917.rds")

#function to assign wep type, when gaps are across several area types
getWEPType <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap_id, wep) %>% #count pixels per mtype per gap
    summarize(count = n())
  #identify dominating management type in gap area
  xx <- data_frame()
  for (i in unique(x$gap_id)) {
    a <- x[x$gap_id == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one management type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several mtypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no management type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other mtype info assign that one to ID
    }
  }
  #xx<- xx %>% mutate(wep = as.factor(recode(wep,
  #                                                `2`="buffer zone",
  #                                                 `4`="core zone")))
  return(xx)
}

wep_clo_mechanism.df.clear <- getWEPType(wep_clo_mechanism.df)
names(wep_clo_mechanism.df.clear) <- c("gap_id", "wep.clear", "count")

#replace wep id with unique wep id from loop

target_df <- wep_clo_mechanism.df %>% 
  left_join(wep_clo_mechanism.df.clear,  
            by = "gap_id") %>% 
  mutate(wep.clear = ifelse(is.na(wep.clear), wep, wep.clear)) %>% 
  select(wep = wep.clear, gap_id, closure_mechanism)


gap_clo_per_id <-  target_df %>% group_by(gap_id,wep) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

# aggregate closure and gap information - prepare df for plotting
# gap_clo_per_id <-  wep_clo_mechanism.df %>% group_by(gap_id,wep) %>%
#   count(closure_mechanism) %>% 
#   mutate(gap_area = sum(n))

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
         closure_mechanism = as.factor(recode(closure_mechanism,
                                              `0`="no closure", 
                                              `1`="horizontal closure",
                                              `2`="vertical closure")),
         wep = as.factor(recode(wep,
                                `1`="no regeneration measures", 
                                `2`="regeneration safeguarding",
                                `3`="planting after natural disturbance",
                                `4`="regeneration preparation", 
                                `5`="medium-term reg. preparation",
                                `6`="core zone")))

#delete gaps without wep information
gap_clo_per_id_nona <- gap_clo_per_id_nona %>% drop_na(wep)

##### if I want to include the no closure pixels, I have to disable following line: !!!!!!!!!!!!!!
#exclude no closure shares
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "horizontal closure"))

gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share))

#aggregate closure share > df without dropping no closure areas
# gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
#   mutate(closure_share_sum = sum(closure_share))

#gap_clo_per_id_nona$gap_area <- as.numeric(gap_clo_per_id_nona$gap_area)

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona<-gap_clo_per_id_nona %>% 
  mutate(gap_area_bins = (cut(gap_area, breaks = c(399,500,1000,5000,10000,50000,450000))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(399,500]`="400-500",
                                     `(500,1e+03]`="500-1000",
                                     `(1e+03,5e+03]`="1000-5000",
                                     `(5e+03,1e+04]`="5000-10,000",
                                     `(1e+04,5e+04]`="10,000-50,000",
                                     `(5e+04,4.5e+05]`="50,000-450,000")))


# prepare plotting
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 18),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=16),
  strip.text.x = element_text(size = 16),
  legend.position="bottom")

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/plots_gap_changes/"
setwd(wd)


#calculate average closure share per wep
clo_per_wep <- gap_clo_per_id_nona %>% group_by(wep, closure_mechanism) %>%
  summarise(avg_closure_share = round(mean(closure_share_sum),2),
            share_mechanism = round(mean(closure_share),2))

# closure share per management type 
tiff("gap_closure_wep-no_clo_917.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_per_id_nona, aes(x=gap.size , y=closure_share_sum)) + geom_boxplot() + facet_wrap(~wep) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [m2]", y= "% of gap area closing")
dev.off()