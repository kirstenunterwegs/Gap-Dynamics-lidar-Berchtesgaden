###############################################
# Identify closure mechanism
##############################################

library(dplyr)
library(tidyr)
library(terra)
library(lidR)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)

# library(sp)
#library(lidR)
# library(ggExtra)
#library(sf)
# library(ForestTools)
# library(ForestGapR)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load layers ----

gaps2009 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

grid <- vect("F:/temp/cornelius_to_kirsten/grid_1km_bdg.gpkg")

###################################################------create forest edge mask for lateral growth classification

# need to loop through a subset of gaps as R crashes when whole gap raster is processed

# e <- ext(797000, 798000 , 5270000 ,5272000 )
# gaps9.sub <- crop(gaps2009, e)

# --- 2009 ---

#create loop where gap layer is cropped to grid cell and processed
grid_id <- grid$id
#i=169
for (i in grid_id) {
  print(i)
  g.sub <- grid[grid$id==i]
  gaps <- gaps2009
  gaps.sub <- crop(gaps2009,g.sub)
  print("finished crop")
  bo <- boundaries(gaps.sub, directions=8, inner=FALSE)
  print("boundaries created")
  writeRaster(bo, paste0("/mnt/public/temp/cornelius_to_kirsten/gap_boundaries9/gap_boundaries_9.", i, ".tif"), overwrite=T)
}

#create SpatRasterCollection from boundary tiles
# bo2009.list <- list.files("F:/temp/cornelius_to_kirsten/gap_boundaries9", full.names = TRUE)
# rsrc <- lapply(bo2009.list,rast) #terra::src
# bo2009.spatcol <- sprc(rsrc) #does not work
# boundaries.2019<- mosaic(bo2009.spatcol)

####### create gap boundary layer for whole NP
bo2009.list <- list.files("F:/temp/cornelius_to_kirsten/gap_boundaries9", full.names = TRUE)
r <- lapply(bo2009.list, raster)

r$overwrite <- TRUE
boundaries.2009 <- do.call(merge, r) # merge all tiles into one rasterlayer

writeRaster(boundaries.2009, "data/gap_boundaries/gap_boundaries_2009.tif")
#reprojecting the boundary layer in QGis to crs 25832


# --- 2017 ---

for (i in grid_id) {
  print(i)
  g.sub <- grid[grid$id==i]
  gaps <- gaps2017
  gaps.sub <- crop(gaps2017,g.sub)
  print("finished crop")
  bo <- boundaries(gaps.sub, directions=8, inner=FALSE)
  print("boundaries created")
  writeRaster(bo, paste0("/mnt/public/temp/cornelius_to_kirsten/gap_boundaries17/gap_boundaries_17.", i, ".tif"), overwrite=T)
}

#create SpatRasterCollection from boundary tiles
# bo2009.list <- list.files("F:/temp/cornelius_to_kirsten/gap_boundaries9", full.names = TRUE)
# rsrc <- lapply(bo2009.list,rast) #terra::src
# bo2009.spatcol <- sprc(rsrc) #does not work
# boundaries.2019<- mosaic(bo2009.spatcol)

####### create gap boundary layer for whole NP
bo2017.list <- list.files("F:/temp/cornelius_to_kirsten/gap_boundaries17", full.names = TRUE)
r <- lapply(bo2017.list, raster)

r$overwrite <- TRUE
boundaries.2017 <- do.call(merge, r) # merge all tiles into one rasterlayer

writeRaster(boundaries.2017, "data/data/gap_boundaries/gap_boundaries_2017.tif")
#reprojecting the boundary layer in QGis to crs 25832



####################################################### classify vertical and horizontal closure ##################################

#load closure areas with growth information
clo_growth_917 <- rast("data/data/closure_area_growth_917.tif")
clo_growth_1721 <- rast("data/data/closure_area_growth_1721.tif")

#load gap boundaries
boundaries.2009 <- rast("data/data/gap_boundaries/gap_boundaries_2009_25832.tif")
boundaries.2017 <- rast("data/data/gap_boundaries/gap_boundaries_2017_25832.tif")

#adjust extents
clo_growth_917 <- crop(clo_growth_917, boundaries.2009)
clo_growth_1721 <- crop(clo_growth_1721, boundaries.2017)

# e <- ext(797000, 798000 , 5270000 ,5272000 )
# clo_growth_917.sub <- crop(clo_growth_917, e)
# boundaries.2009.sub <- crop(boundaries.2009, e)

# need to differ between both timesteps due to differnet growing periods

diff_closure_layer <- clo_growth_917
boundary_layer <- boundaries.2009

gap_closure_mechanism917 <- function(diff_closure_layer, boundary_layer){       # 0.5 m * 8 = 4m height gain

  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 4  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 4  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 4 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 4 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism1721 <- function(diff_closure_layer, boundary_layer){      # 0.5 m * 4 = 2m height gain

  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 2  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 2  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 2 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 2 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism917 <- gap_closure_mechanism917(clo_growth_917, boundaries.2009)
gap_closure_mechanism1721 <- gap_closure_mechanism1721(clo_growth_1721, boundaries.2017)


terra::writeRaster(gap_closure_mechanism917, "data/data/gap_closure_mechanism917_new.tif", overwrite=TRUE) 
terra::writeRaster(gap_closure_mechanism1721, "data/data/gap_closure_mechanism1721_new.tif", overwrite=TRUE)


###----------- analyze closure per gap (size), elevation, aspect, management and forest type ----------- ###


####################################################################################################################
# prepare closure mechanism dfs per timestep
###################################################################################################################

#prepare layers to limit analysis to core zone, below 1800m and with forest type information

# --- load NP information 

foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")
management <- vect("F:/Projects/CanopyDynamicsBDG/data/NP_data/npb_zonierung_22_epsg25832.shp")
aspect<-  rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")

# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))


# 2009 - 2017
###################################################################################################################

#merge closure mechanism with gaps
# gap_closure_mechanism917 <- rast( "processed/closure/gap_closure_mechanism917_new.tif")
# gaps2009 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
# 
# gaps2009 <- crop(gaps2009, gap_closure_mechanism917)
# gap_closure_mechanism_stack <- c(gap_closure_mechanism917, gaps2009)
# 
# #crop to same extent
# elevation.below1800 <- crop(elevation.below1800,gap_closure_mechanism917 )
# foresttype <- crop(foresttype, gap_closure_mechanism917 )
# 
# #mask down to reserach area
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, elevation.below1800)
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, foresttype)
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, core.zone)
# 
# gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
# #exclude pixels without gap (and hence closure):
# gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
# names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")
# 
# saveRDS(gap_closure_mechanism_stack.df,"processed/closure/updated/gap_closure_mechanism_pergap_917.rds" )

gap_closure_mechanism_stack.df <- readRDS("processed/closure/updated/gap_closure_mechanism_pergap_917.rds")

# aggregate closure and gap information 
gap_clo_per_id <-  gap_closure_mechanism_stack.df %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 5 gaps do not experience any closure from 2009-2017 (previously 7, maybe the masking cut off some gap areas)

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


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

#gap_clo_per_id_nona$gap_area <- as.numeric(gap_clo_per_id_nona$gap_area)

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_917<-gap_clo_per_id_nona %>% 
  mutate(gap_area.ha = gap_area/10000,
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.039,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.039,0.1]`="0.04-0.1",
                                     `(0.1,0.2]`="0.1-0.2",
                                     `(0.2,0.3]`="0.2-0.3",
                                     `(0.3,0.4]`="0.3-0.4",
                                     `(0.4,0.5]`="0.4-0.5",
                                     `(0.5,0.6]`="0.5-0.6",
                                     `(0.6,0.7]`="0.6-0.7",
                                     `(0.7,0.8]`="0.7-0.8",
                                     `(0.8,0.9]`="0.8-0.9",
                                     `(0.9,1]`="0.9-1",
                                     `(1,45]`=">1")))

#rough gap size bin classification

# gap_area_bins = (cut(gap_area.ha, breaks = c(0.039,0.05,0.1,0.5,1,45))))%>% 
#   mutate(gap.size = as.factor(recode(gap_area_bins,
#                                      `(0.039,0.05]`="0.04-0.05",
#                                      `(0.05,0.1]`="0.05-0.1",
#                                      `(0.1,0.5]`="0.1-0.5",
#                                      `(0.5,1]`="0.5-1",
#                                      `(1,45]`=">1")))

# mutate(gap.size = as.factor(recode(gap_area_bins,
#                                    `(399,500]`="400-500",
#                                    `(500,1e+03]`="500-1000",
#                                    `(1e+03,5e+03]`="1000-5000",
#                                    `(5e+03,1e+04]`="5000-10,000",
#                                    `(1e+04,5e+04]`="10,000-50,000",
#                                    `(5e+04,4.5e+05]`="50,000-450,000")))

# 2017- 2021
###################################################################################################################

# #merge closure mechanism with gaps
# gap_closure_mechanism1721 <- rast("processed/closure/gap_closure_mechanism1721_new.tif")
# gaps2017 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
# 
# gaps2017 <- crop(gaps2017, gap_closure_mechanism1721)
# gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721, gaps2017)
# 
# #mask down to reserach area
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, elevation.below1800)
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, foresttype)
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, core.zone)
# 
# gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)
# #exclude pixels without gap (and hence closure):
# gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
# names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")
# 
# saveRDS(gap_closure_mechanism_stack.df_1721,"processed/closure/updated/gap_closure_mechanism_pergap_1721.rds" )

gap_closure_mechanism_stack.df_1721 <- readRDS("processed/closure/updated/gap_closure_mechanism_pergap_1721.rds" )

# aggregate closure and gap information - prepare df for plotting
gap_clo_per_id <-  gap_closure_mechanism_stack.df_1721 %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify number of gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 106 gaps do not experience any closure from 2009-2017 (without area reduction 162 gaps)

gaps_no_closing <- subset(gap_clo_per_id, contraction %in% 1)
range(gaps_no_closing$gap_area) # range 401-5678


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

#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))



# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_1721<-gap_clo_per_id_nona %>% 
  mutate(gap_area.ha = gap_area/10000,
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.039,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.039,0.1]`="0.04-0.1",
                                     `(0.1,0.2]`="0.1-0.2",
                                     `(0.2,0.3]`="0.2-0.3",
                                     `(0.3,0.4]`="0.3-0.4",
                                     `(0.4,0.5]`="0.4-0.5",
                                     `(0.5,0.6]`="0.5-0.6",
                                     `(0.6,0.7]`="0.6-0.7",
                                     `(0.7,0.8]`="0.7-0.8",
                                     `(0.8,0.9]`="0.8-0.9",
                                     `(0.9,1]`="0.9-1",
                                     `(1,45]`=">1")))

# mutate(gap.size = as.factor(recode(gap_area_bins,
#                                    `(399,500]`="400-500",
#                                    `(500,1e+03]`="500-1000",
#                                    `(1e+03,5e+03]`="1000-5000",
#                                    `(5e+03,1e+04]`="5000-10,000",
#                                    `(1e+04,5e+04]`="10,000-50,000",
#                                    `(5e+04,4.5e+05]`="50,000-450,000")))


### -------------------------------merge 9-17 and 17-21 dfs for comparison ----------------------- ###

gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
#rearrange timestep labels
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))

# ---add environmental feature information

stats_aspect <- readRDS("processed/gap_features/stats_all_aspect.rds")
# recode years to timesteps for merges
stats_aspect <- stats_aspect %>% mutate( timestep = as.factor(recode(year,
                                                                     `2009`="9-17", 
                                                                     `2017`="17-21")))

#join with forest information
stats_ftype <- readRDS("processed/gap_features/stats_all_ftype.rds")
# recode years to timesteps for merges
stats_ftype <- stats_ftype %>% mutate( timestep = as.factor(recode(year,
                                                                   `2009`="9-17", 
                                                                   `2017`="17-21")))

#join with management information
stats_elevation <- readRDS("processed/gap_features/stats_all_elevation.rds")
# recode years to timesteps for merges
stats_elevation <- stats_elevation %>% mutate( timestep = as.factor(recode(year,
                                                                           `2009`="9-17", 
                                                                           `2017`="17-21")))


gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_aspect[ , c("gap_id","aspect", "pa_ratio", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)
gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_ftype[ , c("gap_id","forest_type", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)
gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_elevation[ , c("gap_id","elevation", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)

# exclude gaps in elevation band above 1800 m
gap_clo_NP_91721 <- gap_clo_NP_91721[gap_clo_NP_91721$elevation != "1800-2000",]

# calculate annual closure rates

gap_clo_NP_91721<- gap_clo_NP_91721 %>% mutate(time = as.numeric(recode(timestep,
                                          `9-17`= 8, 
                                          `17-21`= 4)),
                                          clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                          clo_share_annual = round(closure_share/time,4)*100,
                                          clo_area_sum_annual = round(closure_area_sum/time,4),
                                          clo_area_annual = round(closure_area/time,4))


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
  strip.text.y = element_text(size = 16),
  legend.position="bottom")

require(scales)
library(viridis)

wd <- "F:/Projects/CanopyDynamicsBDG/results/gap_closure/"
setwd(wd)


# ----- closure across NP ------------
library(Hmisc)


clo_NP <- gap_clo_NP_91721 %>% group_by(closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),2),
            share_mechanism_annual = round(weighted.mean(clo_share_annual, time),2),
            percent_share_mechanism = share_mechanism_annual/ avg_clo_share_annual,
            avg_clo_area_annual = round(weighted.mean(clo_area_sum_annual, time),2),
            area_mechanism_annual = round(weighted.mean(clo_area_annual, time),2))

clo_NP.sizebins <- gap_clo_NP_91721 %>% group_by(closure_mechanism, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual,time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            avg_clo_area_annual = round(weighted.mean(clo_area_sum_annual, time),2),
            area_mechanism_annual = round(weighted.mean(clo_area_annual, time),2),
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            n.obs = n())

#plot closure according to gap size bins
tiff("gap_closure_NP_box_violin.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_violin() + 
  geom_boxplot(width=0.1) +
  theme_minimal()+ coord_flip() +
  theme(legend.position="bottom") +My_Theme  +    scale_fill_viridis(discrete = TRUE) +
  labs( x = "gap size [ha]", y= "% of gap area closing annually")
dev.off()

tiff("gap_closure_NP_boxplot.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_boxplot() +
  theme_minimal()+ coord_flip() +
  theme(legend.position="bottom") +My_Theme  +    scale_fill_viridis(discrete = TRUE) +
  labs(x = "gap size [ha]", y= "% of gap area closing annually")
dev.off()

tiff("gap_closure-mechanism_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) + geom_boxplot() +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2") + My_Theme +
  labs( x = "gap size [ha]", y= "% of gap area closing annually", fill= "closure mechanism") 
dev.off()


# plot gap size vs. closure share per mechanism
# 
# pm <- ggplot(gap_clo_NP_91721, aes(x=gap_area, y=clo_share_annual, col=closure_mechanism)) + geom_point()+
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   theme_classic()+  scale_color_brewer(palette="Dark2") + My_Theme +
#   labs(x = "gap size [m2]", y= "% of gap area closing annually", colour= "closure mechanism")
# tiff("gap_closure-no_clo_mechanism.tiff", units="in", width=12, height=8, res=300)
# ggMarginal(pm, type="boxplot", groupColour = TRUE, groupFill = TRUE)
# dev.off()
# 
# 
# # complete closure share vs. gap size
# 
# p<- ggplot(gap_clo_NP_91721, aes(x=gap_area, y=clo_share_sum_annual)) + geom_point() +
#                 theme_classic() +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) + My_Theme+
#   labs(x = "gap size [ha]", y= "% of gap area closing annually")
# 
# 
# tiff("gap_closure-no_clo_all_917.tiff", units="in", width=12, height=8, res=300)
# ggMarginal(p, type="histogram")
# #ggMarginal(p, type="boxplot")
# dev.off()
# 
# # complete closure share vs. timestep
# 
# p<- ggplot(gap_clo_NP_91721, aes(x=gap_area, y=clo_share_sum_annual, col=timestep)) + geom_point() +
#   theme_classic() +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) + My_Theme+
#   labs(x = "gap size [ha]", y= "% of gap area closing annually")
# 
# 
# tiff("gap_closure-timestep.tiff", units="in", width=12, height=8, res=300)
# ggMarginal(p,type="boxplot", groupColour = TRUE, groupFill = TRUE)
# #ggMarginal(p, type="boxplot")
# dev.off()

ggplot(gap_clo_NP_91721  %>% filter(gap_area.ha <0.5) 
       , aes(x=gap_area.ha, y=clo_area_sum_annual)) + geom_point()+
  geom_smooth(colour="blue") +
geom_smooth(colour="red", method = "lm") + My_Theme

ggplot(gap_clo_NP_91721  %>% filter(gap_area.ha <0.5) 
       , aes(x=gap_area.ha, y=clo_share_sum_annual)) + geom_point()+
  geom_smooth(colour="blue") +
  geom_smooth(colour="red", method = "lm") +
  facet_wrap(~timestep) + My_Theme

ggplot(gap_clo_NP_91721  %>% filter(gap_area.ha <0.5) 
       , aes(x=forest_type, y=clo_share_sum_annual)) + geom_point()+
  geom_smooth(colour="blue") +
  geom_smooth(colour="red", method = "lm") + My_Theme # +
  facet_wrap(~elevation)


# mean closure shares

tiff("gap_closure_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + #avg_clo_area_annual
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

# average closure area annual
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_area_annual)) + #avg_clo_area_annual
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average area [m2] per gap closing annually", colour= "forest type")

# point with error bar
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all))

# --- mechanism ---


tiff("gap_closure-mechanism_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
 # facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "closure mechanism")
dev.off()

# point with error bar
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3))


# --- closure by aspect---

#calculate average closure share per aspect
clo_per_aspect <- gap_clo_NP_91721%>% group_by(aspect, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual,time),2),
            share_mechanism_annual = round(weighted.mean(clo_share_annual, time),2),
            percent_share_mechanism = share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            n.obs = n())

clo_per_aspect.sizebins <- gap_clo_NP_91721%>% group_by(aspect, closure_mechanism, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual,time),2),
            share_mechanism_annual = round(weighted.mean(clo_share_annual, time),2),
            percent_share_mechanism = share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            n.obs = n())


# closure share per aspect 
tiff("gap_closure_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_boxplot() + 
   facet_wrap(~aspect) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
dev.off()

# closure mechanism according to aspect 
tiff("gap_closure_mechanism_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) + geom_boxplot() + facet_wrap(~aspect) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
dev.off()

# # - mean sizebin
# 
# # closure share per aspect 
# tiff("gap_closure_aspect_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_aspect.sizebins, aes(x=gap.size , y=avg_clo_share_annual )) + geom_bar(stat = "identity", color="black", position=position_dodge()) + 
#  # facet_wrap(~aspect) +
#   facet_grid(aspect~.) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "gap size [ha]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
# dev.off()
# 
# 
# # closure mechanism according to aspect 
# tiff("gap_closure_mechanism_aspect_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_aspect.sizebins, aes(x=gap.size , y=share_mechanism_annual , fill=closure_mechanism))+ geom_bar(stat = "identity", color="black", position=position_dodge()) + 
#   #facet_wrap(~aspect) +
#   facet_grid(aspect~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
# dev.off()

#-- mean

tiff("gap_closure_aspect_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_aspect, aes(x=aspect , y=avg_clo_share_annual )) + geom_bar(stat = "identity") + 
  # facet_wrap(~aspect) +
 # facet_grid(aspect~.) +
  theme_minimal()  + My_Theme + coord_flip() +
  labs(x = "aspect", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
dev.off()

tiff("gap_closure_mechanism_aspect_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_aspect, aes(x=aspect , y=share_mechanism_annual , fill=closure_mechanism )) + geom_bar(stat = "identity", position="dodge") + 
  # facet_wrap(~aspect) +
  # facet_grid(aspect~.) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "aspect", y= "% of gap area closing annually") +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") 
dev.off()


# --- closure by forest type ---

#calculate average closure share per forest type
clo_per_ftype <- gap_clo_NP_91721%>% group_by(forest_type, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual)

clo_per_ftype.sizebins.mechanism <- gap_clo_NP_91721 %>% group_by(forest_type, closure_mechanism, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            n.obs = n())

clo_per_ftype.sizebins <- gap_clo_NP_91721 %>% group_by(forest_type, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sqrt(wtd.var(clo_share_sum_annual, time)),
            sd_mechanism= sqrt(wtd.var(clo_share_annual, time)),
            n.obs = n())


# closure share per forest type 

tiff("gap_closure_ftype_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_boxplot() + 
   facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually") +   scale_fill_viridis(discrete = TRUE) 
dev.off()


# closure mechanism according to forest type 

tiff("gap_closure_mechanism_ftype_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) + 
  geom_boxplot() + facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
dev.off()

# --- mean closure rates per forest type
tiff("gap_closure_ftype_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_clo_share_annual)) + 
  geom_bar(stat = "identity") +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x="forest type", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

# point with error bar
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all))

tiff("gap_closure_mechanism_ftype_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x="forest type", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

# point with error bar
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3))


# --- mean closure rates per forest type and gap size ---

# closure mechanism mean according to forest type and gap size
# ggplot(clo_per_ftype.sizebins, aes(x=gap.size , y=avg_clo_share_annual, fill=forest_type)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   #facet_wrap(~forest_type) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")

tiff("gap_closure_ftype_size_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  facet_grid(forest_type~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+
  geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
dev.off()

tiff("gap_closure_mechanism_ftype_size_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype.sizebins.mechanism, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  facet_grid(forest_type~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()


# --- closure by elevation ---

ggplot(gap_clo_NP_91721, aes(x=gap_area,fill=factor(elevation))) + geom_density(alpha=0.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))
# gap size is similar distributed across elevation
ggplot(gap_clo_NP_91721, aes(x=gap_area,fill=factor(elevation))) + geom_histogram(bins=100,alpha=0.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))

#calculate average closure share per elevation
clo_per_elevation<- gap_clo_NP_91721 %>% group_by(elevation) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

clo_per_elevation.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

# clo_per_elevation <- gap_clo_NP_91721 %>% group_by(elevation) %>%
#   summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
#             avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             n.obs = n())
#------

clo_per_elevation.sizebins.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, closure_mechanism, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

clo_per_elevation.sizebins <- gap_clo_NP_91721 %>% group_by(elevation, gap.size) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

gap_clo_NP_91721$elevation <- ordered(gap_clo_NP_91721$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
gap_clo_NP_91721$elevation <- factor(gap_clo_NP_91721$elevation,levels=rev(levels(gap_clo_NP_91721$elevation)))

# closure share per elevation 

# tiff("gap_closure_elevation.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) +geom_boxplot() + 
#   # facet_wrap(~elevation) +
#   facet_grid(elevation~.) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "gap size [ha]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
# dev.off()

tiff("gap_closure_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=elevation ,y=clo_share_sum_annual)) +geom_boxplot() + 
  # facet_wrap(~elevation) +
 # facet_grid(elevation~.) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "elevation [m]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
dev.off()


# closure mechanism per elevation 

tiff("gap_closure_mechanism_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=elevation , y=clo_share_annual, fill=closure_mechanism)) + 
  geom_boxplot() +# facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "% of gap area closing annually", colour= "closure mechanism")
dev.off()

# --- mean closure rates per elevation and gap size ---

#order elevation labels
clo_per_elevation.sizebins.mechanism$elevation <- ordered(clo_per_elevation.sizebins.mechanism$elevation,
                                      levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
clo_per_elevation.sizebins.mechanism$elevation <- factor(clo_per_elevation.sizebins.mechanism$elevation,levels=rev(levels(clo_per_elevation.sizebins.mechanism$elevation)))

tiff("gap_closure_elevation_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation, aes(x=elevation , y=avg_clo_share_annual )) + 
  geom_bar(stat = "identity") +
 # facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

tiff("gap_closure_mechanism_elevation_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation.mechanism, aes(x=elevation , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  #facet_wrap(~elevation) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

#reversing elevation labels
clo_per_elevation.sizebins$elevation <- ordered(clo_per_elevation.sizebins$elevation,
                                      levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
clo_per_elevation.sizebins.mechanism$elevation <- ordered(clo_per_elevation.sizebins.mechanism$elevation,
                                                levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))

tiff("gap_closure_elevation_avg_sizebins.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
  geom_bar(stat = "identity") +
  facet_grid(elevation~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
dev.off()

tiff("gap_closure_mechanism_elevation_avg_sizebins.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation.sizebins.mechanism, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  facet_grid(elevation~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()

# ---- closure per elevation and forest type 

clo_per_ftype.elev <- gap_clo_NP_91721 %>% group_by(elevation, forest_type) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

clo_per_ftype.elev.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, forest_type,closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(weighted.mean(clo_share_sum_annual, time),4),
            avg_share_mechanism_annual = round(weighted.mean(clo_share_annual, time),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            n.obs = n())

#clo_per_ftype.elev$elev.ftype <- paste0(clo_per_ftype.elev$elevation ,"-" ,clo_per_ftype.elev$forest_type)

tiff("gap_closure_elev_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype.elev, aes(x=elevation , y=avg_clo_share_annual)) + 
  geom_bar(stat = "identity") +
  facet_grid(forest_type~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
dev.off()

tiff("gap_closure_mechanism_elev_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype.elev.mechanism, aes(x=elevation , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
  geom_bar(stat = "identity", position="dodge") +
  facet_grid(forest_type~.) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
dev.off()


# --- closure by shape complexity ---

# closure share vs. complexity 

tiff("gap_closure_complexity.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_sum_annual)) + geom_point(size = 5, alpha = 0.5) +
  theme_minimal()+ coord_flip()  + My_Theme +
  labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually")  +
  geom_smooth(method=lm) 
dev.off()

# closure mechanism vs. complexity vs. closure share

tiff("gap_closure_mechanism_complexity.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_annual, color=closure_mechanism)) + 
  geom_point(size = 5, alpha = 0.5) +  scale_color_brewer(palette="Dark2", name = "closure mechanism")+
  theme_minimal()+ coord_flip()   + My_Theme +
  labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually", colour= "closure mechanism") +
  geom_smooth(method=lm) 
dev.off()


# closure mechanism vs. complexity vs. time step

tiff("gap_closure_mechanism_complexity.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_annual, color=timestep)) + 
  geom_point(size = 5, alpha = 0.5) +  scale_color_brewer(palette="Dark2", name = "closure mechanism")+
  theme_minimal()+ coord_flip()   + My_Theme +
  labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually", colour= "closure mechanism") 
dev.off()

# p:a vs. aspect
tiff("pa_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=aspect)) + 
  geom_boxplot()+ coord_flip() 
dev.off()

# p:a vs. elevation
tiff("pa_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=elevation)) + 
  geom_boxplot()+ coord_flip() 
dev.off()

# p:a vs. forest type
tiff("pa_forest.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=forest_type)) + 
  geom_boxplot()+ coord_flip() 
dev.off()

#-------------------------------- model check of influence of gap complexity on closure share

gap.closure.lm<-lm(clo_share_annual ~ pa_ratio + aspect  + forest_type, data = gap_clo_NP_91721)

summary(gap.closure.lm)

par(mfrow=c(2,2))
plot(gap.closure.lm)
par(mfrow=c(1,1))

    