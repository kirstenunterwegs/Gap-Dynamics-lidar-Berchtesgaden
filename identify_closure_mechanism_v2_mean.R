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

#saveRDS(gap_clo_NP_91721, "processed/closure/updated/gap_closure_elevation.rds")


# calculate annual closure rates

gap_clo_NP_91721<- gap_clo_NP_91721 %>% mutate(time = as.numeric(recode(timestep,
                                          `9-17`= 8, 
                                          `17-21`= 4)),
                                          clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                          clo_share_annual = round(closure_share/time,4)*100,
                                          clo_area_sum_annual = round(closure_area_sum/time,4),
                                          clo_area_annual = round(closure_area/time,4))

# ---- relabel Larch-Swiss stone pine and Dwarf mountain pine to one forest type class ---

gap_clo_NP_91721 <- gap_clo_NP_91721 %>% mutate(forest_type = as.factor(recode(forest_type,
                                                             `Larch-Swiss stone pine`= "Larch-Pine",
                                                             `Dwarf mountain pine`= "Larch-Pine")))

#label ordering
gap_clo_NP_91721$forest_type <- ordered(gap_clo_NP_91721$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-Pine"))
gap_clo_NP_91721$forest_type <- factor(gap_clo_NP_91721$forest_type,levels=rev(levels(gap_clo_NP_91721$forest_type)))



# --- append lateral + vertical closure info to main df for distribution display ----
gap_clo_NP_91721$id <- as.numeric(paste0(gap_clo_NP_91721$gap_id,gap_clo_NP_91721$time))

gap_clo <- as.data.frame(gap_clo_NP_91721[,c("id", "closure_mechanism",  "gap.size","aspect", "forest_type", "elevation", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)
gap_clo$aspect <- as.character(gap_clo$aspect)
gap_clo$forest_type <- as.character(gap_clo$forest_type)
gap_clo$elevation <- as.character(gap_clo$elevation)

id <- as.character(unique(gap_clo$id))
#i=74
for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  aspect <- unique(sub$aspect)
  ftype <- unique(sub$forest_type)
  elev <- unique(sub$elevation)
  k <- c(i,"lateral + vertical", size, aspect, ftype, elev, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)
gap_clo$aspect <- as.factor(gap_clo$aspect)
gap_clo$forest_type <- as.factor(gap_clo$forest_type)
gap_clo$elevation <- as.factor(gap_clo$elevation)

#order labels
gap_clo$forest_type <- ordered(gap_clo$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-Pine"))
gap_clo$forest_type <- factor(gap_clo$forest_type,levels=rev(levels(gap_clo$forest_type)))

gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

saveRDS(gap_clo, "processed/closure/updated/clo_analysis_ready.rds")
gap_clo <- readRDS("processed/closure/updated/clo_analysis_ready.rds")

# -----

# prepare plotting
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 28),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=18),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16),
  legend.position="top") #bottom

require(scales)
library(viridis)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/results/gap_closure"
setwd(wd)


# ----- closure across NP ------------
#library(Hmisc)


clo_NP <- gap_clo_NP_91721 %>% group_by(closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),2),
            share_mechanism_annual = round(mean(clo_share_annual),2),
            percent_share_mechanism = share_mechanism_annual/ avg_clo_share_annual,
            #avg_clo_area_annual = round(mean(clo_area_sum_annual),2),
            #area_mechanism_annual = round(mean(clo_area_annual),2),
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual))

# clo_NP.agg <- as.data.frame(clo_NP)
# clo_NP.agg <- clo_NP.agg[,c("closure_mechanism", "share_mechanism_annual","sd_mechanism") ]
# 
# #for(i in aspct) {
# # sub <- subset(aspect_gap, aspect %in% i)
#   clo_NP.agg$closure_mechanism <- as.character(clo_NP.agg$closure_mechanism)
#   k <- c( "lateral + vertical", sum(clo_NP.agg$share_mechanism_annual), sum(clo_NP.agg$sd_mechanism))
#   clo_NP.agg <- rbind(clo_NP.agg, k)
#   clo_NP.agg$share_mechanism_annual <- as.numeric(clo_NP.agg$share_mechanism_annual)
#   clo_NP.agg$sd_mechanism <- as.numeric(clo_NP.agg$sd_mechanism)
# #}
#   clo_NP.agg$closure_mechanism <- as.factor(clo_NP.agg$closure_mechanism)


clo_NP.sizebins <- gap_clo_NP_91721 %>% group_by(closure_mechanism, gap.size) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
            avg_share_mechanism_annual = round(mean(clo_share_annual),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            avg_clo_area_annual = round(mean(clo_area_sum_annual),2),
            area_mechanism_annual = round(mean(clo_area_annual),2),
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            sd_all.a = sd(clo_area_sum_annual),
            n.obs = n())

# aggregate closure mechanism information for plotting
clo.size <- as.data.frame(clo_NP.sizebins)
clo.size <- clo.size[,c("closure_mechanism", "gap.size",  "avg_share_mechanism_annual","sd_mechanism") ]

size <- as.character(unique(clo.size$gap.size))
clo.size$closure_mechanism <- as.character(clo.size$closure_mechanism)

for(i in size) {
  sub <- subset(clo.size, gap.size %in% i)
  k <- c("lateral + vertical",i, sum(sub$avg_share_mechanism_annual), sum(sub$sd_mechanism))
  clo.size <- rbind(clo.size, k)
  clo.size$avg_share_mechanism_annual <- as.numeric(clo.size$avg_share_mechanism_annual)
  clo.size$sd_mechanism <- as.numeric(clo.size$sd_mechanism)
}

clo.size$closure_mechanism <- as.factor(clo.size$closure_mechanism)
clo.size$closure_mechanism <-  ordered(clo.size$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  


#plot closure according to gap size bins
# tiff("gap_closure_NP_box_violin.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_violin() + 
#   geom_boxplot(width=0.1) +
#   theme_minimal()+ coord_flip() +
#   theme(legend.position="bottom") +My_Theme  +    scale_fill_viridis(discrete = TRUE) +
#   labs( x = "gap size [ha]", y= "% of gap area closing annually")
# dev.off()


tiff("gap_closure_mechanism_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
dev.off()


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 35),
  axis.text.x = element_text(size = 28),
  axis.text.y = element_text(size = 28),
  axis.title.y = element_text(size = 35),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=18),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16),
  legend.position="top") #bottom

tiff("gap_closure_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo, closure_mechanism %in% "lateral + vertical"), aes(x=gap.size , y=clo_share_annual, fill="green")) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")+ guides(fill = FALSE)  
dev.off()

#--- create panel view plot:

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 35),
  axis.text.x = element_text(size = 26,angle = 45, hjust=1),
  axis.text.y = element_text(size = 26),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=22),
  strip.text.x = element_text(size = 25),
  strip.text.y = element_text(size = 25),
  legend.position="top") #bottom

forest_data <- select(gap_clo, -aspect) # Drop the "aspect" column

new_df <- forest_data %>% 
  gather(category, feature, gap.size:elevation)


new_df$feature <-  ordered(new_df$feature, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1",
                                                      "Beech", "Spruce-fir-beech","Spruce","Larch-Pine",
                                                      "600-800", "800-1000","1000-1200","1200-1400", "1400-1600",  "1600-1800" ))

tiff("gap_closure_mechanism_panel_box.tiff", units="in", width=12, height=8, res=300)
ggplot(new_df, aes(x=feature , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "", y= "% of gap area closing annually", colour= "closure mechanism") +
  facet_grid(~category, scales="free_x", labeller = labeller(category = c("gap.size" = "Gap Size [ha]", "forest_type" = "Forest Type", "elevation" = "Elevation [m]")))
dev.off()

#--- panel plot with gap size and forest type

forest_data <- select(gap_clo, -aspect,-elevation) # Drop the "aspect" and "elevation" column


new_df <- forest_data %>% 
  gather(category, feature, gap.size:forest_type)


new_df$feature <-  ordered(new_df$feature, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1",
                                                      "Beech", "Spruce-fir-beech","Spruce","Larch-Pine"
                                                      ))

tiff("gap_closure_mechanism_panel_box_v2.tiff", units="in", width=12, height=8, res=300)
ggplot(new_df, aes(x=feature , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "", y= "% of gap area closing annually", colour= "closure mechanism") +
  facet_grid(~category, scales="free_x", labeller = labeller(category = c("gap.size" = "Gap Size [ha]", "forest_type" = "Forest Type")))
dev.off()


ggplot(new_df, aes(x=feature , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "", y= "% of gap area closing annually", colour= "closure mechanism") +
  coord_flip()+
  facet_grid(~category, scales="free_x", labeller = labeller(category = c("gap.size" = "Gap Size [ha]", "forest_type" = "Forest Type")))

ggplot(new_df, aes(x=feature , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "", y= "% of gap area closing annually", colour= "closure mechanism") +
  coord_flip()+
  facet_grid(~category, scales="free_y", labeller = labeller(category = c("gap.size" = "Gap Size [ha]", "forest_type" = "Forest Type"))) +
  scale_x_discrete(breaks = c("small", "medium", "large"), labels = c("Small", "Medium", "Large"))


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

#----- full plots closure vs. area scatter
# 
# tiff("gap_closure-abs_area_scatter_below1.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721  %>% filter(gap_area.ha <1) 
#        , aes(x=gap_area.ha, y=clo_area_sum_annual)) + geom_point()+
#   geom_smooth(colour="blue") +
# geom_smooth(colour="red", method = "lm") + My_Theme
# dev.off()
# 
# tiff("gap_closure-abs_area_scatter.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721 # %>% filter(gap_area.ha <1) 
#        , aes(x=gap_area.ha, y=clo_area_sum_annual)) + geom_point()+
#   geom_smooth(colour="blue") +
#   geom_smooth(colour="red", method = "lm") + My_Theme
# dev.off()
# 
# tiff("gap_closure-share_area_scatter_below1_time.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721  %>% filter(gap_area.ha <1) 
#        , aes(x=gap_area.ha, y=clo_share_sum_annual)) + geom_point()+
#   geom_smooth(colour="blue") +
#   geom_smooth(colour="red", method = "lm") +
#   facet_wrap(~timestep) + My_Theme
# dev.off()
# 
# tiff("gap_closure-abs_share_scatter.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721  #%>% filter(gap_area.ha <1) 
#        , aes(x=gap_area.ha, y=clo_share_sum_annual)) + geom_point()+
#   geom_smooth(colour="blue") +
#   geom_smooth(colour="red", method = "lm") + My_Theme #+
#   facet_wrap(~timestep) 
# dev.off()
# 
# 
# tiff("gap_closure-abs_share_scatter_timestep.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721  #%>% filter(gap_area.ha <1) 
#        , aes(x=gap_area.ha, y=clo_share_sum_annual)) + geom_point()+
#   geom_smooth(colour="blue") +
#   geom_smooth(colour="red", method = "lm") + My_Theme +
# facet_wrap(~timestep) 
# dev.off()


# average closure area annual
tiff("gap_closure_mean_are.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_area_annual)) + #avg_clo_area_annual
 # geom_bar(stat = "identity", color="black", position=position_dodge()) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "gap area [m2] closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_area_annual-sd_all.a, ymax=avg_clo_area_annual+sd_all.a))
dev.off()

# mean closure shares


# tiff("gap_closure_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + #avg_clo_area_annual
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   #facet_wrap(~forest_type) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()


# point with error bar

tiff("gap_closure_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all))
dev.off()
# --- mechanism ---


# tiff("gap_closure-mechanism_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#  # facet_wrap(~forest_type) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "closure mechanism")
# dev.off()

# point with error bar
tiff("gap_closure-mechanism_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_NP.sizebins, aes(x=gap.size , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=1)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=1))
dev.off()

tiff("gap_closure-mechanism.tiff", units="in", width=12, height=8, res=300)
ggplot(clo.size, aes(x=gap.size , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=1)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=1))
dev.off()

# --- closure by aspect---

#calculate average closure share per aspect
clo_per_aspect <- gap_clo_NP_91721%>% group_by(aspect, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),2),
            avg_share_mechanism_annual = round(mean(clo_share_annual),2),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            n.obs = n())

# aggregate closure mechanism information for plotting
clo.aspect <- as.data.frame(clo_per_aspect)
clo.aspect <- clo.aspect[,c("closure_mechanism", "aspect",  "avg_share_mechanism_annual","sd_mechanism") ]

aspct <- as.character(unique(clo.aspect$aspect))
clo.aspect$closure_mechanism <- as.character(clo.aspect$closure_mechanism)

for(i in aspct) {
  sub <- subset(clo.aspect, aspect %in% i)
  k <- c("lateral + vertical",i, sum(sub$avg_share_mechanism_annual), sum(sub$sd_mechanism))
  clo.aspect <- rbind(clo.aspect, k)
  clo.aspect$avg_share_mechanism_annual <- as.numeric(clo.aspect$avg_share_mechanism_annual)
  clo.aspect$sd_mechanism <- as.numeric(clo.aspect$sd_mechanism)
}

clo.aspect$closure_mechanism <- as.factor(clo.aspect$closure_mechanism)
clo.aspect$closure_mechanism <-  ordered(clo.aspect$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  


# clo_per_aspect.sizebins <- gap_clo_NP_91721%>% group_by(aspect, closure_mechanism, gap.size) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),2),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),2),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             sd_all = sd(clo_share_sum_annual),
#             sd_mechanism= sd(clo_share_annual),
#             n.obs = n())


# # closure share per aspect  - boxplot
# tiff("gap_closure_aspect.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_sum_annual)) + geom_boxplot() + 
#    facet_wrap(~aspect) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "gap size [ha]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
# dev.off()
# 
# # closure mechanism according to aspect - boxplot
# tiff("gap_closure_mechanism_aspect.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) + geom_boxplot() + facet_wrap(~aspect) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
# dev.off()


#-- mean -bar

# tiff("gap_closure_aspect_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_aspect, aes(x=aspect , y=avg_clo_share_annual )) + geom_bar(stat = "identity") + 
#   # facet_wrap(~aspect) +
#  # facet_grid(aspect~.) +
#   theme_minimal()  + My_Theme + coord_flip() +
#   labs(x = "aspect", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
# dev.off()
# 
# tiff("gap_closure_mechanism_aspect_mean.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_aspect, aes(x=aspect , y=share_mechanism_annual , fill=closure_mechanism )) + geom_bar(stat = "identity", position="dodge") + 
#   # facet_wrap(~aspect) +
#   # facet_grid(aspect~.) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "aspect", y= "% of gap area closing annually") +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") 
# dev.off()

# --mean point - poinrange

tiff("gap_closure_aspect_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_aspect, aes(x=aspect , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "aspect", y= "average % of gap area closing annually")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all))
dev.off()

tiff("gap_closure_mechanism_aspect_mean.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_aspect, aes(x=aspect , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "aspect", y= "average % of gap area closing annually")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3))
dev.off()

tiff("gap_closure_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(clo.aspect, aes(x=aspect , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "aspect", y= "average % of gap area closing annually")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3))
dev.off()


# --- closure by elevation ---


#calculate average closure share per elevation
clo_per_elevation<- gap_clo_NP_91721 %>% group_by(elevation) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
            avg_share_mechanism_annual = round(mean(clo_share_annual),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            n.obs = n())

clo_per_elevation.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
            avg_share_mechanism_annual = round(mean(clo_share_annual),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            n.obs = n())

# aggregate closure mechanism information for plotting
clo.elevation <- as.data.frame(clo_per_elevation.mechanism)
clo.elevation <- clo.elevation[,c("closure_mechanism", "elevation",  "avg_share_mechanism_annual","sd_mechanism") ]

elev <- as.character(unique(clo.elevation$elevation))
clo.elevation$closure_mechanism <- as.character(clo.elevation$closure_mechanism)

for(i in elev) {
  sub <- subset(clo.elevation, elevation %in% i)
  k <- c("lateral + vertical",i, sum(sub$avg_share_mechanism_annual), sum(sub$sd_mechanism))
  clo.elevation <- rbind(clo.elevation, k)
  clo.elevation$avg_share_mechanism_annual <- as.numeric(clo.elevation$avg_share_mechanism_annual)
  clo.elevation$sd_mechanism <- as.numeric(clo.elevation$sd_mechanism)
}

clo.elevation$closure_mechanism <- as.factor(clo.elevation$closure_mechanism)
clo.elevation$closure_mechanism <-  ordered(clo.elevation$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  


# clo_per_elevation <- gap_clo_NP_91721 %>% group_by(elevation) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             n.obs = n())
#------

# clo_per_elevation.sizebins.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, closure_mechanism, gap.size) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             n.obs = n())
# 
# clo_per_elevation.sizebins <- gap_clo_NP_91721 %>% group_by(elevation, gap.size) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             n.obs = n())

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

# tiff("gap_closure_elevation.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=elevation ,y=clo_share_sum_annual)) +geom_boxplot() + 
#   # facet_wrap(~elevation) +
#  # facet_grid(elevation~.) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "elevation [m]", y= "% of gap area closing annually") +    scale_fill_viridis(discrete = TRUE) 
# dev.off()
# 
# 
# # closure mechanism per elevation 
# 
# tiff("gap_closure_mechanism_elevation.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=elevation , y=clo_share_annual, fill=closure_mechanism)) + 
#   geom_boxplot() +# facet_wrap(~elevation) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "elevation [m]", y= "% of gap area closing annually", colour= "closure mechanism")
# dev.off()

# --- mean closure rates per elevation and gap size ---

#order elevation labels
# clo_per_elevation.sizebins.mechanism$elevation <- ordered(clo_per_elevation.sizebins.mechanism$elevation,
#                                       levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
# clo_per_elevation.sizebins.mechanism$elevation <- factor(clo_per_elevation.sizebins.mechanism$elevation,levels=rev(levels(clo_per_elevation.sizebins.mechanism$elevation)))

# tiff("gap_closure_elevation_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_elevation, aes(x=elevation , y=avg_clo_share_annual )) + 
#   geom_bar(stat = "identity") +
#  # facet_wrap(~elevation) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()
# 
# tiff("gap_closure_mechanism_elevation_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_elevation.mechanism, aes(x=elevation , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", position=position_dodge()) +
#   #facet_wrap(~elevation) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()

# --- mean point +pointrange ----

tiff("gap_closure_elevation_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation, aes(x=elevation , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all))
dev.off()

tiff("gap_closure_mechanism_elevation_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_elevation.mechanism, aes(x=elevation , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3))
dev.off()

tiff("gap_closure_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(clo.elevation, aes(x=elevation , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.5)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "elevation [m]", y= "average % of gap area closing annually")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.5))
dev.off()


#reversing elevation labels
# clo_per_elevation.sizebins$elevation <- ordered(clo_per_elevation.sizebins$elevation,
#                                       levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
# clo_per_elevation.sizebins.mechanism$elevation <- ordered(clo_per_elevation.sizebins.mechanism$elevation,
#                                                 levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
# 
# tiff("gap_closure_elevation_avg_sizebins.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_elevation.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
#   geom_bar(stat = "identity") +
#   facet_grid(elevation~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+ 
#   geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
# dev.off()
# 
# tiff("gap_closure_mechanism_elevation_avg_sizebins.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_elevation.sizebins.mechanism, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   facet_grid(elevation~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()


# --- closure by forest type ---

#calculate average closure share per forest type
clo_per_ftype <- gap_clo_NP_91721%>% group_by(forest_type, closure_mechanism) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
            avg_share_mechanism_annual = round(mean(clo_share_annual),4),
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual)

# aggregate closure mechanism information for plotting
clo.forest <- as.data.frame(clo_per_ftype)
clo.forest <- clo.forest[,c("closure_mechanism", "forest_type",  "avg_share_mechanism_annual","sd_mechanism") ]

ftype <- as.character(unique(clo.forest$forest_type))
clo.forest$closure_mechanism <- as.character(clo.forest$closure_mechanism)

for(i in ftype) {
  sub <- subset(clo.forest, forest_type %in% i)
  k <- c("lateral + vertical",i, sum(sub$avg_share_mechanism_annual), sum(sub$sd_mechanism))
  clo.forest <- rbind(clo.forest, k)
  clo.forest$avg_share_mechanism_annual <- as.numeric(clo.forest$avg_share_mechanism_annual)
  clo.forest$sd_mechanism <- as.numeric(clo.forest$sd_mechanism)
}

clo.forest$closure_mechanism <- as.factor(clo.forest$closure_mechanism)
clo.forest$closure_mechanism <-  ordered(clo.forest$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  


# clo_per_ftype.sizebins.mechanism <- gap_clo_NP_91721 %>% group_by(forest_type, closure_mechanism, gap.size) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             sd_all = sd(clo_share_sum_annual),
#             sd_mechanism= sd(clo_share_annual),
#             n.obs = n())
# 
# clo_per_ftype.sizebins <- gap_clo_NP_91721 %>% group_by(forest_type, gap.size) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             sd_all = sd(clo_share_sum_annual),
#             sd_mechanism= sd(clo_share_annual),
#             n.obs = n())


# closure share per forest type 
# 
# tiff("gap_closure_ftype_box.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=forest_type , y=clo_share_sum_annual)) + geom_boxplot() +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "forest type", y= "% of gap area closing annually") +   scale_fill_viridis(discrete = TRUE)
# dev.off()
# 
# 
 # closure mechanism according to forest type 

# boxplots

 tiff("gap_closure_mechanism_ftype_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=forest_type , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot() +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "% of gap area closing annually", colour= "closure mechanism")
 dev.off()

 #order labels
 gap_clo$forest_type <- ordered(gap_clo$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-Pine"))
 gap_clo$forest_type <- factor(gap_clo$forest_type,levels=rev(levels(gap_clo$forest_type)))
 
 tiff("gap_closure_ftype_box.tiff", units="in", width=12, height=8, res=300)
 ggplot(subset(gap_clo, closure_mechanism %in% "lateral + vertical"), aes(x=forest_type , y=clo_share_annual, fill="green")) +
   geom_boxplot() +
   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
   labs(x = "forest type", y= "% of gap area closing annually", colour= "closure mechanism")+ guides(fill = FALSE)   
 dev.off()
 
 
 
# --- mean closure rates per forest type
# tiff("gap_closure_ftype_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype, aes(x=forest_type , y=avg_clo_share_annual)) + 
#   geom_bar(stat = "identity") +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="forest type", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()

# point with error bar
tiff("gap_closure_ftype_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_clo_share_annual)) + 
  geom_point(shape = 21, fill = "black",color = "black", size = 7) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all),linewidth = 2)
dev.off()

# tiff("gap_closure_mechanism_ftype_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype, aes(x=forest_type , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="forest type", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()

# point with error bar
tiff("gap_closure_mechanism_ftype_avg.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype, aes(x=forest_type , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 7, position=position_dodge(width=0.3)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.3),linewidth = 2)
dev.off()


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 24),
  axis.text.x = element_text(size = 23),
  axis.text.y = element_text(size = 23),
  axis.title.y = element_text(size = 24),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=24),
  legend.text = element_text(size=24),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16),
  legend.position="bottom")


tiff("gap_closure_ftype2.tiff", units="in", width=12, height=8, res=300)
ggplot(clo.forest, aes(x=forest_type , y=avg_share_mechanism_annual, colour= closure_mechanism, group= closure_mechanism, fill=closure_mechanism)) + 
  geom_point(shape = 21, size = 9, position=position_dodge(width=0.5)) +
  #facet_wrap(~forest_type) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + 
  scale_colour_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "forest type", y= "average % of gap area closing annually", colour= "forest type")+ 
  geom_pointrange(aes(ymin=avg_share_mechanism_annual-sd_mechanism, ymax=avg_share_mechanism_annual+sd_mechanism), position = position_dodge(width=0.5),linewidth = 2)+
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
dev.off()

# --- mean closure rates per forest type and gap size ---

# closure mechanism mean according to forest type and gap size
# ggplot(clo_per_ftype.sizebins, aes(x=gap.size , y=avg_clo_share_annual, fill=forest_type)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   #facet_wrap(~forest_type) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
# 
# tiff("gap_closure_ftype_size_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype.sizebins, aes(x=gap.size , y=avg_clo_share_annual)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   facet_grid(forest_type~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")+
#   geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
# dev.off()
# 
# tiff("gap_closure_mechanism_ftype_size_avg.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype.sizebins.mechanism, aes(x=gap.size , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", color="black", position=position_dodge()) +
#   facet_grid(forest_type~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x="gap size [ha]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()



# ---- closure per elevation and forest type 

clo_per_ftype.elev <- gap_clo_NP_91721 %>% group_by(elevation, forest_type) %>%
  summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
            avg_share_mechanism_annual = round(mean(clo_share_annual),4),
            percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
            sd_all = sd(clo_share_sum_annual),
            sd_mechanism= sd(clo_share_annual),
            n.obs = n())

# clo_per_ftype.elev.mechanism <- gap_clo_NP_91721 %>% group_by(elevation, forest_type,closure_mechanism) %>%
#   summarise(avg_clo_share_annual = round(mean(clo_share_sum_annual),4),
#             avg_share_mechanism_annual = round(mean(clo_share_annual),4),
#             percent_share_mechanism = avg_share_mechanism_annual/ avg_clo_share_annual,
#             sd_all = sd(clo_share_sum_annual),
#             sd_mechanism= sd(clo_share_annual),
#             n.obs = n())

#clo_per_ftype.elev$elev.ftype <- paste0(clo_per_ftype.elev$elevation ,"-" ,clo_per_ftype.elev$forest_type)

# tiff("gap_closure_elev_ftype.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype.elev, aes(x=elevation , y=avg_clo_share_annual)) + 
#   geom_bar(stat = "identity") +
#   facet_grid(forest_type~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")+ 
#   geom_text(aes(label = paste0("n:",n.obs), hjust = -0.8))
# dev.off()

# ---- mean point + pointrange

tiff("gap_closure_elev_ftype.tiff", units="in", width=12, height=8, res=300)
ggplot(clo_per_ftype.elev, aes(x=elevation , y=avg_clo_share_annual, color= forest_type, fill= forest_type)) + 
  geom_point(shape = 21, size = 6, position=position_dodge(width=1) ) +
 # facet_grid(forest_type~.) +
  theme_minimal()+  My_Theme +  coord_flip() +
  labs(x = "elevation [m]", y= "average % of gap area closing annually")+ 
  # scale_fill_brewer(palette="Set2", name = "forest type")+
  # scale_colour_brewer(palette="Set2", name = "forest type")+
  scale_fill_brewer(palette="Dark2", name = "forest type")+
  scale_colour_brewer(palette="Dark2", name = "forest type")+
  geom_pointrange(aes(ymin=avg_clo_share_annual-sd_all, ymax=avg_clo_share_annual+sd_all), position = position_dodge(width=1))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
dev.off()

# tiff("gap_closure_mechanism_elev_ftype.tiff", units="in", width=12, height=8, res=300)
# ggplot(clo_per_ftype.elev.mechanism, aes(x=elevation , y=avg_share_mechanism_annual, fill=closure_mechanism)) + 
#   geom_bar(stat = "identity", position="dodge") +
#   facet_grid(forest_type~.) +
#   theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
#   labs(x = "elevation [m]", y= "average % of gap area closing annually", colour= "forest type")
# dev.off()


# --- closure by shape complexity ---

# # closure share vs. complexity 
# 
# tiff("gap_closure_complexity.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_sum_annual)) + geom_point(size = 5, alpha = 0.5) +
#   theme_minimal()+ coord_flip()  + My_Theme +
#   labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually")  +
#   geom_smooth(method=lm) 
# dev.off()
# 
# # closure mechanism vs. complexity vs. closure share
# 
# tiff("gap_closure_mechanism_complexity.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_annual, color=closure_mechanism)) + 
#   geom_point(size = 5, alpha = 0.5) +  scale_color_brewer(palette="Dark2", name = "closure mechanism")+
#   theme_minimal()+ coord_flip()   + My_Theme +
#   labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually", colour= "closure mechanism") +
#   geom_smooth(method=lm) 
# dev.off()
# 
# 
# # closure mechanism vs. complexity vs. time step
# 
# tiff("gap_closure_mechanism_complexity.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , y=clo_share_annual, color=timestep)) + 
#   geom_point(size = 5, alpha = 0.5) +  scale_color_brewer(palette="Dark2", name = "closure mechanism")+
#   theme_minimal()+ coord_flip()   + My_Theme +
#   labs(x = "gap complexity (p:a ratio)", y= "% of gap area closing annually", colour= "closure mechanism") 
# dev.off()
# 
# # p:a vs. aspect
# tiff("pa_aspect.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=aspect)) + 
#   geom_boxplot()+ coord_flip() 
# dev.off()
# 
# # p:a vs. elevation
# tiff("pa_elevation.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=elevation)) + 
#   geom_boxplot()+ coord_flip() 
# dev.off()
# 
# # p:a vs. forest type
# tiff("pa_forest.tiff", units="in", width=12, height=8, res=300)
# ggplot(gap_clo_NP_91721, aes(x=pa_ratio , color=forest_type)) + 
#   geom_boxplot()+ coord_flip() 
# dev.off()
# 
# #-------------------------------- model check of influence of gap complexity on closure share
# 
# gap.closure.lm<-lm(clo_share_annual ~ pa_ratio + aspect  + forest_type, data = gap_clo_NP_91721)
# 
# summary(gap.closure.lm)
# 
# par(mfrow=c(2,2))
# plot(gap.closure.lm)
# par(mfrow=c(1,1))

    