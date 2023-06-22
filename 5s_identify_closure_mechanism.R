###############################################
# Identify closure mechanism
##############################################

library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)
library(stringr)  


wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load Gap layers ----

#gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu400.sensitivity.tif") # layer have been cropped previously to the research area
gaps2009<- gap_stack[[1]]
gaps2017<- gap_stack[[2]]



###################################################------create forest edge mask for lateral growth classification


boundaries9 <- boundaries(gaps2009, directions=8, inner=TRUE)
#writeRaster(boundaries9, "processed/sensitivity/gap_boundaries9.tif")
writeRaster(boundaries9, "processed/sensitivity/version.mmu400/gap_boundaries9.tif")

boundaries17 <- boundaries(gaps2017, directions=8, inner=TRUE)
#writeRaster(boundaries17, "processed/sensitivity/gap_boundaries17.tif")
writeRaster(boundaries9, "processed/sensitivity/version.mmu400/gap_boundaries17.tif")

####################################################### classify vertical and horizontal closure ##################################

#load closure areas with growth information
# clo_growth_917 <- rast("processed/sensitivity/closure_area_growth_917.tif")   # mmu100
# clo_growth_1721 <- rast("processed/sensitivity/closure_area_growth_1721.tif")

clo_growth_917 <- rast("processed/sensitivity/version.mmu400/closure_area_growth_917.tif") # mmu400
clo_growth_1721 <- rast("processed/sensitivity/version.mmu400/closure_area_growth_1721.tif")

#load gap boundaries
# boundaries.2009 <- rast("processed/sensitivity/gap_boundaries9.tif")  # mmu100
# boundaries.2017 <- rast("processed/sensitivity/gap_boundaries17.tif")

boundaries.2009 <- rast("processed/sensitivity/version.mmu400/gap_boundaries9.tif") # mmu400
boundaries.2017 <- rast("processed/sensitivity/version.mmu400/gap_boundaries17.tif")

#adjust extents
clo_growth_917 <- crop(clo_growth_917, boundaries.2009)
clo_growth_1721 <- crop(clo_growth_1721, boundaries.2017)

# define function to differentiate between regeneration (vertical) and crown plasticity (horizontal)

# need to differ between both timesteps due to differnet growing periods


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


# terra::writeRaster(gap_closure_mechanism917, "processed/sensitivity/gap_closure_mechanism917.tif", overwrite=TRUE) # mmu100
# terra::writeRaster(gap_closure_mechanism1721, "processed/sensitivity/gap_closure_mechanism1721.tif", overwrite=TRUE)

terra::writeRaster(gap_closure_mechanism917, "processed/sensitivity/version.mmu400/gap_closure_mechanism917.tif", overwrite=TRUE) # mmu400
terra::writeRaster(gap_closure_mechanism1721, "processed/sensitivity/version.mmu400/gap_closure_mechanism1721.tif", overwrite=TRUE)

###----------- analyze closure per gap (size), elevation, aspect, management and forest type ----------- ###


####################################################################################################################
# prepare closure mechanism dfs per timestep
###################################################################################################################

# #prepare layers to limit analysis to core zone, below 1800m and with forest type information
# 
# # --- load NP information 
# 
# foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")
# management <- vect("F:/Projects/CanopyDynamicsBDG/data/NP_data/npb_zonierung_22_epsg25832.shp")
# aspect<-  rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/aspect_2021_classified_1m.tif")
# elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")
# 
# # exclude management zone
# core.zone <- subset(management, management$zone_id == 4, c(1:2))


# 2009 - 2017
###################################################################################################################

#merge closure mechanism with gaps
gap_closure_mechanism917 <- rast( "processed/sensitivity/gap_closure_mechanism917.tif")
gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gaps2009<- gap_stack[[1]]

gaps2009 <- crop(gaps2009, gap_closure_mechanism917)
gap_closure_mechanism_stack <- c(gap_closure_mechanism917, gaps2009)

#crop to same extent
# elevation.below1800 <- crop(elevation.below1800,gap_closure_mechanism917 )
# foresttype <- crop(foresttype, gap_closure_mechanism917 )

#mask down to reserach area
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, elevation.below1800)
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, foresttype)
# gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, core.zone)

gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

# saveRDS(gap_closure_mechanism_stack.df,"processed/closure/updated/gap_closure_mechanism_pergap_917.rds" )
# 
# gap_closure_mechanism_stack.df <- readRDS("processed/closure/updated/gap_closure_mechanism_pergap_917.rds")

# aggregate closure and gap information 
gap_clo_per_id <-  gap_closure_mechanism_stack.df %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#drop gaps < 100 m2 / 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut) #change to 100 for sensitivity!

#gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 100,]
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 2 gaps do not experience any closure from 2009-2017 (previously 7, maybe the masking cut off some gap areas)


gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing

gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
        closure_mechanism = as.factor(recode(closure_mechanism,
                                             `0`="no closure", 
                                             `1`="lateral closure",
                                             `2`="vertical closure")))

##### if I want to include the no closure pixels, I have to disable following line: !!!!
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
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.009, 0.04,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.009,0.04]`="0.01-0.04",
                                     `(0.04,0.1]`="0.04-0.1",
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


# 2017- 2021
###################################################################################################################

# #merge closure mechanism with gaps
gap_closure_mechanism1721 <- rast("processed/sensitivity/gap_closure_mechanism1721.tif")
gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gaps2017<- gap_stack[[2]]


gaps2017 <- crop(gaps2017, gap_closure_mechanism1721)
gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721, gaps2017)

#mask down to reserach area
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, elevation.below1800)
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, foresttype)
# gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, core.zone)

gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")

# saveRDS(gap_closure_mechanism_stack.df_1721,"processed/closure/updated/gap_closure_mechanism_pergap_1721.rds" )
# 
# gap_closure_mechanism_stack.df_1721 <- readRDS("processed/closure/updated/gap_closure_mechanism_pergap_1721.rds" )

# aggregate closure and gap information - prepare df for plotting
gap_clo_per_id <-  gap_closure_mechanism_stack.df_1721 %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)

#gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 100,]
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify number of gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 2 gaps do not experience any closure from 2017-2021 

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
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.009,0.04,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.009,0.04]`="0.01-0.04",
                                     `(0.04,0.1]`="0.04-0.1",
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



### -------------------------------merge 9-17 and 17-21 dfs for comparison ----------------------- ###

gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
#rearrange timestep labels
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))

# ---add environmental feature information

#join with environmental information
#stats <- readRDS("processed/sensitivity/stats_sensitivity.rds") # mmu100
stats <- readRDS("processed/sensitivity/version.mmu400/stats_sensitivity.rds") # mmu400

# recode years to timesteps for merges
stats <- stats %>% mutate( timestep = as.factor(recode(year,
                                                          `2009`="9-17", 
                                                          `2017`="17-21")))


gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats, by = c("gap_id", "timestep"), all.x=TRUE) #[ , c("gap_id",, "forest_type", "elevation", "aspect", "timestep")]


# exclude gaps in elevation band above 1800 m
gap_clo_NP_91721 <- gap_clo_NP_91721[gap_clo_NP_91721$elevation != "1800-2000",]
gap_clo_NP_91721 <- gap_clo_NP_91721[!is.na(gap_clo_NP_91721$elevation),] # only optional if there are NAs in the df

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

gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.01-0.04","0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

# saveRDS(gap_clo, "processed/sensitivity/clo_analysis_ready.rds") # mmu100
gap_clo <- readRDS("processed/sensitivity/clo_analysis_ready.rds")

# saveRDS(gap_clo, "processed/sensitivity/version.mmu400/clo_analysis_ready.rds") # mmu400
# gap_clo <- readRDS("processed/sensitivity/version.mmu400/clo_analysis_ready.rds")

#-------
# reverse gap_id and year/time pasting (to append lateral + vertical closure)

# gap_clo$year <-  as.numeric(str_sub(gap_clo$id, - 1))  
# gap_clo$id <- as.numeric(substring(gap_clo$id, 1, nchar(gap_clo$id) - 1))
# 
# gap_clo <-  select(gap_clo,  -"year")


# ---- load gap closure information of original mmu400 gap layer and merge

gap_clo400 <-  readRDS("processed/sensitivity/version.mmu400/clo_analysis_ready.rds")

# sensitivity.ids <- readRDS("processed/sensitivity/origID_mmu400_sensitivityAoi.rds")
# 
# ids <- subset(sensitivity.ids, year %in% c("9", "17"))
# ids <- ids %>% mutate(year = (recode(year,  # from 2009-2017 = 8 years, from 2017-2021 = 4 years
#                                  `9`= 8, 
#                                  `17`=4)))
# 
# # reverse gap_id and year/time pasting (to append lateral + vertical closure)
# 
# gap_clo400$year <-  as.numeric(str_sub(gap_clo400$id, - 1))  
# gap_clo400$id <- as.numeric(substring(gap_clo400$id, 1, nchar(gap_clo400$id) - 1))
# 
# gap_clo400 <- subset(gap_clo400, id %in% ids$ids & year %in% ids$year)
# 
# gap_clo400 <-  select(gap_clo400,  -"year")


# --- assign mmus and prepare for plotting 

gap_clo400$mmu <- as.factor(400) # indicate mmu gap size
gap_clo$mmu <- as.factor(100)

gap_clo<- rbind(gap_clo, gap_clo400)

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

wd <- "processed/sensitivity/results/"
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


clo_NP.sizebins <- gap_clo_NP_91721 %>% group_by(closure_mechanism, gap.size, mmu) %>%
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


# --- closure rates as boxplots

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

tiff("gap_closure_mechanism_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism") +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)
dev.off()




My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 26),
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=26),
  legend.text = element_text(size=26),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16),
  legend.position="top") #bottom



tiff("gap_closure_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo, closure_mechanism %in% "lateral + vertical"), aes(x=gap.size , y=clo_share_annual, fill=mmu)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  
  scale_fill_brewer(palette="Dark2", name = "mmu") + My_Theme +
  labs(x = "gap size ( ha )", y= "% of gap area closing annually", colour= "mmu")# +
  geom_label(data = subset(gap_clo, closure_mechanism %in% "lateral + vertical") %>% dplyr::group_by(gap.size, mmu) %>% dplyr::summarise(N = n(), clo_share_annual = 22),
             aes(label=paste("n = ", N), fill=mmu), position=position_dodge(width = 0.7))
dev.off()


