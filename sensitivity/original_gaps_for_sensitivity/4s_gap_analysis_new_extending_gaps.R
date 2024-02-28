######################################
# Analysing Gaps
#####################################


library(plyr)
library(dplyr)
library(terra)
library(ggplot2)
require(scales)
library(tidyverse)


# --- load classified Gap layers of new and expanding gaps ----

gap_stack <- rast("data/processed/gaps_sensitivity/gap.stack.mmu400.sensitivity.tif") # layer have been cropped previously to the research area
gaps2017.id<- gap_stack[[2]]
gaps2021.id<- gap_stack[[3]]

# new-expanded classification

gaps2017 <- rast("data/processed/sensitivity/mmu400_height5/gaps2017_new_extended_stable_sensitivity.tif")
gaps2021 <- rast("data/processed/sensitivity/mmu400_height5/gaps2021_new_extended_stable_sensitivity.tif")

# crop gaps ID
gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 

#load expansion and closure layer and extract only expansion areas
exp_clo <- rast("data/processed/sensitivity/mmu400_height5/exp_clo_917_cn2cr2_mmu400n8_filtered.tif")
exp <- classify(exp_clo, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("data/processed/sensitivity/mmu400_height5/exp_clo_1721_cn2cr2_mmu400n8_filtered.tif")
exp1721 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas


# ---- crop layers to research sites:

# crop all layers to same extent

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))


foresttype <- crop(foresttype, gaps2017.id)
aspect<-  crop(aspect, gaps2017.id)
core.zone <- crop(core.zone, gaps2017.id)
elevation.below1800 <- crop(elevation, gaps2017.id)

# --- stack gap information and crop it

#2017

stack2017 <- c(gaps2017.id, gaps2017, exp, foresttype , elevation.below1800, aspect)
names(stack2017) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017 <- mask(stack2017, foresttype)

writeRaster(stack2017, "data/processed/sensitivity/mmu400_height5/stack.2017.all.gap.information.expansion_sensitivity.tif")
gap_stack_2017 <- rast("data/processed/sensitivity/mmu400_height5/stack.2017.all.gap.information.expansion_sensitivity.tif")


df <- as.data.frame(gap_stack_2017, na.rm = FALSE) 

df1 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df1, "data/processed/sensitivity/mmu400_height5/stack_2017_new_exp_df.rds")

#2021

stack21 <- c(gaps2021.id, gaps2021, exp1721, foresttype, elevation.below1800, aspect)
names(stack21) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack21 <- mask(stack21, foresttype)

writeRaster(stack21, "data/processed/sensitivity/mmu400_height5/stack.2021.all.gap.information.expansion_sensitivity.tif")
gap_stack_2021 <- rast("data/processed/sensitivity/mmu400_height5/stack.2021.all.gap.information.expansion_sensitivity.tif")

df <- as.data.frame(gap_stack_2021, na.rm = FALSE) 

df2 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion


write_rds(df2, "data/processed/sensitivity/mmu400_height5/stack_2021_new_exp_df.rds")


# --------calculate features per gap.id

df1 <- readRDS( "data/processed/sensitivity/mmu400_height5/stack_2017_new_exp_df.rds")
df2<- readRDS("data/processed/sensitivity/mmu400_height5/stack_2021_new_exp_df.rds")


gap_features_917 <- df1 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

gap_features_1721 <- df2 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))



saveRDS(gap_features_917,"data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_917.rds")
saveRDS(gap_features_1721,"data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_1721.rds")

