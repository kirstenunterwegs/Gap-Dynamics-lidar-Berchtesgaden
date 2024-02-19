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

# --- load  Gap layers  ----

gaps2009.3 <- rast("processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2017.3 <- rast("processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")

gaps2009.10 <- rast("processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2017.10 <- rast("processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")


###################################################------create forest edge mask for lateral growth classification

# height thres 3m

boundaries9_h3 <- boundaries(gaps2009.3, directions=8, inner=TRUE)
writeRaster(boundaries9_h3, "processed/sensitivity/height_sensitivity/gap_boundaries9_h3.tif")

boundaries17_h3 <- boundaries(gaps2017.3, directions=8, inner=TRUE)
writeRaster(boundaries17_h3, "processed/sensitivity/height_sensitivity/gap_boundaries17_h3.tif")

# height thres 10m

boundaries9_h10 <- boundaries(gaps2009.10, directions=8, inner=TRUE)
writeRaster(boundaries9_h10, "processed/sensitivity/height_sensitivity/gap_boundaries9_h10.tif")

boundaries17_h10 <- boundaries(gaps2017.10, directions=8, inner=TRUE)
writeRaster(boundaries17_h10, "processed/sensitivity/height_sensitivity/gap_boundaries17_h10.tif")


####################################################### classify vertical and horizontal closure ##################################


#load closure areas with growth information

clo_growth_917_h3 <- rast("processed/sensitivity/height_sensitivity/closure_area_growth_917_h3.tif") 
clo_growth_1721_h3 <- rast("processed/sensitivity/height_sensitivity/closure_area_growth_1721_h3.tif")

clo_growth_917_h10 <- rast("processed/sensitivity/height_sensitivity/closure_area_growth_917_h10.tif") 
clo_growth_1721_h10 <- rast("processed/sensitivity/height_sensitivity/closure_area_growth_1721_h10.tif")


#adjust extents
clo_growth_917_h3 <- crop(clo_growth_917_h3, boundaries9_h3)
clo_growth_1721_h3 <- crop(clo_growth_1721_h3, boundaries17_h3)

clo_growth_917_h10 <- crop(clo_growth_917_h10, boundaries9_h10)
clo_growth_1721_h10 <- crop(clo_growth_1721_h10, boundaries17_h10)

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

# height thres 3m 

gap_closure_mechanism917_h3 <- gap_closure_mechanism917(clo_growth_917_h3, boundaries9_h3)
gap_closure_mechanism1721_h3 <- gap_closure_mechanism1721(clo_growth_1721_h3, boundaries17_h3)

# height thres 10m 

gap_closure_mechanism917_h10 <- gap_closure_mechanism917(clo_growth_917_h10, boundaries9_h10)
gap_closure_mechanism1721_h10 <- gap_closure_mechanism1721(clo_growth_1721_h10, boundaries17_h10)


terra::writeRaster(gap_closure_mechanism917_h3, "processed/sensitivity/height_sensitivity/gap_closure_mechanism917_h3.tif") 
terra::writeRaster(gap_closure_mechanism1721_h3, "processed/sensitivity/height_sensitivity/gap_closure_mechanism1721_h3.tif")

terra::writeRaster(gap_closure_mechanism917_h10, "processed/sensitivity/height_sensitivity/gap_closure_mechanism917_h10.tif") 
terra::writeRaster(gap_closure_mechanism1721_h10, "processed/sensitivity/height_sensitivity/gap_closure_mechanism1721_h10.tif")



####################################################################################################################
# prepare closure mechanism dfs per timestep
###################################################################################################################

# --- 2009 - 2017 --- heigh thres 3m

# merge closure mechanism with gaps
gap_closure_mechanism_stack <- c(gap_closure_mechanism917_h3, gaps2009.3)

# convert to df
gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

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
sum(gap_clo_per_id$contraction) # 0 gaps do not experience any closure from 2009-2017

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
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_917_h3<-gap_clo_per_id_nona %>% 
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

# --- 2009 - 2017 --- heigh thres 10 m

# merge closure mechanism with gaps
gap_closure_mechanism_stack <- c(gap_closure_mechanism917_h10, gaps2009.10)

# convert to df
gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

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
sum(gap_clo_per_id$contraction) # 0 gaps do not experience any closure from 2009-2017

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
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_917_h10<-gap_clo_per_id_nona %>% 
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


# --- 2017 - 2021 ---  heigh thres 3m

# merge closure mechanism with gaps
gap_closure_mechanism_stack <- c(gap_closure_mechanism1721_h3, gaps2017.3)

# convert to df
gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

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
sum(gap_clo_per_id$contraction) # 2 gaps do not experience any closure from 2009-2017

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
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_1721_h3<-gap_clo_per_id_nona %>% 
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


# --- 2017 - 2021 ---  heigh thres 10 m

# merge closure mechanism with gaps
gap_closure_mechanism_stack <- c(gap_closure_mechanism1721_h10, gaps2017.10)

# convert to df
gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

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
sum(gap_clo_per_id$contraction) # 2 gaps do not experience any closure from 2009-2017

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
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_1721_h10<-gap_clo_per_id_nona %>% 
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



### -------------------------------merge 9-17 and 17-21 dfs for comparison ----------------------- ###


# --- heigh thres 3m ---

gap_clo_per_id_nona_1721_h3$timestep <- "17-21"
gap_clo_per_id_nona_917_h3$timestep <- "9-17"

gap_clo_NP_91721_h3 <- rbind(gap_clo_per_id_nona_917_h3, gap_clo_per_id_nona_1721_h3)
#rearrange timestep labels
gap_clo_NP_91721_h3$timestep <- factor(gap_clo_NP_91721_h3$timestep , levels=c("9-17", "17-21"))

# calculate annual closure rates

gap_clo_NP_91721_h3<- gap_clo_NP_91721_h3 %>% mutate(time = as.numeric(recode(timestep,
                                          `9-17`= 8, 
                                          `17-21`= 4)),
                                          clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                          clo_share_annual = round(closure_share/time,4)*100,
                                          clo_area_sum_annual = round(closure_area_sum/time,4),
                                          clo_area_annual = round(closure_area/time,4))



# --- append lateral + vertical closure info to main df for distribution display ----

# create unique ID
gap_clo_NP_91721_h3$id <- as.numeric(paste0(gap_clo_NP_91721_h3$gap_id,gap_clo_NP_91721_h3$time))

gap_clo <- as.data.frame(gap_clo_NP_91721_h3[,c("id", "closure_mechanism",  "gap.size", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)

id <- as.character(unique(gap_clo$id))
#i=74
for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  k <- c(i,"lateral + vertical", size, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)

#order labels
gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.01-0.04","0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

gap_clo_h3 <- gap_clo

saveRDS(gap_clo_h3, "processed/sensitivity/height_sensitivity/clo_analysis_ready_h3.rds") 




# --- heigh thres 10m ---

gap_clo_per_id_nona_1721_h10$timestep <- "17-21"
gap_clo_per_id_nona_917_h10$timestep <- "9-17"

gap_clo_NP_91721_h10 <- rbind(gap_clo_per_id_nona_917_h10, gap_clo_per_id_nona_1721_h10)
#rearrange timestep labels
gap_clo_NP_91721_h10$timestep <- factor(gap_clo_NP_91721_h10$timestep , levels=c("9-17", "17-21"))

# calculate annual closure rates

gap_clo_NP_91721_h10<- gap_clo_NP_91721_h10 %>% mutate(time = as.numeric(recode(timestep,
                                                                              `9-17`= 8, 
                                                                              `17-21`= 4)),
                                                     clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                                     clo_share_annual = round(closure_share/time,4)*100,
                                                     clo_area_sum_annual = round(closure_area_sum/time,4),
                                                     clo_area_annual = round(closure_area/time,4))



# --- append lateral + vertical closure info to main df for distribution display ----

# create unique ID
gap_clo_NP_91721_h10$id <- as.numeric(paste0(gap_clo_NP_91721_h10$gap_id,gap_clo_NP_91721_h10$time))

gap_clo <- as.data.frame(gap_clo_NP_91721_h10[,c("id", "closure_mechanism",  "gap.size", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)

id <- as.character(unique(gap_clo$id))
#i=74
for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  k <- c(i,"lateral + vertical", size, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "lateral + vertical"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)

#order labels
gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.01-0.04","0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

gap_clo_h10 <- gap_clo

saveRDS(gap_clo_h10, "processed/sensitivity/height_sensitivity/clo_analysis_ready_h10.rds") 




# ---- load gap closure information of original mmu400 heigh threshold 5m gap layer and merge

gap_clo_h5<-  readRDS("processed/sensitivity/mmu400_height5/clo_analysis_ready.rds")

# reduce to common columns 

columns_to_keep <- colnames(gap_clo_h3)
# Reduce gap_features921_h5 to only those columns in gap_features921_h10
gap_clo_h5 <- gap_clo_h5 %>% select(all_of(columns_to_keep))


# --- assign thresholds and prepare for plotting 

gap_clo_h3$h_thres <- as.factor(3) 
gap_clo_h5$h_thres <- as.factor(5)
gap_clo_h10$h_thres <- as.factor(10)

gap_clo<- rbind(gap_clo_h3, gap_clo_h5, gap_clo_h10)

# -----

# prepare plotting

require(scales)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/results/sensitivity_analysis/height_threshold/"
setwd(wd)


# --- closure rates as boxplots



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
ggplot(subset(gap_clo, closure_mechanism %in% "lateral + vertical"), aes(x=gap.size , y=clo_share_annual, fill=h_thres)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  
  scale_fill_brewer(palette="Set1", name = "Height threshold") + My_Theme +
  labs(x = "gap size ( ha )", y= "% of gap area closing annually", colour= "Height threshold") #+
  geom_label(data = subset(gap_clo, closure_mechanism %in% "lateral + vertical") %>% dplyr::group_by(gap.size, h_thres) %>% dplyr::summarise(N = n(), clo_share_annual = 22),
             aes(label=paste("n = ", N), fill=h_thres), position=position_dodge(width = 0.7))
dev.off()


