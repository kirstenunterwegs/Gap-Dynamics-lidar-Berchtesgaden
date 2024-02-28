###############################################
# Identify closure mechanism
##############################################

library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)


# --- load Gap layers ----

gap_stack <- rast("data/processed/gaps_sensitivity/gap.stack.mmu400.sensitivity.tif") # layer have been cropped previously to the research area
gaps2009<- gap_stack[[1]]
gaps2017<- gap_stack[[2]]



#------create forest edge mask for lateral growth classification


boundaries9 <- boundaries(gaps2009, directions=8, inner=TRUE)
writeRaster(boundaries9, "data/processed/sensitivity/mmu400_height5/gap_boundaries9.tif")

boundaries17 <- boundaries(gaps2017, directions=8, inner=TRUE)
writeRaster(boundaries17, "data/processed/sensitivity/mmu400_height5/gap_boundaries17.tif")


# ---- classify vertical and horizontal closure ----

#load closure areas with growth information
clo_growth_917 <- rast("data/processed/sensitivity/mmu400_height5/closure_area_growth_917.tif")
clo_growth_1721 <- rast("data/processed/sensitivity/mmu400_height5/closure_area_growth_1721.tif")

#load gap boundaries
boundaries.2009 <- rast("data/processed/sensitivity/mmu400_height5/gap_boundaries9.tif")
boundaries.2017 <- rast("data/processed/sensitivity/mmu400_height5/gap_boundaries17.tif")

#adjust extents
clo_growth_917 <- crop(clo_growth_917, boundaries.2009)
clo_growth_1721 <- crop(clo_growth_1721, boundaries.2017)

# define function to differentiate between regeneration (vertical) and crown plasticity (horizontal)

# need to differ between both time steps due to different growing periods

# --- add and deduct 20% to max height gain threshold for sensitivity analysis ---

# lower bound: 0.5 m * 0.8 = 0.4 m
# upper bound: 0.5 m * 1.2 = 0.6 m


# --- lower threshold

gap_closure_mechanism917_low <- function(diff_closure_layer, boundary_layer){       # 0.4 m * 8 yr = 3.2 m height gain

  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 3.2  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 3.2  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 3.2 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 3.2 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism1721_low <- function(diff_closure_layer, boundary_layer){      # 0.4 m * 4 yr = 1.6 m height gain

  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 1.6  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 1.6  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 1.6 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 1.6 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 


# --- higher threshold

gap_closure_mechanism917_high <- function(diff_closure_layer, boundary_layer){       # 0.6 m * 8 yr = 4.8 m height gain
  
  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 4.8  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 4.8  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 4.8 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 4.8 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 

gap_closure_mechanism1721_high <- function(diff_closure_layer, boundary_layer){      # 0.6 m * 4 yr = 2.4 m height gain
  
  closure_mechanism <- diff_closure_layer
  
  # classify change group
  closure_mechanism[diff_closure_layer > 2.4  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity) 
  closure_mechanism[diff_closure_layer >= 2.4  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 2.4 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 2.4 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
} 



gap_closure_mechanism917_low <- gap_closure_mechanism917_low(clo_growth_917, boundaries.2009)
gap_closure_mechanism1721_low <- gap_closure_mechanism1721_low(clo_growth_1721, boundaries.2017)

gap_closure_mechanism917_high <- gap_closure_mechanism917_high(clo_growth_917, boundaries.2009)
gap_closure_mechanism1721_high <- gap_closure_mechanism1721_high(clo_growth_1721, boundaries.2017)


writeRaster(gap_closure_mechanism917_low, "data/processed/sensitivity/growth_thres_sensitivity/gap_closure_mechanism917_low.tif")
writeRaster(gap_closure_mechanism1721_low, "data/processed/sensitivity/growth_thres_sensitivity/gap_closure_mechanism1721_low.tif")

writeRaster(gap_closure_mechanism917_high, "data/processed/sensitivity/growth_thres_sensitivity/gap_closure_mechanism917_high.tif")
writeRaster(gap_closure_mechanism17_high, "data/processed/sensitivity/growth_thres_sensitivity/gap_closure_mechanism1721_high.tif")

# --- prepare dataframes for analysis ----

# --- 2009 - 2017 --- low threshold

gaps2009 <- crop(gaps2009, gap_closure_mechanism917_low)
gap_closure_mechanism_stack <- c(gap_closure_mechanism917_low, gaps2009)


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
sum(gap_clo_per_id$contraction) #0

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor

#recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
        closure_mechanism = as.factor(recode(closure_mechanism,
                                             `0`="no closure", 
                                             `1`="lateral closure",
                                             `2`="vertical closure")))

# if I want to include the no closure pixels, I have to disable following line: !!!!
# exclude no closure shares
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


# --- 2017 - 2021 --- low threshold

gaps2017 <- crop(gaps2017, gap_closure_mechanism1721_low)
gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721_low, gaps2017)


gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")


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

# exclude no closure shares
# if I want to include the no closure pixels, I have to disable following line: !!!!
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


# ----merge 9-17 and 17-21 dfs for comparison ----

gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
#rearrange timestep labels
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))


# calculate annual closure rates

gap_clo_NP_91721<- gap_clo_NP_91721 %>% mutate(time = as.numeric(recode(timestep,
                                          `9-17`= 8, 
                                          `17-21`= 4)),
                                          clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                          clo_share_annual = round(closure_share/time,4)*100,
                                          clo_area_sum_annual = round(closure_area_sum/time,4),
                                          clo_area_annual = round(closure_area/time,4))



# --- append lateral + vertical closure info to main df for distribution display ---- (or total)
gap_clo_NP_91721$id <- as.numeric(paste0(gap_clo_NP_91721$gap_id,gap_clo_NP_91721$time))

gap_clo <- as.data.frame(gap_clo_NP_91721[,c("id", "closure_mechanism",  "gap.size", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)

id <- as.character(unique(gap_clo$id))

for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  aspect <- unique(sub$aspect)
  ftype <- unique(sub$forest_type)
  elev <- unique(sub$elevation)
  k <- c(i,"total", size, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "total"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)

#order labels
gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

gap_clo_low <- gap_clo

saveRDS(gap_clo_low, "data/processed/sensitivity/growth_thres_sensitivity/clo_analysis_ready_low.rds")

##################################################################################

# --- 2009 - 2017 --- high threshold

gaps2009 <- crop(gaps2009, gap_closure_mechanism917_high)
gap_closure_mechanism_stack <- c(gap_closure_mechanism917_high, gaps2009)


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
sum(gap_clo_per_id$contraction) #0

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


# --- 2017 - 2021 --- high threshold

gaps2017 <- crop(gaps2017, gap_closure_mechanism1721_high)
gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721_high, gaps2017)


gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)
#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")


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
##### if I want to include the no closure pixels, I have to disable following line: !!!!
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


# ----merge 9-17 and 17-21 dfs for comparison ----

gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
#rearrange timestep labels
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))


# calculate annual closure rates

gap_clo_NP_91721<- gap_clo_NP_91721 %>% mutate(time = as.numeric(recode(timestep,
                                          `9-17`= 8, 
                                          `17-21`= 4)),
                                          clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                          clo_share_annual = round(closure_share/time,4)*100,
                                          clo_area_sum_annual = round(closure_area_sum/time,4),
                                          clo_area_annual = round(closure_area/time,4))



# --- append lateral + vertical closure info to main df for distribution display ---- (or total)
gap_clo_NP_91721$id <- as.numeric(paste0(gap_clo_NP_91721$gap_id,gap_clo_NP_91721$time))

gap_clo <- as.data.frame(gap_clo_NP_91721[,c("id", "closure_mechanism",  "gap.size", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)

id <- as.character(unique(gap_clo$id))

for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  aspect <- unique(sub$aspect)
  ftype <- unique(sub$forest_type)
  elev <- unique(sub$elevation)
  k <- c(i,"total", size, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "total"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)

# order labels
gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

gap_clo_high <- gap_clo

saveRDS(gap_clo_high, "data/processed/sensitivity/growth_thres_sensitivity/clo_analysis_ready_high.rds")

################################################################################

# ---- load gap closure information of original mmu400 height threshold 5m gap layer and merge

gap_clo_original<-  readRDS("data/processed/sensitivity/mmu400_height5/clo_analysis_ready.rds")


# --- assign thresholds and prepare for plotting 

gap_clo_low$growth_thres <- as.factor(0.4) 
gap_clo_original$growth_thres <- as.factor(0.5)
gap_clo_high$growth_thres <- as.factor(0.6)

gap_clo<- rbind(gap_clo_low, gap_clo_original, gap_clo_high)



# prepare plotting


# --- closure rates as boxplots

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 26),
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=24),
  legend.text = element_text(size=24),
  strip.text.x = element_text(size = 24),
  strip.text.y = element_text(size = 24),
  legend.position="top") #bottom

# append "m" to the factor levels of growth_thres
gap_clo$growth_thres <- paste0(gap_clo$growth_thres, " m")

tiff("data/results/sensitivity_analysis/height_growth_thres/gap_closure_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  
  scale_fill_brewer(palette="Dark2", name = "Closure mechanism") + My_Theme +
  labs(x = "gap size ( ha )", y= "% of gap area closing annually", colour= "Height threshold") +
  facet_wrap(~growth_thres)
dev.off()


# Calculate the share of lateral closure on the sum of lateral and vertical closure per growth threshold

summary <- gap_clo %>%
  group_by(growth_thres) %>%
  summarise(avg_lateral_closure_share = mean(clo_share_annual [closure_mechanism == "lateral closure"]),
            avg_vertical_closure_share = mean(clo_share_annual [closure_mechanism == "vertical closure"]),
            share_lateral_closure = round(avg_lateral_closure_share / (avg_lateral_closure_share + avg_vertical_closure_share),2))

