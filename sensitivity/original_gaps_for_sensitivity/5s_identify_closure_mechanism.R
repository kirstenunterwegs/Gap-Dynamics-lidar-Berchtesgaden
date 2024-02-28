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


# ----- classify vertical and horizontal closure ----

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


terra::writeRaster(gap_closure_mechanism917, "data/processed/sensitivity/mmu400_height5/gap_closure_mechanism917.tif", overwrite=TRUE) 
terra::writeRaster(gap_closure_mechanism1721, "data/processed/sensitivity/mmu400_height5/gap_closure_mechanism1721.tif", overwrite=TRUE)



# --- 2009 - 2017 ---


#merge closure mechanism with gaps
gap_closure_mechanism917 <- rast( "data/processed/sensitivity/mmu400_height5/gap_closure_mechanism917.tif")

gaps2009 <- crop(gaps2009, gap_closure_mechanism917)
gap_closure_mechanism_stack <- c(gap_closure_mechanism917, gaps2009)


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
sum(gap_clo_per_id$contraction) # 0 gaps do not experience any closure from 2009-2017 (previously 7, maybe the masking cut off some gap areas)

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


# --- 2017- 2021 ---


# #merge closure mechanism with gaps
gap_closure_mechanism1721 <- rast("data/processed/sensitivity/mmu400_height5/gap_closure_mechanism1721.tif")

gaps2017 <- crop(gaps2017, gap_closure_mechanism1721)
gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721, gaps2017)


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
#if I want to include the no closure pixels, I have to disable following line: !!!!
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


### -------------------------------merge 9-17 and 17-21 dfs for comparison ----------------------- ###

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



# --- append Total closure info to main df for distribution display 
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
  k <- c(i,"Total", size,sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "Total"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)


#order labels
gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

saveRDS(gap_clo, "data/processed/sensitivity/mmu400_height5/clo_analysis_ready.rds")

