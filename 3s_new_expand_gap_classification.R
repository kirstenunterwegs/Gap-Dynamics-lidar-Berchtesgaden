###########################################
# code to classify gaps as new or expanding gaps
# author: Kirsten Krüger
# affiliation: EDFM - Technische Universität München
####
# code strcuture: assign number 1 to every gap pixel 
#                 add up gap rasters for each observation year
#                 reclassify gaps as new or expanding with focal or boundary functions
###########################################

library(terra)
library(raster)
library(doParallel)
library(dplyr)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/processed/"
setwd(wd)

# --- load Gap layers ----

gap_stack <- rast("gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gaps2009 <- gap_stack[[1]]
gaps2017<- gap_stack[[2]]
gaps2021<- gap_stack[[3]]


# stack gaps and identify amount of gap Ids in new gap ------------------------------

#NP whole
gaps.df <- c(gaps2009, gaps2017, gaps2021)
gaps.df <- as.data.frame(gaps.df, na.rm=FALSE)
names(gaps.df) <- c("gaps2009", "gaps2017", "gaps2021")
gaps.df <- gaps.df[rowSums(is.na(gaps.df)) != ncol(gaps.df), ] # delete pixels without any gap at any moment in time
gaps.df[gaps.df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
gaps.df[is.na(gaps.df)] <- 0

#for subset -------------------------------------------------------------------------

# #crop gap layer to sub area
# subarea <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/new_expanding_gap_classification/new_expanding_gao_classification.gpkg")
# gaps2009<- mask(crop(gaps2009, subarea), subarea)
# gaps2017<- mask(crop(gaps2017, subarea), subarea)
# gaps2021<- mask(crop(gaps2021, subarea), subarea)
# #NP sub
# gaps.df <- c(gaps2009, gaps2017, gaps2021)
# gaps.df <- as.data.frame(gaps.df, na.rm=FALSE)
# names(gaps.df) <- c("gaps2009", "gaps2017", "gaps2021")
# gaps.df <- gaps.df[rowSums(is.na(gaps.df)) != ncol(gaps.df), ] # delete pixels without any gap at any moment in time
# gaps.df[gaps.df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
# gaps.df[is.na(gaps.df)] <- 0 
# 
# #i<- 39606
# class_df <- as.data.frame(ncol(2), nrow(0))
# for (i in unique(gaps.df$gaps2017)) {
#   
#   df <- subset(gaps.df, gaps2017 == i) 
#   gapID_t1 <- unique(df$gaps2009)
#   
#   if (length(gapID_t1) == 1 & 0 %in% gapID_t1) {gap_info <- c(i,0)} # new gap
#   if (length(gapID_t1) == 1 & !(0 %in% gapID_t1)) {gap_info <- c(i,1)} # gap same or sub of existing gap
#   if (length(gapID_t1) == 2 & 0 %in% gapID_t1) {gap_info <- c(i,2)} # extending gap
#   if (length(gapID_t1) == 2 & !(0 %in% gapID_t1)) {gap_info <- c(i,3)} # gap sub of several previously existing gaps
#   if (length(gapID_t1) >2 & !(0 %in% gapID_t1) ) {gap_info <- c(i,4)} # extension of only existing gaps (connecting several gaps)
#   if (length(gapID_t1) >2 &  0 %in% gapID_t1 ) {gap_info <- c(i,5)} # extension of more than 1 gap (connecting several gaps)
#   class_df <- rbind(class_df, gap_info)
#   names(class_df) <- c("gap_id", "class")
# }
# 
# gap_class<- class_df
# gap_class$timestep <- as.factor("9-17")
# 
# gap_class %>% group_by(class) %>% 
#   summarise(n = n())
# 
# #subset raster layer to only new gaps
# ID_vector_extendgap <- gap_class$gap_id[gap_class$class > 0] # get IDs of only extending/stable gaps
# ID_vector_replace <- rep(NA, length(ID_vector_newgap)) # create replace vector for extending gap
# 
# 
# rclmat <- cbind(ID_vector_newgap,ID_vector_replace) # amount of 0 as in ID_vector_newgap
# gaps2017_new <- classify(gaps2017, rclmat, include.lowest=TRUE)



### identify new and extending gaps
#cluster approach for whole NP

# 2009-2017 ---------------------------------------------------------------------

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df$gaps2017), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df, gaps2017 == i) 
  gapID_t1 <- unique(df$gaps2009)
  
  if (length(gapID_t1) == 1 & 0 %in% gapID_t1) {gap_info <- c(i,0)} # new gap
  if (length(gapID_t1) == 1 & !(0 %in% gapID_t1)) {gap_info <- c(i,1)} # gap same or sub of existing gap
  if (length(gapID_t1) == 2 & 0 %in% gapID_t1) {gap_info <- c(i,2)} # expanding gap
 # if (length(gapID_t1) == 2 & !(0 %in% gapID_t1)) {gap_info <- c(i,3)} # gap sub of several previously existing gaps (connecting several gaps without creating new gap area) > impossible
 # if (length(gapID_t1) >2 & !(0 %in% gapID_t1) ) {gap_info <- c(i,4)} # extension of more than 1 gap (connecting several gaps without creating new gap area) > impossible
  if (length(gapID_t1) >2 &  0 %in% gapID_t1 ) {gap_info <- c(i,5)} # extension of more than 1 gap (connecting several gaps)
  gap_info
  
}
stopCluster(cl)
class_df <- as.data.frame(class_df)

saveRDS(class_df, "sensitivity/new_exp_gap_class_917_sensitivity.rds")
class_df<- readRDS("sensitivity/new_exp_gap_class_917_sensitivity.rds")

class_df_917 <- class_df
names(class_df) <- c("gap_id", "class") #rename columns
gap_class_summary <- class_df %>% group_by(class) %>%  #get number of new , stable and expanding gaps
                        summarise(n = n()) %>%
                        mutate(perc = round(n/sum(n),2))


#subset raster layer to only new gaps
ID_vector_stablegap <- class_df$gap_id[class_df$class == 1] # get IDs of only stable/shrinking gaps
ID_vector_extendgap <- class_df$gap_id[class_df$class > 1] # get IDs of only extending gaps
ID_vector_newgap <- class_df$gap_id[class_df$class == 0] # get IDs of only new gaps

ID_vector_replace_stable <- rep(2, length(ID_vector_stablegap)) # create replace vector for extending gap
ID_vector_replace_extended <- rep(1, length(ID_vector_extendgap)) # create replace vector for extending gap
ID_vector_replace_new <- rep(0, length(ID_vector_newgap)) # create replace vector for extending gap

rclmat1 <- cbind(ID_vector_extendgap,ID_vector_replace_extended) # create reclassification matrix 
rclamat2 <- cbind(ID_vector_newgap,ID_vector_replace_new )
rclamat3 <- cbind(ID_vector_stablegap, ID_vector_replace_stable)
rclmat <- rbind(rclmat1, rclamat2, rclamat3)

gaps2017_class<- classify(gaps2017, rclmat, include.lowest=TRUE)
writeRaster(gaps2017_class, "sensitivity/gaps2017_new_extended_stable_sensitivity.tif")


# 2017-2021 ---------------------------------------------------------------------

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df$gaps2021), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df, gaps2021 == i) 
  gapID_t1 <- unique(df$gaps2017)
  
  if (length(gapID_t1) == 1 & 0 %in% gapID_t1) {gap_info <- c(i,0)} # new gap
  if (length(gapID_t1) == 1 & !(0 %in% gapID_t1)) {gap_info <- c(i,1)} # gap same or sub of existing gap
  if (length(gapID_t1) == 2 & 0 %in% gapID_t1) {gap_info <- c(i,2)} # extending gap
  if (length(gapID_t1) == 2 & !(0 %in% gapID_t1)) {gap_info <- c(i,3)} # gap sub of several previously existing gaps
  # if (length(gapID_t1) >2 & !(0 %in% gapID_t1) ) {gap_info <- c(i,4)} # extension of more than 1 gap (connecting several gaps without creating new gap area) > impossible
  if (length(gapID_t1) >2 &  0 %in% gapID_t1 ) {gap_info <- c(i,5)} # extension of more than 1 gap (connecting several gaps)
  gap_info
  
}
stopCluster(cl)
class_df <- as.data.frame(class_df)

saveRDS(class_df, "sensitivity/new_exp_gap_class_1721_sensitivity.rds")
class_df <- readRDS( "sensitivity/new_exp_gap_class_1721_sensitivity.rds")

class_df_1721 <- class_df
names(class_df) <- c("gap_id", "class") #rename columns
gap_class_summary <- class_df %>% group_by(class) %>%  #get number of new , stable and expanding gaps
  summarise(n = n()) %>%
  mutate(perc = round(n/sum(n),2))


#subset raster layer to only new gaps
ID_vector_stablegap <- class_df$gap_id[class_df$class == 1] # get IDs of only stable/shrinking gaps
ID_vector_extendgap <- class_df$gap_id[class_df$class > 1] # get IDs of only extending gaps
ID_vector_newgap <- class_df$gap_id[class_df$class == 0] # get IDs of only new gaps

ID_vector_replace_stable <- rep(2, length(ID_vector_stablegap)) # create replace vector for extending gap
ID_vector_replace_extended <- rep(1, length(ID_vector_extendgap)) # create replace vector for extending gap
ID_vector_replace_new <- rep(0, length(ID_vector_newgap)) # create replace vector for extending gap

rclmat1 <- cbind(ID_vector_extendgap,ID_vector_replace_extended) # create reclassification matrix 
rclamat2 <- cbind(ID_vector_newgap,ID_vector_replace_new )
rclamat3 <- cbind(ID_vector_stablegap, ID_vector_replace_stable)
rclmat <- rbind(rclmat1, rclamat2, rclamat3)

gaps2021_class<- classify(gaps2021, rclmat, include.lowest=TRUE)
writeRaster(gaps2021_class, "sensitivity/gaps2021_new_extended_stable_sensitivity.tif")



