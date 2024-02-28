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


# --- load gap layers ----

gaps2009.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2017.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2021.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")

gaps2009.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2017.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2021.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")



# stack gaps 

# 3 m height thres
gaps.df.3 <- c(gaps2009.3, gaps2017.3, gaps2021.3)
gaps.df.3 <- as.data.frame(gaps.df.3, na.rm=FALSE)
names(gaps.df.3) <- c("gaps2009", "gaps2017", "gaps2021")
gaps.df.3 <- gaps.df.3[rowSums(is.na(gaps.df.3)) != ncol(gaps.df.3), ] # delete pixels without any gap at any moment in time
gaps.df.3[gaps.df.3 == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
gaps.df.3[is.na(gaps.df.3)] <- 0

# 10 m height thres
gaps.df.10 <- c(gaps2009.10, gaps2017.10, gaps2021.10)
gaps.df.10 <- as.data.frame(gaps.df.10, na.rm=FALSE)
names(gaps.df.10) <- c("gaps2009", "gaps2017", "gaps2021")
gaps.df.10 <- gaps.df.10[rowSums(is.na(gaps.df.10)) != ncol(gaps.df.10), ] # delete pixels without any gap at any moment in time
gaps.df.10[gaps.df.10 == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
gaps.df.10[is.na(gaps.df.10)] <- 0


# --- identify new and extending gaps ---

# cluster approach for whole NP

# --- 2009-2017 --- height threshold 3m ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df.3$gaps2017), .combine = rbind) %dopar% { 
  
  df <- subset(gaps.df.3, gaps2017 == i) 
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

gaps2017.3_class<- classify(gaps2017.3, rclmat, include.lowest=TRUE)
writeRaster(gaps2017.3_class, "data/processed/sensitivity/height_sensitivity/gaps2017_new_extended_stable_h3.tif")


# --- 2009-2017 --- height threshold 10m ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df.10$gaps2017), .combine = rbind) %dopar% { 
  
  df <- subset(gaps.df.10, gaps2017 == i) 
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

gaps2017.10_class<- classify(gaps2017.10, rclmat, include.lowest=TRUE)
writeRaster(gaps2017.10_class, "data/processed/sensitivity/height_sensitivity/gaps2017_new_extended_stable_h10.tif")




# --- 2017-2021 --- height threshold 3 m ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df.3$gaps2021), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df.3, gaps2021 == i) 
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

gaps2021.3_class<- classify(gaps2021.3, rclmat, include.lowest=TRUE)
writeRaster(gaps2021.3_class, "data/processed/sensitivity/height_sensitivity/gaps2021_new_extended_stable_h3.tif")



# --- 2017-2021 --- height threshold 10 m ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df.10$gaps2021), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df.10, gaps2021 == i) 
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

gaps2021.10_class<- classify(gaps2021.10, rclmat, include.lowest=TRUE)
writeRaster(gaps2021.10_class, "data/processed/sensitivity/height_sensitivity/gaps2021_new_extended_stable_h10.tif")


