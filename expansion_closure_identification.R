#######################################
# identifying expansion and closure uncertainties
######################################


library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(ForestGapR)
library(ForestTools)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

#--- load layers ---

gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif")
gaps2009 <- gap_stack[[1]]
gaps2017<- gap_stack[[2]]
gaps2021<- gap_stack[[3]]

#--- define functions ---

# simplified function to identify gap expansion and gap closure

gap_change_class <- function(gap_layer1, gap_layer2){
  exp_clo <- rast() #create empty raster to classify
  ext(exp_clo) <- ext(gap_layer1)
  res(exp_clo) <- res(gap_layer1)
  crs(exp_clo) <- crs(gap_layer1)
  # classify change group
  exp_clo[gap_layer1 >0  & is.na(gap_layer2) ] <- 1 #gap closure
  exp_clo[is.na(gap_layer1) & gap_layer2 >0] <- 2 #gap expansion
  return(exp_clo)
} 


#gaps2017 <- crop(gaps2017, gaps2021, snap="near",mask=TRUE)

exp_clo_917 <- gap_change_class(gaps2009, gaps2017)
exp_clo_1721 <- gap_change_class(gaps2017, gaps2021)

terra::writeRaster(exp_clo_917, "processed/sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
terra::writeRaster(exp_clo_1721, "processed/sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")


# --- check on classification uncertainties

# chm_stack <- rast("chm_berchtesgaden_stack_1m.tif")
# chm9 <- chm_stack[[1]]
# chm17<- chm_stack[[2]]
# chm21<- chm_stack[[3]]

#calculate differences

# diff_1721 <- chm21-chm17

# exp_clo_1721 <- rast("exp_clo_1721_cn2cr2_mmu400n8_filtered.tif")
# 
# # subset expansion and closure areas
# exp_1721 <- subst(exp_clo_1721, 1, NA)
# clo_1721 <- subst(exp_clo_1721, 2, NA)
# 
# 
# # clump them (assign ID)
# exp_1721_patches <- patches(exp_1721, directions = 8, allowGaps=FALSE)
# clo_1721_patches <- patches(clo_1721, directions = 8, allowGaps=FALSE)
# 
# terra::writeRaster(exp_1721_patches, "exp_1721_patches_cn2cr2_mmu400n8_filtered.tif")
# terra::writeRaster(clo_1721_patches, "clo_1721_patches_cn2cr2_mmu400n8_filtered.tif")
# 
# exp_1721_patches <- rast("exp_1721_patches_cn2cr2_mmu400n8_filtered.tif")
# clo_1721_patches <- rast("clo_1721_patches_cn2cr2_mmu400n8_filtered.tif")
# 
# #--- define function to mark unclear change signals basing on difference between CHMS ---
# 
# check_exp<- function(change_layer, chm_diff) { 
#   t <- Sys.time()
#   exp_polygon <- as.polygons(change_layer)
#   print("convert gaps to polygon: "); print(Sys.time()-t); t <- Sys.time()
#   change_diff <-  terra::extract(chm_diff, exp_polygon)
#   names(change_diff)[2] <- "chm_diff"
#   print("extract chm in buffer: "); print(Sys.time()-t); t <- Sys.time()
#   change_avg <- change_diff %>%                   #extract average change per expansion area
#     group_by(ID) %>% 
#     summarize(avg_change = mean(chm_diff, na.rm=TRUE))
#   print("extract average change: "); print(Sys.time()-t); t <- Sys.time()
#   change_avg$replace <- ifelse(change_avg$avg_change > 0, 21, 2) #if avg change is growth, than 21 for unclear signal (gap exp, but dom. veg growth), else exp
#   change_layer <- subst(change_layer, from=change_avg$ID, to=change_avg$replace) #replace IDs with either unclear or exp signal
#   print("mark potential false expanison: "); print(Sys.time()-t); t <- Sys.time()
#   return(change_layer)
# }
# 
# # change_layer <- clo_1721_patches
# # chm_diff <- diff_1721
# check_clo<- function(change_layer, chm_diff) { 
#   t <- Sys.time()
#   clo_polygon <- as.polygons(change_layer)
#   print("convert gaps to polygon: "); print(Sys.time()-t); t <- Sys.time()
#   change_diff <-  terra::extract(chm_diff, clo_polygon)
#   names(change_diff)[2] <- "chm_diff"
#   print("extract chm in buffer: "); print(Sys.time()-t); t <- Sys.time()
#   change_avg <- change_diff %>%                   #extract average change per expansion area
#     group_by(ID) %>% 
#     summarize(avg_change = mean(chm_diff, na.rm=TRUE))
#   print("extract average change: "); print(Sys.time()-t); t <- Sys.time()
#   change_avg$replace <- ifelse(change_avg$avg_change < 0, 12, 1) #if avg change is decline, than 12 for unclear signal (gap clo, but dom. veg declone), else clo
#   change_layer <- subst(change_layer, from=change_avg$ID, to=change_avg$replace) #replace IDs with either unclear or clo signal
#   print("mark potential false expanison: "); print(Sys.time()-t); t <- Sys.time()
#   return(change_layer)
# }
# 
# # --- mark unclear expanison and closure areas ---
# 
# clo_1721_check <- check_clo(clo_1721_patches, diff_1721)
# terra::writeRaster(clo_1721_check, "clo_1721_unclearcheck_cn2cr2_mmu400n8_filtered.tif")
# 
# exp_1721_check <- check_exp(exp_1721_patches, diff_1721)
# terra::writeRaster(exp_1721_check, "exp_1721_unclearcheck_cn2cr2_mmu400n8_filtered.tif")

