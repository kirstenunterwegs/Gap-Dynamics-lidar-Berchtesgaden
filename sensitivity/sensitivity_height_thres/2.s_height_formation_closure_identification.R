#######################################
# identifying expansion and closure uncertainties
######################################

library(dplyr)
library(tidyr)
library(terra)
library(ForestGapR)
library(ForestTools)


# --- load gap layers ----

gaps2009.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2017.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2021.3 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")

gaps2009.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2017.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2021.10 <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")


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

gap_change_917.3 <- gap_change_class(gaps2009.3, gaps2017.3)
gap_change_1721.3 <- gap_change_class(gaps2017.3, gaps2021.3)


gap_change_917.10 <- gap_change_class(gaps2009.10, gaps2017.10)
gap_change_1721.10 <- gap_change_class(gaps2017.10, gaps2021.10)


terra::writeRaster(gap_change_917.3, "data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_917.tif")
terra::writeRaster(gap_change_1721.3, "data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_1721.tif")

terra::writeRaster(gap_change_917.10, "data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_917.tif")
terra::writeRaster(gap_change_1721.10, "data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_1721.tif")


# --- extract vegetation growth in gap closure areas per time step ---

exp_clo_917_h3 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_917.tif")
exp_clo_1721_h3 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_1721.tif")

exp_clo_917_h10 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_917.tif")
exp_clo_1721_h10 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_1721.tif")


chm9 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm9_sub_sensitivity.tif")
chm17 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm17_sub_sensitivity.tif")
chm21 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm21_sub_sensitivity.tif")

# get vegetation changes
diff917 <- chm17-chm9
diff1721 <- chm21 - chm17

# extract only closure areas

clo_917_h3 <- classify(exp_clo_917_h3, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas
clo_1721_h3 <- classify(exp_clo_1721_h3, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas

clo_917_h10 <- classify(exp_clo_917_h10, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas
clo_1721_h10 <- classify(exp_clo_1721_h10, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas

# --- extract vegetation growth in closure areas

# height thres 3m

diff917_h3 <- crop(diff917,clo_917_h3 )
clo_growth_917_h3 <-mask(diff917_h3, clo_917_h3) 

diff1721_h3 <- crop(diff1721, clo_1721_h3)
clo_growth_1721_h3 <-mask(diff1721_h3, clo_1721_h3) 

# height thres 10m

diff917_h10 <- crop(diff917,clo_917_h10 )
clo_growth_917_h10 <-mask(diff917_h10, clo_917_h10) 

diff1721_h10 <- crop(diff1721, clo_1721_h10)
clo_growth_1721_h10 <-mask(diff1721_h10, clo_1721_h10) 


writeRaster(clo_growth_917_h3 , "data/processed/sensitivity/height_sensitivity/closure_area_growth_917_h3.tif")
writeRaster(clo_growth_1721_h3 , "data/processed/sensitivity/height_sensitivity/closure_area_growth_1721_h3.tif")

writeRaster(clo_growth_917_h10 , "data/processed/sensitivity/height_sensitivity/closure_area_growth_917_h10.tif")
writeRaster(clo_growth_1721_h10 , "data/processed/sensitivity/height_sensitivity/closure_area_growth_1721_h10.tif")
