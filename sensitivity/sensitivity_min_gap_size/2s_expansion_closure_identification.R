#######################################
# identifying expansion and closure uncertainties
######################################

# --- load libaries

library(dplyr)
library(tidyr)
library(terra)


# --- load gap layers with MMU 100 m^2 ---


gaps2009 <- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")
gaps2017<- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")
gaps2021<- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")

# --- define functions ---

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

terra::writeRaster(exp_clo_917, "data/processed/sensitivity/mmu_sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
terra::writeRaster(exp_clo_1721, "data/processed/sensitivity/mmu_sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")


# --- extract vegeation growth in gap closure areas per time step ---

exp_clo_917 <- rast("data/processed/sensitivity/mmu_sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
exp_clo_1721 <- rast("data/processed/sensitivity/mmu_sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")

chm9 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm9_sub_sensitivity.tif")
chm17 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm17_sub_sensitivity.tif")
chm21 <- rast("data/processed/gaps_sensitivity/CHM_sensitivity_area/chm21_sub_sensitivity.tif")

# get vegetation changes
diff917 <- chm17-chm9
diff1721 <- chm21 - chm17

# extract only closure areas
clo_917 <- classify(exp_clo_917, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas
clo_1721 <- classify(exp_clo_1721, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas

# extract vegetation growth in closure areas

diff917 <- crop(diff917,clo_917 )
clo_growth_917 <-mask(diff917, clo_917) 

diff1721 <- crop(diff1721, clo_1721)
clo_growth_1721 <-mask(diff1721, clo_1721) 

writeRaster(clo_growth_917 , "data/processed/sensitivity/mmu_sensitivity/closure_area_growth_917.tif")
writeRaster(clo_growth_1721 , "data/processed/sensitivity/mmu_sensitivity/closure_area_growth_1721.tif")
