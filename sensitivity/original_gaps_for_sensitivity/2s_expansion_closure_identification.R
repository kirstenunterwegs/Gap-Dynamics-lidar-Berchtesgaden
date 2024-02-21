#######################################
# identifying expansion and closure 
######################################

library(dplyr)
library(tidyr)
library(terra)
library(ForestGapR)
library(ForestTools)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

#--- load layers ---

gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu400.sensitivity.tif")
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

terra::writeRaster(exp_clo_917, "processed/sensitivity/mmu400_height5/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
terra::writeRaster(exp_clo_1721, "processed/sensitivity/mmu400_height5/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")


# --- extract vegeation growth in gap closure areas per time step ---

chm9 <- rast("processed/gaps_sensitivity/CHM_sensitivity_area/chm9_sub_sensitivity.tif")
chm17 <- rast("processed/gaps_sensitivity/CHM_sensitivity_area/chm17_sub_sensitivity.tif")
chm21 <- rast("processed/gaps_sensitivity/CHM_sensitivity_area/chm21_sub_sensitivity.tif")

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

writeRaster(clo_growth_917 , "processed/sensitivity/mmu400_height5/closure_area_growth_917.tif")
writeRaster(clo_growth_1721 , "processed/sensitivity/mmu400_height5/closure_area_growth_1721.tif")

