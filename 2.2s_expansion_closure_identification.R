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

#gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif")
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

# terra::writeRaster(exp_clo_917, "processed/sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
# terra::writeRaster(exp_clo_1721, "processed/sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")

terra::writeRaster(exp_clo_917, "processed/sensitivity/version.mmu400/exp_clo_917_cn2cr2_mmu400n8_filtered.tif")
terra::writeRaster(exp_clo_1721, "processed/sensitivity/version.mmu400/exp_clo_1721_cn2cr2_mmu400n8_filtered.tif")


# --- extract vegeation growth in gap closure areas per time step ---

#exp_clo_917 <- rast("processed/sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
#exp_clo_1721 <- rast("processed/sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")

chm9 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/chm9_artifacts_masked.tif")
chm17 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/chm21_artifacts_masked.tif")

chm9 <- crop(chm9, chm21)
chm17 <- crop(chm17, chm21)

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

# writeRaster(clo_growth_917 , "processed/sensitivity/closure_area_growth_917.tif")
# writeRaster(clo_growth_1721 , "processed/sensitivity/closure_area_growth_1721.tif")

writeRaster(clo_growth_917 , "processed/sensitivity/version.mmu400/closure_area_growth_917.tif")
writeRaster(clo_growth_1721 , "processed/sensitivity/version.mmu400/closure_area_growth_1721.tif")
