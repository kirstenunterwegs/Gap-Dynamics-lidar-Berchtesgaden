library(terra)

#load extent of semsitivity analysis area
sensitivity.a <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/raw/sensitivity_analysis_area.gpkg")


# load and crop original CHM for creating new gap layers
wd <- "F:/Projects/CanopyDynamicsBDG/data/CHM_data/"
setwd(wd)

chm9 <- rast("chm9_artifacts_masked.tif")
chm17<- rast("chm17_artifacts_masked.tif")
chm21 <- rast("chm21_artifacts_masked.tif")

chm9.crop <-  crop(chm9, sensitivity.a, snap="near", mask=T)
chm17.crop <- crop(chm17, sensitivity.a, snap="near", mask=T)
chm21.crop <-  crop(chm21, sensitivity.a, snap="near", mask=T)

writeRaster(chm9.crop,"F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm9_sub_sensitivity.tif")
writeRaster(chm17.crop,"F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm17_sub_sensitivity.tif")
writeRaster(chm21.crop,"F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm21_sub_sensitivity.tif")

# load and crop sensitivity gap layers
wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/processed/"
setwd(wd)

# mmu 100m^2

gaps9.100 <- rast("gaps_sensitivity/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")
gaps17.100 <- rast("gaps_sensitivity/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")
gaps21.100 <- rast("gaps_sensitivity/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")

gaps9.100.crop <- crop(gaps9.100, sensitivity.a, snap="near", mask=T)
gaps17.100.crop <- crop(gaps17.100, sensitivity.a, snap="near", mask=T)
gaps21.100.crop <- crop(gaps21.100, sensitivity.a, snap="near", mask=T)

gaps.100.crop <- c(gaps9.100.crop, gaps17.100.crop,  gaps21.100.crop)

#load forest time and elevation to crop to comparable research area
foresttype <- rast("environment_features/forest_type2020_1m.tif")
foresttype<- crop(foresttype, sensitivity.a, snap="near", mask=T)
elevation.below1800 <- rast("environment_features/elevation_below1800_200steps.tif")
elevation.below1800 <- crop(elevation.below1800, sensitivity.a, snap="near", mask=T)

gaps.100.crop <- mask(gaps.100.crop, foresttype)
gaps.100.crop <- mask(gaps.100.crop, elevation.below1800)

plot(gaps.100.crop)
names(gaps.100.crop) <- c("gaps.9", "gaps.17", "gaps.21")

writeRaster(gaps.100.crop, "gaps_sensitivity/gap.stack.mmu100.sensitivity.tif")

#t <- rast("gaps_sensitivity/gap.stack.mmu100.sensitivity.tif")
