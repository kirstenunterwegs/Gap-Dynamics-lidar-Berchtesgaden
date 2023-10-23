library(terra)

#--- crop artifacts masked CHM to subarea for sensitivity analysis

chm9 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm9_artifacts_masked.tif")
chm17<- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("F:/Projects/CanopyDynamicsBDG/data/CHM_data/CHM_data/chm21_artifacts_masked.tif")

#subarea <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/sensitivity_analysis_subarea.gpkg")
subarea <- vect("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/variability_analysis_subarea.shp")
closed.forest <- vect("F:/Projects/CanopyDynamicsBDG/data/data/Waldmaske/closed_forest.shp")
closed.forest <- project(closed.forest, subarea)

chm9.sub <- crop(mask(chm9, subarea), subarea)
chm9.sub.f <- crop(chm9.sub, closed.forest, mask =T)

chm17.sub <- crop(mask(chm17, subarea), subarea)
chm17.sub <- crop(chm17.sub, closed.forest, mask =T)

chm21.sub <- crop(mask(chm21, subarea), subarea)
chm21.sub <- crop(chm21.sub, closed.forest, mask =T)

writeRaster(chm9.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm9_sub_sensitivity.tif")
writeRaster(chm17.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm17_sub_sensitivity.tif")
writeRaster(chm21.sub, "F:/Projects/CanopyDynamicsBDG/data/sensitivity/CHM_subarea/chm21_sub_sensitivity.tif")

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
