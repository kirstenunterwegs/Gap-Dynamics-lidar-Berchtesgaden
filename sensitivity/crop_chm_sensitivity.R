#################################################
#
# Crop CHMS to the area for sensitivity analysis
#
#################################################


library(terra)

#--- crop artifacts masked CHM to subarea for sensitivity analysis

chm9 <- rast("data/processed/CHM_data/chm9_artifacts_masked.tif")
chm17<- rast("data/processed/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("data/processed/CHM_data/chm21_artifacts_masked.tif")

subarea <- vect("data/raw/sensitivity_analysis_area.gpkg")
closed.forest <- vect("data/raw/closed_forest_corezone.gpkg")


chm9.sub <- crop(mask(chm9, subarea), subarea)
chm9.sub.f <- crop(chm9.sub, closed.forest, mask =T)

chm17.sub <- crop(mask(chm17, subarea), subarea)
chm17.sub <- crop(chm17.sub, closed.forest, mask =T)

chm21.sub <- crop(mask(chm21, subarea), subarea)
chm21.sub <- crop(chm21.sub, closed.forest, mask =T)

writeRaster(chm9.sub, "data/processed/gaps_sensitivity/CHM_sensitivity_area/chm9_sub_sensitivity.tif")
writeRaster(chm17.sub, "data/processed/gaps_sensitivity/CHM_sensitivity_area/chm17_sub_sensitivity.tif")
writeRaster(chm21.sub, "data/processed/gaps_sensitivity/CHM_sensitivity_area/chm21_sub_sensitivity.tif")


