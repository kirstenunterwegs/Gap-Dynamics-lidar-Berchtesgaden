#########################################################
# Retrievig gap change areas from consecutive gap layers
#########################################################

#libaries
library(terra)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)


# --- load gap layers ----

#gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif")
gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu400.sensitivity.tif")
gaps2009 <- gap_stack[[1]]
gaps2017<- gap_stack[[2]]
gaps2021<- gap_stack[[3]]

# gaps2009 <- rast("processed/gaps_sensitivity/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")
# gaps2017 <- rast("processed/gaps_sensitivity/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")
# gaps2021 <- rast("processed/gaps_sensitivity/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu100n8_filtered_woheight.tif")

# --- function to retrieve gap changes ---

GapChangeDecTerra <- function (gap_layer1, gap_layer2)  # adapted from ForestRGap GapChangeDec function
{
  gap_layer1[!is.na(gap_layer1)] <- 1
  gap_layer2[!is.na(gap_layer2)] <- 2
  gap_layer1[is.na(gap_layer1)] <- 0
  gap_layer2[is.na(gap_layer2)] <- 0
  gap_diff <- gap_layer2 - gap_layer1
  gap_diff[gap_diff != 2] <- NA
  return(gap_diff)
}

getGapChanges <- function(gap1, gap2){ #gap 1 = t1 gaps, gap2 = t2 gaps
  change1 <- GapChangeDecTerra(gap_layer1 = gap1, gap_layer2 = gap2)
  change2 <- GapChangeDecTerra(gap_layer1 = gap2, gap_layer2 = gap1) #change detection twice to get gains and losses
  gap_changes <- terra::merge(change1, change2)
  return(gap_changes)
}

# --- get gap changes ----
#ensure same extent 
#ext(gaps2009) <- ext(gaps2021)
#ext(gaps2017) <- ext(gaps2021)

gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE) 

gap_change_917 <- getGapChanges(gaps2009, gaps2017)
gap_change_1721 <- getGapChanges(gaps2017, gaps2021)


# terra::writeRaster(gap_change_917, "processed/sensitivity/gap_change_917_cn2cr2_mmu100n8.tif",overwrite=TRUE)
# terra::writeRaster(gap_change_1721, "processed/sensitivity/gap_change_1721_cn2cr2_mmu100n8.tif",overwrite=TRUE)

terra::writeRaster(gap_change_917, "processed/sensitivity/version.mmu400/gap_change_917_cn2cr2_mmu400n8.tif",overwrite=TRUE)
terra::writeRaster(gap_change_1721, "processed/sensitivity/version.mmu400/gap_change_1721_cn2cr2_mmu400n8.tif",overwrite=TRUE)

