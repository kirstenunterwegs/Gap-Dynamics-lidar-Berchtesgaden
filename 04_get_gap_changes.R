#########################################################
#
# Retrieving gap change areas from consecutive gap layers
#
#########################################################

# --- libaries ---

library(terra)

# --- set working directory ---

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)


# --- load gap layers ----

gaps2009 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")


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

getGapChanges <- function(gap_1, gap2){ #gap 1 = t1 gaps, gap2 = t2 gaps
  change1 <- GapChangeDecTerra(gap_layer1 = gap_1, gap_layer2 = gap2)
  change2 <- GapChangeDecTerra(gap_layer1 = gap2, gap_layer2 = gap_1) #change detection twice to get gains and losses
  gap_changes <- terra::merge(change1, change2)
  return(gap_changes)
}

# --- get gap changes ----

gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE) 

gap_change_917 <- getGapChanges(gaps2009, gaps2017)
gap_change_1721 <- getGapChanges(gaps2017, gaps2021)

# --- write gap changes to file ---

terra::writeRaster(gap_change_917, "processed/gap_change/gap_change_917_cn2cr2_mmu400n8.tif",overwrite=TRUE)
terra::writeRaster(gap_change_1721, "processed/gap_change/gap_change_1721_cn2cr2_mmu400n8.tif",overwrite=TRUE)
