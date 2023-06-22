library(terra)

setwd("C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/")

#load extent of semsitivity analysis area
sensitivity.a <- vect("raw/sensitivity_analysis_area.gpkg")


# load and subset original gap layer to sensitivity area

# mmu 400m^2

gaps9 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps17 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps21 <- rast("processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop gaps to sensitivity area for stacking
gaps9 <- crop(gaps9, sensitivity.a, snap="near", mask=T)
gaps17 <- crop(gaps17, sensitivity.a, snap="near", mask=T)
gaps21 <- crop(gaps21, sensitivity.a, snap="near", mask=T)

gap.stack <- c(gaps9, gaps17,  gaps21)
names(gap.stack) <- c(9,17,21) # observation years

# prepare loop to extract IDs for sensitivity area

sensitivity.ids <- data.frame(matrix(ncol = 2, nrow = 0)) 
names(sensitivity.ids) <- c("id", "year")

numbers <- c(1,2,3)

for (r in numbers) {
  
  print(r)  
  gaps <- gap.stack[[r]] # subset to gap layer of respective year
    
    ids <- extract(gaps, sensitivity.a) # extract raster values
    ids <- na.omit(ids)
    names(ids) <- c("x", "id")
    ids <-unique(ids$id) # get unique ID values
    
    year <- rep(as.numeric(names(gaps)), length(ids)) # create array of respective obs year
    ids <- cbind(ids, year)

  sensitivity.ids<- rbind(sensitivity.ids, ids) # fuse with main df

}

saveRDS(sensitivity.ids, "processed/sensitivity/origID_mmu400_sensitivityAoi.rds")

