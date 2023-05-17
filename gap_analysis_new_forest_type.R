######################################
# Analysing Gaps
#####################################

library(sf)
library(plyr)
library(dplyr)
library(dplyr)
library(terra)
library(ForestGapR)
library(ggplot2)
library(ForestTools)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)
library(landscapemetrics)
require(scales)
library(tidyverse)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load CHM and Gap layers ----


gaps2009 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# --- load NP information 

foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")

# crop layer to same extent

gaps2009 <- crop(gaps2009, foresttype) 
gaps2017 <- crop(gaps2017, foresttype) 
gaps2021 <-crop(gaps2021, foresttype) 


# --- define functions -----

# ---- forest type



#function to assign forest type
getForestType <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, ftype) %>% #count pixels per ftype per gap
    dplyr::summarize(count = n())
  #identify dominating forest type in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(ftype = as.factor(recode(ftype,
                                                    `1`="Beech",
                                                    `2`="Spruce-fir-beech",
                                                    `4`="Spruce",
                                                    `5`="Larch-Swiss stone pine",
                                                    `6`="Dwarf mountain pine")))
  return(xx)
}


Gap_Stats_ftype <- function (gap_layer, foresttype, year) 
{  t <- Sys.time()
  gap_list <- data.frame(terra::freq(gap_layer))
  gap_list <- gap_list[!is.na(gap_list[, 1]), ]

print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, foresttype), na.rm = FALSE)
names(gap_chm) <- c("ID", "ftype")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()

# identify forest type

gap_list$ftype <- as.data.frame(getForestType(gap_chm))[,"ftype"]

print("get forest type: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "count", "forest_type", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#------- calculate gap stats per forest type
stats_ftype_2009<- Gap_Stats_ftype(gaps2009, foresttype, "2009")
stats_ftype_2017 <- Gap_Stats_ftype(gaps2017, foresttype, "2017")
stats_ftype_2021 <- Gap_Stats_ftype(gaps2021, foresttype, "2021")

stats_all_ftype <- rbind(stats_ftype_2009, stats_ftype_2017, stats_ftype_2021)

saveRDS(stats_all_ftype, "processed/gap_features/stats_all_ftype.rds")

