######################################
# Analysing Gaps
#####################################

library(terra)
library(sf)
library(plyr)
library(dplyr)


wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/processed/"
setwd(wd)

# --- load Gap layers ----

gap_stack <- rast("gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gaps2009 <- gap_stack[[1]]
gaps2017<- gap_stack[[2]]
gaps2021<- gap_stack[[3]]


# --- load NP information 

foresttype <- rast("environment_features/forest_type2020_1m.tif")
aspect<-  rast("environment_features/aspect_2021_classified_1m.tif")
elevation <- rast("environment_features/elevation_below1800_200steps.tif")

#crop environmental features to sensitivity gap area
foresttype <- crop(foresttype, gaps2009)
aspect <- crop(aspect, gaps2009)
elevation <- crop(elevation, gaps2009)


# --- define functions -----


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

#function to assign elevation class
getElevation <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, elevation) %>% #count pixels per elevation class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
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
  xx<- xx %>% mutate(elevation = as.factor(recode(elevation,
                                                  `1`="600-800",
                                                  `2`="800-1000",
                                                  `3`="1000-1200",
                                                  `4`="1200-1400",
                                                  `5`="1400-1600",
                                                  `6`="1600-1800",
                                                  `7`="1800-2000")))
  return(xx)
}

#function to assign aspect class
getAspect <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, aspect) %>% #count pixels per aspect class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one aspect type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several aspect classes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no aspect info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other aspect info assign that one to ID
    }
  }
  xx<- xx %>% mutate(aspect = as.factor(recode(aspect,
                                               `1`="North",
                                               `2`="East",
                                               `3`="South",
                                               `4`="West")))
  return(xx)
}



Gap_Stats_ftype <- function (gap_layer, foresttype, elevation, aspect, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list <- gap_list[!is.na(gap_list[, 1]), ]

print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, foresttype, elevation, aspect), na.rm = FALSE)
names(gap_chm) <- c("ID", "ftype", "elevation", "aspect")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()

# identify forest type
gap_list$ftype <- as.data.frame(getForestType(gap_chm))[,"ftype"]

print("get forest type: "); print(Sys.time()-t); t <- Sys.time()

# identify elevation
gap_list$elevation <- as.data.frame(getElevation(gap_chm))[,"elevation"]

print("get elevation: "); print(Sys.time()-t); t <- Sys.time()

# identify aspect
gap_list$aspect <- as.data.frame(getAspect(gap_chm))[,"aspect"]

print("get aspect: "); print(Sys.time()-t); t <- Sys.time()


gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "count", "forest_type", "elevation", "aspect", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#------- calculate gap stats per forest type
stats_2009<- Gap_Stats_ftype(gaps2009, foresttype, elevation, aspect, "2009")
stats_2017 <- Gap_Stats_ftype(gaps2017, foresttype, elevation, aspect, "2017")
stats_2021 <- Gap_Stats_ftype(gaps2021, foresttype, elevation, aspect, "2021")

stats_all <- rbind(stats_2009, stats_2017, stats_2021)

saveRDS(stats_all, "sensitivity/stats_sensitivity.rds")



