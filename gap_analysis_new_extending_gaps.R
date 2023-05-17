######################################
# Analysing Gaps
#####################################


library(plyr)
library(dplyr)
library(terra)
library(ForestGapR)
library(ggplot2)
require(scales)
library(tidyverse)
#library(sf)
# library(ForestTools)
# library(lattice)
# library(latticeExtra)
# library(rasterVis)
# library(RColorBrewer)
# library(sp)
# library(landscapemetrics)


wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load classified Gap layers of new and expanding gaps ----

# chm_stack <- rast("chm_berchtesgaden_stack_1m_maskedartifacts.tif")
# chm9 <- chm_stack[[1]]
# chm17<- chm_stack[[2]]
# chm21<- chm_stack[[3]]

#chm17 <- rast("i:/Fonda/workspace/berchtesgaden/lidar/2017/berchtesgaden_2017_chm_1m.tif")
# 
gaps2017.id <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021.id <- rast("processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

gaps2017 <- rast("processed/creation/gaps2017_new_extended_stable.tif")
gaps2021 <- rast("processed/creation/gaps2021_new_extended_stable.tif")

# crop gaps ID
gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 

#load expansion and closure layer and extract only expansion areas
exp_clo <- rast("processed/creation//expansion_closure_917_cn2cr2_mmu800.tif")
exp <- classify(exp_clo, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("processed/creation//expansion_closure_1721_cn2cr2_mmu800.tif")
exp1721 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

# --- load NP information 

foresttype <- rast("raw/forest_types2020.tif")
management <- vect("raw/npb_zonierung_22_epsg25832.shp")
aspect<-  rast("processed/environment_features/aspect_2021_classified_1m.tif")
elevation <- rast("processed/environment_features/berchtesgaden_2021_classified_200steps_dtm_1m.tif")
closed.forest <- vect("raw/closed_forest_epsg25832.shp")

# ---- crop layers to research sites:

# exclude sites > 1800 m
m <- c(0,6,1, 6,8,NA) #7 and 8 are areas > 1800 m elevation
rclmat <- matrix(m, ncol=3, byrow=TRUE)
elevation.1800 <- classify(elevation, rclmat, include.lowest=TRUE)
elevation.1800.poly <- as.polygons(elevation.1800, trunc=TRUE, dissolve=TRUE, values=TRUE, #prepare masking vector for elevation with area < 1800m
            na.rm=TRUE, na.all=FALSE, extent=FALSE)

elevation.below1800 <- mask(elevation, elevation.1800.poly)

writeRaster(elevation.below1800, "processed/environment_features/elevation_below1800_200steps.tif")

elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")

#freq(elevation.below1800) #check if I really excluded all pixels >1800 m

# extract core zone to exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

# -- rasterize NP information ---NOT needed here

# #crate empty raster for rasterization
# r_1 <- rast()
# ext(r_1) <- ext(gaps2017)
# terra::res(r_1) <- terra::res(gaps2017)  
# terra::crs(r_1) <- terra::crs(gaps2017)
# 
# # rasterize forest type Information
# ftype <- terra::rasterize(foresttype, r_1, field="type")

# --- reclassify forest types to combine spruce-fir-beech (2) and spruce-fir (3)
foresttype <- subst(foresttype, 3, 2) 
# resample to match resolution
foresttype <- resample(foresttype, elevation.below1800, method="near")

writeRaster(foresttype, "processed/environment_features/forest_type2020_1m.tif")
foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")

# crop all layers to same extent

gaps2017.id <- crop(gaps2017.id, elevation.below1800)
gaps2017<- crop(gaps2017, elevation.below1800)

# --- stack gap information and crop it

#2017

stack2017 <- c(gaps2017.id, gaps2017, exp, foresttype , elevation.below1800, aspect)
names(stack2017) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017 <- mask(stack2017, foresttype)
stack2017 <- mask(stack2017, core.zone)
stack2017 <- mask(stack2017, closed.forest)
stack2017 <- mask(stack2017, elevation.below1800)

writeRaster(stack2017, "processed/creation/updated/stack.2017.all.gap.information.expansion.tif")
gap_stack_2017 <- rast("processed/creation/updated/stack.2017.all.gap.information.expansion.tif")


df <- as.data.frame(gap_stack_2017, na.rm = FALSE) 

df1 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df1 <- df1[!is.na(df1$forest_type),] # exclude all areas with no forest type information
# df1 <- df1[!is.na(df1$elevation),]# exclude all areas > 1800 m (NA in this case, as it was re-coded above) 

write_rds(df1, "processed/creation/updated/stack_2017_new_exp_df.rds")
df1 <- readRDS( "processed/creation/updated/stack_2017_new_exp_df.rds")

#2021

stack21 <- c(gaps2021.id, gaps2021, exp1721, foresttype, elevation.below1800, aspect)
names(stack21) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack21 <- mask(stack21, foresttype)
stack21 <- mask(stack21, core.zone)
stack21 <- mask(stack21, closed.forest)
stack21 <- mask(stack21, elevation.below1800)

writeRaster(stack21, "processed/creation/updated/stack.2021.all.gap.information.expansion.tif", overwrite=T)
gap_stack_2021 <- rast("processed/creation/updated/stack.2021.all.gap.information.expansion.tif")

df <- as.data.frame(gap_stack_2021, na.rm = FALSE) 

df2 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df2 <- df2[!is.na(df2$forest_type),] # exclude all areas with no forest type information
# df2 <- df2[!is.na(df2$elevation),]# exclude all areas > 1800 m (na in this case, as it was recoded above) 

write_rds(df2, "processed/creation/updated/stack_2021_new_exp_df.rds")
df2<- readRDS("processed/creation/updated/stack_2021_new_exp_df.rds")


#-------calculate area shares per category 

gap_stack_2017.closed <- rast("processed/creation/updated/stack.2017.all.gap.information.expansion.tif")

df.area <- as.data.frame(gap_stack_2017.closed, na.rm = FALSE) 

#recode forest type to drop NAs
# df.area <- df.area %>% mutate( ftype = as.numeric(recode(forest_type,
#                                                         'broadleaved' = 1,
#                                                         'coniferous' = 2,
#                                                         'larch-dominant' = 3,
#                                                         'mixed stands' = 4)))

#delete all pixels with no forest type information, as we do not consider these areas + all areas >1800 m 
df.area.nona <- df.area %>% drop_na(elevation)
df.area.nona <- df.area.nona %>% drop_na(forest_type)
df.area.nona <- df.area.nona[df.area.nona$elevation != 7,] #pixels which have been left at boundaries of crop area

head(df.area.nona)  

#keep only necessary environmental feature columns
keeps <- c("elevation","aspect", "forest_type")
df.area.nona = df.area.nona[keeps] 

saveRDS(df.area.nona, "processed/creation/updated/df.area.nona.rds")
df.area.nona <- readRDS("processed/creation/updated/df.area.nona.rds") 


df.area.long <- gather(df.area.nona, category, class)

head(df.area.long)

area_share_class <- df.area.long %>% group_by(category,class) %>% 
  summarize(total_area =round((sum(!is.na(class)))/10000, 2)) # divided by 10.000 to get ha (res 1m)

area_share_class <- area_share_class %>% group_by(category) %>% 
  mutate(total_area_category =round((sum(total_area)), 2),
            class_area_perc = round(total_area/ total_area_category, 4)) 

#exclude NAs
# area_share_class<- area_share_class[!is.na(area_share_class$class),]

# recode category labels
class.name <- c("North", "East", "South", "West", 
                "600-800", "800-1000", "1000-1200", "1200-1400", "1400-1600", "1600-1800", 
               # "broadleaved", "coniferous", "larch-dominant", "mixed stands") #old forest type classification
                "Beech", "Spruce-fir-beech", "Spruce", "Larch-Swiss stone pine", "Dwarf mountain pine")
area_share_class$class.name <- class.name

write_rds(area_share_class, "processed/creation/updated/area_share_per_class.rds")


# --------calculate area shares per elevation - forest type combination !!!!! neu berechnen wenn notwendig!!!
df.area.nona$elev.ftype = paste0(df.area.nona$elevation, df.area.nona$forest_type)

area_share_class.elev.ftype <- df.area.nona %>% group_by(elev.ftype) %>% 
  summarize(total_area =round((sum(!is.na(elev.ftype)))/10000, 2)) # divided by 10.000 to get ha (res 1m)

area_share_class.elev.ftype  <- area_share_class.elev.ftype %>% 
  mutate(total_area_category =round((sum(total_area)), 2),
         class_area_perc = round(total_area/ total_area_category, 4)*100) 

elev.ftype <- c("600-800-Beech", "600-800-Spruce-fir-beech", "600-800-Spruce", "600-800-Larch-Swiss stone pine","600-800-Dwarf mountain pine",
                "800-1000-Beech","800-1000-Spruce-fir-beech","800-1000-Spruce", "800-1000-Larch-Swiss stone pine","800-1000-Dwarf mountain pine",
                "1000-1200-Beech", "1000-1200-Spruce-fir-beech","1000-1200-Spruce","1000-1200-Larch-Swiss stone pine", "1000-1200-Dwarf mountain pine",
                "1200-1400-Beech","1200-1400-Spruce-fir-beech","1200-1400-Spruce","1200-1400-Larch-Swiss stone pine", "1200-1400-Dwarf mountain pine",
                "1400-1600-Beech", "1400-1600-Spruce-fir-beech","1400-1600-Spruce","1400-1600-Larch-Swiss stone pine", "1400-1600-Dwarf mountain pine",
                "1600-1800-Spruce-fir-beech","1600-1800-Spruce","1600-1800-Larch-Swiss stone pine","1600-1800-Dwarf mountain pine" )
area_share_class.elev.ftype$elevation.ftype <- elev.ftype

write_rds(area_share_class.elev.ftype, "processed/creation/updated/area_share_per_class.elev.ftype.rds")

area_share_class.elev.ftype.sub40 <- area_share_class.elev.ftype %>% filter(total_area >40)

write_rds(area_share_class.elev.ftype.sub40, "processed/creation/updated/area_share_per_class.elev.ftype.sub40.rds")

# --------calculate features per gap.id

gap_features_917 <- df1 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

gap_features_1721 <- df2 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))


# functions to identify major forest type, elevation and aspect per gap and hence expansion ------------------

getForestType <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, forest_type) %>% #count pixels per ftype per gap
    summarize(count = n())
  #identify dominating forest type in gap area
  xx <- data_frame()
  for (i in unique(gap_chm$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(forest_type = as.factor(recode(forest_type,
                                                  `1`="Beech",
                                                  `2`="Spruce-fir-beech",
                                                  `4`="Spruce",
                                                  `5`="Larch-Swiss stone pine",
                                                  `6`="Dwarf mountain pine")))
  return(xx)
}

#function to assign elevation class
getElevation <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, elevation) %>% #count pixels per elevation class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
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
                                                  `7`="1800-2000",
                                                  `8`="2000-2800")))
  return(xx)
}


#function to assign elevation class
getAspect <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, aspect) %>% #count pixels per aspect class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
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


#2017

ftype <- getForestType(df1)
elevation <- getElevation(df1)
aspect <- getAspect(df1)

gap_features_917 <- merge(gap_features_917, ftype[,c("gap.id","forest_type")], by = "gap.id", all.x = TRUE)
gap_features_917 <- merge(gap_features_917, elevation[,c("gap.id","elevation")], by = "gap.id", all.x = TRUE)
gap_features_917 <- merge(gap_features_917, aspect[,c("gap.id","aspect")], by = "gap.id", all.x = TRUE)

saveRDS(gap_features_917,"processed/creation/updated/gap_features_new_expanding_917.rds")
gap_features_917 <- readRDS("processed/creation/updated/gap_features_new_expanding_917.rds")


#2021

ftype <- getForestType(df2)
elevation <- getElevation(df2)
aspect <- getAspect(df2)

gap_features_1721 <- merge(gap_features_1721, ftype[,c("gap.id","forest_type")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, elevation[,c("gap.id","elevation")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, aspect[,c("gap.id","aspect")], by = "gap.id", all.x = TRUE)

saveRDS(gap_features_1721,"processed/creation/updated/gap_features_new_expanding_1721.rds")
gap_features_1721 <- readRDS("processed/creation/updated/gap_features_new_expanding_1721.rds")


# analyse new and expanding gaps ----------------------------------------------------------

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 19,angle = 45, hjust=1),
  axis.text.y = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=20),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"))
# legend.position="bottom")

  
#calculate amount of expansion (check with MA) 2017
sum(gap_features_917$exp.area.ha)/8 # 9.48  fits values of MA (9.89 ha/yr, as I excluded areas >1800m)
sum(gap_features_917$new.exp == "new") # 213
sum(gap_features_917$new.exp == "expanding") #3750
sum(gap_features_917$new.exp == "stable") #122

#calculate amount of expansion (check with MA) 2021
sum(gap_features_1721$exp.area.ha)/4 # 51.38197 fits values of MA (48.5ha/yr, as I excluded areas >1800m)
sum(gap_features_1721$new.exp == "new") # 686
sum(gap_features_1721$new.exp == "expanding") #3465
sum(gap_features_1721$new.exp == "stable") #0

gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921 <- rbind(gap_features_917, gap_features_1721)
gap_features921 <- subset(gap_features921, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis
gap_features921 <- gap_features921[gap_features921$elevation != "1800-2000",]
gap_features921 <- gap_features921[gap_features921$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the croping of the reserach area

#calculate annual gap creation rate

# all

gap.creation <- gap_features921 %>% group_by(new.exp, year) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp)%>%
  mutate(avg.gap.creation.annual = round(weighted.mean(gap.creation.annual, time),2))

#--------- area scaling
gap.creation$area.scaling.factor <- 100/4000 #scaling to gap creation per 100 ha (percent!), total reserach area is 4000 ha
gap.creation$avg.gap.creation.annual.scaled <- gap.creation$avg.gap.creation.annual*gap.creation$area.scaling.factor


# forest type and elevation


gap.creation.ftype.elev <- gap_features921 %>% group_by(year, forest_type, elevation) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by( forest_type, elevation)%>%
  mutate(avg.gap.creation.annual = round(weighted.mean(gap.creation.annual, time),2))
gap.creation.ftype.elev$elev.ftype <- paste0(gap.creation.ftype.elev$elevation ,"-" ,gap.creation.ftype.elev$forest_type)

#----------------- area scaling

area_share_class.elev.ftype <- readRDS("processed/creation/updated/area_share_per_class.elev.ftype.rds")

gap.creation.ftype.elev.scaled <- gap.creation.ftype.elev[,c( "forest_type", "elevation", "elev.ftype", "avg.gap.creation.annual")]
gap.creation.ftype.elev.scaled <- gap.creation.ftype.elev.scaled[!duplicated(gap.creation.ftype.elev.scaled), ]

#merge with area share information
gap.creation.ftype.elev.scaled <- merge(gap.creation.ftype.elev.scaled, area_share_class.elev.ftype[,c("elevation.ftype", "class_area_perc", "total_area", "total_area_category" )],
                                   by.x = "elev.ftype", by.y = "elevation.ftype", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.ftype.elev.scaled$area.scaling.factor <- 100/gap.creation.ftype.elev.scaled$total_area

gap.creation.ftype.elev.scaled$avg.gap.creation.annual_ascaled <- round(gap.creation.ftype.elev.scaled$avg.gap.creation.annual * gap.creation.ftype.elev.scaled$area.scaling.factor,2)

gap.creation.ftype.elev.scaled$elevation <- ordered(gap.creation.ftype.elev.scaled$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
gap.creation.ftype.elev.scaled$elevation <- factor(gap.creation.ftype.elev.scaled$elevation,levels=rev(levels(gap.creation.ftype.elev.scaled$elevation)))



# --- forest type

gap.creation.ftype <- gap_features921 %>% group_by(new.exp, year, forest_type) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, forest_type)%>%
  mutate(avg.gap.creation.annual = round(weighted.mean(gap.creation.annual, time),2),
         sd = sd(gap.creation.annual))


area_share_class <- readRDS("processed/creation/updated/area_share_per_class.rds")
#----------------------------------------------------- area scaling 

gap.creation.ftype.scaled <- gap.creation.ftype[,c("new.exp", "forest_type", "avg.gap.creation.annual", "sd")]
gap.creation.ftype.scaled <- gap.creation.ftype.scaled[!duplicated(gap.creation.ftype.scaled), ]

#merge with area share information
gap.creation.ftype.scaled <- merge(gap.creation.ftype.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                by.x = "forest_type", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.ftype.scaled$area.scaling.factor <- 100/gap.creation.ftype.scaled$total_area

gap.creation.ftype.scaled$avg.gap.creation.annual_ascaled <- round(gap.creation.ftype.scaled$avg.gap.creation.annual * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled$sd_ascaled <- round(gap.creation.ftype.scaled$sd * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled <- gap.creation.ftype.scaled %>% group_by(forest_type) %>%
  mutate(avg.gap.creation.annual_ascaled.both = sum(avg.gap.creation.annual_ascaled))


# --- elevation

gap.creation.elevation<- gap_features921 %>% group_by(new.exp, year, elevation) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, elevation)%>%
  mutate(avg.gap.creation.annual = round(weighted.mean(gap.creation.annual, time),2))

gap.creation.elevation$elevation <- ordered(gap.creation.elevation$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
gap.creation.elevation$elevation <- factor(gap.creation.elevation$elevation,levels=rev(levels(gap.creation.elevation$elevation)))

#----------------------------------------------------- area scaling 

gap.creation.elevation.scaled <- gap.creation.elevation[,c("new.exp", "elevation", "avg.gap.creation.annual")]
gap.creation.elevation.scaled <- gap.creation.elevation.scaled[!duplicated(gap.creation.elevation.scaled), ]

#merge with area share information
gap.creation.elevation.scaled <- merge(gap.creation.elevation.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                   by.x = "elevation", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.elevation.scaled$area.scaling.factor <- 100/gap.creation.elevation.scaled$total_area

gap.creation.elevation.scaled$avg.gap.creation.annual_ascaled <- round(gap.creation.elevation.scaled$avg.gap.creation.annual * gap.creation.elevation.scaled$area.scaling.factor,2)
gap.creation.elevation.scaled <- gap.creation.elevation.scaled %>% group_by(elevation) %>%
  mutate(avg.gap.creation.annual_ascaled.both = sum(avg.gap.creation.annual_ascaled))


# --- aspect

gap.creation.aspect<- gap_features921 %>% group_by(new.exp, year, aspect) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, aspect)%>%
  mutate(avg.gap.creation.annual = round(weighted.mean(gap.creation.annual, time),2))

#----------------------------------------------------- area scaling 

gap.creation.aspect.scaled <- gap.creation.aspect[,c("new.exp", "aspect", "avg.gap.creation.annual")]
gap.creation.aspect.scaled <- gap.creation.aspect.scaled[!duplicated(gap.creation.aspect.scaled), ]

#merge with area share information
gap.creation.aspect.scaled <- merge(gap.creation.aspect.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                       by.x = "aspect", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.aspect.scaled$area.scaling.factor <- 100/gap.creation.aspect.scaled$total_area

gap.creation.aspect.scaled$avg.gap.creation.annual_ascaled <- round(gap.creation.aspect.scaled$avg.gap.creation.annual * gap.creation.aspect.scaled$area.scaling.factor,2)
gap.creation.aspect.scaled <- gap.creation.aspect.scaled %>% group_by(aspect) %>%
  mutate(avg.gap.creation.annual_ascaled.both = sum(avg.gap.creation.annual_ascaled))
#-----------------------------------------------------


#bring elevation labels in right order

gap_features921 <-  gap_features921 %>% mutate(elevation = as.factor(recode(elevation,
                                                               `1`="600-800",
                                                               `2`="800-1000",
                                                               `3`="1000-1200",
                                                               `4`="1200-1400",
                                                               `5`="1400-1600",
                                                               `6`="1600-1800")))
gap_features921$elevation <- ordered(gap_features921$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
gap_features921$elevation <- factor(gap_features921$elevation,levels=rev(levels(gap_features921$elevation)))

# --------------------------------------------------- graphs

wd <- "F:/Projects/CanopyDynamicsBDG/results/graphs/gap_creation/"
setwd(wd)


#new gaps size vs. expansions size - Density plots

tiff("new_exp_density.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=exp.area.ha, fill=factor(new.exp))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  guides(fill=guide_legend(title="type")) +
   My_Theme +  scale_fill_brewer(palette = "Dark2")+ labs( x="gap creation in log [ha]") +
  theme_classic()
dev.off()

ggplot(gap_features921, aes(x=area.ha, fill=factor(new.exp))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  guides(fill=guide_legend(title="type")) +
  My_Theme +  scale_fill_brewer(palette = "Dark2")+ labs( x="gap creation in log [ha]") +
  theme_classic()

# range of new and expanding gaps

q <- c(0.01, 0.25, 0.5, 0.75, 0.99)

gap_features921 %>% group_by(new.exp) %>% 
  summarise(range = max(area.ha)-min(area.ha),
            mean.ha = round(mean(area.ha),2),
            quant1 = round(quantile(area.ha, probs = q[1]),3),
            quant25 = round(quantile(area.ha, probs = q[2]),3), 
            quant50 = round(quantile(area.ha, probs = q[3]),3),
            quant75 = round(quantile(area.ha, probs = q[4]),3),
            quant99 = round(quantile(area.ha, probs = q[5]),3))



tiff("new_exp_density_forest.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=area.ha, fill=factor(forest_type))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  guides(fill=guide_legend(title="forest type")) +
  facet_wrap(~new.exp)+ My_Theme +
   scale_fill_brewer(palette = "Dark2")+ labs( x="gap creation in log [ha]")+
  theme_classic()
dev.off()


tiff("new_exp_density_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=area.ha, fill=factor(elevation))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  guides(fill=guide_legend(title="elevation"))+
  facet_wrap(~new.exp)+ My_Theme+
  scale_fill_brewer(palette = "Blues")+ labs( x="gap creation in log [ha]")+
  theme_classic()
dev.off()


tiff("new_exp_density_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=area.ha, fill=factor(aspect))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  guides(fill=guide_legend(title="type"))+
  facet_wrap(~new.exp)+ My_Theme +
  scale_fill_brewer(palette = "Dark2")+ labs( x="gap creation in log [ha]")+
  theme_classic()
dev.off()

# gap-size-frequency distributions

#calculate nunmber per gap size
# gap.freq <- gap_features921 %>% group_by(area.ha, new.exp) %>%
#   summarize(freq = n())
# 
# ggplot(gap.freq) + geom_jitter(aes(x=area.ha, y=freq, col=new.exp), size=2, alpha=0.5) +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))
# 
# gaps_stats <- gap_features921
# gaps_stats$area.m2 <- gaps_stats$area.ha * 10000
# gaps_stats <- gaps_stats %>% rename( gap_area = area.m2)
# 
# gaps_stats.new <- subset(gaps_stats, new.exp %in% "new")
# gaps_stats.exp <- subset(gaps_stats, new.exp %in% "expanding")
# 
# #all
# GapSizeFDist(
#   gaps_stats = gaps_stats, method = "Hanel_2017", col = "forestgreen", pch = 16, cex = 1,
#   axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
# )
# axis(1)
# axis(2)
# grid(4, 4)
# 
# #new
# GapSizeFDist(
#   gaps_stats = gaps_stats.new, method = "Hanel_2017", col = "red", pch = 16, cex = 1,
#   axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
# )
# axis(1)
# axis(2)
# grid(4, 4)
# title("new gaps")
# 
# #expading
# GapSizeFDist(
#   gaps_stats = gaps_stats.exp, method = "Hanel_2017", col = "blue", pch = 16, cex = 1,
#   axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
# )
# axis(1)
# axis(2)
# grid(4, 4)
# title("expanding gaps")


#new vs. expanding 

tiff("area_ new_exp-timesteps.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation, aes(x=new.exp, y=gap.creation.ha, fill=year))+ 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="gap creation", y="area of gap creation [ha]") +
  guides(fill=guide_legend(title="obs. period")) 
dev.off()


tiff("area_ new_exp.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation, aes(x=new.exp, y=avg.gap.creation.annual.scaled))+ 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="gap creation mechanism", y="area of annual gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="obs. period")) 
dev.off()

# per environmental feature

# forest type ------------------------------------------------------------------


# forest boxplots

ggplot(gap_features921, aes(y=exp.area.ha,  fill= factor(forest_type)))+ 
  geom_boxplot() + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="forest type", y="gap creation area") +
  guides(fill=guide_legend(title="forest type")) +
  facet_wrap(~new.exp)



# scaled

tiff("area_ new_exp_ftype_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.ftype.scaled, aes(x=new.exp, y=avg.gap.creation.annual_ascaled,  fill= factor(forest_type)))+ 
  geom_bar(stat = "identity",position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="gap creation mechanism", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="forest type")) 
dev.off()

tiff("area_ new_exp_ftype_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.ftype.scaled, aes(x=forest_type, y=avg.gap.creation.annual_ascaled,  fill= factor(new.exp)))+ 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="forest type", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="gap creation mechanism")) #+
  # geom_errorbar(aes(ymin=avg.gap.creation.annual_ascaled-sd_ascaled, ymax=avg.gap.creation.annual_ascaled+sd_ascaled), width=.2,
  #               position=position_dodge(.9)) 
dev.off()

tiff("area_ new_exp_ftype_combined.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.ftype.scaled, aes(x=forest_type, y=avg.gap.creation.annual_ascaled.both))+ 
  geom_bar(stat = "identity") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="forest type", y="annual area of gap creation [ha/yr/100 ha]") 
dev.off()


# elevation --------------------------------------------------------------------

# scaled

tiff("area_ new_exp_elevation_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.elevation.scaled, aes(x=new.exp, y=avg.gap.creation.annual_ascaled,  fill= factor(elevation)))+ 
  geom_bar(stat = "identity",position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Blues")+ My_Theme+ labs( x="gap creation mechanism", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="elevation")) 
dev.off()

tiff("area_ new_exp_elevation_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.elevation.scaled, aes(x=elevation, y=avg.gap.creation.annual_ascaled,  fill= factor(new.exp)))+ 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="elevation", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="gap creation mechanism")) 
dev.off()


tiff("area_ new_exp_elevation_combined.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.elevation.scaled, aes(x=elevation, y=avg.gap.creation.annual_ascaled.both))+ 
  geom_bar(stat = "identity") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="elevation", y="annual area of gap creation [ha/yr/100 ha]")
dev.off()


# forest type AND elevation ----------------------------------------------------


# not scaled
ggplot(gap_features921, aes(x=elevation, y=exp.area.ha, fill=forest_type))+ 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="elevation [m]", y="area of gap creation [ha]") +
  guides(fill=guide_legend(title="forest type")) +
  coord_flip() 

# scaled
tiff("area_ftype_elevation_combined.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.ftype.elev.scaled, aes(x=elevation, y=avg.gap.creation.annual_ascaled, fill=forest_type))+ 
  geom_bar(stat = "identity", position =  position_dodge(width = .9)) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="elevation [m]", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="forest type")) +
  coord_flip()+
  geom_text(aes(label = paste0("ha:",total_area), group = forest_type, hjust = 1.1), position=position_dodge(width = .9)) #
dev.off()


# aspect -----------------------------------------------------------------------

# scaled


tiff("area_ new_exp_aspect_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.aspect.scaled, aes(x=new.exp, y=avg.gap.creation.annual_ascaled,  fill= factor(aspect)))+ 
  geom_bar(stat = "identity",position = position_dodge()) + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="gap creation mechanism", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="aspect")) 
dev.off()

tiff("area_ new_exp_aspect_seperated.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.aspect.scaled, aes(x=aspect, y=avg.gap.creation.annual_ascaled,  fill= factor(new.exp)))+ 
  geom_bar(stat = "identity", position="dodge") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="aspect", y="annual area of gap creation [ha/yr/100 ha]") +
  guides(fill=guide_legend(title="gap creation mechanism")) 
dev.off()

tiff("area_ new_exp_aspect_combined.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation.aspect.scaled, aes(x=aspect, y=avg.gap.creation.annual_ascaled.both))+ 
  geom_bar(stat = "identity") + 
  theme_classic()+  scale_fill_brewer(palette = "Dark2")+ My_Theme+ labs( x="aspect", y="annual area of gap creation [ha/yr/100 ha]") 
dev.off()



