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
library(ggthemes)



wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load classified Gap layers of new and expanding gaps ----

gap_stack <- rast("processed/gaps_sensitivity/gap.stack.mmu100.sensitivity.tif") # layer have been cropped previously to the research area
gaps2017.id<- gap_stack[[2]]
gaps2021.id<- gap_stack[[3]]

# new-expanded classification

gaps2017 <- rast("processed/sensitivity/gaps2017_new_extended_stable_sensitivity.tif")
gaps2021 <- rast("processed/sensitivity/gaps2021_new_extended_stable_sensitivity.tif")

# crop gaps ID
gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 

#load expansion and closure layer and extract only expansion areas
exp_clo <- rast("processed/sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
exp <- classify(exp_clo, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("processed/sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")
exp1721 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

# --- load NP information 

foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")
aspect<-  rast("processed/environment_features/aspect_2021_classified_1m.tif")
elevation <- rast("processed/environment_features/berchtesgaden_2021_classified_200steps_dtm_1m.tif")
closed.forest.core <- vect("raw/closed_forest_corezone.gpkg")
elevation <- rast("processed/environment_features/elevation_below1800_200steps.tif")

# ---- crop layers to research sites:

# crop all layers to same extent

foresttype <- rast("processed/environment_features/forest_type2020_1m.tif")
aspect<-  rast("processed/environment_features/aspect_2021_classified_1m.tif")
closed.forest.core <- vect("raw/closed_forest_corezone.gpkg")
elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")


foresttype <- crop(foresttype, gaps2017.id)
aspect<-  crop(aspect, gaps2017.id)
closed.forest.core <- crop(closed.forest.core, gaps2017.id)
elevation.below1800 <- crop(elevation, gaps2017.id)

# --- stack gap information and crop it

#2017

stack2017 <- c(gaps2017.id, gaps2017, exp, foresttype , elevation.below1800, aspect)
names(stack2017) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017 <- mask(stack2017, foresttype)
# stack2017 <- mask(stack2017, core.zone)
# stack2017 <- mask(stack2017, closed.forest)
# stack2017 <- mask(stack2017, elevation.below1800)

writeRaster(stack2017, "processed/sensitivity/stack.2017.all.gap.information.expansion_sensitivity.tif")
gap_stack_2017 <- rast("processed/sensitivity/stack.2017.all.gap.information.expansion_sensitivity.tif")


df <- as.data.frame(gap_stack_2017, na.rm = FALSE) 

df1 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df1 <- df1[!is.na(df1$forest_type),] # exclude all areas with no forest type information
# df1 <- df1[!is.na(df1$elevation),]# exclude all areas > 1800 m (NA in this case, as it was re-coded above) 

write_rds(df1, "processed/sensitivity/stack_2017_new_exp_df.rds")

#2021

stack21 <- c(gaps2021.id, gaps2021, exp1721, foresttype, elevation.below1800, aspect)
names(stack21) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack21 <- mask(stack21, foresttype)
# stack21 <- mask(stack21, core.zone)
# stack21 <- mask(stack21, closed.forest)
# stack21 <- mask(stack21, elevation.below1800)

writeRaster(stack21, "processed/sensitivity/stack.2021.all.gap.information.expansion_sensitivity.tif", overwrite=T)
gap_stack_2021 <- rast("processed/sensitivity/stack.2021.all.gap.information.expansion_sensitivity.tif")

df <- as.data.frame(gap_stack_2021, na.rm = FALSE) 

df2 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df2 <- df2[!is.na(df2$forest_type),] # exclude all areas with no forest type information
# df2 <- df2[!is.na(df2$elevation),]# exclude all areas > 1800 m (na in this case, as it was recoded above) 

write_rds(df2, "processed/sensitivity/stack_2021_new_exp_df.rds")



#-------calculate area shares per category - not done for sensitivity
# 
# gap_stack_2017.closed <- rast("processed/creation/updated/stack.2017.all.gap.information.expansion.tif")
# 
# df.area <- as.data.frame(gap_stack_2017.closed, na.rm = FALSE) 
# 
# #recode forest type to drop NAs
# # df.area <- df.area %>% mutate( ftype = as.numeric(recode(forest_type,
# #                                                         'broadleaved' = 1,
# #                                                         'coniferous' = 2,
# #                                                         'larch-dominant' = 3,
# #                                                         'mixed stands' = 4)))
# 
# #delete all pixels with no forest type information, as we do not consider these areas + all areas >1800 m 
# df.area.nona <- df.area %>% drop_na(elevation)
# df.area.nona <- df.area.nona %>% drop_na(forest_type)
# df.area.nona <- df.area.nona[df.area.nona$elevation != 7,] #pixels which have been left at boundaries of crop area
# 
# head(df.area.nona)  
# 
# #keep only necessary environmental feature columns
# keeps <- c("elevation","aspect", "forest_type")
# df.area.nona = df.area.nona[keeps] 
# 
# saveRDS(df.area.nona, "processed/creation/updated/df.area.nona.rds")
# df.area.nona <- readRDS("processed/creation/updated/df.area.nona.rds") 
# 
# 
# df.area.long <- gather(df.area.nona, category, class)
# 
# head(df.area.long)
# 
# area_share_class <- df.area.long %>% group_by(category,class) %>% 
#   summarize(total_area =round((sum(!is.na(class)))/10000, 2)) # divided by 10.000 to get ha (res 1m)
# 
# area_share_class <- area_share_class %>% group_by(category) %>% 
#   mutate(total_area_category =round((sum(total_area)), 2),
#             class_area_perc = round(total_area/ total_area_category, 4)) 
# 
# #exclude NAs
# # area_share_class<- area_share_class[!is.na(area_share_class$class),]
# 
# # recode category labels
# class.name <- c("North", "East", "South", "West", 
#                 "600-800", "800-1000", "1000-1200", "1200-1400", "1400-1600", "1600-1800", 
#                # "broadleaved", "coniferous", "larch-dominant", "mixed stands") #old forest type classification
#                 "Beech", "Spruce-fir-beech", "Spruce", "Larch-Swiss stone pine", "Dwarf mountain pine")
# area_share_class$class.name <- class.name
# 
# write_rds(area_share_class, "processed/creation/updated/area_share_per_class.rds")
# 
# 
# # --------calculate area shares per elevation - forest type combination !!!!! neu berechnen wenn notwendig!!!
# df.area.nona$elev.ftype = paste0(df.area.nona$elevation, df.area.nona$forest_type)
# 
# area_share_class.elev.ftype <- df.area.nona %>% group_by(elev.ftype) %>% 
#   summarize(total_area =round((sum(!is.na(elev.ftype)))/10000, 2)) # divided by 10.000 to get ha (res 1m)
# 
# area_share_class.elev.ftype  <- area_share_class.elev.ftype %>% 
#   mutate(total_area_category =round((sum(total_area)), 2),
#          class_area_perc = round(total_area/ total_area_category, 4)*100) 
# 
# elev.ftype <- c("600-800-Beech", "600-800-Spruce-fir-beech", "600-800-Spruce", "600-800-Larch-Swiss stone pine","600-800-Dwarf mountain pine",
#                 "800-1000-Beech","800-1000-Spruce-fir-beech","800-1000-Spruce", "800-1000-Larch-Swiss stone pine","800-1000-Dwarf mountain pine",
#                 "1000-1200-Beech", "1000-1200-Spruce-fir-beech","1000-1200-Spruce","1000-1200-Larch-Swiss stone pine", "1000-1200-Dwarf mountain pine",
#                 "1200-1400-Beech","1200-1400-Spruce-fir-beech","1200-1400-Spruce","1200-1400-Larch-Swiss stone pine", "1200-1400-Dwarf mountain pine",
#                 "1400-1600-Beech", "1400-1600-Spruce-fir-beech","1400-1600-Spruce","1400-1600-Larch-Swiss stone pine", "1400-1600-Dwarf mountain pine",
#                 "1600-1800-Spruce-fir-beech","1600-1800-Spruce","1600-1800-Larch-Swiss stone pine","1600-1800-Dwarf mountain pine" )
# area_share_class.elev.ftype$elevation.ftype <- elev.ftype
# 
# write_rds(area_share_class.elev.ftype, "processed/creation/updated/area_share_per_class.elev.ftype.rds")
# 
# area_share_class.elev.ftype.sub40 <- area_share_class.elev.ftype %>% filter(total_area >40)
# 
# write_rds(area_share_class.elev.ftype.sub40, "processed/creation/updated/area_share_per_class.elev.ftype.sub40.rds")

# --------calculate features per gap.id

df1 <- readRDS( "processed/sensitivity/stack_2017_new_exp_df.rds")
df2<- readRDS("processed/sensitivity/stack_2021_new_exp_df.rds")


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

saveRDS(gap_features_917,"processed/sensitivity/gap_features_new_expanding_917.rds")


#2021

ftype <- getForestType(df2)
elevation <- getElevation(df2)
aspect <- getAspect(df2)

gap_features_1721 <- merge(gap_features_1721, ftype[,c("gap.id","forest_type")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, elevation[,c("gap.id","elevation")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, aspect[,c("gap.id","aspect")], by = "gap.id", all.x = TRUE)

saveRDS(gap_features_1721,"processed/sensitivity/gap_features_new_expanding_1721.rds")
gap_features_1721 <- readRDS("processed/sensitivity/gap_features_new_expanding_1721.rds")


# analyse new and expanding gaps ----------------------------------------------------------


gap_features_917 <- readRDS("processed/sensitivity/gap_features_new_expanding_917.rds")
gap_features_1721 <- readRDS("processed/sensitivity/gap_features_new_expanding_1721.rds")


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
#----!!! change from 400 to 100 !!!
gap_features921 <- gap_features921[gap_features921$area.ha >= 0.01,] #delete gaps smaller than 100m2, as they emerged out of the cropping of the reserach area
gap_features921$mmu <- as.factor("100")

# ----- load mmu 400 gap data

gap_features_917 <- readRDS("processed/sensitivity/version.mmu400/gap_features_new_expanding_917.rds")
gap_features_1721 <- readRDS("processed/sensitivity/version.mmu400/gap_features_new_expanding_1721.rds")

gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921.400 <- rbind(gap_features_917, gap_features_1721)
gap_features921.400 <- subset(gap_features921.400, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis
gap_features921.400 <- gap_features921.400[gap_features921.400$elevation != "1800-2000",]

gap_features921.400 <- gap_features921.400[gap_features921.400$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the cropping of the reserach area
gap_features921.400$mmu <- as.factor("400")


# combine both mmu datasets

gap_features921 <- rbind(gap_features921, gap_features921.400)
gap_features921 <- na.omit(gap_features921)

#calculate annual gap creation rate

# all

gap.creation <- gap_features921 %>% group_by(new.exp, year, mmu) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, mmu)%>%
  mutate(avg.gap.creation.annual = round(mean(gap.creation.annual),2),
         sd = sd(gap.creation.annual,2),
         median = round(median(gap.creation.annual),2),
         q5 = quantile(gap.creation.annual, 0.05),
         q95 = quantile(gap.creation.annual, 0.95))


#--------- area scaling for sensitivity area (214 ha)
gap.creation$area.scaling.factor <- 100/214 #scaling to gap creation per 100 ha (percent!), sensitivity reserach area is 214 ha
gap.creation$avg.gap.creation.annual.scaled <- gap.creation$avg.gap.creation.annual*gap.creation$area.scaling.factor
gap.creation$sd_ascaled <- round(gap.creation$sd * gap.creation$area.scaling.factor,2)

#with median and quantiles
gap.creation$median.scaled <- gap.creation$median*gap.creation$area.scaling.factor
gap.creation$q5_ascaled <- round(gap.creation$q5 * gap.creation$area.scaling.factor,2)
gap.creation$q95_ascaled <- round(gap.creation$q95 * gap.creation$area.scaling.factor,2)




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

gap.creation_unique <- gap.creation %>%
  distinct(median, .keep_all = TRUE)


# --------------------------------------------------- graphs

wd <- "processed/Sensitivity/results/"
setwd(wd)

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position="top")


#new gaps size vs. expansions size - Density plots

tiff("new_exp_density.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=area.ha, fill=factor(new.exp))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = label_number(accuracy = 0.1))+
  guides(fill=guide_legend(title="mechanism of gap creation")) +
  theme_classic() +
  My_Theme +  
  #scale_fill_brewer(palette = "Dark2")+ 
  scale_fill_colorblind()+
  labs( x="gap size log10 [ha]") 
dev.off()


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.3, 0.8))



tiff("new_exp_density_creation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=exp.area.ha, fill=factor(mmu))) + geom_density(alpha=.5) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = label_number(accuracy = 0.1))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  guides(fill=guide_legend(title="mmu")) +
  theme_classic() +
  My_Theme +  
  #scale_fill_brewer(palette = "Dark2")+ 
  scale_fill_colorblind()+
  labs( x="Size of gap formation area in log10 (ha)", y ="Density") +
  facet_grid(~new.exp)
dev.off()



My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(4, "lines"),
  legend.position = c(0.8, 0.8))


#new vs. expanding


tiff("area_new_exp.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation_unique, aes(x=new.exp, y=median.scaled, colour=mmu, group=mmu, fill=mmu)) + 
  geom_point(shape = 21, size = 10, position = position_dodge(width = 0.5)) +
  theme_classic() + coord_flip() +
  scale_color_colorblind(name = "mmu") +
  scale_fill_colorblind(name = "mmu") +
  My_Theme +
  labs(x = "Formation mechanism", y = expression("Annual rate of gap formation (ha per 100 ha)")) + 
  geom_pointrange(aes(ymin=q5_ascaled, ymax=q95_ascaled), linewidth = 3, position = position_dodge(width = 0.5))

dev.off()



