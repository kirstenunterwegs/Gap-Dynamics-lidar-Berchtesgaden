######################################
# Analysing Gaps
#####################################

library(sf)
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

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

# --- load CHM and Gap layers ----

chm_names <- list("chm9_fs1", "chm17_fs1", "chm21_fs1", "chm9_fs2", "chm17_fs2", "chm21_fs2",
                  "chm9_fs3", "chm17_fs3", "chm21_fs3","chm9_fs4", "chm17_fs4", "chm21_fs4")

# load CHM croped to focus sites + 500m Buffer
chm_fs1_crop <- rast("chm_focus_site1_large_stack.tif")
chm9_fs1 <- chm_fs1_crop[[1]]
chm17_fs1<- chm_fs1_crop[[2]]
chm21_fs1<- chm_fs1_crop[[3]]
chm_fs2_crop <- rast("chm_focus_site2_large_stack.tif")
chm9_fs2 <- chm_fs2_crop[[1]]
chm17_fs2<- chm_fs2_crop[[2]]
chm21_fs2<- chm_fs2_crop[[3]]
chm_fs3_crop <- rast("chm_focus_site3_large_stack.tif")
chm9_fs3 <- chm_fs3_crop[[1]]
chm17_fs3<- chm_fs3_crop[[2]]
chm21_fs3<- chm_fs3_crop[[3]]
chm_fs4_crop <- rast("chm_focus_site4_large_stack.tif")
chm9_fs4 <- chm_fs4_crop[[1]]
chm17_fs4<- chm_fs4_crop[[2]]
chm21_fs4<- chm_fs4_crop[[3]]

chm_list <- list(chm9_fs1, chm17_fs1, chm21_fs1,
                 chm9_fs2, chm17_fs2, chm21_fs2,
                 chm9_fs3, chm17_fs3, chm21_fs3,
                 chm9_fs4, chm17_fs4, chm21_fs4)
names(chm_list) <- chm_names

# load gap layer stacks

gap_stack_fs1 <- rast("gap_layers_fs1_erosion.tif")
gap_stack_fs2 <- rast("gap_layers_fs2_erosion.tif")
gap_stack_fs3 <- rast("gap_layers_fs3_erosion.tif")
gap_stack_fs4 <- rast("gap_layers_fs4_erosion.tif")


gap_list <- list(gap_stack_fs1$gaps_9_fs1,gap_stack_fs1$gaps_17_fs1,gap_stack_fs1$gaps_21_fs1,
                 gap_stack_fs2$gaps_9_fs2,gap_stack_fs2$gaps_17_fs2,gap_stack_fs2$gaps_21_fs2,
                 gap_stack_fs3$gaps_9_fs3,gap_stack_fs3$gaps_17_fs3,gap_stack_fs3$gaps_21_fs3,
                 gap_stack_fs4$gaps_9_fs4,gap_stack_fs4$gaps_17_fs4,gap_stack_fs4$gaps_21_fs4 )
names(gap_list) <- chm_names


# load gap polygons
polygons_erosion <- list()
polygons_erosion <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_erosion_", n ,"/","gaps_polygons_erosion_", n ,".shp", sep=""))
})
names(polygons_erosion) <- chm_names

polygons_no_erosion <- list()
polygons_no_erosion <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_", n ,"/","gaps_polygons_", n ,".shp", sep=""))
})
names(polygons_no_erosion) <- chm_names

# --- Gap statistics ---

stats <- readRDS("gap_stats.RData")

stats_9_fs1 <- as.data.frame(stats$chm9_fs1)
stats_9_fs1$year <- as.factor(2009)
stats_9_fs1$site <- as.factor(1)
stats_17_fs1 <- as.data.frame(stats$chm17_fs1)
stats_17_fs1$year <- as.factor(2017)
stats_17_fs1$site <- as.factor(1)
stats_21_fs1 <- as.data.frame(stats$chm21_fs1)
stats_21_fs1$year <- as.factor(2021)
stats_21_fs1$site <- as.factor(1)

stats_fs1 <- rbind(stats_9_fs1, stats_17_fs1, stats_21_fs1)

stats_9_fs2 <- as.data.frame(stats$chm9_fs2)
stats_9_fs2$year <- as.factor(2009)
stats_9_fs2$site <- as.factor(2)
stats_17_fs2 <- as.data.frame(stats$chm17_fs2)
stats_17_fs2$year <- as.factor(2017)
stats_17_fs2$site <- as.factor(2)
stats_21_fs2 <- as.data.frame(stats$chm21_fs2)
stats_21_fs2$year <- as.factor(2021)
stats_21_fs2$site <- as.factor(2)

stats_fs2 <- rbind(stats_9_fs2, stats_17_fs2, stats_21_fs2)

stats_9_fs3 <- as.data.frame(stats$chm9_fs3)
stats_9_fs3$year <- as.factor(2009)
stats_9_fs3$site <- as.factor(3)
stats_17_fs3 <- as.data.frame(stats$chm17_fs3)
stats_17_fs3$year <- as.factor(2017)
stats_17_fs3$site <- as.factor(3)
stats_21_fs3 <- as.data.frame(stats$chm21_fs3)
stats_21_fs3$year <- as.factor(2021)
stats_21_fs3$site <- as.factor(3)

stats_fs3 <- rbind(stats_9_fs3, stats_17_fs3, stats_21_fs3)

stats_9_fs4 <- as.data.frame(stats$chm9_fs4)
stats_9_fs4$year <- as.factor(2009)
stats_9_fs4$site <- as.factor(4)
stats_17_fs4 <- as.data.frame(stats$chm17_fs4)
stats_17_fs4$year <- as.factor(2017)
stats_17_fs4$site <- as.factor(4)
stats_21_fs4 <- as.data.frame(stats$chm21_fs4)
stats_21_fs4$year <- as.factor(2021)
stats_21_fs4$site <- as.factor(4)

stats_fs4 <- rbind(stats_9_fs4, stats_17_fs4, stats_21_fs4)


stats_all <- rbind(stats_fs1, stats_fs2, stats_fs3, stats_fs4)
stats_all <- stats_all[stats_all$gap_id != 0, ]



area <- ggplot(stats_all, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() 

range <- ggplot(stats_all, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() 

mean <- ggplot(stats_all, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() 

max <- ggplot(stats_all, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() 

# ---calculate perimeter to area ratio

#perimeter to area ratio (P:A) :circular gap will have the lowest P:A and as P:A increases the shape of gaps becomes
# more complex.

perimeter <- lapply(polygons_no_erosion, function(n){
  perim(n)
})
perimeter_unlist <- unlist(perimeter)

stats_all$perimeter <- perimeter_unlist 
stats_all$pa_ratio <- stats_all$perimeter/stats_all$gap_area  

pa <- ggplot(stats_all, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() 


# ------ calculate gap stats -----

Gap_Stats_apply <- function (gap_layer, chm_layer, gap_polygon) 
{ gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * raster::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
# extract raster values per gap
gap_chm <- terra::extract(chm_layer, gap_polygon)
names(gap_chm)[2] <- "chm_values"
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second colum
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
gap_list$perimeter <- perim(gap_polygon)

gap_list$year <- sub("_.*", "", names(chm_layer))
gap_list$site <- sub(".*_", "", names(chm_layer))

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]
# add perimeter to calculate perimeter/area ratio
colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "year", "site")
return(gap_list)
}

#stats_list <- mapply(function(x,y,z){
#  Gap_Stats_apply(gap_layer=x, chm_layer=y, gap_polygon=z)}, x= gap_list , y= chm_list, z= polygons_erosion) 

stats9_1 <- Gap_Stats_apply(gap_list$chm9_fs1, chm_list$chm9_fs1, polygons_erosion$chm9_fs1)
stats_17_1 <- Gap_Stats_apply(gap_list$chm17_fs1, chm_list$chm17_fs1, polygons_erosion$chm17_fs1)
stats_21_1 <- Gap_Stats_apply(gap_list$chm21_fs1, chm_list$chm21_fs1, polygons_erosion$chm21_fs1)

stats9_2 <- Gap_Stats_apply(gap_list$chm9_fs2, chm_list$chm9_fs2, polygons_erosion$chm9_fs2)
stats_17_2 <- Gap_Stats_apply(gap_list$chm17_fs2, chm_list$chm17_fs2, polygons_erosion$chm17_fs2)
stats_21_2 <- Gap_Stats_apply(gap_list$chm21_fs2, chm_list$chm21_fs2, polygons_erosion$chm21_fs2)

stats9_3 <- Gap_Stats_apply(gap_list$chm9_fs3, chm_list$chm9_fs3, polygons_erosion$chm9_fs3)
stats_17_3 <- Gap_Stats_apply(gap_list$chm17_fs3, chm_list$chm17_fs3, polygons_erosion$chm17_fs3)
stats_21_3 <- Gap_Stats_apply(gap_list$chm21_fs3, chm_list$chm21_fs3, polygons_erosion$chm21_fs3)

stats9_4 <- Gap_Stats_apply(gap_list$chm9_fs4, chm_list$chm9_fs4, polygons_erosion$chm9_fs4)
stats_17_4 <- Gap_Stats_apply(gap_list$chm17_fs4, chm_list$chm17_fs4, polygons_erosion$chm17_fs4)
stats_21_4 <- Gap_Stats_apply(gap_list$chm21_fs4, chm_list$chm21_fs4, polygons_erosion$chm21_fs4)

stats_all_erosion <- rbind(stats9_1, stats9_2, stats9_3, stats9_4, 
                           stats_17_1, stats_17_2, stats_17_3, stats_17_4, 
                           stats_21_1, stats_21_2, stats_21_3, stats_21_4)


#stats_trial

area_e <- ggplot(stats_all_erosion, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter")

max_e <- ggplot(stats_all_erosion, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter")

range_e <- ggplot(stats_all_erosion, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter")

mean_e <- ggplot(stats_all_erosion, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter")

stats_all_erosion$pa_ratio <- stats_all_erosion$perimeter/stats_all_erosion$gap_area
pa_e <- ggplot(stats_all_erosion, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter")


Gap_Stats_apply5 <- function (gap_layer, chm_layer, gap_polygon) 
{ gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * raster::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
# extract raster values per gap
gap_chm <- terra::extract(chm_layer, gap_polygon)
names(gap_chm)[2] <- "chm_values"
gap_chm<- gap_chm %>%
  mutate(chm_values = ifelse(chm_values > 5, NA, chm_values)) #filter out values above 5
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second colum
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
gap_list$perimeter <- perim(gap_polygon)

gap_list$year <- sub("_.*", "", names(chm_layer))
gap_list$site <- sub(".*_", "", names(chm_layer))

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]
# add perimeter to calculate perimeter/area ratio
colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "perimeter", "year", "site")
return(gap_list)
}


stats9_1 <- Gap_Stats_apply5(gap_list$chm9_fs1, chm_list$chm9_fs1, polygons_erosion$chm9_fs1)
stats_17_1 <- Gap_Stats_apply5(gap_list$chm17_fs1, chm_list$chm17_fs1, polygons_erosion$chm17_fs1)
stats_21_1 <- Gap_Stats_apply5(gap_list$chm21_fs1, chm_list$chm21_fs1, polygons_erosion$chm21_fs1)

stats9_2 <- Gap_Stats_apply5(gap_list$chm9_fs2, chm_list$chm9_fs2, polygons_erosion$chm9_fs2)
stats_17_2 <- Gap_Stats_apply5(gap_list$chm17_fs2, chm_list$chm17_fs2, polygons_erosion$chm17_fs2)
stats_21_2 <- Gap_Stats_apply5(gap_list$chm21_fs2, chm_list$chm21_fs2, polygons_erosion$chm21_fs2)

stats9_3 <- Gap_Stats_apply5(gap_list$chm9_fs3, chm_list$chm9_fs3, polygons_erosion$chm9_fs3)
stats_17_3 <- Gap_Stats_apply5(gap_list$chm17_fs3, chm_list$chm17_fs3, polygons_erosion$chm17_fs3)
stats_21_3 <- Gap_Stats_apply5(gap_list$chm21_fs3, chm_list$chm21_fs3, polygons_erosion$chm21_fs3)

stats9_4 <- Gap_Stats_apply5(gap_list$chm9_fs4, chm_list$chm9_fs4, polygons_erosion$chm9_fs4)
stats_17_4 <- Gap_Stats_apply5(gap_list$chm17_fs4, chm_list$chm17_fs4, polygons_erosion$chm17_fs4)
stats_21_4 <- Gap_Stats_apply5(gap_list$chm21_fs4, chm_list$chm21_fs4, polygons_erosion$chm21_fs4)

stats_all_erosion <- rbind(stats9_1, stats9_2, stats9_3, stats9_4, 
                           stats_17_1, stats_17_2, stats_17_3, stats_17_4, 
                           stats_21_1, stats_21_2, stats_21_3, stats_21_4)



area_e5 <- ggplot(stats_all_erosion, aes(x=site, y=gap_area, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter and delte val >5")

max_e5 <- ggplot(stats_all_erosion, aes(x=site, y=chm_max, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter and delte val >5")

range_e5 <- ggplot(stats_all_erosion, aes(x=site, y=chm_range, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter and delte val >5")

mean_e5 <- ggplot(stats_all_erosion, aes(x=site, y=chm_mean, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter and delte val >5")

stats_all_erosion$pa_ratio <- stats_all_erosion$perimeter/stats_all_erosion$gap_area
pa_e5 <- ggplot(stats_all_erosion, aes(x=site, y=pa_ratio, color=year)) + geom_boxplot() + theme_minimal() + labs(title="with Erosion filter and delte val >5")

######################################################################################
# Subset gaps per elevation level, management, other regional criteria
######################################################################################