######################################
# Analysing Gaps
#####################################


library(plyr)
library(dplyr)
library(terra)
library(ggplot2)
require(scales)
library(tidyverse)


# --- load  Gap layers  ----

gaps2009.3id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2017.3id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")
gaps2021.3id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height3_mmu400n8.tif")

gaps2009.10id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2017.10id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")
gaps2021.10id <- rast("data/processed/gaps_sensitivity/height_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_height10_mmu400n8.tif")

# # crop gaps ID
# gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
# gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 

# new-expanded classification

gaps2017_h3_class <- rast("data/processed/sensitivity/height_sensitivity/gaps2017_new_extended_stable_h3.tif")
gaps2021_h3_class <- rast("data/processed/sensitivity/height_sensitivity/gaps2021_new_extended_stable_h3.tif")

gaps2017_h10_class <- rast("data/processed/sensitivity/height_sensitivity/gaps2017_new_extended_stable_h10.tif")
gaps2021_h10_class <- rast("data/processed/sensitivity/height_sensitivity/gaps2021_new_extended_stable_h10.tif")


#load expansion and closure layer and extract only expansion areas

# height thres 3m
exp_clo917 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_917.tif")
exp917_h3 <- classify(exp_clo917, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height3_1721.tif")
exp1721_h3 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

# height thres 10m
exp_clo917 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_917.tif")
exp917_h10 <- classify(exp_clo917, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("data/processed/sensitivity/height_sensitivity/formation_closure_mmu400n8_height10_1721.tif")
exp1721_h10 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas


# --- load NP information 

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

# ---- crop layers to research sites:

foresttype <- crop(foresttype, gaps2009.3id)
aspect<-  crop(aspect, gaps2009.3id)
elevation.below1800 <- crop(elevation.below1800, gaps2009.3id)
core.zone <- crop(core.zone, gaps2009.3id)

# --- stack gap information and crop it

# 2017 - height thres 3

stack2017_h3 <- c(gaps2017.3id, gaps2017_h3_class, exp917_h3, foresttype , elevation.below1800, aspect)
names(stack2017_h3) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017_h3 <- mask(stack2017_h3, foresttype)
stack2017_h3 <- mask(stack2017_h3, core.zone)
stack2017_h3 <- mask(stack2017_h3, elevation.below1800)

df <- as.data.frame(stack2017_h3, na.rm = FALSE) 

df <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df, "data/processed/sensitivity/height_sensitivity/stack_2017_new_exp_df_h3.rds")

# 2017 - height thres 10

stack2017_h10 <- c(gaps2017.10id, gaps2017_h10_class, exp917_h10, foresttype , elevation.below1800, aspect)
names(stack2017_h10) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017_h10 <- mask(stack2017_h10, foresttype)
stack2017_h10 <- mask(stack2017_h10, core.zone)
stack2017_h10 <- mask(stack2017_h10, elevation.below1800)

df <- as.data.frame(stack2017_h10, na.rm = FALSE) 

df <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df, "data/processed/sensitivity/height_sensitivity/stack_2017_new_exp_df_h10.rds")

# 2021 - height thres 3

stack2021_h3 <- c(gaps2021.3id, gaps2021_h3_class, exp917_h3, foresttype , elevation.below1800, aspect)
names(stack2021_h3) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2021_h3 <- mask(stack2021_h3, foresttype)
stack2021_h3 <- mask(stack2021_h3, core.zone)
stack2021_h3 <- mask(stack2021_h3, elevation.below1800)

df <- as.data.frame(stack2021_h3, na.rm = FALSE) 

df <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df, "data/processed/sensitivity/height_sensitivity/stack_2021_new_exp_df_h3.rds")

# 2017 - height thres 10

stack2021_h10 <- c(gaps2021.10id, gaps2021_h10_class, exp917_h10, foresttype , elevation.below1800, aspect)
names(stack2021_h10) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2021_h10 <- mask(stack2021_h10, foresttype)
stack2021_h10 <- mask(stack2021_h10, core.zone)
stack2021_h10 <- mask(stack2021_h10, elevation.below1800)

df <- as.data.frame(stack2021_h10, na.rm = FALSE) 

df <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df, "data/processed/sensitivity/height_sensitivity/stack_2021_new_exp_df_h10.rds")


rm(list = ls())

# --------calculate features per gap.id

df_917_h3 <- readRDS( "data/processed/sensitivity/height_sensitivity/stack_2017_new_exp_df_h3.rds")
df_1721_h3<- readRDS("data/processed/sensitivity/height_sensitivity/stack_2021_new_exp_df_h3.rds")

df_917_h10 <- readRDS( "data/processed/sensitivity/height_sensitivity/stack_2017_new_exp_df_h10.rds")
df_1721_h10<- readRDS("data/processed/sensitivity/height_sensitivity/stack_2021_new_exp_df_h10.rds")


# height thres 3m 

gap_features_917_h3 <- df_917_h3 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

gap_features_1721_h3 <- df_1721_h3 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

# height thres 10m 

gap_features_917_h10 <- df_917_h10 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

gap_features_1721_h10 <- df_1721_h10 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

# analyse new and expanding gaps ----------------------------------------------------------

# h thres 3m

gap_features_1721_h3$year <- as.factor("17-21")
gap_features_917_h3$year <- as.factor("9-17")

gap_features921_h3 <- rbind(gap_features_917_h3, gap_features_1721_h3)
gap_features921_h3 <- subset(gap_features921_h3, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis
gap_features921_h3 <- gap_features921_h3[gap_features921_h3$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the cropping of the reserach area
gap_features921_h3$h_thres <- as.factor("3") # add indicator for height threshold

# h thres 10m

gap_features_1721_h10$year <- as.factor("17-21")
gap_features_917_h10$year <- as.factor("9-17")

gap_features921_h10 <- rbind(gap_features_917_h10, gap_features_1721_h10)
gap_features921_h10 <- subset(gap_features921_h10, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis
gap_features921_h10 <- gap_features921_h10[gap_features921_h10$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the cropping of the reserach area
gap_features921_h10$h_thres <- as.factor("10") # add indicator for height threshold

# --- load mmu 400 height threshold 5 m gap data --- (used in paper)

gap_features_917 <- readRDS("data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_917.rds")
gap_features_1721 <- readRDS("data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_1721.rds")

gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921_h5 <- rbind(gap_features_917, gap_features_1721)
gap_features921_h5 <- subset(gap_features921_h5, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis

gap_features921_h5 <- gap_features921_h5[gap_features921_h5$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the cropping of the reserach area
gap_features921_h5$h_thres <- as.factor("5")

# reduce to common columns 

columns_to_keep <- colnames(gap_features921_h10)
# Reduce gap_features921_h5 to only those columns in gap_features921_h10
gap_features921_h5 <- gap_features921_h5 %>% select(all_of(columns_to_keep))


# combine all thresholds

gap_features921 <- rbind(gap_features921_h3, gap_features921_h5, gap_features921_h10)
gap_features921 <- na.omit(gap_features921) # 2 rows with all NAs deleted

#calculate annual gap creation rate

# all

gap.creation <- gap_features921 %>% group_by(new.exp, year, h_thres) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, h_thres)%>%
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


gap.creation_unique <- gap.creation %>%
  distinct(median, .keep_all = TRUE)


# --------------------------------------------------- graphs

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


tiff("data/results/sensitivity_analysis/height_threshold/new_exp_density_creation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=exp.area.ha, fill=factor(h_thres))) + geom_density(alpha=.5) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = label_number(accuracy = 0.1))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  guides(fill=guide_legend(title="Height threshold (m)")) +
  theme_classic() +
  My_Theme +  
  scale_fill_brewer(palette = "Set1")+ 
  #scale_fill_colorblind()+
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


tiff("data/results/sensitivity_analysis/height_threshold/area_new_exp.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation_unique, aes(x=new.exp, y=median.scaled, colour=h_thres, group=h_thres, fill=h_thres)) + 
  geom_point(shape = 21, size = 10, position = position_dodge(width = 0.5)) +
  theme_classic() + coord_flip() +
  scale_fill_brewer(palette = "Set1", name = "Height threshold")+ 
  scale_color_brewer(palette = "Set1", name = "Height threshold")+ 
  #scale_fill_colorblind(name = "Height threshold") +
  My_Theme +
  labs(x = "Formation mechanism", y = expression("Annual rate of gap formation (ha per 100 ha)")) + 
  geom_pointrange(aes(ymin=q5_ascaled, ymax=q95_ascaled), linewidth = 3, position = position_dodge(width = 0.5))

dev.off()



