######################################
# Analysing Gaps
#####################################


library(plyr)
library(dplyr)
library(terra)
library(ggplot2)
require(scales)
library(tidyverse)


# --- load classified Gap layers of new and expanding gaps ----

gaps2009 <- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm9_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")
gaps2017<- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm17_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")
gaps2021<- rast("data/processed/gaps_sensitivity/min_size_sensitivity/chm21_sub_sensitivity_patchid_cn2cr2_mmu100n8.tif")

# new-expanded classification

gaps2017 <- rast("data/processed/sensitivity/mmu_sensitivity/gaps2017_new_extended_stable_sensitivity.tif")
gaps2021 <- rast("data/processed/sensitivity/mmu_sensitivity/gaps2021_new_extended_stable_sensitivity.tif")


# crop gaps ID

gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 


#load expansion and closure layer and extract only expansion areas

exp_clo <- rast("data/processed/sensitivity/mmu_sensitivity/exp_clo_917_cn2cr2_mmu100n8_filtered.tif")
exp <- classify(exp_clo, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("data/processed/sensitivity/mmu_sensitivity/exp_clo_1721_cn2cr2_mmu100n8_filtered.tif")
exp1721 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas


# --- load NP information 

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

# ---- crop layers to research sites:

# crop all layers to same extent

foresttype <- crop(foresttype, gaps2017.id)
aspect<-  crop(aspect, gaps2017.id)
closed.forest.core <- crop(closed.forest.core, gaps2017.id)
elevation.below1800 <- crop(elevation, gaps2017.id)

# --- stack gap information and crop it

# 2017

stack2017 <- c(gaps2017.id, gaps2017, exp, foresttype , elevation.below1800, aspect)
names(stack2017) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017 <- mask(stack2017, foresttype)
stack2017 <- mask(stack2017, core.zone)
stack2017 <- mask(stack2017, elevation.below1800)

writeRaster(stack2017, "data/processed/sensitivity/mmu_sensitivity/stack.2017.all.gap.information.expansion_sensitivity.tif")
gap_stack_2017 <- rast("data/processed/sensitivity/mmu_sensitivity/stack.2017.all.gap.information.expansion_sensitivity.tif")


df <- as.data.frame(gap_stack_2017, na.rm = FALSE) 

df1 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df1, "data/processed/sensitivity/mmu_sensitivity/stack_2017_new_exp_df.rds")

# 2021

stack21 <- c(gaps2021.id, gaps2021, exp1721, foresttype, elevation.below1800, aspect)
names(stack21) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack21 <- mask(stack21, foresttype)
stack21 <- mask(stack21, core.zone)
stack21 <- mask(stack21, elevation.below1800)

writeRaster(stack21, "data/processed/sensitivity/mmu_sensitivity/stack.2021.all.gap.information.expansion_sensitivity.tif")
gap_stack_2021 <- rast("data/processed/sensitivity/mmu_sensitivity/stack.2021.all.gap.information.expansion_sensitivity.tif")

df <- as.data.frame(gap_stack_2021, na.rm = FALSE) 

df2 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion

write_rds(df2, "data/processed/sensitivity/mmu_sensitivity/stack_2021_new_exp_df.rds")



# --------calculate features per gap.id

df1 <- readRDS( "data/processed/sensitivity/mmu_sensitivity/stack_2017_new_exp_df.rds")
df2<- readRDS("data/processed/sensitivity/mmu_sensitivity/stack_2021_new_exp_df.rds")


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


# analyse new and expanding gaps ----------------------------------------------------------


gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921 <- rbind(gap_features_917, gap_features_1721)
gap_features921 <- subset(gap_features921, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis

#----!!! change from 400 to 100 !!!
gap_features921 <- gap_features921[gap_features921$area.ha >= 0.01,] #delete gaps smaller than 100m2, as they emerged out of the cropping of the reserach area
gap_features921$mmu <- as.factor("100")

# ----- load mmu 400 gap data

gap_features_917 <- readRDS("data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_917.rds")
gap_features_1721 <- readRDS("data/processed/sensitivity/mmu400_height5/gap_features_new_expanding_1721.rds")

gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921.400 <- rbind(gap_features_917, gap_features_1721)
gap_features921.400 <- subset(gap_features921.400, new.exp %in% c("new", "expanding")) #exclude stable gaps for the analysis
gap_features921.400 <- gap_features921.400[gap_features921.400$elevation != "1800-2000",]

gap_features921.400 <- gap_features921.400[gap_features921.400$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the cropping of the reserach area
gap_features921.400$mmu <- as.factor("400")

# combine both mmu datasets

# reduce to common columns 
columns_to_keep <- colnames(gap_features921)
# Reduce gap_features921_h5 to only those columns in gap_features921.400
gap_features921.400 <- gap_features921.400 %>% select(all_of(columns_to_keep))

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



tiff("data/results/sensitivity_analysis/mmu/new_exp_density_creation.tiff", units="in", width=12, height=8, res=300)
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


tiff("data/results/sensitivity_analysis/mmu/area_new_exp.tiff", units="in", width=12, height=8, res=300)
ggplot(gap.creation_unique, aes(x=new.exp, y=median.scaled, colour=mmu, group=mmu, fill=mmu)) + 
  geom_point(shape = 21, size = 10, position = position_dodge(width = 0.5)) +
  theme_classic() + coord_flip() +
  scale_color_colorblind(name = "mmu") +
  scale_fill_colorblind(name = "mmu") +
  My_Theme +
  labs(x = "Formation mechanism", y = expression("Annual rate of gap formation (ha per 100 ha)")) + 
  geom_pointrange(aes(ymin=q5_ascaled, ymax=q95_ascaled), linewidth = 3, position = position_dodge(width = 0.5))

dev.off()



