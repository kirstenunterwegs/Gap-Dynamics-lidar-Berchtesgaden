#############################
# Extract environmental features of Reserach area
#############################

# load libaries

library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)

# set working directory

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)


# --- load NP information 

foresttype <- rast("raw/forest_final/forest_types2020.tif")
management <- vect("raw/npb_zonierung_22_epsg25832.shp") # where did I create this elevation layer?
aspect<-  rast("processed/environment_features/aspect_2021_classified_1m.tif")
elevation <- rast("processed/environment_features/berchtesgaden_2021_classified_200steps_dtm_1m.tif")
closed.forest <- vect("raw/closed_forest_epsg25832.shp")


# --- extract environmental information for research area
# criteria:
# in the core zone
# < 1800 m
# within closed forest area
# where we have forest type information

# ---- crop environmental layers to research sites:

# --- exclude sites > 1800 m
m <- c(0,6,1, 6,8,NA) #7 and 8 are areas > 1800 m elevation
rclmat <- matrix(m, ncol=3, byrow=TRUE)
elevation.1800 <- classify(elevation, rclmat, include.lowest=TRUE)
elevation.1800.poly <- as.polygons(elevation.1800, trunc=TRUE, dissolve=TRUE, values=TRUE, #prepare masking vector for elevation with area < 1800m
                                   na.rm=TRUE, na.all=FALSE, extent=FALSE)

elevation.below1800 <- mask(elevation, elevation.1800.poly)

writeRaster(elevation.below1800, "processed/environment_features/elevation_below1800_200steps.tif")

elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")

#freq(elevation.below1800) #check if I really excluded all pixels >1800 m


# --- extract core zone to exclude management zone

core.zone <- subset(management, management$zone_id == 4, c(1:2))


# --- reclassify forest types 

# combine spruce-fir-beech (2) and spruce-fir (3)
foresttype <- subst(foresttype, 3, 2)

# combine Larch-Swiss stone pine (5) and Dwarf mountain pine (6)
foresttype <- subst(foresttype, 6, 5)


# resample to match resolution
foresttype <- resample(foresttype, elevation.below1800, method="near")

writeRaster(foresttype, "processed/environment_features/forest_type2020_reclass_1m.tif")
foresttype <- rast("processed/environment_features/forest_type2020_reclass_1m.tif")


#-------calculate area shares per category of study area


environ_stack <- c(foresttype , elevation.below1800, aspect)
names(environ_stack) <- c("forest_type", "elevation", "aspect")

# crop all environmentral layers to the study area
environ_stack.study <- mask(environ_stack, core.zone)
environ_stack.study <- mask(environ_stack.study, closed.forest)
environ_stack.study <- mask(environ_stack.study, foresttype)
environ_stack.study <- mask(environ_stack.study, elevation.below1800)

writeRaster(environ_stack.study, "processed/environment_features/stack_environment_studyarea.tif")

rm(list = ls()) # clear workspace

environ_stack.study <- rast("processed/environment_features/stack_environment_studyarea.tif")

df.area <- as.data.frame(environ_stack.study, na.rm = FALSE) 


#delete all pixels with no forest type information, as we do not consider these areas + all areas >1800 m 
df.area.nona <- df.area %>% drop_na(elevation)
df.area.nona <- df.area.nona %>% drop_na(forest_type)
df.area.nona <- df.area.nona[df.area.nona$elevation != 7,] #pixels which have been left at boundaries of crop area

head(df.area.nona)  


saveRDS(df.area.nona, "processed/environment_features/df.area.nona.rds") # old: processed/creation/updated/df.area.nona.rds

rm(list = ls()) # clear workspace

df.area.nona <- readRDS("processed/environment_features/df.area.nona.rds") 


# bring into long format for analysis

df.area.long <- gather(df.area.nona, category, class)

head(df.area.long)

area_share_class <- df.area.long %>% group_by(category,class) %>% 
  summarize(total_area =round((sum(!is.na(class)))/10000, 2)) # divided by 10.000 to get ha (res 1m)

area_share_class <- area_share_class %>% group_by(category) %>% 
  mutate(total_area_category =round((sum(total_area)), 2),
         class_area_perc = round(total_area/ total_area_category, 4)) 



# recode category labels
class.name <- c("North", "East", "South", "West", 
                "600-800", "800-1000", "1000-1200", "1200-1400", "1400-1600", "1600-1800", 
                "Beech", "Spruce-fir-beech", "Spruce", "Larch-Pine")
area_share_class$class.name <- class.name

saveRDS(area_share_class, "processed/environment_features/area_share_per_class_studyarea.rds")

area_share_class %>%
  group_by(category) %>%
  summarise(sum_class_area_perc = sum(class_area_perc))



# --- load area shares

area_share_class <- readRDS("processed/environment_features/area_share_per_class_studyarea.rds")


#--- prepare df for area share plotting


areas_categories <- area_share_class[,c("category", "total_area", "class.name")]
areas_categories$id <- seq(1, nrow(areas_categories))

data <- as.data.frame(areas_categories)
names(data) <- c("group", "value", "individual", "id")
data$group <- as.factor(data$group) 


### --- plotting area share ----

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 13),
  axis.title.y = element_text(size = 18),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=16),
  strip.text.x = element_text(size = 20))


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

base_data <- base_data %>% mutate(group = recode(group,
                                    `forest_type`="forest type",
                                    `elevation`="elevation [m]"))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 500, xend = start, yend = 500), colour = "grey", alpha=1, size=0.7 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1000, xend = start, yend = 1000), colour = "grey", alpha=1, size=0.7 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1500, xend = start, yend = 1500), colour = "grey", alpha=1, size=0.7 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2000, xend = start, yend = 2000), colour = "grey", alpha=1, size=0.7 , inherit.aes = FALSE ) +

  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(500, 1000, 1500, 2000), label = c("500 ha", "1000 ha", "1500 ha", "2000 ha") , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +

  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-1800,5500) + #change according to value range
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -30, xend = end, yend = -30), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -28, label=group), hjust=c(0.9,0.9,-0.1),vjust=c(1.7,-1.5,0) ,colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Dark2") 
p

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/results/"
setwd(wd)
tiff("area_shares_NP.png", units="in", width=12, height=8, res=300)
p
dev.off()

