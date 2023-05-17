#############################
# analyze gap changes
#############################

library(plyr)
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)
library(sp)
#library(ggsankey)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load area shares

area_share_class <- readRDS("processed/creation/updated/area_share_per_class.rds")

# 
# ################################# --- identify gap share on total area per year --- ##################################
# 
# # function to convert chm and gap information into a df
# 
# get_gap_chm_df <- function(chm, gap, ftype, mtype, elevation, aspect, yr) {
#   stack <- c(chm, gap, ftype, mtype, elevation, aspect)
#   df <- as.data.frame(stack, na.rm=FALSE)
#   names(df)<- c("chm", "gaps", "forest_type", "management", "elevation", "aspect")
#   #df <- df[rowSums(is.na(df)) != ncol(df), ] # delete pixels not in research area
#   df <- df[!is.na(df$chm), ] # delete pixels not in research area
#   df$year <- as.factor(yr)
#   return(df)
# }
# 
# df.2009 <- get_gap_chm_df(chm9, gaps2009, ftype, mtype, elevation, aspect, "2009")
# df.2017 <- get_gap_chm_df(chm17, gaps2017, ftype, mtype, elevation, aspect, "2017")
# df.2021 <- get_gap_chm_df(chm21, gaps2021, ftype, mtype, elevation, aspect, "2021")
# 
# df_all <- rbind(df.2009, df.2017, df.2021)
# 
# saveRDS(df_all,"i:/Fonda/workspace/berchtesgaden/gaps/gap_share_npinfo.rds" )
# 
# df_all<- readRDS("i:/Fonda/workspace/berchtesgaden/gaps/gap_share_npinfo.rds")

#relabel and combine Larch-Swiss stone pine + Dwarf mountain pine



area_share_class <- area_share_class %>% mutate(class.name = recode(class.name,
                                                             `Larch-Swiss stone pine`= "Larch-Pine",
                                                             `Dwarf mountain pine`= "Larch-Pine"))


  sub <- subset(area_share_class, class.name %in% "Larch-Pine")
  sub <- sub %>% mutate(total_area = sum(total_area),
                        class= 3)
  k <- sub[1,]
  area_share_class <- rbind(area_share_class, k)
  area_share_class<-area_share_class[!(area_share_class$category=="forest_type" & area_share_class$class== 6 |area_share_class$category=="forest_type" & area_share_class$class== 5),]


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



########### ---- research area per category (elevation, aspect, management, ftype) ----- ##################################

#--- prepare df for area share plotting

# areas_categories <- df_all %>% group_by(year, forest_type, management, elevation, aspect) %>% 
#   summarize(total_area =round((sum(!is.na(chm)))/10000,2), # divided by 10 000 to get ha
#             area_gaps = round((sum(!is.na(gaps)))/10000,2),
#             area_nogaps = total_area - area_gaps)

areas_categories <- area_share_class[,c("category", "total_area", "class.name")]
areas_categories$id <- seq(1, nrow(areas_categories))

data <- as.data.frame(areas_categories)
names(data) <- c("group", "value", "individual", "id")
data$group <- as.factor(data$group) 
#--------------------------------------------------------------------------
# long_categories <- areas_categories %>% gather(condition, area_ha, -c( total_area))
# long_categories$year <- as.factor(long_categories$year)
# 
# areas_categories$elevation <- as.factor(areas_categories$elevation)
# areas_categories$management <- as.factor(areas_categories$management)
# areas_categories$aspect <- as.factor(areas_categories$aspect)
# 
# 
# areas_categories9 <- subset(areas_categories, year== 2009)
# areas_categories9 <- areas_categories9[-c(1,7:8)]
# 
# areas_categories9 <- areas_categories9 %>% mutate(management = as.factor(recode(management,
#                                        `2`="buffer zone",
#                                        `4`="core zone")) )%>%
#   mutate(aspect = as.factor(recode(aspect,
#                                        `1`="North",
#                                        `2`="East",
#                                        `3`="South",
#                                        `4`="West"))) %>%
#   mutate(elevation = as.factor(recode(elevation,
#                                     `1`="600-800",
#                                     `2`="800-1000",
#                                     `3`="1000-1200",
#                                     `4`="1200-1400",
#                                     `5`="1400-1600",
#                                     `6`="1600-1800",
#                                     `7`="1800-2000",
#                                     `8`="2000-2800")))
#  
# data <- gather(areas_categories9, key = category, value = subcategory, -c("total_area") )
# names(data) <- c("value", "group", "individual")
# data <- data %>% group_by(group,individual)%>% #summarize area information per category
#  summarise(value = sum(value))
# 
# data[is.na(data)] <- "no information" # exchange NAs
# data <- data[-13, ] #remove row with 3.24 ha of no information in forest type
# data[17,3] <- data[17,3] + 3.24 # add 3.24 ha of no entry in forest type to no information group
# #remove no management information 0.37ha
# data <- data[! data$individual== NaN,] 
# 
# data <- as.data.frame(data) #convert from tibble to df
# data$group <- as.factor(data$group) 
# data <- data %>% mutate(group = recode(group,
#                                `forest_type`="forest type",
#                                `elevation`="elevation [m]"))
# 


# Create dataset
# data2 <- data.frame(
#   individual=paste( "Mister ", seq(1,60), sep=""), #individual category
#   group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) , # elevation/aspect/management/ftype
#   value=sample( seq(10,100), 60, replace=T) #area share value
# )

### --- plotting area share ----

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

