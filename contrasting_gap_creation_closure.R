library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# load data
area_share_class <- readRDS("processed/creation/updated/area_share_per_class.rds")
creation <- readRDS("processed/creation/updated/gap_creation_elevation.rds")
closure <- readRDS("processed/closure/updated/gap_closure_elevation.rds")

#--- process data for comparison

# gap creation
creation <- as.data.frame(creation[,c("elevation", "new.exp", "gap.creation.ha") ])

elev <- as.character(unique(creation$elevation))

for(i in elev) {
  sub <- subset(creation, elevation %in% i)
  k <- c(i, "gap.creation.ha", sum(sub$gap.creation.ha))
  creation <- rbind(creation, k)
  creation$gap.creation.ha <- as.numeric(creation$gap.creation.ha)
}

creation$new.exp <- as.factor(creation$new.exp)

creation <- as.data.frame(creation[25:30,c("elevation", "gap.creation.ha")])

# gap closure
closure<- closure %>% 
  mutate(closure_area.ha = round(closure_area/10000,4))%>%
  group_by(elevation) %>%
  summarise(gap.closure.ha = sum(closure_area.ha))

creation.closure <- merge(creation, closure, by = "elevation")


#--- area scaling 


#merge with area share information
creation.closure.scaled <- merge(creation.closure, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                       by.x = "elevation", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
creation.closure.scaled$area.scaling.factor <- 100/creation.closure.scaled$total_area

creation.closure.scaled$gap.creation.ha.scaled <- round(creation.closure.scaled$gap.creation.ha * creation.closure.scaled$area.scaling.factor,2)
creation.closure.scaled$gap.closure.ha.scaled <- round(creation.closure.scaled$gap.closure.ha * creation.closure.scaled$area.scaling.factor,2)

creation.closure.scaled$elevation <- ordered(creation.closure.scaled$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
#creation.closure.scaled$elevation <- factor(creation.closure.scaled$elevation,levels=rev(levels(creation.closure.scaled$elevation)))

sum_creation_closure <- creation.closure.scaled %>% summarize(creation = sum(gap.creation.ha),
                                      closure = sum(gap.closure.ha))

#creation closure
# 278.72  351.52

# re-label to ecological elevation classes

creation.closure.elev.class <- creation.closure.scaled %>% mutate(elev.class = as.factor(recode(elevation,
                                                                 '600-800' = "submontan",
                                                                 '800-1000' = "submontan",
                                                                 '1000-1200' = "montan",
                                                                 '1200-1400' = "montan",
                                                                 '1400-1600' = "montan",
                                                                 '1600-1800' = "sub-alpine"))) %>%
                           group_by(elev.class) %>%
                           summarise(gap.creation.ha.scaled = sum(gap.creation.ha.scaled),
                                     gap.closure.ha.scaled = sum(gap.closure.ha.scaled))

#--- plotting

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 28),
  axis.text.y = element_text(size = 28),
  axis.title.y = element_text(size = 28),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=28),
  legend.text = element_text(size=24),
  strip.text.x = element_text(size = 24),
  legend.position = c(0.87, 0.25))

wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/results/"
setwd(wd)

tiff("creation_closure_elevation.tiff", units="in", width=12, height=8, res=300)
ggplot(creation.closure.scaled, aes(x=gap.creation.ha.scaled, y=gap.closure.ha.scaled, col=elevation)) + 
  geom_point(size=10)+
  geom_abline() +
  xlim(0,12)+ ylim(0,12) +
  theme_classic()+My_Theme +
  scale_color_brewer(palette = "BrBG")+ labs( x="gap creation [ha/100ha]", y= "gap closure [ha/100ha]", color="elevation [m]" )
dev.off()

tiff("creation_closure_elevation.class.tiff", units="in", width=12, height=8, res=300)
ggplot(creation.closure.elev.class, aes(x=gap.creation.ha.scaled, y=gap.closure.ha.scaled, col=elev.class)) + 
  geom_point(size=8)+
  geom_abline() +
  xlim(0,24.5)+ ylim(0,24.5) +
  theme_classic()+My_Theme +
  #scale_color_brewer(palette = "Set1")+ 
  scale_color_manual(values = rev(wes_palette("Chevalier1", n = 3)))+ #GrandBudapest1  Moonrise1
  labs( x="gap creation [ha/100ha]", y= "gap closure [ha/100ha]", color="elevation zone" )
dev.off()

