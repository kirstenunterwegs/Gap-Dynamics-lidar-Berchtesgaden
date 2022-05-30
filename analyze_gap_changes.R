#############################
# analyze gap changes
#############################


library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(ForestGapR)
library(ggplot2)
library(ForestTools)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)
library(ggsankey)

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

# min 400
gap_stack_fs1 <- rast("gap_layers_fs1_400.tif")
gap_stack_fs2 <- rast("gap_layers_fs2_400.tif")
gap_stack_fs3 <- rast("gap_layers_fs3_400.tif")
gap_stack_fs4 <- rast("gap_layers_fs4_400.tif")


gap_list_400 <- list(gap_stack_fs1$gaps_9_fs1,gap_stack_fs1$gaps_17_fs1,gap_stack_fs1$gaps_21_fs1,
                     gap_stack_fs2$gaps_9_fs2,gap_stack_fs2$gaps_17_fs2,gap_stack_fs2$gaps_21_fs2,
                     gap_stack_fs3$gaps_9_fs3,gap_stack_fs3$gaps_17_fs3,gap_stack_fs3$gaps_21_fs3,
                     gap_stack_fs4$gaps_9_fs4,gap_stack_fs4$gaps_17_fs4,gap_stack_fs4$gaps_21_fs4 )
names(gap_list_400) <- chm_names

polygons_400 <- list()
polygons_400 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_400_", n ,"/","gaps_polygons_400_", n ,".shp", sep=""))
})
names(polygons_400) <- chm_names

# --- identify gap share on total area per year ---

# function to convert chm and gap information into a df
get_gap_chm_df <- function(chm, gap, fs, yr) {
  stack <- c(chm, gap)
  df <- as.data.frame(stack, na.rm=FALSE)
  names(df)<- c("chm", "gaps")
  df <- df[rowSums(is.na(df)) != ncol(df), ] # delete pixels not in research area
  df$site <- fs
  df$year <- yr
  return(df)
}


fs1_9 <- get_gap_chm_df(chm9_fs1, gap_stack_fs1$gaps_9_fs1, 1, 2009) #fs1
fs1_17<- get_gap_chm_df(chm17_fs1, gap_stack_fs1$gaps_17_fs1, 1, 2017)
fs1_21<- get_gap_chm_df(chm21_fs1, gap_stack_fs1$gaps_21_fs1, 1, 2021)
fs2_9 <- get_gap_chm_df(chm9_fs2, gap_stack_fs2$gaps_9_fs2, 2, 2009) #fs2
fs2_17<- get_gap_chm_df(chm17_fs2, gap_stack_fs2$gaps_17_fs2, 2, 2017)
fs2_21<- get_gap_chm_df(chm21_fs2, gap_stack_fs2$gaps_21_fs2, 2, 2021)
fs3_9 <- get_gap_chm_df(chm9_fs3, gap_stack_fs3$gaps_9_fs3, 3, 2009) #fs3
fs3_17<- get_gap_chm_df(chm17_fs3, gap_stack_fs3$gaps_17_fs3, 3, 2017)
fs3_21<- get_gap_chm_df(chm21_fs3, gap_stack_fs3$gaps_21_fs3, 3, 2021)
fs4_9 <- get_gap_chm_df(chm9_fs4, gap_stack_fs4$gaps_9_fs4, 4, 2009) #fs4
fs4_17<- get_gap_chm_df(chm17_fs4, gap_stack_fs4$gaps_17_fs4, 4, 2017)
fs4_21<- get_gap_chm_df(chm21_fs4, gap_stack_fs4$gaps_21_fs4, 4, 2021)

df_all <- rbind(fs1_9, fs1_17, fs1_21, fs2_9, fs2_17, fs2_21, fs3_9, fs3_17, fs3_21, fs4_9, fs4_17, fs4_21)


areas <- df_all %>% group_by(site, year) %>% 
  summarize(total_area =round((sum(!is.na(chm)))/40000, 2), # divided by 40000 to get ha !!!! /10 000 when res 1m
            area_gaps = round((sum(!is.na(gaps)))/40000,2),
            area_nogaps = total_area - area_gaps)

long <- areas %>% gather(condition, area_ha, -c(year, site, total_area))
long$year <- as.factor(long$year)

ggplot(long, aes(fill=condition, y=area_ha, x=year, label = area_ha)) + 
  geom_bar(position="stack", stat="identity") + theme_minimal()  +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +facet_wrap(~site)





# ---- classify gap changes in gap expansion and gap closure ---

# function to identify gap expansion and closure

gap_change_class <- function(chm9, chm17, gap_change){
  exp_clo <- rast() #create empty raster to classify
  ext(exp_clo) <- ext(chm9)
  res(exp_clo) <- res(chm9)
  crs(exp_clo) <- crs(chm9)
  # classify change group
  exp_clo[gap_change == 2 & chm9 <= 5 & chm17 <=5] <- 1 #steady gaps which fall out of gap detetction or are now included (prev not)
  exp_clo[gap_change == 2 & chm9 > 5 & chm17 <=5] <- 2 #gap expansion
  exp_clo[gap_change == 2 & chm9 <= 5 & chm17 >5] <- 3 #gap closure
  exp_clo[gap_change == 2 & chm9 >5 & chm17 >5] <- 4 #steady vefetation
  return(exp_clo)
}   

# simplified function to identify gap expansion and gap closure

gap_change_class <- function(gap_layer1, gap_layer2){
  exp_clo <- rast() #create empty raster to classify
  ext(exp_clo) <- ext(gap_layer1)
  res(exp_clo) <- res(gap_layer1)
  crs(exp_clo) <- crs(gap_layer1)
  # classify change group
  exp_clo[gap_layer1 >0  & is.na(gap_layer2) ] <- 1 #gap closure
  exp_clo[is.na(gap_layer1) & gap_layer2 >0] <- 2 #gap expansion
  return(exp_clo)
} 

exp_clo_fs1_917 <- gap_change_class(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1) #fs1

# apply reclassification function
exp_clo_fs1_917 <- gap_change_class(chm9_fs1, chm17_fs1, gap_stack_fs1$gap_changes_fs1_917) #fs1
exp_clo_fs1_1721 <- gap_change_class(chm17_fs1, chm21_fs1, gap_stack_fs1$gap_changes_fs1_1721)

exp_clo_fs2_917 <- gap_change_class(chm9_fs2, chm17_fs2, gap_stack_fs2$gap_changes_fs2_917) #fs2
exp_clo_fs2_1721 <- gap_change_class(chm17_fs2, chm21_fs2, gap_stack_fs2$gap_changes_fs2_1721)

exp_clo_fs3_917 <- gap_change_class(chm9_fs3, chm17_fs3, gap_stack_fs3$gap_changes_fs3_917) #fs3
exp_clo_fs3_1721 <- gap_change_class(chm17_fs3, chm21_fs3, gap_stack_fs3$gap_changes_fs3_1721)

exp_clo_fs4_917 <- gap_change_class(chm9_fs4, chm17_fs4, gap_stack_fs4$gap_changes_fs4_917) #fs4
exp_clo_fs4_1721 <- gap_change_class(chm17_fs4, chm21_fs4, gap_stack_fs4$gap_changes_fs4_1721)

#exp_clo[gap_stack_fs1$gap_changes_fs1_917 == 2 & chm9_fs1 <= 5 & chm17_fs1 <=5] <- 1 #steady gap
#exp_clo[gap_stack_fs1$gap_changes_fs1_917 == 2 & chm9_fs1 > 5 & chm17_fs1 <=5] <- 2 #gap expansion
#exp_clo[gap_stack_fs1$gap_changes_fs1_917 == 2 & chm9_fs1 <= 5 & chm17_fs1 >5] <- 3 #gap closure
#exp_clo[gap_stack_fs1$gap_changes_fs1_917 == 2 & chm9_fs1 >5 & chm17_fs1 >5] <- 4 #steady vefetation

# --------------  apply simplified reclassification function
exp_clo_fs1_917 <- gap_change_class(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1) #fs1
exp_clo_fs1_1721 <- gap_change_class(gap_stack_fs1$gaps_17_fs1, gap_stack_fs1$gaps_21_fs1)

exp_clo_fs2_917 <- gap_change_class(gap_stack_fs2$gaps_9_fs2, gap_stack_fs2$gaps_17_fs2) #fs2
exp_clo_fs2_1721 <- gap_change_class(gap_stack_fs2$gaps_17_fs2, gap_stack_fs2$gaps_21_fs2)

exp_clo_fs3_917 <- gap_change_class(gap_stack_fs3$gaps_9_fs3, gap_stack_fs3$gaps_17_fs3)#fs3
exp_clo_fs3_1721 <- gap_change_class(gap_stack_fs3$gaps_17_fs3, gap_stack_fs3$gaps_21_fs3)

exp_clo_fs4_917 <- gap_change_class(gap_stack_fs4$gaps_9_fs4, gap_stack_fs4$gaps_17_fs4) #fs4
exp_clo_fs4_1721 <- gap_change_class(gap_stack_fs4$gaps_17_fs4, gap_stack_fs4$gaps_21_fs4)


#plot gap changes
myTheme <- viridisTheme()
myTheme <- BuRdTheme()
#myTheme$panel.background$col = 'gray' 

#fs1
levelplot(exp_clo_fs1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change fs1 9-17"),
          colorkey=list(at=seq(0, 2, 1),labels=list(at=c(1, 2), 
                                                    labels=c( "gap closure", "gap expansion"))), scales=list(draw=FALSE)) + 
  layer_(sp.polygons(as(polygons_400$chm9_fs1, "Spatial"), fill='grey', alpha=0.5))

levelplot(exp_clo_fs1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change fs1 17-21"),
          colorkey=list(at=seq(0, 2, 1),labels=list(at=c(1, 2), 
                                                    labels=c( "gap closure", "gap expansion"))), scales=list(draw=FALSE)) +
  layer_(sp.polygons(as(polygons_400$chm17_fs1, "Spatial"), fill='grey', alpha=0.5))

#fs2
levelplot(exp_clo_fs2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change fs2 9-17"),
          colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), 
                                                    labels=c("gap area newly or not anymore detetcted", "gap expansion", "gap closure"))), scales=list(draw=FALSE))
levelplot(exp_clo_fs2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change fs2 17-21"),
          colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), 
                                                    labels=c("gap area newly or not anymore detetcted", "gap expansion", "gap closure"))), scales=list(draw=FALSE))
#fs3
levelplot(exp_clo_fs3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change fs3 9-17"),
          colorkey=list(at=seq(0, 2, 1),labels=list(at=c(1, 2), 
                                                    labels=c( "gap closure", "gap expansion"))), scales=list(draw=FALSE)) + 
  layer_(sp.polygons(as(polygons_400$chm9_fs3, "Spatial"), fill='grey', alpha=0.5))

levelplot(exp_clo_fs3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change fs3 17-21"),
          colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), 
                                                    labels=c("gap area newly or not anymore detetcted", "gap expansion", "gap closure"))), scales=list(draw=FALSE))
#fs4
levelplot(exp_clo_fs4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change fs4 9-17"),
          colorkey=list(at=seq(0, 2, 1),labels=list(at=c(1, 2), 
                                                    labels=c( "gap closure", "gap expansion"))), scales=list(draw=FALSE)) + 
  layer_(sp.polygons(as(polygons_400$chm9_fs4, "Spatial"), fill='grey', alpha=0.5))

levelplot(exp_clo_fs4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change fs4 17-21"),
          colorkey=list(at=seq(0, 2, 1),labels=list(at=c(1, 2), 
                                                    labels=c( "gap closure", "gap expansion"))), scales=list(draw=FALSE)) + 
  layer_(sp.polygons(as(polygons_400$chm17_fs4, "Spatial"), fill='grey', alpha=0.5))

# --- analyze gap closure and gap expansion 

# function to convert chm and gap information into a df
get_cloexp_df <- function(chm, change, fs, ts) {
  stack <- c(chm, change)
  df <- as.data.frame(stack, na.rm=FALSE)
  names(df)<- c("chm", "gap_change_class")
  df <- df[rowSums(is.na(df)) != ncol(df), ] # delete pixels not in research area
  df[df == "NaN"] <- 4 # DAs muss ich eigentlich löschen!!!!!
  df$site <- as.factor(fs)
  df$timestep <- as.factor(ts)
  return(df)
}

clo_exp_df1917 <- get_cloexp_df(chm9_fs1, exp_clo_fs1_917, 1, "9-17") #fs1
clo_exp_df11721 <- get_cloexp_df(chm17_fs1, exp_clo_fs1_1721, 1, "17-21")
clo_exp_df2917 <- get_cloexp_df(chm9_fs2, exp_clo_fs2_917, 2, "9-17") #fs2
clo_exp_df21721 <- get_cloexp_df(chm17_fs2, exp_clo_fs2_1721, 2, "17-21")
clo_exp_df3917 <- get_cloexp_df(chm9_fs3, exp_clo_fs3_917, 3, "9-17") #fs3
clo_exp_df31721 <- get_cloexp_df(chm17_fs3, exp_clo_fs3_1721, 3, "17-21")
clo_exp_df4917 <- get_cloexp_df(chm9_fs4, exp_clo_fs4_917, 4, "9-17") #fs4
clo_exp_df41721 <- get_cloexp_df(chm17_fs4, exp_clo_fs4_1721, 4, "17-21")

clo_exp_df_all <- rbind(clo_exp_df1917, clo_exp_df11721,
                        clo_exp_df2917, clo_exp_df21721,
                        clo_exp_df3917, clo_exp_df31721,
                        clo_exp_df4917, clo_exp_df41721)

# ich sollte an diesem Punkt die share der change area am gesamten Gebiet berechnen
# dann pixel welche keine change area sind (also mit NA) rausschmeißen und die Flächen der einzelnen change classes berechnen

clo_exp_df_summary <- clo_exp_df_all %>%
  group_by(gap_change_class, site, timestep) %>%
  summarise(area = n()/40000, 0) %>% # /40000 to get ha, change to /10000 when res 1m
  mutate(area_share = area/sum(area),
         gap_change_class = as.factor(recode(gap_change_class,
                                             `1`="gap area newly or not anymore detetcted", #through vegetation growth or loss
                                             `2`="gap expansion through vegetation loss",
                                             `3`="gap closure through vegetation growth",
                                             `4`="steady vegetation")))

clo_exp_df_summary <- clo_exp_df_all %>%
  group_by(gap_change_class, site, timestep) %>%
  summarise(area = n()/40000, 0) %>% # /40000 to get ha, change to /10000 when res 1m
  mutate(area_share = area/sum(area),
         gap_change_class = as.factor(recode(gap_change_class,
                                             `1`="gap closure through vegetation growth", #through vegetation growth or loss
                                             `2`="gap expansion through vegetation loss")),
         years_ts = recode(timestep,
                           "9-17" = 8,
                           "17-21" = 4),
         area_per_yr = area/years_ts)



# delete vegetation to not distort results 

clo_exp_df_summary2 <- subset(clo_exp_df_summary, gap_change_class != "steady vegetation")

#absolute values
ggplot(clo_exp_df_summary2, aes(x= timestep, y= area, fill=gap_change_class, label=area )) +  
  geom_bar(position="stack", stat="identity") + ylab("dynamic area in ha") +
  theme_minimal() +facet_wrap(~site)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +facet_wrap(~site)
#relative values
ggplot(clo_exp_df_summary2, aes(x= timestep, y= area_per_yr, fill=gap_change_class, label=area )) +  
  geom_bar(position="stack", stat="identity") + ylab("dynamic area in ha/yr") +
  theme_minimal() +facet_wrap(~site)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +facet_wrap(~site)

# --- classify trajectories ---

#gaps9 <- gap_stack_fs1$gaps_9_fs1
#gaps17 <-gap_stack_fs1$gaps_17_fs1
#gaps21 <- gap_stack_fs1$gaps_21_fs1

classify_trajectories <- function(gaps9, gaps17, gaps21){
  trajectories <- rast() #create empty raster to classify
  ext(trajectories) <- ext(gaps9)
  res(trajectories) <- res(gaps9)
  crs(trajectories) <- crs(gaps9)
  # classify trajectories
  trajectories[is.nan(gaps9) & is.nan(gaps17) & is.nan(gaps21) ] <- 1 #steady vegetation
  trajectories[is.nan(gaps9) & is.nan(gaps17) & gaps21 >0 ] <- 2 #steady vegetation - gap expansion
  trajectories[is.nan(gaps9) & gaps17 >0 & is.nan(gaps21)  ] <- 3 #gap expansion - gap closure
  trajectories[gaps9 >0 & is.nan(gaps17) & is.nan(gaps21)  ] <- 4 #gap closure - steady vegetation
  trajectories[is.nan(gaps9) & gaps17 >0  & gaps21 >0   ] <- 5 #gap expansion - steady gap
  trajectories[ gaps9 >0 &  is.nan(gaps17) & gaps21 >0   ] <- 6 # gap closure - gap expansion 
  trajectories[ gaps9 >0 & gaps17 >0 & is.nan(gaps21)  ] <- 7 # steady gap - gap closure 
  trajectories[ gaps9 >0 & gaps17 >0 & gaps21 >0 ] <- 8 # steady gap 
  return(trajectories)
} 

trajectories_fs1 <- classify_trajectories(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1, gap_stack_fs1$gaps_21_fs1)
trajectories_fs2 <- classify_trajectories(gap_stack_fs2$gaps_9_fs2, gap_stack_fs2$gaps_17_fs2, gap_stack_fs2$gaps_21_fs2)
trajectories_fs3 <- classify_trajectories(gap_stack_fs3$gaps_9_fs3, gap_stack_fs3$gaps_17_fs3, gap_stack_fs3$gaps_21_fs3)
trajectories_fs4 <- classify_trajectories(gap_stack_fs4$gaps_9_fs4, gap_stack_fs4$gaps_17_fs4, gap_stack_fs4$gaps_21_fs4)

library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)


levelplot(trajectories_fs1, margin = FALSE, main= list("Trajectories focus site 1"),
          colorkey=list(at=seq(0,8,1),labels=list(at=c(1, 2, 3, 4, 5, 6, 7, 8), 
                                                  labels=c("steady vegetation", "steady vegetation - gap expansion", "gap expansion - gap closure",
                                                           "gap closure - steady vegetation", "gap expansion - steady gap", "gap closure - gap expansion",
                                                           "steady gap - gap closure ", "steady gap "))), scales=list(draw=FALSE))

levelplot(trajectories_fs2, margin = FALSE, main= list("Trajectories focus site 2"),
          colorkey=list(at=seq(0,8,1),labels=list(at=c(1, 2, 3, 4, 5, 6, 7, 8), 
                                                  labels=c("steady vegetation", "steady vegetation - gap expansion", "gap expansion - gap closure",
                                                           "gap closure - steady vegetation", "gap expansion - steady gap", "gap closure - gap expansion",
                                                           "steady gap - gap closure ", "steady gap "))), scales=list(draw=FALSE))

levelplot(trajectories_fs3, margin = FALSE,  main= list("Trajectories focus site 3"),
          colorkey=list(at=seq(0,8,1),labels=list(at=c(1, 2, 3, 4, 5, 6, 7, 8), 
                                                  labels=c("steady vegetation", "steady vegetation - gap expansion", "gap expansion - gap closure",
                                                           "gap closure - steady vegetation", "gap expansion - steady gap", "gap closure - gap expansion",
                                                           "steady gap - gap closure ", "steady gap "))), scales=list(draw=FALSE))

levelplot(trajectories_fs4, margin = FALSE, main= list("Trajectories focus site 4"),
          colorkey=list(at=seq(0,8,1),labels=list(at=c(1, 2, 3, 4, 5, 6, 7, 8), 
                                                  labels=c("steady vegetation", "steady vegetation - gap expansion", "gap expansion - gap closure",
                                                           "gap closure - steady vegetation", "gap expansion - steady gap", "gap closure - gap expansion",
                                                           "steady gap - gap closure ", "steady gap "))), scales=list(draw=FALSE))


# -- create sankey diagram for gap trajectories
stack_gaps_fs1 <- c(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1, gap_stack_fs1$gaps_21_fs1)
trajectories_fs1_df<- as.data.frame(stack_gaps_fs1, na.rm=FALSE)

change_vertical_bins["Change_class"][change_vertical_bins["Change_class"] == 1] <- "height loss"

#trajectories_fs1_df["gaps_9_fs1"][trajectories_fs1_df["gaps_9_fs1"] == "NaN"] <- "vegetation"
trajectories_fs1_df[trajectories_fs1_df == "NaN"] <- "vegetation"
trajectories_fs1_df[trajectories_fs1_df != "NaN"] <- "gap"
colnames(trajectories_fs1_df) <- c("2009", "2017", "2021")

df <- trajectories_fs1_df %>%
  make_long("2009", "2017", "2021")

# Chart 1
pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = FALSE)
pl <- pl + geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl <- pl +  theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(title = "Sankey diagram for gaps")
pl <- pl + labs(subtitle = "flow unit is pixels per gap")
#pl <- pl + labs(caption = "@techanswers88")
pl <- pl + labs(fill = 'Nodes')
pl

#vertical barplot with share on different trajectories
trajectories_fs1_df1 <- as.data.frame(trajectories_fs1) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area = (n * res(trajectories_fs1)[1] * res(trajectories_fs1)[2])/10000,
site = as.factor(1),
trajectory = as.factor(recode(lyr.1,
                              `1`="steady vegetation",
                              `2`="steady vegetation - gap expansion",
                              `3`="gap expansion - gap closure",
                              `4`= "gap closure - steady vegetation",
                              `5`="gap expansion - steady gap",
                              `6`="gap closure - gap expansion",
                              `7`="steady gap - gap closure ",
                              `8`="steady gap ")))

ggplot(trajectories_fs1_df1, aes(x=trajectory, y=area_ha)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_minimal()


trajectories_fs2_df2 <- as.data.frame(trajectories_fs2) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area_ha = round((n * res(trajectories_fs2)[1] * res(trajectories_fs2)[2])/10000, 2), 
         site = as.factor(2),
         trajectory = as.factor(recode(lyr.1,
                                       `1`="steady vegetation",
                                       `2`="steady vegetation - gap expansion",
                                       `3`="gap expansion - gap closure",
                                       `4`= "gap closure - steady vegetation",
                                       `5`="gap expansion - steady gap",
                                       `6`="gap closure - gap expansion",
                                       `7`="steady gap - gap closure ",
                                       `8`="steady gap ")))

ggplot(trajectories_fs2_df2, aes(x=trajectory, y=area_ha)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_minimal()

trajectories_fs3_df1 <- as.data.frame(trajectories_fs3) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area_ha = round((n * res(trajectories_fs3)[1] * res(trajectories_fs3)[2])/10000, 2), 
         site = as.factor(3),
         trajectory = as.factor(recode(lyr.1,
                                       `1`="steady vegetation",
                                       `2`="steady vegetation - gap expansion",
                                       `3`="gap expansion - gap closure",
                                       `4`= "gap closure - steady vegetation",
                                       `5`="gap expansion - steady gap",
                                       `6`="gap closure - gap expansion",
                                       `7`="steady gap - gap closure ",
                                       `8`="steady gap ")))

ggplot(trajectories_fs3_df1, aes(x=trajectory, y=area_ha)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_minimal()

trajectories_fs4_df1 <- as.data.frame(trajectories_fs4) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area_ha = round((n * res(trajectories_fs4)[1] * res(trajectories_fs4)[2])/10000, 2), 
         site = as.factor(4),
         trajectory = as.factor(recode(lyr.1,
                                       `1`="steady vegetation",
                                       `2`="steady vegetation - gap expansion",
                                       `3`="gap expansion - gap closure",
                                       `4`= "gap closure - steady vegetation",
                                       `5`="gap expansion - steady gap",
                                       `6`="gap closure - gap expansion",
                                       `7`="steady gap - gap closure ",
                                       `8`="steady gap ")))

ggplot(trajectories_fs4_df1, aes(x=trajectory, y=area_ha)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_minimal()

trajectories_all_df <- rbind(trajectories_fs1_df1, trajectories_fs2_df2, trajectories_fs3_df1, trajectories_fs4_df1)

ggplot(trajectories_all_df, aes(x=trajectory, y=area_ha, fill=site)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) + coord_flip() + theme_minimal() #+facet_wrap(~site)

# --- aggregate trajectories 

#trajectories <- trajectories_fs1

aggregate_trajectories <- function(trajectories){
  trajectories_agg <- rast() #create empty raster to classify
  ext(trajectories_agg ) <- ext(trajectories)
  res(trajectories_agg ) <- res(trajectories)
  crs(trajectories_agg ) <- crs(trajectories)
  # classify trajectories
  trajectories_agg [trajectories == 1 ] <- 1                    #steady vegetation
  trajectories_agg [trajectories == 2 | trajectories == 5 ] <- 2 # gap expansion
  trajectories_agg [trajectories == 4 | trajectories == 7 ] <- 3 # gap closure
  trajectories_agg [trajectories == 8 ] <- 4                    #steady gap
  trajectories_agg [trajectories == 3 | trajectories == 6 ] <- 5 #twice state change
  
  return(trajectories_agg )
} 

trajectories_agg_fs1 <- aggregate_trajectories(trajectories_fs1)
trajectories_agg_fs2 <- aggregate_trajectories(trajectories_fs2)
trajectories_agg_fs3 <- aggregate_trajectories(trajectories_fs3)
trajectories_agg_fs4 <- aggregate_trajectories(trajectories_fs4)

rainbTheme5 <- rasterTheme(region = rainbow(n = 5))
# add par.setting =rainbTheme5

levelplot(trajectories_agg_fs1, margin = FALSE, main= list("Trajectories aggregated focus site 1"),
          colorkey=list(at=seq(0,5,1),labels=list(at=c(1, 2, 3, 4, 5), 
                                                  labels=c("steady vegetation", "gap expansion", "gap closure",
                                                           "steady gap ", "twice state change"))), scales=list(draw=FALSE))

levelplot(trajectories_agg_fs2, margin = FALSE, main= list("Trajectories aggregated focus site 2"),
          colorkey=list(at=seq(0,5,1),labels=list(at=c(1, 2, 3, 4, 5), 
                                                  labels=c("steady vegetation", "gap expansion", "gap closure",
                                                           "steady gap ", "twice state change"))), scales=list(draw=FALSE))

levelplot(trajectories_agg_fs3, margin = FALSE, main= list("Trajectories aggregated focus site 3"),
          colorkey=list(at=seq(0,5,1),labels=list(at=c(1, 2, 3, 4, 5), 
                                                  labels=c("steady vegetation", "gap expansion", "gap closure",
                                                           "steady gap ", "twice state change"))), scales=list(draw=FALSE))

levelplot(trajectories_agg_fs4, margin = FALSE,  main= list("Trajectories aggregated focus site 4"),
          colorkey=list(at=seq(0,5,1),labels=list(at=c(1, 2, 3, 4, 5), 
                                                  labels=c("steady vegetation", "gap expansion", "gap closure",
                                                           "steady gap ", "twice state change"))), scales=list(draw=FALSE))

# area per trajectory
as.data.frame(trajectories_agg_fs1) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area = n * res(trajectories_agg_fs1)[1] * res(trajectories_agg_fs1)[2])

as.data.frame(trajectories_agg_fs2) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area = n * res(trajectories_agg_fs2)[1] * res(trajectories_agg_fs2)[2])

as.data.frame(trajectories_agg_fs3) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area = n * res(trajectories_agg_fs3)[1] * res(trajectories_agg_fs3)[2])

as.data.frame(trajectories_agg_fs4) %>%
  group_by(lyr.1) %>%
  tally() %>%
  mutate(area = n * res(trajectories_agg_fs4)[1] * res(trajectories_agg_fs4)[2])

# --- calculate splitting and merging factor ---

# get dataframe of gap development trhough time
#focus site 1
gaps_fs1 <- c(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1, gap_stack_fs1$gaps_21_fs1)
gaps_fs1_df <- as.data.frame(gaps_fs1, na.rm=FALSE)
gaps_fs1_df <- gaps_fs1_df[rowSums(is.na(gaps_fs1_df)) != ncol(gaps_fs1_df), ] # delete pixels without any gap at any moment in time
gaps_fs1_df[gaps_fs1_df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
#gaps_fs1_df[is.na(gaps_fs1_df)] <- 0 

#focus site 3
gaps_fs2 <- c(gap_stack_fs2$gaps_9_fs2, gap_stack_fs2$gaps_17_fs2, gap_stack_fs2$gaps_21_fs2)
gaps_fs2_df <- as.data.frame(gaps_fs2, na.rm=FALSE)
gaps_fs2_df <- gaps_fs2_df[rowSums(is.na(gaps_fs2_df)) != ncol(gaps_fs2_df), ] # delete pixels without any gap at any moment in time
gaps_fs2_df[gaps_fs2_df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel

#focus site 3
gaps_fs3 <- c(gap_stack_fs3$gaps_9_fs3, gap_stack_fs3$gaps_17_fs3, gap_stack_fs3$gaps_21_fs3)
gaps_fs3_df <- as.data.frame(gaps_fs3, na.rm=FALSE)
gaps_fs3_df <- gaps_fs3_df[rowSums(is.na(gaps_fs3_df)) != ncol(gaps_fs3_df), ] # delete pixels without any gap at any moment in time
gaps_fs3_df[gaps_fs3_df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel

#focus site 4
gaps_fs4 <- c(gap_stack_fs4$gaps_9_fs4, gap_stack_fs4$gaps_17_fs4, gap_stack_fs4$gaps_21_fs4)
gaps_fs4_df <- as.data.frame(gaps_fs4, na.rm=FALSE)
gaps_fs4_df <- gaps_fs4_df[rowSums(is.na(gaps_fs4_df)) != ncol(gaps_fs4_df), ] # delete pixels without any gap at any moment in time
gaps_fs4_df[gaps_fs4_df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel


#Loop to detetct splits

#fs1 -9-17
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs1_df$gaps_9_fs1)) {
  
  df <- subset(gaps_fs1_df, gaps_9_fs1 == i) 
  splits <- unique(df$gaps_17_fs1)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits1_917 <- splits_df
splits1_917$site <- as.factor(1)
splits1_917$timestep <- as.factor("9-17")

#fs1 17_19
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs1_df$gaps_17_fs1)) {
  
  df <- subset(gaps_fs1_df, gaps_17_fs1 == i) 
  splits <- unique(df$gaps_21_fs1)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits1_1721 <- splits_df
splits1_1721$site <- as.factor(1)
splits1_1721$timestep <- as.factor("17-21")

#fs2 -9-17
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs2_df$gaps_9_fs2)) {
  
  df <- subset(gaps_fs2_df, gaps_9_fs2 == i) 
  splits <- unique(df$gaps_17_fs2)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits2_917 <- splits_df
splits2_917$site <- as.factor(2)
splits2_917$timestep <- as.factor("9-17")

#fs1 17_19
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs2_df$gaps_17_fs2)) {
  
  df <- subset(gaps_fs2_df, gaps_17_fs2 == i) 
  splits <- unique(df$gaps_21_fs2)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits2_1721 <- splits_df
splits2_1721$site <- as.factor(2)
splits2_1721$timestep <- as.factor("17-21")

#fs3 -9-17
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs3_df$gaps_9_fs3)) {
  
  df <- subset(gaps_fs3_df, gaps_9_fs3 == i) 
  splits <- unique(df$gaps_17_fs3)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits3_917 <- splits_df
splits3_917$site <- as.factor(3)
splits3_917$timestep <- as.factor("9-17")

#fs3 17_19
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs3_df$gaps_17_fs3)) {
  
  df <- subset(gaps_fs3_df, gaps_17_fs3 == i) 
  splits <- unique(df$gaps_21_fs3)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits3_1721 <- splits_df
splits3_1721$site <- as.factor(3)
splits3_1721$timestep <- as.factor("17-21")

#fs4 -9-17
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs4_df$gaps_9_fs4)) {
  
  df <- subset(gaps_fs4_df, gaps_9_fs4 == i) 
  splits <- unique(df$gaps_17_fs4)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits4_917 <- splits_df
splits4_917$site <- as.factor(4)
splits4_917$timestep <- as.factor("9-17")

#fs1 17_19
splits_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs4_df$gaps_17_fs4)) {
  
  df <- subset(gaps_fs4_df, gaps_17_fs4 == i) 
  splits <- unique(df$gaps_21_fs4)
  
  if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
  if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and nothing into veg pool
  if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
  splits_df <- rbind(splits_df, split_info)
  names(splits_df) <- c("gap_id", "no_splits", "split_yes_no")
}

splits4_1721 <- splits_df
splits4_1721$site <- as.factor(4)
splits4_1721$timestep <- as.factor("17-21")

#-----------------
# Test
i=14
df <- subset(gaps_fs1_df, gaps_9_fs1 == 14) 
splits <- unique(df$gaps_17_fs1)

if (length(splits) == 1 & splits == 0) {split_info <- c(i,0,0)} # gap disappears (veg growth)
if (length(splits) == 1 & splits != 0) {split_info <- c(i,1,0)} # gap stays as one gap
if (length(splits) == 2 & 0 %in% splits) {split_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
if (length(splits) == 2 & !(0 %in% splits)) {split_info <- c(i,2,1)} # gap splits once and not into veg pool
if (length(splits) >2 & !(0 %in% splits) ) {split_info <- c(i, length(splits),1)} # gap splits more than twice and not into veg pool
if (length(splits) >2 &  0 %in% splits ) {split_info <- c(i, (length(splits)-1),1)} # gap splits more than twice
splits_df <- rbind(splits_df, split_info)
#------------------

splits_all <- rbind(splits1_917, splits1_1721, splits2_917, splits2_1721, splits3_917, splits3_1721, splits4_917, splits4_1721)

# -- calculate split factors

splits_all_noVeg <- splits_all[splits_all$gap_id != 0, ] # exclude the vegetation pool

split_factor_all <- splits_all_noVeg %>%
  group_by(site, timestep) %>%
  summarise(n_gaps_split = sum(split_yes_no[split_yes_no ==1]),
            sum_gap_splits = sum(no_splits[split_yes_no ==1]),
            n_gaps = n(),
            split_factor = sum_gap_splits/n_gaps_split,
            split_share = n_gaps_split/ n_gaps,
            global_split_factor = sum_gap_splits / n_gaps)

ggplot(split_factor_all, aes(x=timestep, y=global_split_factor, color=site, group=site)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "gobal split factor")

ggplot(split_factor_all, aes(x=timestep, y=split_share, color=site, group=site)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "split_share")

ggplot(split_factor_all, aes(x=timestep, y=split_factor , color=site, group=site)) + 
  geom_point(size=3) + geom_line() + theme_minimal() +  labs(title = "split factor") +geom_jitter(width = 0.1, height = 0.1, size=3)

#splits_df_gaps[splits_df_gaps$split_yes_no ==1, ] # df of gap, which split
#sum(splits_df_gaps$no_splits) # sum of splits (but of all gaps, also those which don't split)
#sum(splits_df_gaps$split_yes_no == 1) # amount of gaps which split

# calculate split factor (sum of gap splits of gaps which split / number of gaps which split)
split_factor <- (sum(splits_df_gaps[splits_df_gaps$split_yes_no ==1, ] $no_splits)) / (sum(splits_df_gaps$split_yes_no == 1))

# calculate share of gaps which split (number of gaps which split / number of gaps at t0)
split_share <- (sum(splits_df_gaps$split_yes_no == 1)) / nrow(splits_df_gaps)

# calculate global split factor (sum of gap splits of gaps which split / number of gaps at t0 )
global_split_factor <- (sum(splits_df_gaps[splits_df_gaps$split_yes_no ==1, ] $no_splits)) / nrow(splits_df_gaps)


# --- Loop to detect merges

#fs1 -9-17
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs1_df$gaps_17_fs1)) {
  
  df <- subset(gaps_fs1_df, gaps_17_fs1 == i) 
  merges <- unique(df$gaps_9_fs1)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges1_917 <- merges_df
merges1_917$site <- as.factor(1)
merges1_917$timestep <- as.factor("9-17")

#fs1 17-21
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs1_df$gaps_21_fs1)) {
  
  df <- subset(gaps_fs1_df, gaps_21_fs1 == i) 
  merges <- unique(df$gaps_17_fs1)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges1_1721 <- merges_df
merges1_1721$site <- as.factor(1)
merges1_1721$timestep <- as.factor("17-21")


#fs2 -9-17
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs2_df$gaps_17_fs2)) {
  
  df <- subset(gaps_fs2_df, gaps_17_fs2 == i) 
  merges <- unique(df$gaps_9_fs2)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges2_917 <- merges_df
merges2_917$site <- as.factor(2)
merges2_917$timestep <- as.factor("9-17")

#fs2 17-21
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs2_df$gaps_21_fs2)) {
  
  df <- subset(gaps_fs2_df, gaps_21_fs2 == i) 
  merges <- unique(df$gaps_17_fs2)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges2_1721 <- merges_df
merges2_1721$site <- as.factor(2)
merges2_1721$timestep <- as.factor("17-21")


#fs3 -9-17
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs3_df$gaps_17_fs3)) {
  
  df <- subset(gaps_fs3_df, gaps_17_fs3 == i) 
  merges <- unique(df$gaps_9_fs3)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges3_917 <- merges_df
merges3_917$site <- as.factor(3)
merges3_917$timestep <- as.factor("9-17")

#fs1 17-21
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs3_df$gaps_21_fs3)) {
  
  df <- subset(gaps_fs3_df, gaps_21_fs3 == i) 
  merges <- unique(df$gaps_17_fs3)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges3_1721 <- merges_df
merges3_1721$site <- as.factor(3)
merges3_1721$timestep <- as.factor("17-21")



#fs4 -9-17
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs4_df$gaps_17_fs4)) {
  
  df <- subset(gaps_fs4_df, gaps_17_fs4 == i) 
  merges <- unique(df$gaps_9_fs4)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges4_917 <- merges_df
merges4_917$site <- as.factor(4)
merges4_917$timestep <- as.factor("9-17")

#fs4 17-21
merges_df <- as.data.frame(ncol(2), nrow(0))
for (i in unique(gaps_fs4_df$gaps_21_fs4)) {
  
  df <- subset(gaps_fs4_df, gaps_21_fs4 == i) 
  merges <- unique(df$gaps_17_fs4)
  
  if (length(merges) == 1 & merges == 0) {merge_info <- c(i,0,0)} # gap disappears (veg growth)
  if (length(merges) == 1 & merges != 0) {merge_info <- c(i,1,0)} # gap stays as one gap
  if (length(merges) == 2 & 0 %in% merges) {merge_info <- c(i,1,0)} # gap stays as one and one part goes to vegetation
  if (length(merges) == 2 & !(0 %in% merges)) {merge_info <- c(i,2,1)} # gap splits once (without parts ging to vegetation pool)
  if (length(merges) >2 & !(0 %in% merges) ) {merge_info <- c(i, length(merges),1)} # gap splits more than twice and nothing into veg pool
  if (length(merges) >2 &  0 %in% merges ) {merge_info <- c(i, (length(merges)-1),1)} # gap splits more than twice
  merges_df <- rbind(merges_df, merge_info)
  names(merges_df) <- c("gap_id", "no_merges", "merge_yes_no")
}

merges4_1721 <- merges_df
merges4_1721$site <- as.factor(4)
merges4_1721$timestep <- as.factor("17-21")

# -- calculate merge factors

merges_all <- rbind(merges1_917, merges1_1721, merges2_917, merges2_1721, merges3_917, merges3_1721, merges4_917, merges4_1721)

# -- calculate split factors

merges_all_noVeg <- merges_all[merges_all$gap_id != 0, ] # exclude the vegetation pool

merge_factor_all <- merges_all_noVeg %>%
  group_by(site, timestep) %>%
  summarise(n_gaps_merge = sum(merge_yes_no[merge_yes_no ==1]),
            sum_gap_merges = sum(no_merges[merge_yes_no ==1]),
            n_gaps = n(),
            merge_factor = sum_gap_merges/n_gaps_merge,
            merge_share = n_gaps_merge/ n_gaps,
            global_merge_factor = sum_gap_merges / n_gaps)


ggplot(merge_factor_all, aes(x=timestep, y=global_merge_factor, color=site, group=site)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "gobal merge factor")

ggplot(merge_factor_all, aes(x=timestep, y=merge_share, color=site, group=site)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "merge share")

ggplot(merge_factor_all, aes(x=timestep, y=merge_factor ,  group=site)) + 
  geom_point(size=3) + geom_line() + theme_minimal() +  labs(title = "merge factor") + facet_wrap(~site)


merges_df_gaps <- merges_df[merges_df$gap_id != 0, ] # exclude the vegetation pool
#splits_df_gaps[splits_df_gaps$split_yes_no ==1, ] # df of gap, which split
#sum(splits_df_gaps$no_splits) # sum of splits (but of all gaps, also those which don't split)
#sum(splits_df_gaps$split_yes_no == 1) # amount of gaps which split

# calculate split factor (sum of gap splits of gaps which split / number of gaps which split)
merge_factor <- (sum(merges_df_gaps[merges_df_gaps$merge_yes_no ==1, ] $no_merges)) / (sum(merges_df_gaps$merge_yes_no == 1))

# calculate share of gaps which split (number of gaps which split / number of gaps at t0)
merge_share <- (sum(merges_df_gaps$merge_yes_no == 1)) / nrow(merges_df_gaps)

# calculate global split factor (sum of gap splits of gaps which split / number of gaps at t0 )
global_merge_factor <- (sum(merges_df_gaps[merges_df_gaps$merge_yes_no ==1, ] $no_merges)) / nrow(merges_df_gaps)


# plot factors and compare between focus sites
# I can use the dataframes from the loops to identify the share of gaps, which go into the vegetation pool
# same for share of gaps, emerging from vegetation pool

# --- compare merge and split factor ---

merge_factor_all$split_or_merge <- as.factor("merge")
split_factor_all$split_or_merge <- as.factor("split")
names(merge_factor_all) <- c("site", "timestep", "n_merge_splits", "sum_merges_splits", "n_gaps", "factor", "share", "global_factor", "split_or_merge")
names(split_factor_all) <- c("site", "timestep", "n_merge_splits", "sum_merges_splits", "n_gaps", "factor", "share", "global_factor", "split_or_merge")

split_merges_all <- rbind(merge_factor_all, split_factor_all)

ggplot(split_merges_all, aes(x=timestep, y=global_factor, color=split_or_merge, group=split_or_merge)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "gobal merge and split factor") + facet_wrap(~site)

ggplot(split_merges_all, aes(x=timestep, y=share, color=split_or_merge, group=split_or_merge)) + 
  geom_point() + geom_line() + theme_minimal() +  labs(title = "share of gaps merging and splitting") + facet_wrap(~site)

ggplot(split_merges_all, aes(x=timestep, y=factor , color=split_or_merge, group=split_or_merge)) + 
  geom_point(size=3) + geom_line() + theme_minimal() +  labs(title = "merge and split factor") + facet_wrap(~site)

# calculate share of vegs which go to veg pool

#--- extract vegetation change per gap change class and plot

# --- extracting CHM DIFF for gap changes ---

GetChangeCHM <- function(chm1, chm2, gapchange) { #chm1 ist newer, chm2 is older
  diff <- chm1 - chm2  # create simple difference
  # gap_change_dir <- raster::mask(raster::raster(diff_class), gapchange)
  gap_change_dir <- mask(diff, gapchange)
  return(gap_change_dir)
}

diff_fs1_917 <- GetChangeCHM(chm17_fs1, chm9_fs1, gap_stack_fs1$gap_changes_fs1_917)
diff_fs1_1721 <- GetChangeCHM(chm21_fs1, chm17_fs1, gap_stack_fs1$gap_changes_fs1_1721)
# für alle sites und Zeitschritte berechnen

# 2009-2017
change_vertical_bins1 <- c(gap_stack_fs1$gap_change_dir1_917, diff_fs1_917, chm9_fs1 )
names(change_vertical_bins1) <- c("Change_class", "height_change", "height")
change_vertical_bins1 <- as.data.frame(change_vertical_bins1)
change_vertical_bins1<-change_vertical_bins1 %>% 
  mutate(Canopy_height_bins = cut(height, breaks = c(-1,1,2,3,4,5,10,20,30, 40)),
         time_step = as.factor("2009-2017"))

# 2017-2021
change_vertical_bins2 <- c(gap_stack_fs1$gap_change_dir1_1721, diff_fs1_1721, chm17_fs1 )
names(change_vertical_bins2) <- c("Change_class", "height_change", "height")
change_vertical_bins2 <- as.data.frame(change_vertical_bins2)
change_vertical_bins2<-change_vertical_bins2 %>% 
  mutate(Canopy_height_bins = cut(height, breaks = c(-1,1,2,3,4,5,10,20,30, 40)),
         time_step = as.factor("2017-2021"))

change_vertical_bins <- rbind(change_vertical_bins1, change_vertical_bins2)


change_vertical_bins["Change_class"][change_vertical_bins["Change_class"] == 1] <- "height loss"
change_vertical_bins["Change_class"][change_vertical_bins["Change_class"] == 2] <- "stable height"
change_vertical_bins["Change_class"][change_vertical_bins["Change_class"] == 3] <- "height gain"
change_vertical_bins$Change_class <- as.factor(change_vertical_bins$Change_class)

ggplot(change_vertical_bins, aes(x= height_change, y=Canopy_height_bins, col=Change_class)) + 
  geom_boxplot() + theme_minimal()# +facet_wrap(~Change_class)
ggplot(change_vertical_bins, aes(x= height_change, y=Canopy_height_bins)) +
  geom_boxplot() + theme_minimal()+facet_wrap(~Change_class)
# all three change classes
ggplot(change_vertical_bins, aes(x= height_change, y=Canopy_height_bins, col=time_step)) + 
  geom_boxplot() + theme_minimal()+facet_wrap(~Change_class)
# height gain only
ggplot(change_vertical_bins[ which(change_vertical_bins$Change_class=='height gain'), ], aes(x= height_change, y=Canopy_height_bins, col=time_step)) + 
  geom_boxplot() + theme_minimal()


# für beide Zeitschritte in einen DF
# für alle sites und Zeitschritte (innerhalb des plots zwischen Zeitschritten unterscheiden und mit facet-wrap zwischen sites)

# --- trace individual gaps ----

#focus site 1
gaps_fs1 <- c(gap_stack_fs1$gaps_9_fs1, gap_stack_fs1$gaps_17_fs1, gap_stack_fs1$gaps_21_fs1)
gaps_fs1_df <- as.data.frame(gaps_fs1, na.rm=FALSE)
gaps_fs1_df <- gaps_fs1_df[rowSums(is.na(gaps_fs1_df)) != ncol(gaps_fs1_df), ] # delete pixels without any gap at any moment in time

df <- gaps_fs1_df %>%
  make_long(gaps_9_fs1, gaps_17_fs1, gaps_21_fs1)
df[df == "NaN"] <- 999

# Chart 1
pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = FALSE)
pl <- pl + geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl <- pl +  theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(title = "Sankey diagram for gaps")
pl <- pl + labs(subtitle = "flow unit is pixels per gap")
#pl <- pl + labs(caption = "@techanswers88")
pl <- pl + labs(fill = 'Nodes')
pl

# focus site 2

gaps_fs2 <- c(gap_stack_fs2$gaps_9_fs2, gap_stack_fs2$gaps_17_fs2, gap_stack_fs2$gaps_21_fs2)
gaps_fs2_df <- as.data.frame(gaps_fs2, na.rm=FALSE)
gaps_fs2_df <- gaps_fs2_df[rowSums(is.na(gaps_fs2_df)) != ncol(gaps_fs2_df), ] # delete pixels without any gap at any moment in time

df2 <- gaps_fs2_df %>%
  make_long(gaps_9_fs2, gaps_17_fs2, gaps_21_fs2)
df2[df2 == "NaN"] <- 999

# Chart 1
pl2 <- ggplot(df2, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
)
pl2 <- pl2 +geom_sankey(flow.alpha = 0.5
                        , node.color = "black"
                        ,show.legend = FALSE)
pl2 <- pl2 +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl2 <- pl2 +  theme_bw()
pl2 <- pl2 + theme(legend.position = "none")
pl2 <- pl2 +  theme(axis.title = element_blank()
                    , axis.text.y = element_blank()
                    , axis.ticks = element_blank()  
                    , panel.grid = element_blank())
pl2 <- pl2 + scale_fill_viridis_d(option = "inferno")
pl2 <- pl2 + labs(title = "Sankey diagram for gaps")
pl2 <- pl2 + labs(subtitle = "flow unit is pixels per gap")
#pl <- pl + labs(caption = "@techanswers88")
pl2 <- pl2 + labs(fill = 'Nodes')
pl2

# focus site 3

gaps_fs3 <- c(gap_stack_fs3$gaps_9_fs3, gap_stack_fs3$gaps_17_fs3, gap_stack_fs3$gaps_21_fs3)
gaps_fs3_df <- as.data.frame(gaps_fs3,na.rm=FALSE)
gaps_fs3_df <- gaps_fs3_df[rowSums(is.na(gaps_fs3_df)) != ncol(gaps_fs3_df), ] # delete pixels without any gap at any moment in time

df3 <- gaps_fs3_df %>%
  make_long(gaps_9_fs3, gaps_17_fs3, gaps_21_fs3)

df3[df3 == "NaN"] <- 999

# Chart 1
pl3 <- ggplot(df3, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
)
pl3 <- pl3 +geom_sankey(flow.alpha = 0.5
                        , node.color = "black"
                        ,show.legend = FALSE)
pl3 <- pl3 +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl3 <- pl3 +  theme_bw()
pl3 <- pl3 + theme(legend.position = "none")
pl3 <- pl3 +  theme(axis.title = element_blank()
                    , axis.text.y = element_blank()
                    , axis.ticks = element_blank()  
                    , panel.grid = element_blank())
pl3 <- pl3 + scale_fill_viridis_d(option = "inferno")
pl3 <- pl3 + labs(title = "Sankey diagram for gaps")
pl3 <- pl3 + labs(subtitle = "flow unit is pixels per gap")
#pl <- pl + labs(caption = "@techanswers88")
pl3 <- pl3 + labs(fill = 'Nodes')
pl3

# focus site 4

gaps_fs4 <- c(gap_stack_fs4$gaps_9_fs4, gap_stack_fs4$gaps_17_fs4, gap_stack_fs4$gaps_21_fs4)
gaps_fs4_df <- as.data.frame(gaps_fs4,na.rm=FALSE)
gaps_fs4_df <- gaps_fs4_df[rowSums(is.na(gaps_fs4_df)) != ncol(gaps_fs4_df), ] 

df4 <- gaps_fs4_df %>%
  make_long(gaps_9_fs4, gaps_17_fs4, gaps_21_fs4)

df4[df4 == "NaN"] <- 999

# Chart 1
pl4 <- ggplot(df4, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
)
pl4 <- pl4 +geom_sankey(flow.alpha = 0.5
                        , node.color = "black"
                        ,show.legend = FALSE)
pl4 <- pl4 +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl4 <- pl4 +  theme_bw()
pl4 <- pl4 + theme(legend.position = "none")
pl4 <- pl4 +  theme(axis.title = element_blank()
                    , axis.text.y = element_blank()
                    , axis.ticks = element_blank()  
                    , panel.grid = element_blank())
pl4 <- pl4 + scale_fill_viridis_d(option = "inferno")
pl4 <- pl4 + labs(title = "Sankey diagram for gaps")
pl4 <- pl4 + labs(subtitle = "flow unit is pixels per gap")
#pl <- pl + labs(caption = "@techanswers88")
pl4 <- pl4 + labs(fill = 'Nodes')
pl4
