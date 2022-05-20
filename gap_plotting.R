####################################################################
# Plotting gap and gap change detetction
# author: Kirsten Krüger
###################################################################

#---- libaries 

library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(sp)

wd <- "i:/Fonda/workspace/berchtesgaden/gaps/"
setwd(wd)

#---- load CHMs, gap layers and polygons ------

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

# load gap layer stacks

gap_stack1_100 <- rast("gap_layers_fs1_100.tif")
gap_stack2_100 <- rast("gap_layers_fs2_100.tif")
gap_stack3_100 <- rast("gap_layers_fs3_100.tif")
gap_stack4_100 <- rast("gap_layers_fs4_100.tif")

gap_stack1_250 <- rast("gap_layers_fs1_250.tif")
gap_stack2_250 <- rast("gap_layers_fs2_250.tif")
gap_stack3_250 <- rast("gap_layers_fs3_250.tif")
gap_stack4_250 <- rast("gap_layers_fs4_250.tif")

gap_stack1_400 <- rast("gap_layers_fs1_400.tif")
gap_stack2_400 <- rast("gap_layers_fs2_400.tif")
gap_stack3_400 <- rast("gap_layers_fs3_400.tif")
gap_stack4_400 <- rast("gap_layers_fs4_400.tif")

#load gap polygons 

chm_names <- list("chm9_fs1", "chm17_fs1", "chm21_fs1", "chm9_fs2", "chm17_fs2", "chm21_fs2",
                  "chm9_fs3", "chm17_fs3", "chm21_fs3","chm9_fs4", "chm17_fs4", "chm21_fs4")

polygons_100 <- list()
polygons_100 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_100_", n ,"/","gaps_polygons_100_", n ,".shp", sep=""))
})
names(polygons_100) <- chm_names


polygons_250 <- list()
polygons_250 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_250_", n ,"/","gaps_polygons_250_", n ,".shp", sep=""))
})
names(polygons_250) <- chm_names


polygons_400 <- list()
polygons_400 <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_400_", n ,"/","gaps_polygons_400_", n ,".shp", sep=""))
})
names(polygons_400) <- chm_names


#--- theme preperation

myTheme <- viridisTheme()
myTheme$panel.background$col = 'gray' 


# ---- Focus site 1

#plotting

site1_9_100 <-levelplot(chm9_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 100"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_100$chm9_fs1, "Spatial"), fill='red'))
site1_17_100 <- levelplot(chm17_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm17_fs1, "Spatial"), fill='red'))
site1_21_100 <- levelplot(chm21_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm21_fs1, "Spatial"), fill='red'))
site1_change_917_100 <- levelplot(gap_stack1_100$gap_change_dir1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 100"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site1_change_1721_100 <- levelplot(gap_stack1_100$gap_change_dir1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 100"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site1_9_250 <-levelplot(chm9_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 250"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_250$chm9_fs1, "Spatial"), fill='red'))
site1_17_250 <- levelplot(chm17_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm17_fs1, "Spatial"), fill='red'))
site1_21_250 <- levelplot(chm21_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm21_fs1, "Spatial"), fill='red'))
site1_change_917_250 <- levelplot(gap_stack1_250$gap_change_dir1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 250"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site1_change_1721_250 <- levelplot(gap_stack1_250$gap_change_dir1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 250"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))             

site1_9_400 <-levelplot(chm9_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 400"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_400$chm9_fs1, "Spatial"), fill='red'))
site1_17_400 <- levelplot(chm17_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm17_fs1, "Spatial"), fill='red'))
site1_21_400 <- levelplot(chm21_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm21_fs1, "Spatial"), fill='red'))
site1_change_917_400 <- levelplot(gap_stack1_400$gap_change_dir1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 400"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site1_change_1721_400 <- levelplot(gap_stack1_400$gap_change_dir1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 400"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))             


# ---- Focus site 2

#plotting

site2_9_100 <-levelplot(chm9_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 100"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_100$chm9_fs2, "Spatial"), fill='red'))
site2_17_100 <- levelplot(chm17_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm17_fs2, "Spatial"), fill='red'))
site2_21_100 <- levelplot(chm21_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm21_fs2, "Spatial"), fill='red'))
site2_change_917_100 <- levelplot(gap_stack2_100$gap_change_dir2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 100"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site2_change_1721_100 <- levelplot(gap_stack2_100$gap_change_dir2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 100"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site2_9_250 <-levelplot(chm9_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 250"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_250$chm9_fs2, "Spatial"), fill='red'))
site2_17_250 <- levelplot(chm17_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm17_fs2, "Spatial"), fill='red'))
site2_21_250 <- levelplot(chm21_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm21_fs2, "Spatial"), fill='red'))
site2_change_917_250 <- levelplot(gap_stack2_250$gap_change_dir2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 250"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site2_change_1721_250 <- levelplot(gap_stack2_250$gap_change_dir2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 250"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

site2_9_400<-levelplot(chm9_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 400"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_400$chm9_fs2, "Spatial"), fill='red'))
site2_17_400 <- levelplot(chm17_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm17_fs2, "Spatial"), fill='red'))
site2_21_400 <- levelplot(chm21_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm21_fs2, "Spatial"), fill='red'))
site2_change_917_400 <- levelplot(gap_stack2_400$gap_change_dir2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 400"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site2_change_1721_400 <- levelplot(gap_stack2_400$gap_change_dir2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 400"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

# --- Focus Site 3 ------

#plotting
site3_9_100 <-levelplot(chm9_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 100"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_100$chm9_fs3, "Spatial"), fill='red'))
site3_17_100 <- levelplot(chm17_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm17_fs3, "Spatial"), fill='red'))
site3_21_100 <- levelplot(chm21_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm21_fs3, "Spatial"), fill='red'))
site3_change_917_100 <- levelplot(gap_stack3_100$gap_change_dir3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 100"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site3_change_1721_100 <- levelplot(gap_stack3_100$gap_change_dir3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 100"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site3_9_250 <-levelplot(chm9_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 250"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_250$chm9_fs3, "Spatial"), fill='red'))
site3_17_250 <- levelplot(chm17_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm17_fs3, "Spatial"), fill='red'))
site3_21_250 <- levelplot(chm21_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm21_fs3, "Spatial"), fill='red'))
site3_change_917_250 <- levelplot(gap_stack3_250$gap_change_dir3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 250 "),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site3_change_1721_250 <- levelplot(gap_stack3_250$gap_change_dir3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 250 "),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

site3_9_400 <-levelplot(chm9_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 400"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_400$chm9_fs3, "Spatial"), fill='red'))
site3_17_400 <- levelplot(chm17_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm17_fs3, "Spatial"), fill='red'))
site3_21_400 <- levelplot(chm21_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm21_fs3, "Spatial"), fill='red'))
site3_change_917_400 <- levelplot(gap_stack3_400$gap_change_dir3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 400 "),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site3_change_1721_400 <- levelplot(gap_stack3_400$gap_change_dir3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 400 "),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     
# --- Focus Site 4 ------


site4_9_100 <-levelplot(chm9_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 100"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_100$chm9_fs4, "Spatial"), fill='red'))
site4_17_100 <- levelplot(chm17_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm17_fs4, "Spatial"), fill='red'))
site4_21_100 <- levelplot(chm21_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 100"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_100$chm21_fs4, "Spatial"), fill='red'))
site4_change_917_100 <- levelplot(gap_stack4_100$gap_change_dir4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 100"),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site4_change_1721_100 <- levelplot(gap_stack4_100$gap_change_dir4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 100"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site4_9_250 <-levelplot(chm9_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 250"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_250$chm9_fs4, "Spatial"), fill='red'))
site4_17_250 <- levelplot(chm17_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm17_fs4, "Spatial"), fill='red'))
site4_21_250 <- levelplot(chm21_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 250"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_250$chm21_fs4, "Spatial"), fill='red'))
site4_change_917_250 <- levelplot(gap_stack4_250$gap_change_dir4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 250 "),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site4_change_1721_250 <- levelplot(gap_stack4_250$gap_change_dir4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 250"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

site4_9_400 <-levelplot(chm9_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 min 400"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_400$chm9_fs4, "Spatial"), fill='red'))
site4_17_400 <- levelplot(chm17_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm17_fs4, "Spatial"), fill='red'))
site4_21_400 <- levelplot(chm21_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 min 400"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_400$chm21_fs4, "Spatial"), fill='red'))
site4_change_917_400 <- levelplot(gap_stack4_400$gap_change_dir4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change 9-17 min 400 "),
                                  colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site4_change_1721_400 <- levelplot(gap_stack4_400$gap_change_dir4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change 17-21 min 400"),
                                   colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

