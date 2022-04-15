####################################################################
# Plotting gap and gap change detetction
# author: Kirsten Krüger
###################################################################

#---- libaries & workingspace

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

gap_stack_fs1 <- rast("gap_layers_fs1.tif")
gap_stack_fs2 <- rast("gap_layers_fs2.tif")
gap_stack_fs3 <- rast("gap_layers_fs3.tif")
gap_stack_fs4 <- rast("gap_layers_fs4.tif")

gap_stack_fs1_erosion <- rast("gap_layers_fs1_erosion.tif")
gap_stack_fs2_erosion <- rast("gap_layers_fs2_erosion.tif")
gap_stack_fs3_erosion <- rast("gap_layers_fs3_erosion.tif")
gap_stack_fs4_erosion <- rast("gap_layers_fs4_erosion.tif")


#load gap polygons 

chm_names <- list("chm9_fs1", "chm17_fs1", "chm21_fs1", "chm9_fs2", "chm17_fs2", "chm21_fs2",
                  "chm9_fs3", "chm17_fs3", "chm21_fs3","chm9_fs4", "chm17_fs4", "chm21_fs4")

polygons <- list()
polygons <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_", n ,"/","gaps_polygons_", n ,".shp", sep=""))
})
names(polygons) <- chm_names


polygons_erosion <- list()
polygons_erosion <- lapply(chm_names, function(n) {
  vect(paste("gaps_polygons_erosion_", n ,"/","gaps_polygons_erosion_", n ,".shp", sep=""))
})
names(polygons_erosion) <- chm_names


# ---- Focus site 1
myTheme <- viridisTheme()
myTheme$panel.background$col = 'gray' 

#plotting
site1_9 <-levelplot(chm9_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 "), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons$chm9_fs1, "Spatial"), fill='red'))
site1_17 <- levelplot(chm17_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm17_fs1, "Spatial"), fill='red'))
site1_21 <- levelplot(chm21_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm21_fs1, "Spatial"), fill='red'))
site1_change_917 <- levelplot(gap_stack_fs1$gap_change_dir1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 "),
                              colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site1_change_1721 <- levelplot(gap_stack_fs1$gap_change_dir1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21"),
                               colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site1_9_e <-levelplot(chm9_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 erosion filter"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_erosion$chm9_fs1, "Spatial"), fill='red'))
site1_17_e <- levelplot(chm17_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm17_fs1, "Spatial"), fill='red'))
site1_21_e <- levelplot(chm21_fs1, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm21_fs1, "Spatial"), fill='red'))
site1_change_917_erosion <- levelplot(gap_stack_fs1_erosion$gap_change_dir1_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 erosion filter "),
                                      colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site1_change_1721_erosion <- levelplot(gap_stack_fs1_erosion$gap_change_dir1_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21 erosion filter "),
                                       colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))             

# ---- Focus site 2

#plotting

site2_9 <-levelplot(chm9_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 "), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons$chm9_fs2, "Spatial"), fill='red'))
site2_17 <- levelplot(chm17_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm17_fs2, "Spatial"), fill='red'))
site2_21 <- levelplot(chm21_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm21_fs2, "Spatial"), fill='red'))
site2_change_917 <- levelplot(gap_stack_fs2$gap_change_dir2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 "),
                              colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site2_change_1721 <- levelplot(gap_stack_fs2$gap_change_dir2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21"),
                               colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site2_9_e <-levelplot(chm9_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 erosion filter"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_erosion$chm9_fs2, "Spatial"), fill='red'))
site2_17_e <- levelplot(chm17_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm17_fs2, "Spatial"), fill='red'))
site2_21_e <- levelplot(chm21_fs2, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm21_fs2, "Spatial"), fill='red'))
site2_change_917_erosion <- levelplot(gap_stack_fs2_erosion$gap_change_dir2_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 erosion filter "),
                                      colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site2_change_1721_erosion <- levelplot(gap_stack_fs2_erosion$gap_change_dir2_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21 erosion filter "),
                                       colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

# --- Focus Site 3 ------

#plotting
site3_9 <-levelplot(chm9_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 "), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons$chm9_fs3, "Spatial"), fill='red'))
site3_17 <- levelplot(chm17_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm17_fs3, "Spatial"), fill='red'))
site3_21 <- levelplot(chm21_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm21_fs3, "Spatial"), fill='red'))
site3_change_917 <- levelplot(gap_stack_fs3$gap_change_dir3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 "),
                              colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site3_change_1721 <- levelplot(gap_stack_fs3$gap_change_dir3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21"),
                               colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site3_9_e <-levelplot(chm9_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 erosion filter"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_erosion$chm9_fs3, "Spatial"), fill='red'))
site3_17_e <- levelplot(chm17_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm17_fs3, "Spatial"), fill='red'))
site3_21_e <- levelplot(chm21_fs3, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm21_fs3, "Spatial"), fill='red'))
site3_change_917_erosion <- levelplot(gap_stack_fs3_erosion$gap_change_dir3_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 erosion filter "),
                                      colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site3_change_1721_erosion <- levelplot(gap_stack_fs3_erosion$gap_change_dir3_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21 erosion filter "),
                                       colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

# --- Focus Site 4 ------


site4_9 <-levelplot(chm9_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 "), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons$chm9_fs4, "Spatial"), fill='red'))
site4_17 <- levelplot(chm17_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm17_fs4, "Spatial"), fill='red'))
site4_21 <- levelplot(chm21_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 "), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons$chm21_fs4, "Spatial"), fill='red'))
site4_change_917 <- levelplot(gap_stack_fs4$gap_change_dir4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 "),
                              colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site4_change_1721 <- levelplot(gap_stack_fs4$gap_change_dir4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21"),
                               colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))

site4_9_e <-levelplot(chm9_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2009 erosion filter"), scales=list(draw=FALSE)) + # CHM and gaps
  layer(sp.polygons(as(polygons_erosion$chm9_fs4, "Spatial"), fill='red'))
site4_17_e <- levelplot(chm17_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2017 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm17_fs4, "Spatial"), fill='red'))
site4_21_e <- levelplot(chm21_fs4, margin = FALSE, par.setting =GrTheme(), main= list("Gaps 2021 erosion filter"), scales=list(draw=FALSE)) + 
  layer(sp.polygons(as(polygons_erosion$chm21_fs4, "Spatial"), fill='red'))
site4_change_917_erosion <- levelplot(gap_stack_fs4_erosion$gap_change_dir4_917, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 9-17 erosion filter "),
                                      colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))
site4_change_1721_erosion <- levelplot(gap_stack_fs4_erosion$gap_change_dir4_1721, margin = FALSE, par.setting =myTheme, main= list("Gap change direction 17-21 erosion filter "),
                                       colorkey=list(at=seq(0, 3, 1),labels=list(at=c(1, 2, 3), labels=c("loss", "steady", "gain"))), scales=list(draw=FALSE))     

