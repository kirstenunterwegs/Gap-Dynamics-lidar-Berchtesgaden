########################################################
#
# Script for preparation of lidar data for CHM creation (2009)
#
# author: Dirk Pflugmacher (Humboldt Universität zu Berlin)
# last revised: 20.10.2023
#
#######################################################

# --- libaries ---

library(lidR)
library(sf)
library(dplyr)

# --- working directory and data loading ---

root <- "i:/Fonda/workspace/berchtesgaden/lidar/2009"

p_reproj_sued <- file.path(root, "las_reproj_sued")
p_reproj_mitte <- file.path(root, "las_reproj_mitte")

p_tiles_sued <- file.path(root, "las_tiles_sued")
p_tiles_mitte <- file.path(root, "las_tiles_mitte")

ascii_sued <- list.files(file.path(root, "ascii_sued"), "fpl$", full.names = T)
ascii_mitte <- list.files(file.path(root, "ascii_mitte"), "fpl$", full.names = T)


#--- functions ---

ascii2las <- function(fpl_file, outpath=NULL, proj=5678, out_proj=25832) {
  
  # proj <- 31468
  # proj <- 5678 # gk4
  # proj <- 25832 # utm 32 with etrs datum
  
  lpl_file <- paste0(substr(fpl_file, 1, nchar(fpl_file)-4), '.lpl')
  las_file <- paste0(substr(fpl_file, 1, nchar(fpl_file)-4), '.las')
  
  if (!is.null(outpath)) las_file <- file.path(outpath, basename(las_file))
  
  if (!file.exists(lpl_file)) {
    message(paste("Missing: ", lpl_file))
    stop()
  }

  if (!file.exists(las_file)) {
    
    message(paste("Process: ", lpl_file))
    
    ds1 <- read.table(fpl_file, header=F)
    names(ds1) <- c('X', 'Y', 'Z', "Classification")
    ds1$ReturnNumber <- as.integer(1)
    ds1$NumberOfReturns <- as.integer(2)
    
    ds2 <- read.table(lpl_file, header=F)
    ds2 <- ds2[, 1:4]
    names(ds2) <- c('X', 'Y', 'Z', "Classification")
    ds2$ReturnNumber <- as.integer(2)
    ds2$NumberOfReturns <- as.integer(2)
    
    ds <- rbind(ds1, ds2)
    rm(ds1, ds2)
    
    if (!is.null(out_proj)) {
      ds <- sf::st_as_sf(ds, coords=c("X", "Y"), crs=proj)
      ds <- sf::st_transform(ds, crs=out_proj) 
      ds <- dplyr::bind_cols(as.data.frame(st_coordinates(ds)), st_drop_geometry(ds))
    } else {
      out_proj <- proj
    }
    
    ascii_class <- ds$Classification
    
    # ds$Classification[ascii_class == 0] <- as.integer(0) # First returns were never classified
    ds$Classification[ascii_class == 1] <- as.integer(2) # Bodenpunkt (sicher)
    ds$Classification[ascii_class == 2] <- as.integer(1) # Bodenpunkt (unsicher)
    ds$Classification[ascii_class == 3] <- as.integer(3) # Objektpunkt, Vegetationspunkt
    ds$Classification[ascii_class == 4] <- as.integer(6) # Laserpunkte innerhalb eines Gebäudes (aus DFK) mit umgebenden Puffer (z.Zt. 1.5 m)

    
    ld <- lidR::LAS(data=ds) #, header = header, check = TRUE)

    lidR::epsg(ld) <- as.numeric(out_proj)
    
    lidR::writeLAS(ld, las_file)
    
    # ds.sf <- st_as_sf(ds, coords=c('X','Y'), crs=out_proj)
    # write_sf(ds.sf, sub("[.]fpl", ".shp", fn_xyz))
    
    message(las_file)
  } else {
    message(paste('Skipped', las_file))
  }
}


#--- compile original tile codes -----

bsued <- sub(".fpl", "", basename(ascii_sued))
bmitt <- sub(".fpl", "", basename(ascii_mitte))

tiles <- data.frame(code=bmitt, region="mitte")
tiles <- rbind(tiles, data.frame(code=bsued, region="sued"))
write.csv(tiles, file.path(root, "berchtesgaden_2009_tiles_ascii.csv"), row.names=F, quote=F)


#--- import and reproject selected ------

require(foreach)

tns <- c("4560_5277", "4561_5277", "4560_5276", "4561_5276")
i <- NULL
for (tn in tns) i <- c(i, grep(tn, ascii_sued))
ascii_sued <- ascii_sued[i]

# ascii2las("z:/Fonda/workspace/berchtesgaden/lidar/2009/ascii_sued/4581_5284.fpl")

foreach(i = ascii_sued, .combine=c, .packages=c("sf", "lidR", "dplyr")) %do% {
  res <- ascii2las(i, outpath=p_reproj_sued)
}


foreach(i = ascii_mitte, .combine=c, .packages=c("sf", "lidR", "dplyr")) %do% {
  ascii2las(i, outpath=p_reproj_mitte)
}



#--- import and reproject (parallel) ------

require(foreach)
require(doParallel)



cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

foreach(i = ascii_sued, .combine=c, .packages=c("sf", "lidR", "dplyr")) %dopar% {
  ascii2las(i, outpath=p_reproj_sued, proj=31468)
}

parallel::stopCluster(cl)


#--- select ----

cat <- read_sf(file.path(root, "berchtesgaden_2009_tiles_ascii_selected.gpkg"))
cat <- basename(cat$filename)

p_all <- file.path(root, "las_reproj_sued")
p_outsider <- file.path(root, "las_reproj_sued_out")

for (fn in cat) {
  fn_in <- file.path(p_all, fn)
  fn_out <- file.path(p_outsider, fn)
  file.rename(fn_out, fn_in)
}


#--- retile --------------


ctg <- readLAScatalog(file.path(root, "las_reproj_sued"))

opt_output_files(ctg) <- paste0(file.path(root, "las"), "/{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 1000 # retile to 250 m
opt_chunk_alignment(ctg) <- c(784000, 5275000)
ctg_rt <- catalog_retile(ctg) # apply retile
# plot(ctg_rt) # some plotting

ctgn <- readLAScatalog(file.path(root, "las"))
g <- st_as_sf(ctgn)
write_sf(g, file.path(root, "berchtesgaden_2009_tiles.gpkg")) 


#--- density --------------

ctg <- readLAScatalog(file.path(root, "las"))

dens <- grid_density(ctg, 10)

raster::writeRaster(dens, file.path(root, "grid_density.tif"), format="GTiff")

