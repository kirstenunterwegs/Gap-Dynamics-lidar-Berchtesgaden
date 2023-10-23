########################################################
#
# Script for preparation of lidar data for CHM creation (2017)
#
# author: Dirk Pflugmacher (Humboldt Universität zu Berlin)
# last revised: 20.10.2023
#
#######################################################

# --- libaries ---

library(lidR)
library(sf)

# --- wd ----

root <- "f:/spring/berchtesgaden/lidar/2017"

#--- select ----

cat <- read_sf(file.path(root, "catalog_utm_selected.gpkg"))
cat <- basename(cat$filename)

p_all <- file.path(root, "las_epsg5678")
p_outsider <- file.path(root, "las_epsg5678_out")

for (fn in cat) {
  fn_in <- file.path(p_all, fn)
  fn_out <- file.path(p_outsider, fn)
  file.rename(fn_out, fn_in)
}

#--- fix false tile ----

fn_corrupt <- "i:/Fonda/workspace/berchtesgaden/lidar/2017/las_epsg5678/4562_5276_all.las"
las_corrupt <- readLAS(fn_corrupt)

las_fix <- las_corrupt[las_corrupt$gpstime>0,]

file.rename(fn_corrupt, paste0(substr(fn_corrupt, 1, nchar(fn_corrupt)-3), "corrupt"))
writeLAS(las_fix, fn_corrupt)

#--- reproject ----

p_all <- file.path(root, "las_epsg5678")
p_las <- file.path(root, "las_reproj")
fns <- list.files(p_all, ".las$")

# install pdal via conda
# pdal translate -i 4559_5276_all.las -o 4559_5276_utm.las -f filters.reprojection --filters.reprojection.in_srs="EPSG:5678" --filters.reprojection.out_srs="EPSG:25832"
# for %f in (*.las) do pdal translate -i %f -o %f_utm.las -f filters.reprojection --filters.reprojection.in_srs="EPSG:5678" --filters.reprojection.out_srs="EPSG:25832"

Sys.setenv(PATH="%PATH%;C:/Users/pflugmad/Miniconda3/Library/bin")
Sys.setenv(GDAL_DATA="C:/Users/pflugmad/Miniconda3/Library/share/gdal")
Sys.setenv(PROJ_LIB="C:/Users/pflugmad/Miniconda3/Library/share/proj")

for (fn in fns) {
  
  in_file <- file.path(p_all, fn)
  out_file <- file.path(p_las, fn)
  
  if (file.exists(out_file)) {
    message(paste("Skipped", out_file))
  } else {
    system2(command = "pdal",
            args = c('translate', '-i', in_file, '-o', out_file, '-f', 'filters.reprojection', 
                     '--filters.reprojection.in_srs="EPSG:5678"', 
                     '--filters.reprojection.out_srs="EPSG:25832"'))
    message(out_file)
  }

}


#--- retile --------------

p_reproj <- file.path(root, "las_reproj")
p_las <- file.path(root, "las")


ctg_all <- readLAScatalog(p_reproj)
write_sf(st_as_sf(ctg_all), file.path(root, "berchtesgaden_2017_reproj.gpkg")) 


ctg <- ctg_all # [ctg_all$Max.X < 792000 & ctg_all$Max.Y > 5277000,]
lidR:::catalog_laxindex(ctg_all)

opt_output_files(ctg) <- paste0(p_las, "/{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 1000 # retile to 250 m
opt_chunk_alignment(ctg) <- c(784000, 5275000)
ctg_rt <- catalog_retile(ctg) # apply retile
# plot(ctg_rt) # some plotting

ctgn <- readLAScatalog(p_las)
write_sf(st_as_sf(ctgn), file.path(root, "berchtesgaden_2017_tiles.gpkg")) 

lidR:::catalog_laxindex(p_las)

#--- testing ----

p <- "i:/Fonda/workspace/berchtesgaden/lidar/2017/las/"
p2009_dk4 <- "i:/Fonda/workspace/berchtesgaden/lidar/2009/las_dk4/"
 
ctg <- readLAScatalog(p)

plot(ctg)

g <- st_as_sf(ctg)
st_crs(g) <- 5678
write_sf(g, file.path(p, "catalog_5678.gpkg"))


p2009_dk4 <- "i:/Fonda/workspace/berchtesgaden/lidar/2009/las_dk4/"
ctg09_dk4 <- readLAScatalog(p2009_dk4)

plot(ctg09_dk4)

sf09_dk4 <- st_as_sf(ctg09_dk4)
st_crs(sf09_dk4) <- 5678
write_sf(sf09_dk4, file.path(p2009_dk4, "catalog_5678.gpkg"))
sf09_utm <- st_transform(sf09_dk4, crs=25832)
write_sf(sf09_utm, file.path(p2009_dk4, "catalog_25832.gpkg"))
