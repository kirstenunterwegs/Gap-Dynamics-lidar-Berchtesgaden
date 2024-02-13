###########################################################################

# Identifying number of and area covered by gaps in reserach area per year

##########################################################################


library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)


wd <- "C:/Users/ge92vuh/Documents/MA_gap_dynamics/data/"
setwd(wd)

# --- load layers ----

gaps2009 <- rast("processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop gaps to CHM 
gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE) 

gap.stack <- c(gaps2009, gaps2017, gaps2021)

#prepare layers to limit analysis to core zone, below 1800m and with forest type information

# --- load NP information 

foresttype <- rast("processed/environment_features/forest_type2020_reclass_1m.tif")
management <- vect("F:/Projects/CanopyDynamicsBDG/data/NP_data/npb_zonierung_22_epsg25832.shp")
elevation.below1800 <- rast("processed/environment_features/elevation_below1800_200steps.tif")


# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

#mask down to reserach area
gap.stack <- mask(gap.stack, core.zone)
gap.stack <- mask(gap.stack, elevation.below1800)
gap.stack <- mask(gap.stack, foresttype)

writeRaster(gap.stack, "processed/gaps_final/gaps_masked_reserach_area_paper.tif")

# convert masked gap stack to data frame & write to disk
gap.stack.df <- as.data.frame(gap.stack)
saveRDS(gap.stack.df,"processed/gaps_final/gaps_masked_reserach_area_paper.df.rds" )

gap.stack.df <- readRDS("processed/gaps_final/gaps_masked_reserach_area_paper.df.rds")
names(gap.stack.df) <- c("gaps9", "gaps17", "gaps21")

# --- summarize number & area of gaps:

# Gather the data into long format
gathered_df <- gap.stack.df %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "gap_id") %>%
  filter(!if_all(.cols = everything(), is.na))%>% # Drop rows with all NAs
  count(year, gap_id) %>%
  filter(n >= 400) %>% # filter out gaps < 400m2, which emerged due to cropping of gaps layers to reserach area
  filter(!is.nan(gap_id))# Remove rows where gap_id is NaN

# overall gap area
sum(gathered_df$n / 10000) # 2730.522 ha

# Count the occurrences of each value in each column
gap_counts <- gathered_df %>%
  group_by(year) %>%
  summarize(number_gaps = length(unique(gap_id)))

# overall number of gaps
sum(gap_counts$number_gaps) # 11331 



# # Extract unique values for gaps9
# unique_values_gaps9 <- unique(gathered_df$value[gathered_df$year == "gaps9"])
# 
# # Extract unique values for gaps17
# unique_values_gaps17 <- unique(gathered_df$value[gathered_df$year == "gaps17"])
# 
# # Extract unique values for gaps21
# unique_values_gaps21 <- unique(gathered_df$value[gathered_df$year == "gaps21"])

