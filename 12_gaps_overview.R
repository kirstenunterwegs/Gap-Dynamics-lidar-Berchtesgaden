
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
# +++ add closed forest for masking!!!

# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

#mask down to reserach area
gap.stack <- mask(gap.stack, core.zone)
gap.stack <- mask(gap.stack, elevation.below1800)
gap.stack <- mask(gap.stack, foresttype)

writeRaster(gap.stack, "processed/gaps_final/gaps_masked_reserach_area_paper.tif")

# filter out gaps < 400m2 ---------------------------------------------------------!!!!!!!!!!!!!!!!!!!!!

n_9 <- unique(gap.stack$berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight) 
n_17 <- unique(gap.stack$berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight)
n_21 <- unique(gap.stack$berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight)

n_gaps <- length(n_9$berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight) + 
          length(n_17$berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight)+ 
          length(n_21$berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight)

gap.stack.df <- as.data.frame(gap.stack)
saveRDS(gap.stack.df,"processed/gaps_final/gaps_masked_reserach_area_paper.df.rds" )

gap.stack.df <- readRDS("processed/gaps_final/gaps_masked_reserach_area_paper.df.rds")
names(gap.stack.df) <- c("gaps9", "gaps17", "gaps21")

# --- summarize the size of gaps:

# Gather the data into long format
gathered_df <- gap.stack.df %>%
  pivot_longer(cols = everything(), names_to = "column_name", values_to = "value") %>%
  drop_na()  # Drop rows with NA values

# Count the occurrences of each value in each column
value_counts <- gathered_df %>%
  group_by(column_name, value) %>%
  summarise(count = n())

sum(value_counts$count / 10000) # 2745.052 ha


# Extract unique values for gaps9
unique_values_gaps9 <- unique(gathered_df$value[gathered_df$column_name == "gaps9"])

# Extract unique values for gaps17
unique_values_gaps17 <- unique(gathered_df$value[gathered_df$column_name == "gaps17"])

# Extract unique values for gaps21
unique_values_gaps21 <- unique(gathered_df$value[gathered_df$column_name == "gaps21"])

