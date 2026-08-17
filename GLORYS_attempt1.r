###Adding GLORYS and MODIS data to mark-recap data##
# -------------------------------------------------------------------------
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# -------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(CopernicusMarine)
library(rerddap)
library(rerddapXtracto)
library(CFtime)

# -------------------------------------------------------------------------
# 2. LOAD AND CLEAN THE DATASET
# -------------------------------------------------------------------------
# Load your dataset
df_mr <- read_csv("Data/NCEAS_Cobia_Recaps_202605-1778080360235.csv")
df_ac <- read_csv('/Users/rodemann/Data/detections_2plus_202606.csv')
head(df_ac)

df_ac <- df_ac %>% dplyr::filter(decimalLongitude >= -80.9)
# Ensure date formatting is correct and remove rows with missing coordinates
df_clean <- df_mr %>%
  mutate(datecollected = as.Date(datecollected)) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(datecollected)) %>% 
  mutate(cat_num = paste0(catalognumber...13, '-', interaction_order))

write.csv(file = 'Data/MR_clean.csv', df_clean)

df_ai <- df_ac %>% dplyr::slice(1:1000)
write.csv(file = 'Data/ac_ai.csv', df_ai)

df_ac_clean <- df_ac %>%
  # 1. Parse the datetime string into a usable datetime object
  mutate(datetime = ymd_hms(dateCollectedUTC)) %>%
  
  # 2. Sort the data by fish identifier and datetime for accurate chronological ordering
  arrange(catalogNumber, datetime) %>%

  dplyr::filter(!is.na(datetime)) %>% 
  
  # 3. Group by individual fish to calculate the interaction order
  group_by(catalogNumber) %>%
  mutate(
    interaction_order = row_number(),
    # Extract date and time into separate columns
    datecollected = as.Date(datetime),
    timeofday = format(datetime, "%H:%M:%S"),
    # Add an empty length column, as length isn't natively in this detections table
    length = NA 
  ) %>%
  ungroup() %>%
  
  # 4. Rename columns to map to the MR_clean standard
  rename(
    project_shortcode = collectionCode,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  ) %>%
  
  # 5. Select only the requested columns
  select(
    project_shortcode,
    catalogNumber, # Represents the individual fish
    length,
    latitude,
    longitude,
    datecollected,
    timeofday,
    interaction_order
  ) 

# Preview the transformed dataset
head(df_ac_clean)

df_mr_subset <- df_mr %>%
  # Rename the identifier column to match the acoustic dataset (adjust if your primary ID is 'cat_num')
  rename(catalogNumber = catalognumber...13) %>% 
  
  # Select only the matching columns
  select(
    project_shortcode,
    catalogNumber,
    length,
    latitude,
    longitude,
    datecollected,
    timeofday,
    interaction_order
  ) %>%
  
  # Add the source identifier
  mutate(data_source = "mark-recapture") %>% 
  mutate(datecollected = as.Date(datecollected))

# 3. Add the source identifier to the clean acoustic dataset
df_ac_ready <- df_ac_clean %>%
  mutate(data_source = "acoustic") %>% 
  mutate(timeofday = hms::as_hms(timeofday))
head(df_mr_subset)
head(df_ac_ready)
# 4. Combine the two datasets together
combined_dataset <- bind_rows(df_mr_subset, df_ac_ready) %>%
  # Optional: arrange by fish ID and date to keep an individual's timeline organized
  arrange(catalogNumber, datecollected, timeofday)

# Preview the final combined dataset
head(combined_dataset)


# -------------------------------------------------------------------------
# 3. EXTRACT MODIS AQUA CHLOROPHYLL (Via NOAA ERDDAP)
# -------------------------------------------------------------------------
print("Starting MODIS Aqua Extraction...")

# MODIS Aqua started providing data around July 2002
df_modis <- combined_dataset %>% 
  filter(datecollected >= as.Date("2003-01-01")) %>% 
  filter(datecollected <= as.Date('2022-07-23')) %>% 
  dplyr::filter(!is.na(longitude)) %>% 
  dplyr::filter(!is.na(latitude))

# Define the NOAA ERDDAP dataset ID for MODIS Aqua Level 3 (Daily)
modis_info <- rerddap::info('erdMH1chla1day')

# Extract the chlorophyll data matching the exact x, y, and t coordinates
# We use a tiny search radius (xlen/ylen = 0.02 degrees) to get the nearest pixel
chl_extract <- rxtracto(
  dataInfo = modis_info,
  parameter = 'chlorophyll',
  tcoord = df_modis$datecollected,
  xcoord = df_modis$longitude,
  ycoord = df_modis$latitude,
  xlen = 0.02, 
  ylen = 0.02
)

#doesn't work for huge dataset. do it in batches
library(rerddapXtracto)
library(future)
library(furrr)
library(purrr)

# Reduce workers to 2 to prevent ERDDAP connection drops
plan(multisession, workers = 8)

# Wrap rxtracto so failed chunks return NULL instead of killing the worker process
safe_rxtracto <- possibly(
  .f = function(chunk) {
    rxtracto(
      dataInfo = modis_info,
      parameter = "chlorophyll",
      tcoord = chunk$datecollected,
      xcoord = chunk$longitude,
      ycoord = chunk$latitude,
      xlen = 0.05,
      ylen = 0.05
    )
  },
  otherwise = NULL
)


# 2. Combine the data frames directly

# 1. Convert each chunk list into a flat N-row data frame
df_list <- lapply(results, function(res) {
  if (is.null(res)) return(NULL)
  as.data.frame(res)
})

df_modis_extract <- do.call(rbind, df_list)
head(df_modis_extract)
dfm <- df_modis_extract %>% mutate(longitude = (requested.lon.min+requested.lon.max)/2, latitude = (requested.lat.min+requested.lat.max)/2) %>% dplyr::select(mean.chlorophyll, requested.date, longitude, latitude) %>% rename(datecollected = requested.date)

dfm_com <- merge(df_modis, dfm, by = c('datecollected', 'latitude', 'longitude'))

dfm_com <- dfm_com %>% distinct()
write.csv(dfm_com, file = 'Data/MODIS_chl_MRAC.csv')


#get zones from NOAA commercial fishing zones

dfm_com_zone <- dfm_com %>%
  mutate(
    # Force negative longitude
    lon_neg = -abs(longitude),
    
    zone = case_when(
      # South Texas Coast
      lon_neg <= -96.5 & latitude >= 26.5 ~ 20,
      lon_neg <= -96.5 & latitude >= 24.0 ~ 21,
      
      # Northern & Western Gulf Strips
      lon_neg > -96.5 & lon_neg <= -95.5 ~ 19,
      lon_neg > -95.5 & lon_neg <= -94.5 ~ 18,
      lon_neg > -94.5 & lon_neg <= -93.5 ~ 17,
      lon_neg > -93.5 & lon_neg <= -92.5 ~ 16,
      lon_neg > -92.5 & lon_neg <= -91.5 ~ 15,
      lon_neg > -91.5 & lon_neg <= -90.5 ~ 14,
      lon_neg > -90.5 & lon_neg <= -89.5 ~ 13,
      lon_neg > -89.5 & lon_neg <= -89.0 ~ 12,
      lon_neg > -89.0 & lon_neg <= -88.0 ~ 11,
      lon_neg > -88.0 & lon_neg <= -87.0 ~ 10,
      lon_neg > -87.0 & lon_neg <= -86.0 ~ 9,
      lon_neg > -86.0 & lon_neg <= -85.0 & latitude >= 29.0 ~ 8,
      
      # Florida Shelf Bands
      lon_neg > -85.0 & latitude >= 29.0 ~ 7,
      lon_neg > -85.0 & latitude >= 28.0 & latitude < 29.0 ~ 6,
      lon_neg > -85.0 & latitude >= 27.0 & latitude < 28.0 ~ 5,
      lon_neg > -85.0 & latitude >= 26.0 & latitude < 27.0 ~ 4,
      lon_neg > -85.0 & latitude >= 25.0 & latitude < 26.0 ~ 3,
      lon_neg > -85.0 & lon_neg <= -82.0 & latitude >= 24.0 & latitude < 25.0 ~ 2,
      lon_neg > -82.0 & latitude < 24.8 ~ 1,
      
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-lon_neg) # Drop temporary coordinate column

write.csv(file = 'Data/MRAC_chl_zone.csv', dfm_com_zone)

library(dplyr)
library(CopernicusMarine)
library(terra)

# Install reticulate if you don't have it
install.packages("reticulate")
library(reticulate)

# Install the official Copernicus Marine python toolbox
# This creates a local python environment with the tool installed
py_install("copernicusmarine", pip = TRUE)

## 1. Load your dataset and get the bounding box
df <- read.csv("Data/MRAC_chl_zone.csv") %>%
  mutate(date = as.Date(datecollected))

# Import the Python copernicusmarine module
cm <- import("copernicusmarine")

# 1. Log in securely (This hits the correct, updated server)
cm$login(username = "jrodemann", password = "Cynoscion!2025")

# 2. Download the data 
# Note: In the python toolbox, this function is called 'subset'
cm$subset(
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
  output_filename = "gulf_ocean_model.nc",
  output_directory = getwd(),
  variables = list("thetao", "bottomT", "so", "uo", "vo"),
  minimum_longitude = min(df$longitude, na.rm = TRUE) - 0.5,
  maximum_longitude = max(df$longitude, na.rm = TRUE) + 0.5,
  minimum_latitude = min(df$latitude, na.rm = TRUE) - 0.5,
  maximum_latitude = max(df$latitude, na.rm = TRUE) + 0.5,
  # Format dates as strings ("YYYY-MM-DD")
  start_datetime = as.character(min(df$date, na.rm = TRUE)),
  end_datetime = as.character(max(df$date, na.rm = TRUE)),
  minimum_depth = 0.49402499198913574,
  maximum_depth = 0.494025
)


library(terra)
library(future)
library(furrr)

# 1. Set up parallel processing
# Limit workers to 4-6 to prevent your RAM from overflowing again
plan(multisession, workers = 8) 

# 2. Get unique dates from your dataset
unique_dates <- as.character(unique(df$date))


# 3. Run the extraction in parallel across the dates
results_list <- future_map(unique_dates, function(d) {
  
  # -- INSIDE THE PARALLEL WORKER --
  library(terra)
  library(dplyr)
  
  # Connect to the NetCDF file LOCALLY inside the worker
  worker_raster <- rast("gulf_ocean_model.nc")
  
  # Filter the raw data frame for this specific date
  df_day <- df %>% filter(date == as.Date(d))
  pts_day <- vect(df_day, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # FIX: Convert the raster's time dimension to a standard Date before matching
  date_match <- as.Date(time(worker_raster)) == as.Date(d)
  
  # Ensure there is at least one matching layer for this date
  if (any(date_match, na.rm = TRUE)) {
    
    # Subset the raster for this specific date
    raster_day <- worker_raster[[date_match]]
    
    # Extract the data
    ext_day <- extract(raster_day, pts_day, method = "simple") 
    df_day_extracted <- bind_cols(df_day, ext_day)
    return(df_day_extracted)
    
  } else {
    return(NULL) # Return NULL if this date isn't in the NetCDF
  }
  
}, .options = furrr_options(seed = TRUE)) 

# 4. Shut down parallel workers
plan(sequential)

# 5. Bind all the daily chunks back together
df_mrac_glor <- data.table::rbindlist(results_list, fill = T) %>%
  select(-ID) 
head(df_mrac_glor)
df_all <- df_mrac_glor %>% rename(temp = `thetao_depth=0.49402499_1`, salinity = `so_depth=0.49402499_1`, east_velo = `uo_depth=0.49402499_1`, north_velo = `vo_depth=0.49402499_1`, bottom_temp = bottomT_1)

write.csv(file = 'Data/mrac_final.csv', df_all)

df <- read_csv('Data/mrac_final.csv')

df_int2 <- df %>% dplyr::filter(data_source == 'mark-recapture') %>% group_by(catalogNumber) %>% dplyr::filter(any(interaction_order == 2) & !any(interaction_order == 1)) %>%
  ungroup()

write.csv(df_int2, file = 'Data/interaction_only2.csv')

#add in 7 regions from other map
head(df)
library(sf)
regions_df <- data.frame(
  region = c("MX_Shelf", "TX_Shelf", "NCentral_GOA_Shelf", 
             "FL_Panhandle", "WFL_Shelf", "SFL_Shelf", "EFL_Shelf"),
  xmin = c(-98.0, -98.0, -93.5, -87.5, -87.5, -84.0, -81.5),
  xmax = c(-87.5, -93.5, -87.5, -84.0, -81.5, -79.0, -79.0),
  ymin = c( 18.0,  26.0,  26.0,  28.5,  25.0,  23.0,  25.0),
  ymax = c( 26.0,  31.0,  31.0,  31.0,  28.5,  25.0,  31.0)
)

# 2. Function to convert min/max coordinates into a closed spatial polygon
make_poly <- function(xmin, xmax, ymin, ymax) {
  st_polygon(list(matrix(c(xmin, ymin,
                           xmax, ymin,
                           xmax, ymax,
                           xmin, ymax,
                           xmin, ymin), ncol = 2, byrow = TRUE)))
}

# 3. Create a spatial (sf) dataframe containing the region polygons
regions_sf <- regions_df %>%
  rowwise() %>%
  mutate(geometry = list(make_poly(xmin, xmax, ymin, ymax))) %>%
  ungroup() %>%
  # Assign standard WGS84 coordinate reference system (EPSG:4326)
  st_as_sf(crs = 4326)

# 4. Convert your main dataset into a spatial object
# (Assuming your dataset is called 'df' and longitudes are negative)
df_spatial <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# 5. Find the nearest region for every point
# This handles both points strictly inside the boxes AND points in the deep Gulf
nearest_indices <- st_nearest_feature(df_spatial, regions_sf)

# 6. Append the matched region names back to your original dataset
df$region <- regions_sf$region[nearest_indices]

df_slice <- df %>% dplyr::slice(1:1000)
write.csv(df_slice, file = 'Data/sample_dat.csv')

df <- df %>%
  mutate(date_parsed = as.POSIXct(datecollected)) %>%
  arrange(catalogNumber, date_parsed)

total_fish_per_source <- df %>%
  group_by(data_source) %>%
  summarise(total_fish = n_distinct(catalogNumber), .groups = "drop")


# ==============================================================================
# TABLE 1a: Count of interaction_order == 1 by REGION
# ==============================================================================
table1a_interaction_1_region <- df %>%
  filter(interaction_order == 1) %>%
  # Group only by data_source and region
  group_by(data_source, region) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(data_source, region)

# ==============================================================================
# TABLE 1b: Count of interaction_order == 1 by ZONE
# ==============================================================================
table1b_interaction_1_zone <- df %>%
  filter(interaction_order == 1) %>%
  # Group only by data_source and zone
  group_by(data_source, zone) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(data_source, zone)

# ==============================================================================
# UPDATED TABLE 2: Number and Percentage of fish that DID NOT LEAVE their zone
# ==============================================================================
table2_stayed_zone <- df %>%
  filter(!is.na(zone)) %>%
  group_by(data_source, catalogNumber) %>%
  summarise(n_zones = n_distinct(zone), .groups = "drop") %>%
  # Filter to fish that only registered in 1 zone
  filter(n_zones == 1) %>%
  group_by(data_source) %>%
  summarise(fish_that_stayed_in_zone = n(), .groups = "drop") %>%
  # Join the totals to calculate the percentage
  left_join(total_fish_per_source, by = "data_source") %>%
  mutate(
    percentage = round((fish_that_stayed_in_zone / total_fish) * 100, 1)
  )

# ==============================================================================
# UPDATED TABLE 3: Number and Percentage of fish that DID NOT LEAVE their region
# ==============================================================================
table3_stayed_region <- df %>%
  filter(!is.na(region)) %>%
  group_by(data_source, catalogNumber) %>%
  summarise(n_regions = n_distinct(region), .groups = "drop") %>%
  # Filter to fish that only registered in 1 region
  filter(n_regions == 1) %>%
  group_by(data_source) %>%
  summarise(fish_that_stayed_in_region = n(), .groups = "drop") %>%
  # Join the totals to calculate the percentage
  left_join(total_fish_per_source, by = "data_source") %>%
  mutate(
    percentage = round((fish_that_stayed_in_region / total_fish) * 100, 1))
# ==============================================================================
# TABLE 4: Number of fish by their HIGHEST interaction_order
# ==============================================================================
table4_max_interaction <- df %>%
  group_by(data_source, catalogNumber) %>%
  # Find the max interaction order for each fish
  summarise(max_interaction = max(interaction_order, na.rm = TRUE), .groups = "drop") %>%
  # Count how many fish reached each max interaction level
  group_by(data_source, max_interaction) %>%
  summarise(number_of_fish = n(), .groups = "drop") %>%
  arrange(data_source, desc(max_interaction))

# ==============================================================================
# TABLE 5: Zone Transitions (Zone X -> Zone Y)
# ==============================================================================
table5_zone_transitions <- df %>%
  filter(!is.na(zone)) %>%
  group_by(catalogNumber) %>%
  # Keep only the rows where a fish ACTUALLY changed zones 
  # (Filters out consecutive pings in the same zone)
  filter(row_number() == 1 | zone != lag(zone)) %>%
  # Look ahead to the next zone the fish visited
  mutate(next_zone = lead(zone)) %>%
  filter(!is.na(next_zone)) %>%
  ungroup() %>%
  # Group by source and the transition combo to count them
  group_by(data_source, zone, next_zone) %>%
  summarise(number_of_fish = n(), .groups = "drop") %>%
  # Format neatly as "X -> Y"
  mutate(transition = paste(zone, "->", next_zone)) %>%
  select(data_source, transition, number_of_fish) %>%
  arrange(data_source, desc(number_of_fish))

# ==============================================================================
# TABLE 6: Region Transitions (Region X -> Region Y)
# ==============================================================================
table6_region_transitions <- df %>%
  filter(!is.na(region)) %>%
  group_by(catalogNumber) %>%
  # Keep only the rows where a fish ACTUALLY changed regions
  filter(row_number() == 1 | region != lag(region)) %>%
  mutate(next_region = lead(region)) %>%
  filter(!is.na(next_region)) %>%
  ungroup() %>%
  group_by(data_source, region, next_region) %>%
  summarise(number_of_fish = n(), .groups = "drop") %>%
  mutate(transition = paste(region, "->", next_region)) %>%
  select(data_source, transition, number_of_fish) %>%
  arrange(data_source, desc(number_of_fish))

# ==============================================================================
# View the generated tables
# ==============================================================================
print(table1a_interaction_1_region)
print(table1b_interaction_1_zone)
print(table2_stayed_zone)
print(table3_stayed_region)
print(table4_max_interaction)
print(table5_zone_transitions)
print(table6_region_transitions)

write.csv(table1a_interaction_1_region, file = 'Data/summary_tables/region_tag_summary.csv')
write.csv(table1b_interaction_1_zone, file = 'Data/summary_tables/zone_tag_summary.csv')
write.csv(table2_stayed_zone, file = 'Data/summary_tables/number_stayed_zone.csv')
write.csv(table3_stayed_region, file = 'Data/summary_tables/number_stayed_region.csv')
write.csv(table4_max_interaction, file = 'Data/summary_tables/max_interactions.csv')
write.csv(table5_zone_transitions, file = 'Data/summary_tables/zone_transitions.csv')
write.csv(table6_region_transitions, file = 'Data/summary_tables/region_transitions.csv')
