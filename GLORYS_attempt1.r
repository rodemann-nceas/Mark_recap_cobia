###Adding GLORYS and MODIS data to mark-recap data##
# -------------------------------------------------------------------------
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# -------------------------------------------------------------------------
# Install packages if you don't have them:
install.packages('CFtime')
remotes::install_github("rmendels/rerddapXtracto")

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
df <- read_csv("Data/NCEAS_Cobia_Recaps_202605-1778080360235.csv")

# Ensure date formatting is correct and remove rows with missing coordinates
df_clean <- df %>%
  mutate(datecollected = as.Date(datecollected)) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(datecollected)) %>% 
  mutate(cat_num = paste0(catalognumber...13, '-', interaction_order))

# -------------------------------------------------------------------------
# 3. EXTRACT MODIS AQUA CHLOROPHYLL (Via NOAA ERDDAP)
# -------------------------------------------------------------------------
print("Starting MODIS Aqua Extraction...")

# MODIS Aqua started providing data around July 2002
df_modis <- df_clean %>% 
  filter(datecollected >= as.Date("2003-01-01")) %>% 
  filter(datecollected <= as.Date('2022-07-26'))

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

# Bind the extracted mean chlorophyll back to the dataset
df_modis$chlorophyll_mg_m3 <- chl_extract$`mean chlorophyll`

names(df_modis)

df_modis <- df_modis %>% dplyr::select(-catalognumber...23) %>% rename(catalognumber = catalognumber...13) 
df_modis <- df_modis %>% mutate(cat_num = paste0(catalognumber, '-', interaction_order))

# Re-join with the master dataset (leaves NA for dates before MODIS)

df_clean <- df_clean %>%
  left_join(df_modis %>% select(cat_num, chlorophyll_mg_m3), by = "cat_num") %>% distinct()

write.csv(df_clean, file = 'Data/MODIS_chl_MR.csv')

# -------------------------------------------------------------------------
# 4. EXTRACT GLORYS DATA (Copernicus Marine)
# -------------------------------------------------------------------------
print("Starting GLORYS Extraction...")

# GLORYS12V1 covers 1993 to ~present. Filter out pre-1993 data.
df_glorys <- df_clean %>%
  filter(datecollected >= as.Date("1993-01-01")) %>% 
  mutate(cat_num = paste0(catalognumber...13, '-', interaction_order))

# Ensure you are logged into Copernicus Marine API
# Run this once in your console if you haven't: cms_login()

# Pre-allocate columns for our new variables
df_glorys <- df_glorys %>%
  mutate(
    glorys_sst = NA_real_,
    glorys_bottom_temp = NA_real_,
    glorys_sss = NA_real_,
    glorys_u_velocity = NA_real_,
    glorys_v_velocity = NA_real_
  )

# GLORYS Product IDs and Layers
product_id <- "GLOBAL_MULTIYEAR_PHY_001_030"
layer_3d <- "cmems_mod_glo_phy_my_0.083deg_P1D-m" # For thetao, so, uo, vo
layer_2d <- "cmems_mod_glo_phy_mybottom_0.083deg_P1D-m" # For bottomT

# Loop through each row to download just the pixel we need, extract it, and delete the file
for (i in 1:nrow(df_glorys)) {
  
  lat <- df_glorys$latitude[i]
  lon <- df_glorys$longitude[i]
  date <- df_glorys$datecollected[i]
  
  cat(sprintf("Processing %d of %d: Date %s at [%.2f, %.2f]...\n", i, nrow(df_glorys), date, lon, lat))
  
  tryCatch({
    # --- 1. Extract 3D Variables (Surface Temp, Salinity, Velocities) ---
    # We set z = c(0, 0.5) to only grab the surface layer (0 meters)
    file_3d <- cms_download_subset(
      product = product_id,
      layer = layer_3d,
      variable = c("thetao", "so", "uo", "vo"),
      x = c(lon - 0.05, lon + 0.05),
      y = c(lat - 0.05, lat + 0.05),
      t = c(as.character(date), as.character(date)),
      z = c(0, 0.5), 
      destination = tempdir(),
      overwrite = TRUE
    )
    
    # Read the downloaded NetCDF and extract the exact point value
    rast_3d <- rast(file_3d)
    pt_val_3d <- terra::extract(rast_3d, matrix(c(lon, lat), ncol = 2))
    
    df_glorys$glorys_sst[i]       <- pt_val_3d$thetao
    df_glorys$glorys_sss[i]       <- pt_val_3d$so
    df_glorys$glorys_u_velocity[i] <- pt_val_3d$uo
    df_glorys$glorys_v_velocity[i] <- pt_val_3d$vo
    
    # --- 2. Extract 2D Variable (Seafloor / Bottom Temperature) ---
    file_2d <- cms_download_subset(
      product = product_id,
      layer = layer_2d,
      variable = "bottomT",
      x = c(lon - 0.05, lon + 0.05),
      y = c(lat - 0.05, lat + 0.05),
      t = c(as.character(date), as.character(date)),
      destination = tempdir(),
      overwrite = TRUE
    )
    
    rast_2d <- rast(file_2d)
    pt_val_2d <- terra::extract(rast_2d, matrix(c(lon, lat), ncol = 2))
    
    df_glorys$glorys_bottom_temp[i] <- pt_val_2d$bottomT
    
    # Cleanup: Delete the temp files to save space
    file.remove(file_3d, file_2d)
    
  }, error = function(e) {
    message(paste("Failed to pull data for row", i, "-", e$message))
  })
}

# -------------------------------------------------------------------------
# 5. FINAL EXPORT
# -------------------------------------------------------------------------
# Re-merge the completed GLORYS data back into the main pipeline
final_df <- df_clean %>%
  left_join(
    df_glorys %>% select(cat_num, glorys_sst, glorys_bottom_temp, glorys_sss, glorys_u_velocity, glorys_v_velocity),
    by = "cat_num"
  )

# Export the fully enriched dataset
write_csv(final_df, "NCEAS_Cobia_Recaps_Enriched.csv")
print("Extraction Complete! File saved as NCEAS_Cobia_Recaps_Enriched.csv")
