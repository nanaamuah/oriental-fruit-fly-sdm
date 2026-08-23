# 01_downloads.R
# Retrieves everything the analysis is built from: the Bactrocera dorsalis
# occurrence records, the Tephritidae records that later become the target-group
# background, the WorldClim predictor stack, and the EPPO distribution table.
# Everything lands in data/raw and is never edited in place.


# Section 1: Load required libraries
library(here)
library(rgbif)
library(terra)
library(geodata)


# Section 2: Translate scientific name into taxon key and download occurrence data from GBIF
# The taxon key is required to download occurrence data from GBIF
taxon_key <- name_backbone(name = "Bactrocera dorsalis")$usageKey

# Download occurrence data from GBIF using the taxon key
# Before downloading, use occ_count() to check the number of records available for the taxon key in specified location against global occurrence data
print(occ_count(taxonKey = taxon_key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE)) #11757
print(occ_count(taxonKey = taxon_key, continent = "Africa", hasCoordinate = TRUE, hasGeospatialIssue = FALSE)) # 891
# A huge difference in the number of records available for the taxon key in specified location against global occurrence data.
# This may be because the species is native to Asia and has been introduced to Africa.
# Proceed with global occurrence, sf can later be used to clip the study area to avoid losing records that are inside the study area.

# Before download, add your GBIF username and password to your .Renviron file in the project root directory. This is required to download occurrence data from GBIF.
# usethis::edit_r_environ()
# # Add the following lines to your .Renviron file:
# GBIF_USER=your_gbif_username
# GBIF_PWD=your_gbif_password
# GBIF_EMAIL=your_gbif_email            # Restart R session after adding your GBIF username and password to your .Renviron file.
# Sys.getenv("GBIF_USER")          # Confirm that your GBIF username is stored in your .Renviron file.

# occ_download() is asynchronous. GBIF queues the request, mints a citable DOI for it,
# and returns a key. This route is used instead of occ_search() precisely because of
# that DOI, which is what makes the download citable and reproducible.
request_data <- occ_download(
  pred_and(
    pred("taxonKey", taxon_key),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE)
  )
)

# Blocks until GBIF reports the download as ready. Record the DOI it prints.
occ_download_wait(request_data)

global_occurrence_data <- occ_download_get(
  key = request_data,
  path = here("data", "raw")
)


# Section 3: Download WorldClim climatic predictors
# 2.5 arc-minutes is roughly 4.5 km at the equator. Finer resolutions buy nothing here
# because the occurrence records are thinned to one per cell anyway.
worldclim_data <- worldclim_global(var = "bio", res = 2.5, path = here("data", "raw"))

# Load the WorldClim data into a raster stack
worldclim_stack <- rast(list.files(
  here("data", "raw", "climate", "wc2.1_2.5m"),
  pattern = "\\.tif$",
  full.names = TRUE
))

# Check all the layers in the raster stack loaded successfully (19 layers are expected)
# as well as other parameters
print(worldclim_stack) # nlyr should be 19

# Check names of all the layers in the raster stack to confirm sequence from 1 to 19
print(names(worldclim_stack)) # also ensures there are no duplicates


# Section 4: Download EPPO data for Bactrocera dorsalis
# The CSV file can be downloaded directly from the EPPO Global Database website and loaded
# The csv url can be copied and pasted as shown below to download the data directly into R
# DACUDO is the EPPO code for B. dorsalis. This table is used only as an independent
# check on the modelled surface, never as model input.
download.file(
  url = "https://gd.eppo.int/taxon/DACUDO/download/distribution_csv",
  destfile = here("data", "raw", "eppo_distribution_data.csv"),
  mode = "wb"
)


# Section 5: Obtain occurrence data on Tephritidae as target group background
# Exclude Bactrocera dorsalis from the target group background data to avoid bias in the model
# The logic: fruit flies are collected by the same traps and the same people, so where
# other Tephritidae were recorded is a map of where anyone was looking. Using that as the
# background absorbs survey effort instead of treating it as climate.
# Obtain the taxon key for Tephritidae
teph_key <- name_backbone(name = "Tephritidae")$usageKey

# Download occurrence data for Tephritidae from GBIF using the taxon key and exclude Bactrocera dorsalis from the target group background data
# pred_not() drops dorsalis at the query stage. Records that share a grid cell with a
# dorsalis presence are removed later, in 03_analysis.R.
teph_data <- occ_download(
  pred_and(
    pred("taxonKey", teph_key),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("occurrenceStatus", "PRESENT"),
    pred_not(pred("taxonKey", taxon_key))
  )
)

occ_download_wait(teph_data)

teph_occurrence_data <- occ_download_get(
  key = teph_data,
  path = here("data", "raw")
)
