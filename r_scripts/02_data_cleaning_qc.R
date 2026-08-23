# 02_data_cleaning_qc.R
# Turns the raw GBIF downloads into two modelling tables: cleaned dorsalis presences
# and cleaned Tephritidae records for the background. Every filter writes a row to a
# cleaning log, so the drop from raw count to final count can be read off rather than
# reconstructed. Ends with the Day-2 kill-switch check on African record spread.


# Section 1: Load required libraries
library(here)
library(readr)
library(dplyr)
library(CoordinateCleaner)
library(terra)
library(sf)


# Section 2: Load and inspect raw Bactrocera dorsalis occurrence data
# guess_max = Inf forces readr to scan every row before assigning column types. The
# Darwin Core archive has 230 columns and many are empty for thousands of rows, so a
# short guess window silently mistypes them.
bd_data <- read_tsv(
  here("data", "raw", "0019221-260806074905277.zip"),
  guess_max = Inf
)

# First look at the data. Coordinate ranges outside -180/180 and -90/90, or any missing
# coordinates, would mean the GBIF filters did not do what was asked of them.
print(
  bd_data %>%
    summarise(
      total_records = n(),
      unique_countries = n_distinct(countryCode, na.rm = TRUE),
      missing_longitude = sum(is.na(decimalLongitude)),
      missing_latitude = sum(is.na(decimalLatitude)),
      min_longitude = min(decimalLongitude, na.rm = TRUE),
      max_longitude = max(decimalLongitude, na.rm = TRUE),
      min_latitude = min(decimalLatitude, na.rm = TRUE),
      max_latitude = max(decimalLatitude, na.rm = TRUE),
      duplicate_gbifid = sum(duplicated(gbifID))
    )
)

# How clustered are the records? A single coordinate carrying hundreds of records is
# normal for trapping studies and is the reason thinning is needed in Section 2.2.
print(
  bd_data %>%
    count(decimalLongitude, decimalLatitude) %>%
    filter(n > 1) %>%
    summarise(
      repeated_coordinates = n(),
      records_at_repeated_coordinates = sum(n),
      max_records_at_one_coordinate = max(n)
    )
)


# Section 2.1: Clean problematic coordinates
# Five tests, chosen because each flags a coordinate that cannot be a real observation:
# administrative centroids, identical lat/lon, the GBIF office in Copenhagen, museum and
# herbarium addresses, and exact zeros. Capital-city and sea tests were run during scoping
# but are not applied here. Capitals flagged 433 records, most of them genuine urban
# observations, and the sea test is highly sensitive to coastline resolution. Neither is a
# specific enough indicator of error to delete on.
bd_qc_flags <- clean_coordinates(
  x = bd_data,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  tests = c(
    "centroids",
    "equal",
    "gbif",
    "institutions",
    "zeros"
  ),
  value = "spatialvalid"
)

# Count what each test caught. The columns are flags, so TRUE means the record passed;
# the negation counts failures.
print(
  bd_qc_flags %>%
    summarise(
      centroid_flags = sum(!.cen),
      equal_coordinate_flags = sum(!.equ),
      gbif_headquarters_flags = sum(!.gbf),
      institution_flags = sum(!.inst),
      zero_coordinate_flags = sum(!.zer)
    )
)

# The log starts at the raw count and gains one row per filter.
cleaning_log <- data.frame(
  stage = "Raw GBIF occurrences",
  records = nrow(bd_data)
)

# Keep only records that passed all five tests. Records flagged by more than one test are
# dropped once, so the loss is smaller than the sum of the individual flag counts.
bd_clean <- bd_qc_flags %>%
  filter(.cen, .inst, .equ, .gbf, .zer)

cleaning_log <- bind_rows(
  cleaning_log,
  data.frame(
    stage = "After coordinate artefact filtering",
    records = nrow(bd_clean)
  )
)

print(cleaning_log)


# Section 2.2: Assign occurrences to WorldClim cells and spatially thin
bioclim_stack <- rast(list.files(
  here("data", "raw", "climate", "wc2.1_2.5m"),
  pattern = "\\.tif$",
  full.names = TRUE
))

# Only bio_1 is needed here. cells = TRUE returns the raster cell index for each point,
# which is what the thinning step keys on, and the bio_1 value doubles as a coverage test.
bd_to_bioclim <- extract(
  x = subset(bioclim_stack, 1),
  y = data.frame(
    decimal_longitude = bd_clean$decimalLongitude,
    decimal_latitude = bd_clean$decimalLatitude
  ),
  method = "simple",
  cells = TRUE
)

# A missing bio_1 means the point fell outside the WorldClim land mask, usually just
# offshore. Those records cannot enter a climate model, so they go here. This is also how
# the sea-flag question from Section 2.1 gets settled: by predictor availability rather
# than by a coastline polygon.
bd_clean <- bd_clean %>%
  mutate(
    cell = bd_to_bioclim$cell,
    bio1 = bd_to_bioclim$wc2.1_2.5m_bio_1
  ) %>%
  filter(!is.na(bio1))

cleaning_log <- bind_rows(
  cleaning_log,
  data.frame(
    stage = "After removing occurrences without WorldClim coverage",
    records = nrow(bd_clean)
  )
)

# One record per grid cell. Exact-coordinate duplicates were deliberately not removed as a
# separate step, because thinning at cell level already handles both exact repeats and
# distinct coordinates that carry identical climate values.
bd_clean <- bd_clean %>%
  distinct(cell, .keep_all = TRUE)

cleaning_log <- bind_rows(
  cleaning_log,
  data.frame(
    stage = "After thinning to one occurrence per WorldClim cell",
    records = nrow(bd_clean)
  )
)

write_csv(
  cleaning_log,
  here("data", "processed", "cleaning_log.csv")
)

# Drop the working climate column and the CoordinateCleaner flag columns, which all begin
# with a dot, before writing the modelling table.
bd_clean <- bd_clean %>%
  select(-bio1, -starts_with("."))

write_csv(
  bd_clean,
  here("data", "processed", "bd_clean.csv")
)


# Section 3: Load and inspect raw Tephritidae occurrence data
# This archive is much larger and messier than the dorsalis one. Reading every column as
# character with quoting disabled avoids parse failures on unescaped quotes in free-text
# fields; the two coordinate columns are converted back to numeric immediately after.
teph_raw <- read_tsv(
  unz(
    here("data", "raw", "0022617-260806074905277.zip"),
    "occurrence.txt"
  ),
  col_types = cols(.default = "c"),
  quote = ""
)

glimpse(teph_raw)

teph_raw <- teph_raw %>%
  mutate(
    decimalLongitude = as.numeric(decimalLongitude),
    decimalLatitude = as.numeric(decimalLatitude)
  )

# Same inspection as the dorsalis table, repeated because this file was parsed differently.
teph_raw %>%
  summarise(
    total_records = n(),
    unique_countries = n_distinct(countryCode, na.rm = TRUE),
    missing_longitude = sum(is.na(decimalLongitude)),
    missing_latitude = sum(is.na(decimalLatitude)),
    min_longitude = min(decimalLongitude, na.rm = TRUE),
    max_longitude = max(decimalLongitude, na.rm = TRUE),
    min_latitude = min(decimalLatitude, na.rm = TRUE),
    max_latitude = max(decimalLatitude, na.rm = TRUE),
    duplicate_gbifid = sum(duplicated(gbifID))
  ) %>%
  glimpse()

# Confirms pred_not() worked in 01_downloads.R. 7930834 is the dorsalis taxon key. Both
# counts must be zero, including the accepted key, which catches records filed under a
# synonym that GBIF later resolved to dorsalis.
print(
  teph_raw %>%
    summarise(
      no_of_dorsalis = sum(taxonKey == "7930834", na.rm = TRUE),
      no_of_accepted_dorsalis = sum(
        acceptedTaxonKey == "7930834",
        na.rm = TRUE
      )
    )
)


# Section 3.1: Clean problematic coordinates
# The background is put through exactly the same five tests as the presences. If the two
# sides were cleaned differently, any difference the model found between them could be an
# artefact of the cleaning rather than of climate.
teph_cleaning_log <- data.frame(
  stage = "Raw Tephritidae GBIF occurrences",
  records = nrow(teph_raw)
)

teph_qc_flags <- clean_coordinates(
  x = teph_raw,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  tests = c(
    "centroids",
    "equal",
    "gbif",
    "institutions",
    "zeros"
  ),
  value = "spatialvalid"
)

teph_clean <- teph_qc_flags %>%
  filter(.cen, .inst, .equ, .gbf, .zer)

nrow(teph_clean)

teph_cleaning_log <- bind_rows(
  teph_cleaning_log,
  data.frame(
    stage = "After coordinate artefact filtering",
    records = nrow(teph_clean)
  )
)


# Section 3.2: Assign occurrences to WorldClim cells and spatially thin
# Identical treatment to Section 2.2, and for the same reason: presences and background
# must be thinned on the same grid or the comparison between them is not like for like.
teph_to_bioclim <- extract(
  x = subset(bioclim_stack, 1),
  y = data.frame(
    decimal_longitude = teph_clean$decimalLongitude,
    decimal_latitude = teph_clean$decimalLatitude
  ),
  method = "simple",
  cells = TRUE
)

teph_clean <- teph_clean %>%
  mutate(
    cell = teph_to_bioclim$cell,
    bio1 = teph_to_bioclim$wc2.1_2.5m_bio_1
  ) %>%
  filter(!is.na(bio1))

teph_cleaning_log <- bind_rows(
  teph_cleaning_log,
  data.frame(
    stage = "After removing occurrences without WorldClim coverage",
    records = nrow(teph_clean)
  )
)

teph_clean <- teph_clean %>%
  distinct(cell, .keep_all = TRUE)

teph_cleaning_log <- bind_rows(
  teph_cleaning_log,
  data.frame(
    stage = "After thinning to one occurrence per WorldClim cell",
    records = nrow(teph_clean)
  )
)

write_csv(
  teph_cleaning_log,
  here("data", "processed", "teph_cleaning_log.csv")
)

teph_clean <- teph_clean %>%
  select(-bio1, -starts_with("."))

write_csv(
  teph_clean,
  here("data", "processed", "teph_clean.csv")
)


# Section 4: Day-2 kill-switch check
# The stated stopping rule: if the cleaned African records occupy fewer than about 25
# distinct 100 km cells, the model would be describing a handful of field campaigns rather
# than a climate envelope, and the species should be switched to Spodoptera frugiperda.
bd_africa <- bd_clean %>%
  filter(continent == "AFRICA")

nrow(bd_africa)

# EPSG:6933 is an equal-area projection. The grid has to be built in metres, not degrees,
# or a 100 km cell would not be 100 km outside the equator.
bd_africa_sf <- bd_africa %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  ) %>%
  st_transform(crs = 6933)

bd_africa_grid <- st_make_grid(
  bd_africa_sf,
  cellsize = 100000
)

# The number printed here is the decision. Count of distinct occupied 100 km cells.
print(
  bd_africa_sf %>%
    st_intersects(bd_africa_grid) %>%
    unlist() %>%
    unique() %>%
    length()
)
