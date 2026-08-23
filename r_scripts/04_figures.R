# 04_figures.R
# Turns the objects written by 03_analysis.R into the finished figures, maps and tables
# in outputs/. This script does no modelling and fits nothing. It reads only from
# data/processed, so a figure can be restyled without refitting anything.


# Section 1: Load required libraries
library(readr)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(terra)
library(sf)


# Section 2: Load required data
cv_results <- read_csv(here("data", "processed", "cv_fold_results.csv"))


# Section 3: Build the headline validation figure
# Long format is needed so the three metrics can be faceted down the rows. The factor
# levels fix the panel order; without them ggplot sorts alphabetically and Boyce lands
# first, which reads oddly when AUC is the metric most readers look for.
cv_long <- cv_results %>%
  pivot_longer(cols = c(AUC, TSS, Boyce), names_to = "metric", values_to = "value")

cv_long <- cv_long %>%
  mutate(
    metric = factor(metric, levels = c("AUC", "TSS", "Boyce")),
    validation = factor(
      validation,
      levels = c("random", "spatial"),
      labels = c("Random", "Spatial")
    )
  )

# Individual folds are plotted as well as the mean, deliberately. The mean alone would
# hide the fact that the spatial folds disagree violently with each other, and that
# spread is as much a part of the finding as the drop in the average.
validation_plot <- ggplot(
  cv_long,
  aes(x = validation, y = value)
) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.7,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.15
  ) +
  facet_grid(
    metric ~ model,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Random cross-validation overestimates model performance",
    subtitle = "Points show individual folds; large points and error bars show mean ± SD",
    x = "Validation scheme",
    y = "Validation score"
  )

validation_plot

ggsave(
  validation_plot,
  filename = here("outputs", "figures", "validation_plot.png"),
  width = 8,
  height = 6,
  dpi = 300
)


# Section 4: Create Validation Metric Table
cv_summary <- read_csv(
  here("data", "processed", "cv_summary.csv")
)

# Grouping by model means each algorithm's drop is computed against its own random-fold
# score, rather than against a pooled baseline that would mix MaxEnt and GLM together.
validation_metric_table <- cv_summary %>%
  group_by(model) %>%
  mutate(
    AUC_drop = mean_AUC[validation == "random"] -
      mean_AUC[validation == "spatial"],
    TSS_drop = mean_TSS[validation == "random"] -
      mean_TSS[validation == "spatial"],
    Boyce_drop = mean_Boyce[validation == "random"] -
      mean_Boyce[validation == "spatial"]
  ) %>%
  mutate(
    AUC_drop = if_else(validation == "spatial", AUC_drop, NA_real_),
    TSS_drop = if_else(validation == "spatial", TSS_drop, NA_real_),
    Boyce_drop = if_else(validation == "spatial", Boyce_drop, NA_real_)
  ) %>%
  ungroup()

print(validation_metric_table)

# The drop is written only on the spatial row. Repeating it on both rows would imply the
# random-fold score has a drop of its own, which it does not; it is the baseline.
validation_table_display <- validation_metric_table %>%
  mutate(
    Algorithm = model,
    Validation = recode(
      validation,
      random = "Random",
      spatial = "Spatial"
    ),
    AUC = sprintf("%.3f ± %.3f", mean_AUC, sd_AUC),
    TSS = sprintf("%.3f ± %.3f", mean_TSS, sd_TSS),
    Boyce = sprintf("%.3f ± %.3f", mean_Boyce, sd_Boyce),
    `AUC drop` = if_else(
      is.na(AUC_drop),
      "—",
      sprintf("%.3f", AUC_drop)
    ),
    `TSS drop` = if_else(
      is.na(TSS_drop),
      "—",
      sprintf("%.3f", TSS_drop)
    ),
    `Boyce drop` = if_else(
      is.na(Boyce_drop),
      "—",
      sprintf("%.3f", Boyce_drop)
    )
  ) %>%
  select(
    Algorithm,
    Validation,
    AUC,
    TSS,
    Boyce,
    `AUC drop`,
    `TSS drop`,
    `Boyce drop`
  ) %>%
  arrange(Algorithm, factor(Validation, levels = c("Random", "Spatial")))

print(validation_table_display)

write_csv(
  validation_table_display,
  here("outputs", "tables", "validation_metric_table.csv")
)


# Section 5: West Africa Suitability Map with EPPO overlay
maxent_wa <- rast(
  here("data", "processed", "maxent_wa_suitability.tif")
)

# The gpkg carries the eppo_category column built in 03_analysis.R, so the boundary layer
# and the EPPO overlay are the same object read once.
west_africa_boundary <- st_read(
  here("data", "processed", "west_africa_eppo.gpkg"),
  quiet = TRUE
)

# Saint Helena is in Natural Earth's Western Africa subregion but sits far out in the
# Atlantic. Dropping it is a plotting decision only; it was never in the model extent.
west_africa_plot_boundary <- west_africa_boundary %>%
  filter(iso_a2 != "SH")

maxent_wa_plot <- crop(
  maxent_wa,
  vect(west_africa_plot_boundary)
)

# ggplot cannot draw a SpatRaster directly, so the raster is converted to a data frame of
# cell centres. na.rm drops the masked ocean cells rather than plotting them as grey.
maxent_wa_df <- as.data.frame(
  maxent_wa_plot,
  xy = TRUE,
  na.rm = TRUE
)

# The overlay works by mapping country outline colour to EPPO status. Fill is already
# spoken for by the suitability surface, so colour is the only channel left, and an
# outline reads cleanly over a continuous fill. The comparison this invites is the point:
# where the model predicts high suitability inside a country EPPO has never recorded, that
# is either a surveillance gap or a false positive, and the map cannot tell you which.
west_africa_suitability_plot <- ggplot() +
  geom_tile(
    data = maxent_wa_df,
    aes(
      x = x,
      y = y,
      fill = MaxEnt_suitability
    ),
    width = res(maxent_wa)[1],
    height = res(maxent_wa)[2]
  ) +
  geom_sf(
    data = west_africa_plot_boundary,
    aes(colour = eppo_category),
    fill = NA,
    linewidth = 0.6
  ) +
  scale_fill_viridis_c(
    name = "Climatic suitability",
    limits = c(0, 1)
  ) +
  scale_colour_manual(
    name = "EPPO recorded distribution",
    values = c(
      "EPPO documented present" = "#d7191c",
      "No EPPO record" = "grey30"
    )
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.2))
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Climatic suitability for Bactrocera dorsalis in West Africa",
    x = "Longitude",
    y = "Latitude",
    caption = "Saint Helena excluded from map extent for legibility. No EPPO record means no report, not confirmed absence."
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

print(west_africa_suitability_plot)

ggsave(
  west_africa_suitability_plot,
  filename = here("outputs", "maps", "west_africa_suitability_plot.png"),
  width = 9,
  height = 6,
  dpi = 300
)


# Section 6: Ghana Suitability Map with EPPO overlay
maxent_ghana <- rast(
  here("data", "processed", "maxent_ghana_suitability.tif")
)

ghana_boundary <- st_read(
  here("data", "processed", "west_africa_eppo.gpkg"),
  quiet = TRUE
) %>%
  filter(admin == "Ghana")

maxent_ghana_df <- as.data.frame(
  maxent_ghana,
  xy = TRUE,
  na.rm = TRUE
)

# Same overlay logic as Section 5. With one country the legend carries a single entry,
# which is the intended reading: it states Ghana's EPPO status on the face of the map.
ghana_suitability_plot <- ggplot() +
  geom_tile(
    data = maxent_ghana_df,
    aes(
      x = x,
      y = y,
      fill = MaxEnt_suitability
    ),
    width = res(maxent_ghana)[1],
    height = res(maxent_ghana)[2]
  ) +
  geom_sf(
    data = ghana_boundary,
    aes(colour = eppo_category),
    fill = NA,
    linewidth = 0.7
  ) +
  scale_fill_viridis_c(
    name = "Climatic suitability",
    limits = c(0, 1)
  ) +
  scale_colour_manual(
    name = "EPPO recorded distribution",
    values = c(
      "EPPO documented present" = "#d7191c",
      "No EPPO record" = "grey30"
    )
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.2))
  ) +
  coord_sf(
    expand = FALSE
  ) +
  labs(
    title = "Climatic suitability for Bactrocera dorsalis in Ghana",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(ghana_suitability_plot)

ggsave(
  here("outputs", "maps", "ghana_maxent_suitability.png"),
  ghana_suitability_plot,
  width = 7,
  height = 7,
  dpi = 300
)


# Section 7: Create Response Curves
response_curves <- read_csv(
  here("data", "processed", "maxent_response_curves.csv")
)

# WorldClim layer names carry no meaning to a reader. Recoding to BIO codes with units
# makes each panel self-explanatory, and the explicit factor levels put the temperature
# variables before the precipitation ones instead of ordering them BIO1, BIO13, BIO14.
response_curves <- response_curves %>%
  mutate(
    predictor = recode(
      predictor,
      "wc2.1_2.5m_bio_1"  = "BIO1: Annual mean temperature (°C)",
      "wc2.1_2.5m_bio_2"  = "BIO2: Mean diurnal range (°C)",
      "wc2.1_2.5m_bio_3"  = "BIO3: Isothermality (%)",
      "wc2.1_2.5m_bio_8"  = "BIO8: Mean temperature of wettest quarter (°C)",
      "wc2.1_2.5m_bio_9"  = "BIO9: Mean temperature of driest quarter (°C)",
      "wc2.1_2.5m_bio_13" = "BIO13: Precipitation of wettest month (mm)",
      "wc2.1_2.5m_bio_14" = "BIO14: Precipitation of driest month (mm)",
      "wc2.1_2.5m_bio_15" = "BIO15: Precipitation seasonality",
      "wc2.1_2.5m_bio_18" = "BIO18: Precipitation of warmest quarter (mm)",
      "wc2.1_2.5m_bio_19" = "BIO19: Precipitation of coldest quarter (mm)"
    ),
    predictor = factor(
      predictor,
      levels = c(
        "BIO1: Annual mean temperature (°C)",
        "BIO2: Mean diurnal range (°C)",
        "BIO3: Isothermality (%)",
        "BIO8: Mean temperature of wettest quarter (°C)",
        "BIO9: Mean temperature of driest quarter (°C)",
        "BIO13: Precipitation of wettest month (mm)",
        "BIO14: Precipitation of driest month (mm)",
        "BIO15: Precipitation seasonality",
        "BIO18: Precipitation of warmest quarter (mm)",
        "BIO19: Precipitation of coldest quarter (mm)"
      )
    )
  )

# Free scales on both axes. A shared x axis would compress the temperature panels against
# precipitation ranges running to thousands of millimetres and flatten every curve.
response_curve_plot <- ggplot(
  response_curves,
  aes(
    x = predictor_value,
    y = suitability
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ predictor,
    scales = "free",
    ncol = 2
  ) +
  labs(
    title = "MaxEnt climatic response curves for Bactrocera dorsalis",
    subtitle = "Axes vary among predictors to emphasise response shape",
    x = "Predictor value",
    y = "Predicted climatic suitability"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(response_curve_plot)

ggsave(
  here("outputs", "figures", "maxent_response_curves.png"),
  response_curve_plot,
  width = 10,
  height = 12,
  dpi = 300
)


# Section 8: Create Variable Contribution Table
variable_contribution <- read_csv(
  here("data", "processed", "maxent_variable_contribution.csv")
)

# The rows arrive already sorted by contribution from 03_analysis.R, so row_number() is
# the rank. Rounding to one decimal is as much precision as a heuristic measure supports.
variable_contribution_table <- variable_contribution %>%
  mutate(
    Predictor = recode(
      predictor,
      "wc2.1_2.5m_bio_1"  = "BIO1: Annual mean temperature",
      "wc2.1_2.5m_bio_2"  = "BIO2: Mean diurnal range",
      "wc2.1_2.5m_bio_3"  = "BIO3: Isothermality",
      "wc2.1_2.5m_bio_8"  = "BIO8: Temperature of wettest quarter",
      "wc2.1_2.5m_bio_9"  = "BIO9: Temperature of driest quarter",
      "wc2.1_2.5m_bio_13" = "BIO13: Wettest-month precipitation",
      "wc2.1_2.5m_bio_14" = "BIO14: Driest-month precipitation",
      "wc2.1_2.5m_bio_15" = "BIO15: Precipitation seasonality",
      "wc2.1_2.5m_bio_18" = "BIO18: Warmest-quarter precipitation",
      "wc2.1_2.5m_bio_19" = "BIO19: Coldest-quarter precipitation"
    ),
    Contribution = round(contribution, 1),
    Rank = row_number()
  ) %>%
  select(
    Rank,
    Predictor,
    Contribution
  )

print(variable_contribution_table)

# Should come to roughly 100. A large shortfall would mean the row-name matching in
# 03_analysis.R missed a predictor.
print(
  paste(
    "Total contribution:",
    sum(variable_contribution_table$Contribution),
    "%"
  )
)

write_csv(
  variable_contribution_table,
  here("outputs", "tables", "maxent_variable_contribution.csv")
)
