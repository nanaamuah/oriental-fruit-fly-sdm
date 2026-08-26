# 03_analysis.R
# The modelling script. Screens predictors for collinearity, builds the presence and
# target-group background table, fits MaxEnt and a ridge GLM, then scores both under
# random and spatially blocked cross-validation. The gap between those two scores is
# the point of the whole project. Everything is written to data/processed for 04 to read.


# Section 1: Load required libraries
library(here)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(usdm)
library(readr)
library(predicts)
library(glmnet)
library(blockCV)
library(lwgeom)
library(modEvA)


# Section 2: Prepare Environmental Predictors
# The full global stack is loaded, not a cropped one. Cropping happens next, but the
# model is fitted on global presences, so the global stack is needed for extraction.
bioclim_data <- rast(list.files(
  here("data", "raw", "climate", "wc2.1_2.5m"),
  pattern = "\\.tif$",
  full.names = TRUE
))
# Refer to 01_downloads.R on verification of bioclim data


# Section 3: Create West Africa boundary
# Natural Earth's "Western Africa" subregion defines the prediction extent. Note this
# includes Saint Helena, which sits far out in the Atlantic and stretches the map extent;
# 04_figures.R drops it for plotting only.
west_africa <- ne_countries(
    scale = "medium",
    returnclass = "sf"
) %>%
  filter(subregion == "Western Africa")

plot(st_geometry(west_africa), main = "West Africa Boundary")

# crop() cuts to the bounding box, mask() sets everything outside the polygons to NA.
# Both are needed; crop alone would leave ocean cells inside the box.
bioclim_wa <- bioclim_data %>%
  crop(vect(west_africa)) %>%
  mask(vect(west_africa))

print(bioclim_wa)


# Section 4: Screen predictors for collinearity
# vifstep drops variables one at a time, recomputing after each removal, until every
# remaining variable has VIF below 10. Screening is done on the West African stack rather
# than globally, because collinearity between bioclim variables is regional: two variables
# that separate cleanly worldwide can be near-identical across West Africa alone.
vif_result <- vifstep(x = bioclim_wa, th = 10)

selected_clim <- vif_result@results$Variables

# The retained set is applied to both stacks: global for fitting, West African for
# prediction. The layers must match in name and order or predict() silently misaligns them.
bioclim_selected <- subset(bioclim_data, selected_clim)

bioclim_wa_selected <- subset(bioclim_wa, selected_clim)

names(bioclim_selected)
names(bioclim_wa_selected)


# Section 5: Load the cleaned occurrence tables from 02_data_cleaning_qc.R
dorsalis_occ <- read_csv(here("data", "processed", "bd_clean.csv"), guess_max = Inf)
teph_occ <- read_csv(here("data", "processed", "teph_clean.csv"), guess_max = Inf)

# Section 6: Extract predictor values at occurrence and background points
# method = "simple" takes the value of the cell the point falls in, with no interpolation.
# That is the right choice here because the records were already thinned to one per cell,
# so the cell is the unit of observation.
bd_env <- extract(
  x = bioclim_selected,
  y = data.frame(decimal_longitude = dorsalis_occ$decimalLongitude,
    decimal_latitude = dorsalis_occ$decimalLatitude),
  method = "simple"
)

print(colSums(is.na(bd_env)))

teph_env <- extract(
  x = bioclim_selected,
  y = data.frame(decimal_longitude = teph_occ$decimalLongitude,
    decimal_latitude = teph_occ$decimalLatitude),
  method = "simple"
)
print(colSums(is.na(teph_env)))

# Expect 0 NAs for all columns in both dorsalis and tephritidae environmental data

# Section 7: Construct the target group background
teph_background <- teph_occ %>%
  filter(!cell %in% dorsalis_occ$cell)

# Any Tephritidae record sharing a cell with a dorsalis presence is dropped. Leaving them
# in would put the same cell on both sides of the response, which is the contaminated
# controls problem: the model would be asked to separate a cell from itself.
nrow(teph_background)
nrow(teph_occ)
print(n_distinct(teph_occ$cell))

# 10,000 background points is the MaxEnt default and is large enough to characterise the
# available climate space without the fit becoming dominated by background density.
set.seed(50)
background_sample <- slice_sample(teph_background, n = 10000)

# connect sampled cells to their climate values

teph_env <- teph_env %>%
  mutate(cell = teph_occ$cell)

background_env <- teph_env %>%
  filter(cell %in% background_sample$cell)

nrow(background_env) # Expect 10,000 rows in the background environmental data


# prepare the presence side
bd_env <- bd_env %>%
  mutate(cell = dorsalis_occ$cell)

# Presences and background are stacked into one table with a 1/0 response. Note what this
# 0 means: not an absence, but a place someone looked and recorded a different fruit fly.
model_data <- bd_env %>%
  select(-ID, -cell) %>%
  mutate(response = 1)

background_cont <- background_env %>%
  select(-ID, -cell) %>%
  mutate(response = 0)

model_data <- bind_rows(model_data, background_cont)

model_data %>%
  nrow()
model_data %>%
  count(response) %>%
  print()

# Expect 10,000 for 0, and 1089 for 1


# Section 8: Fit both models and set up the two validation schemes
# MaxEnt first, then a ridge GLM as a contrast. Two algorithms are used so that the
# validation result cannot be dismissed as a quirk of one fitting procedure.
model_x <- model_data %>%
  select(-response)

model_y <- model_data$response

maxent_fit <- MaxEnt(x = model_x, p = model_y)

show(maxent_fit)

# Ridge (alpha = 0) rather than lasso, because the retained predictors are still
# correlated after VIF screening and ridge shrinks correlated coefficients together
# instead of arbitrarily picking one and zeroing the rest.
glm_fit <- cv.glmnet(
  x = as.matrix(model_x),
  y = model_y,
  family = "binomial",
  alpha = 0 # Ridge regression
)

print(glm_fit$lambda.min)
print(glm_fit$lambda.1se) # picked as the ridge penalty

# blockCV needs coordinates; MaxEnt needs climate values. Rather than carry both in one
# object, a parallel spatial table is built in exactly the same row order, so fold indices
# from one can be used to subset the other. The identical() check below is what guarantees
# that alignment holds.
presence_spatial <- dorsalis_occ %>%
  select(cell, decimalLongitude, decimalLatitude) %>%
  mutate(response = 1)

background_spatial <- background_env %>%
  select(cell) %>%
  left_join(
    select(teph_occ, cell, decimalLongitude, decimalLatitude),
    by = "cell"
  ) %>%
  mutate(response = 0)

model_spatial <- bind_rows(presence_spatial, background_spatial)
glimpse(model_spatial)
count(model_spatial, response) %>%
  print() # Expect 10,000 for 0, and 1089 for 1

print(anyNA(model_spatial)) # Expect FALSE

print(identical(model_spatial$response, model_data$response)) # Expect TRUE

model_spatial_sf <- st_as_sf(
  model_spatial,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

# Reprojected to equal-area metres so that block sizes are real distances.
model_spatial_proj <- st_transform(model_spatial_sf, crs = 6933)

sf_use_s2(FALSE) # Disable s2 geometry for spatial autocorrelation analysis (when turned off, the function will use lwgeom)
# This is the step that decides block size. cv_spatial_autocor fits variograms to the
# data and returns the median range of spatial autocorrelation. Using that number rather
# than a round figure is what makes the blocks defensible; the sensitivity test in
# Section 11 then shows how much the answer depends on it.
auto_cor <- cv_spatial_autocor(
  x = model_spatial_proj,
  column = "response",
  plot = TRUE
)
sf_use_s2(TRUE) # Re-enable s2 geometry for other spatial operations

print(auto_cor[c("range", "range_table")])

# The block size is a number the reader needs in order to judge the blocked design, so it
# is written out rather than left inside an object in the session. cv_spatial_autocor
# returns the range in the units of the projection, which is metres under EPSG:6933.
block_size_record <- tibble::tibble(
  block_size = c("half", "selected", "double"),
  size_m = c(auto_cor$range / 2, auto_cor$range, auto_cor$range * 2)
) %>%
  dplyr::mutate(size_km = round(size_m / 1000, 1))

print(block_size_record)

# presence_bg = TRUE tells blockCV that response 0 is background rather than absence, so
# it balances presences across folds instead of treating the classes symmetrically.
# iteration = 100 lets it try 100 block-to-fold assignments and keep the most balanced.
set.seed(50)
spatial_folds <- cv_spatial(
  model_spatial_proj,
  k = 5,
  size = auto_cor$range,
  selection = "random",
  iteration = 100,
  column = "response",
  presence_bg = TRUE
)

# The conventional scheme, for contrast: folds assigned at random, stratified by response
# so each fold holds roughly the same presence-to-background ratio. Nothing stops a test
# point from sitting next door to a training point, which is exactly the flattery being
# measured.
set.seed(50)
random_folds <- model_data %>%
  group_by(response) %>%
  mutate(fold = sample(rep(1:5, length.out = n()))) %>%
  ungroup()

count(random_folds, fold, response) %>%
  print()


# Section 9: Run MaxEnt cross-validation under both validation schemes
# The same model specification is scored twice, changing only how the data are split.
# Within each fold the classification threshold is taken from the training data
# (max_spec_sens) and then applied to the held-out data. Taking it from the test data
# would let the threshold adapt to the answer and would hide the effect being measured.
random_results <- tibble()

for (i in 1:5) {

  random_train <- random_folds %>%
    filter(fold != i)

  random_test <- random_folds %>%
    filter(fold == i)

  train_x <- random_train %>%
    select(-response, -fold)

  train_y <- random_train$response

  fold_model <- MaxEnt(
    x = train_x,
    p = train_y
  )

  test_x <- random_test %>%
    select(-response, -fold)

  test_pred <- predict(
    object = fold_model,
    x = test_x
  )

  pred_presence <- test_pred[random_test$response == 1]
  pred_background <- test_pred[random_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = random_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_model,
    x = train_x
  )

  train_pred_presence <- train_pred[random_train$response == 1]
  train_pred_background <- train_pred[random_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    validation = "random",
    model = "MaxEnt",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  random_results <- bind_rows(
    random_results,
    fold_result
  )
}

print(random_results)


# Same loop, spatially blocked folds. cv_spatial returns row indices rather than a fold
# column, so the split is done by indexing model_data instead of by filtering.
spatial_results <- tibble()

for (i in 1:5) {

  train_index <- spatial_folds$folds_list[[i]][[1]]
  test_index <- spatial_folds$folds_list[[i]][[2]]

  spatial_train <- model_data[train_index, ]
  spatial_test <- model_data[test_index, ]

  train_x <- spatial_train %>%
    select(-response)

  train_y <- spatial_train$response

  fold_model <- MaxEnt(
    x = train_x,
    p = train_y
  )

  test_x <- spatial_test %>%
    select(-response)

  test_pred <- predict(
    object = fold_model,
    x = test_x
  )

  pred_presence <- test_pred[spatial_test$response == 1]
  pred_background <- test_pred[spatial_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = spatial_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_model,
    x = train_x
  )

  train_pred_presence <- train_pred[spatial_train$response == 1]
  train_pred_background <- train_pred[spatial_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    validation = "spatial",
    model = "MaxEnt",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  spatial_results <- bind_rows(
    spatial_results,
    fold_result
  )
}

print(spatial_results)


# Section 10: Quantify the performance drop
# The headline number. Mean metrics under each scheme, then the difference between them.
# A small drop means the model learned climate. A large drop means it learned where
# entomologists have been working.
maxent_results <- bind_rows(random_results, spatial_results)

maxent_summary <- maxent_results %>%
  group_by(validation) %>%
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    mean_TSS = mean(TSS),
    sd_TSS = sd(TSS),
    mean_Boyce = mean(Boyce),
    sd_Boyce = sd(Boyce)
  )


performance_drop <- maxent_summary %>%
  summarise(
    AUC_drop = mean_AUC[validation == "random"] - mean_AUC[validation == "spatial"],
    TSS_drop = mean_TSS[validation == "random"] - mean_TSS[validation == "spatial"],
    Boyce_drop = mean_Boyce[validation == "random"] - mean_Boyce[validation == "spatial"]
  )


# Section 11: Test sensitivity to block size
# The estimated autocorrelation range is one number with uncertainty around it, so the
# blocked validation is repeated at half and double that size. If the conclusion survives
# all three, it does not rest on the exact block size chosen.

# Half size blocks
set.seed(50)

spatial_folds_half <- cv_spatial(
  x = model_spatial_proj,
  k = 5,
  size = auto_cor$range / 2,
  selection = "random",
  iteration = 100,
  column = "response",
  presence_bg = TRUE
)


spatial_results_half <- tibble()

for (i in 1:5) {

  train_index <- spatial_folds_half$folds_list[[i]][[1]]
  test_index <- spatial_folds_half$folds_list[[i]][[2]]

  spatial_train <- model_data[train_index, ]
  spatial_test <- model_data[test_index, ]

  train_x <- spatial_train %>%
    select(-response)

  train_y <- spatial_train$response

  fold_model <- MaxEnt(
    x = train_x,
    p = train_y
  )

  test_x <- spatial_test %>%
    select(-response)

  test_pred <- predict(
    object = fold_model,
    x = test_x
  )

  pred_presence <- test_pred[spatial_test$response == 1]
  pred_background <- test_pred[spatial_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = spatial_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_model,
    x = train_x
  )

  train_pred_presence <- train_pred[spatial_train$response == 1]
  train_pred_background <- train_pred[spatial_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    block_size = "half",
    validation = "spatial",
    model = "MaxEnt",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  spatial_results_half <- bind_rows(
    spatial_results_half,
    fold_result
  )
}

print(spatial_results_half)


# Double size blocks
set.seed(50)

spatial_folds_double <- cv_spatial(
  x = model_spatial_proj,
  k = 5,
  size = auto_cor$range * 2,
  selection = "random",
  iteration = 100,
  column = "response",
  presence_bg = TRUE
)


spatial_results_double <- tibble()

for (i in 1:5) {

  train_index <- spatial_folds_double$folds_list[[i]][[1]]
  test_index <- spatial_folds_double$folds_list[[i]][[2]]

  spatial_train <- model_data[train_index, ]
  spatial_test <- model_data[test_index, ]

  train_x <- spatial_train %>%
    select(-response)

  train_y <- spatial_train$response

  fold_model <- MaxEnt(
    x = train_x,
    p = train_y
  )

  test_x <- spatial_test %>%
    select(-response)

  test_pred <- predict(
    object = fold_model,
    x = test_x
  )

  pred_presence <- test_pred[spatial_test$response == 1]
  pred_background <- test_pred[spatial_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = spatial_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_model,
    x = train_x
  )

  train_pred_presence <- train_pred[spatial_train$response == 1]
  train_pred_background <- train_pred[spatial_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    block_size = "double",
    validation = "spatial",
    model = "MaxEnt",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  spatial_results_double <- bind_rows(
    spatial_results_double,
    fold_result
  )
}

print(spatial_results_double)


spatial_results_sel <- spatial_results %>%
  mutate(block_size = "selected")

block_sensitivity_results <- bind_rows(
  spatial_results_sel,
  spatial_results_half,
  spatial_results_double
)


block_sensitivity_summary <- block_sensitivity_results %>%
  group_by(block_size) %>%
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    mean_TSS = mean(TSS),
    sd_TSS = sd(TSS),
    mean_Boyce = mean(Boyce),
    sd_Boyce = sd(Boyce)
  )

print(block_sensitivity_summary)


# Section 12: Repeat both validation schemes with the ridge GLM
# glmnet needs a matrix rather than a data frame, and predictions need an explicit lambda
# and type = "response" to come back on the probability scale. The seed is offset per fold
# so glmnet's own internal cross-validation for lambda differs between folds.
glm_random_results <- tibble()

for (i in 1:5) {

  random_train <- random_folds %>%
    filter(fold != i)

  random_test <- random_folds %>%
    filter(fold == i)

  train_x <- random_train %>%
    select(-response, -fold) %>%
    as.matrix()

  train_y <- random_train$response

  set.seed(50 + i)

  fold_glm <- glmnet::cv.glmnet(
    x = train_x,
    y = train_y,
    family = "binomial",
    alpha = 0
  )

  test_x <- random_test %>%
    select(-response, -fold) %>%
    as.matrix()

  test_pred <- predict(
    object = fold_glm,
    newx = test_x,
    s = "lambda.1se",
    type = "response"
  ) %>%
    as.numeric()

  pred_presence <- test_pred[random_test$response == 1]
  pred_background <- test_pred[random_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = random_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_glm,
    newx = train_x,
    s = "lambda.1se",
    type = "response"
  ) %>%
    as.numeric()

  train_pred_presence <- train_pred[random_train$response == 1]
  train_pred_background <- train_pred[random_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    validation = "random",
    model = "GLM",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  glm_random_results <- bind_rows(
    glm_random_results,
    fold_result
  )
}

print(glm_random_results)




glm_spatial_results <- tibble()

for (i in 1:5) {

  train_index <- spatial_folds$folds_list[[i]][[1]]
  test_index <- spatial_folds$folds_list[[i]][[2]]

  spatial_train <- model_data[train_index, ]
  spatial_test <- model_data[test_index, ]

  train_x <- spatial_train %>%
    select(-response) %>%
    as.matrix()

  train_y <- spatial_train$response

  set.seed(50 + i)

  fold_glm <- glmnet::cv.glmnet(
    x = train_x,
    y = train_y,
    family = "binomial",
    alpha = 0
  )

  test_x <- spatial_test %>%
    select(-response) %>%
    as.matrix()

  test_pred <- predict(
    object = fold_glm,
    newx = test_x,
    s = "lambda.1se",
    type = "response"
  ) %>%
    as.numeric()

  pred_presence <- test_pred[spatial_test$response == 1]
  pred_background <- test_pred[spatial_test$response == 0]

  fold_eval <- pa_evaluate(
    p = pred_presence,
    a = pred_background
  )

  auc_value <- fold_eval@stats$auc

  boyce_value <- modEvA::Boyce(
    obs = spatial_test$response,
    pred = test_pred,
    plot = FALSE,
    simplif = TRUE
  )

  train_pred <- predict(
    object = fold_glm,
    newx = train_x,
    s = "lambda.1se",
    type = "response"
  ) %>%
    as.numeric()

  train_pred_presence <- train_pred[spatial_train$response == 1]
  train_pred_background <- train_pred[spatial_train$response == 0]

  train_eval <- pa_evaluate(
    p = train_pred_presence,
    a = train_pred_background
  )

  fold_threshold <- train_eval@thresholds$max_spec_sens

  sensitivity <- mean(pred_presence >= fold_threshold)
  specificity <- mean(pred_background < fold_threshold)

  tss_value <- sensitivity + specificity - 1

  fold_result <- tibble(
    validation = "spatial",
    model = "GLM",
    fold = i,
    n_presence = length(pred_presence),
    n_background = length(pred_background),
    AUC = auc_value,
    threshold = fold_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss_value,
    Boyce = boyce_value
  )

  glm_spatial_results <- bind_rows(
    glm_spatial_results,
    fold_result
  )
}

print(glm_spatial_results)



# Section 13: Combine all folds and summarise
all_cv_results <- bind_rows(
  random_results,
  spatial_results,
  glm_random_results,
  glm_spatial_results
)


cv_summary <- all_cv_results %>%
  group_by(validation, model) %>%
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    mean_TSS = mean(TSS),
    sd_TSS = sd(TSS),
    mean_Boyce = mean(Boyce),
    sd_Boyce = sd(Boyce),
    .groups = "drop"
  )

print(cv_summary)



# Section 14: Predict suitability across West Africa
# Both models are projected onto the West African stack. Values are a relative suitability
# index, not a probability that the fly is present.
maxent_wa_suitability <- predict(
  object = maxent_fit,
  x = bioclim_wa_selected
)

# terra::predict cannot call a cv.glmnet object directly, because glmnet's predict method
# needs newx as a matrix and a lambda choice. This wrapper supplies both.
glm_predict <- function(model, data) {
  as.numeric(
    predict(
      model,
      newx = as.matrix(data),
      s = "lambda.1se",
      type = "response"
    )
  )
}

glm_wa_suitability <- predict(
  bioclim_wa_selected,
  glm_fit,
  fun = glm_predict,
  na.rm = TRUE
)

names(glm_wa_suitability) <- "GLM_suitability"
names(maxent_wa_suitability) <- "MaxEnt_suitability"


print(maxent_wa_suitability)
print(glm_wa_suitability)


# Section 15: Crop the prediction to Ghana
# Cropped from the West African surface rather than refitted, so the Ghana map and the
# regional map show the same model at two scales.

ghana_boundary <- west_africa %>%
  filter(admin == "Ghana")

maxent_ghana_suitability <- maxent_wa_suitability %>%
  crop(vect(ghana_boundary)) %>%
  mask(vect(ghana_boundary))

glm_ghana_suitability <- glm_wa_suitability %>%
  crop(vect(ghana_boundary)) %>%
  mask(vect(ghana_boundary))

print(maxent_ghana_suitability)
print(glm_ghana_suitability)



# Section 16: Response curves and variable contributions
# partialResponse varies one predictor across its range while holding the others at their
# means, which shows the shape of the fitted relationship rather than its importance.
maxent_response_curves <- partialResponse(
  model = maxent_fit,
  data = model_x,
  var = names(model_x),
  nsteps = 50,
  plot = FALSE
)

print(names(maxent_response_curves))
print(head(maxent_response_curves[[1]]))


response_curves_tidy <- bind_rows(
  lapply(names(maxent_response_curves), function(v) {
    tibble(
      predictor = v,
      predictor_value = maxent_response_curves[[v]][[v]],
      suitability = maxent_response_curves[[v]]$p
    )
  })
)

glimpse(response_curves_tidy)


# Percent contribution is pulled out of the MaxEnt results matrix by matching row names
# ending in ".contribution". These are heuristic and depend on the order variables entered
# the fit, so they rank predictors rather than partition explained variance.
maxent_variable_contribution <- tibble(
  predictor = gsub(
    ".contribution",
    "",
    rownames(maxent_fit@results)[grep(
      ".contribution",
      rownames(maxent_fit@results),
      fixed = TRUE
    )]
  ),
  contribution = maxent_fit@results[
    grep(
      ".contribution",
      rownames(maxent_fit@results),
      fixed = TRUE
    ),
    1
  ]
) %>%
  arrange(desc(contribution))

print(maxent_variable_contribution)


# Section 17: Attach EPPO recorded distribution as an independent check
# EPPO is an official phytosanitary record and never entered the model, so it is a genuinely
# external comparison. Countries with no EPPO record are labelled as such rather than as
# absences, since a missing record means nobody reported, not that the pest is absent.
eppo_data <- read_csv(here("data", "raw", "eppo_distribution_data.csv"), guess_max = Inf)

count(eppo_data, Status) %>%
  print()

count(eppo_data, country, sort = TRUE) %>%
  filter(n > 1) %>%
  print()


eppo_wa <- eppo_data %>%
  filter(`country code` %in% west_africa$iso_a2)

west_africa_eppo <- west_africa %>%
  left_join(
    eppo_wa %>%
      select(`country code`, Status),
    by = c("iso_a2" = "country code")
  )

count(west_africa_eppo, Status) %>%
  print()


west_africa_eppo <- west_africa_eppo %>%
  mutate(
    eppo_category = if_else(
      is.na(Status),
      "No EPPO record",
      "EPPO documented present"
    )
  )

west_africa_eppo %>%
  count(eppo_category) %>%
  st_drop_geometry() %>%
  print()



# Section 18: Save analysis outputs
# Everything 04_figures.R needs goes to data/processed. Nothing is written to outputs/
# from here, so modelling stays separate from figure and table generation.

writeRaster(
  maxent_wa_suitability,
  here("data", "processed", "maxent_wa_suitability.tif"),
  overwrite = TRUE
)

writeRaster(
  glm_wa_suitability,
  here("data", "processed", "glm_wa_suitability.tif"),
  overwrite = TRUE
)

writeRaster(
  maxent_ghana_suitability,
  here("data", "processed", "maxent_ghana_suitability.tif"),
  overwrite = TRUE
)

writeRaster(
  glm_ghana_suitability,
  here("data", "processed", "glm_ghana_suitability.tif"),
  overwrite = TRUE
)


write_csv(
  all_cv_results,
  here("data", "processed", "cv_fold_results.csv")
)

write_csv(
  cv_summary,
  here("data", "processed", "cv_summary.csv")
)

write_csv(
  performance_drop,
  here("data", "processed", "performance_drop.csv")
)

write_csv(
  block_size_record,
  here("data", "processed", "block_size_record.csv")
)

write_csv(
  block_sensitivity_results,
  here("data", "processed", "block_sensitivity_results.csv")
)

write_csv(
  block_sensitivity_summary,
  here("data", "processed", "block_sensitivity_summary.csv")
)

write_csv(
  response_curves_tidy,
  here("data", "processed", "maxent_response_curves.csv")
)

write_csv(
  maxent_variable_contribution,
  here("data", "processed", "maxent_variable_contribution.csv")
)


st_write(
  west_africa_eppo,
  here("data", "processed", "west_africa_eppo.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)