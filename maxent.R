library(raster)
library(sf)
library(dismo)
library(tidyverse)
library(terra)

# Read Baseline Data (predictors)
# Linda used 5 predictors:

  # o	Maximum temperature of the warmest month (BIO5)
  # o	Annual precipitation (BIO12)
  # o	Soil pH
  # o	Soil organic carbon
  # o	Available water capacity

# set working directory to the main code directory

# setwd(" ")

drx <- "input/maxent/data/baseline/"

# Change file names to make all layers match across scenarios

file.rename(from = list.files(path = drx, pattern = "*soc", full.names = TRUE),
            to = "input/maxent/data/baseline/soc.tif")
file.rename(from = list.files(path = drx, pattern = "*awc", full.names = TRUE),
            to = "input/maxent/data/baseline/awc.tif")
file.rename(from = list.files(path = drx, pattern = "*ph", full.names = TRUE),
            to = "input/maxent/data/baseline/ph.tif")
file.rename(from = list.files(path = drx, pattern = "*p05", full.names = TRUE),
            to = "input/maxent/data/baseline/p05.tif")
file.rename(from = list.files(path = drx, pattern = "*p12", full.names = TRUE),
            to = "input/maxent/data/baseline/p12.tif")

# Create stack with all rasters

files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

baseline <- stack()

for(i in 1:length(files.baseline)) {
  baseline <- stack(files.baseline[i], baseline)
}

# Read Occurrence Data

occurrence_NSW <- read_sf(dsn = "input/maxent/data", layer = "KoalaRecords_1979_2019_FINAL_40005") %>%
                  select(SF_06) %>%
                  subset(SF_06 == 1)
occ <- st_coordinates(occurrence_NSW)
occ <- as.data.frame(occ)

# Sample Random Background Points,
# selected 50000 following Valavi et al 2021 and Linda's model
# Points include the KMA + coast

# read KMA + coast layer to define extent of background points

kma_coast <- raster("input/maxent/data/extent_background_coast.tif")

bground <- randomPoints(kma_coast, 50000, p = occ,
                        extf=1.1, excludep=TRUE, prob=FALSE,
                        cellnumbers=FALSE, tryf=3, warn=2, lonlatCorrection=TRUE) %>%
                        as_tibble() %>%
                        mutate(occ = 0)

plot(bground)
write.csv(bground,"input/maxent/data/bground.csv")
names(bground) <- c("X", "Y", "occ")

# Witholding a 20% sample for testing, 80% for training

fold <- kfold(occ, k=5)
occtest <- occ[fold == 1, ]
occtrain <- occ[fold != 1, ]

# Combine Presence (training) and Background Points
train_bground <- bind_rows(occtrain %>% mutate(occ = 1), bground)

# Fixing Problems with memory

memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) # Assign memory limit

# Extract environmental data at each occurrence
covars_train <- terra::extract(baseline,
                         train_bground %>% dplyr::select(-occ),
                         method = "bilinear",
                         list = FALSE) %>%
                        as.data.frame()

# Fit a MaxEnt model with the training parameters NO tunING ----
maxent_notuning_baseline <- dismo::maxent(x = covars_train,
                                   p = train_bground$occ,
                                   path = "output/maxent/maxent_baseline_notuning_50000/")

save(maxent_notuning_baseline, file = "/output/maxent/maxent_baseline_notunning_50000/maxent_notuning_baseline.RData")

# Predict to the entire of NSW

predictions_notun <- predict(maxent_notuning_baseline, baseline)
save(predictions_notun, file = "output/maxent/maxent_baseline_notunning_50000/prediction_notunning_baseline.RData")

# Compare observed data to model data

# Evaluate and Find the right threshold
eval_no_tuning <- evaluate(p = occtest,
                 a = bground %>%
                   dplyr::select(X,Y),
                 model = maxent_notuning_baseline,
                 x = baseline)

thresh_no_tuning <- threshold(eval_no_tuning)

# Convert predictions into binary
pred_binary_notun <- overlay(predictions_notun,
                             fun = function(x){
                               ifelse(x <=thresh_no_tuning$spec_sens,NA,1)})

save(pred_binary_notun, file = "output/maxent/maxent_baseline_notunning_50000/pred_binary_notun.RData")


# tunING MaxEnt, using maxent_param function from Valavi et al 2021 ----

nfolds <- ifelse(sum(occtrain$occ) < 10, 2, 5)

param_optim_50000 <- maxent_param(data = train_bground,
                                  k = nfolds,
                                  filepath = "output/maxent/maxent_baseline_model/")

save(param_optim_50000, file = "output/maxent/maxent_baseline_model/param_optim.RData")

# Fit a MaxEnt model with the training parameters tuning

maxent_tuning_baseline <- dismo::maxent(x = covars_train,
                                           p = train_bground$occ,
                                           args = param_optim_50000,
                                           path = "output/maxent/maxent_baseline_model/")

save(maxent_tuning_baseline, file = "output/maxent/maxent_baseline_model/baseline_model.RData")

# Predict to the entire of NSW

predictions_tuned <- predict(maxent_tuning_baseline, baseline)
save(predictions_tuned, file = "output/maxent/predictions/baseline.RData")
writeRaster(predictions_tuned, filename = "output/maxent/predictions/baseline.tif")

# Compare observed data to model data

# Find the right threshold

eval_tuning <- evaluate(p = occtest,
                            a = bground %>%
                              dplyr::select(X,Y),
                            model = maxent_tuning_baseline,
                            x = baseline)

save(eval_tuning, file = "output/maxent/maxent_baseline_model/eval_tuning.RData")

thresh_tuning <- threshold(eval_tuning)

# Convert model data into binary

pred_binary_tun <- overlay(predictions_tuned,
                             fun = function(x){
                               ifelse(x <=thresh_tuning$spec_sens,NA,1)})

writeRaster(pred_binary_tun, filename = "output/maxent/predictions_binary/baseline_bi.tif")
