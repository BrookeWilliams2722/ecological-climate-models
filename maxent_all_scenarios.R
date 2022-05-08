library(raster)
library(sf)
library(dismo)
library(tidyverse)
library(terra)

# set working directory to the main code directory

# setwd(" ")

# Predict to the different climate change scenarios

# load baseline model

load("output/maxent_baseline_model/baseline_model.RData")

# Read data and rename to make all layers consistent

drx <- "input/maxent/data/projections/"

files_old <- list.files(path = drx, pattern = "*p12", full.names = TRUE, recursive = TRUE,
                      include.dirs = TRUE)

files_new <- str_remove_all(files_old, "[_]\\d{4}") #[] means start with \\d digits {} number of digits [0-9] numbers from 0-9 + more than 1
files_new2 <- str_remove_all(files_new, "[_]R") # # https://regenerativetoday.com/a-complete-beginners-guide-to-regular-expressions-in-r/
files_new3 <- str_remove_all(files_new2, "\\CCCMA31[0-9]_")
files_new4 <- str_remove_all(files_new3, "\\ECHAM5[0-9]_")
files_new5 <- str_remove_all(files_new4, "CSIRO[-]MK30[0-9]_")

file.rename(from = files_old, to = files_new5)

# Read all layers again

drx <- "input/maxent/data/projections/"

all_layers <- list.files(path = drx, pattern = "*tif$", full.names = TRUE, recursive = TRUE,
                        include.dirs = TRUE)

# Split data into the different scenarios

files_list <- split(all_layers, dirname(all_layers))

# Obtain scenario names

names_files <- unique(dirname(all_layers)) %>%
                str_remove(., "input/maxent/data/projections/") %>%
                str_replace_all(., "/", "_")

names(files_list) <- names_files

names_CCCMA <- names_files[1:21]
names_CSIRO <- names_files[22:42]
names_ECHAM <- names_files[43:63]

# Divide full list into the different scenarios

files_CCCMA <- files_list[1:21]
files_CSIRO <- files_list[22:42]
files_ECHAM <- files_list[43:63]

# Create lists of raster stacks

rasters_CCCMA <- lapply(files_CCCMA, function(x){stack(x)})
rasters_CSIRO <- lapply(files_CSIRO, function(x){stack(x)})
rasters_ECHAM <- lapply(files_ECHAM, function(x){stack(x)})


# Function to predict to the other scenarios

pred_CCCMA <- lapply(rasters_CCCMA, function(x){raster::predict(x,maxent_tuning_baseline)})
save(pred_CCCMA, file = "output/maxent/predictions/pred_CCCMA.RData")
stack_CCCMA <- stack(pred_CCCMA)
writeRaster(stack_CCCMA, filename="output/maxent/predictions/C.tif", bylayer=TRUE,
            suffix = names_CCCMA)

pred_CSIRO <- lapply(rasters_CSIRO, function(x){raster::predict(x,maxent_tuning_baseline)})
save(pred_CSIRO, file = "output/maxent/predictions/pred_CSIRO.RData")
stack_CSIRO <- stack(pred_CSIRO)
writeRaster(stack_CSIRO, filename="output/maxent/predictions/C.tif", bylayer=TRUE,
            suffix = names_CSIRO)

pred_ECHAM <- lapply(rasters_ECHAM, function(x){raster::predict(x,maxent_tuning_baseline)})
save(pred_ECHAM, file = "output/maxent/predictions/pred_ECHAM.RData")
stack_ECHAM <- stack(pred_ECHAM)
writeRaster(stack_ECHAM, filename="output/maxent/predictions/C.tif", bylayer=TRUE,
            suffix = names_ECHAM)
