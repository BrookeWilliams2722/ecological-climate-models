# load required libraries
library(tidyverse)
library(parallel)
library(gdal)
library(raster)
library(sp)

# location of functions
#source("functions.r")

# load tree data
TreeData <- read_csv(file = "input/tree_data_from_jill.csv")

# load species list
Species <- read_csv(file = "input/species_list_from_jill.csv")
