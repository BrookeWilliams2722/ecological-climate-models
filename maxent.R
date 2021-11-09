library(raster)
library(sf)
library(dismo)
library(tidyverse)


# read baseline data (predictors)  
drx <- "data/baseline"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

baseline <- stack() 

for(i in 1:length(files.baseline)) { 
  baseline <- stack(files.baseline[i], baseline)
}

# read the occurrence data 

occurrence_NSW <- read_sf(dsn = "data", layer = "KoalaRecords_1979_2019_FINAL_40005") %>% 
                  select(SF_06) %>% 
                  subset(SF_06 == 1)

coordinates <- st_coordinates(occurrence_NSW)

coordinates <- as.data.frame(coordinates)


# next -> background points
# next -> run model 