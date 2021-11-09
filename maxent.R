library(raster)

# read baseline data
drx <- "data/baseline"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

baseline <- stack() 

for(i in 1:length(files.baseline)) { 
  baseline <- stack(files.baseline[i], baseline)
}

# testing this 
