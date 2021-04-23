# load required libraries
library(tidyverse)
library(gbm)
library(dismo)

# location of functions
#source("functions.r")

# load tree plot data
TreeData <- read_csv(file = "input/tree_plot_data.csv")

# load species labels
Species_Labels <- read_csv(file = "input/species_labels.csv")

# load list of coastal species
Species <- read_csv(file = "input/tree_species_list.csv")

# join species labels and select
Species <- Species %>% left_join(Species_Labels, by = c("Species" = "SynonymName")) %>% dplyr::select(PATNLabel)

# select columns from plot data and plots to include
# based on the species and 18 covariates

# DYNAMIC CLIMATE VARIABLES
# ce_radseas Radiation of Seasonality: Coefficient of Variation (bio23)
# ct_tempann Annual Mean Temperature (bio1)
# ct_tempiso Isothermality 2/7 (bio3)
# ct_tempmtcp Min Temperature of Coldest Period (bio6)
# cw_precipdp Precipitation of Driest Period (bio14)
# cw_precipwp Precipitation of Wettest Period (bio13)
# DYNAMIC SOIL VARIABLES
# sp_awc000_100 Available water capacity proportionally combined depths from 0 to 100 cm
# so_ph000_100 pH (calcium chloride) proportionally combined depths from 0 to 100 cm
# so_soc000_100 Organic Carbon proportionally combined depths from 0 to 100 cm
# STATIC BIOPHYSICAL VARIABLES
# dl_strmdstge2 Euclidean distance to 2nd order streams and above
# dl_strmdstge6 Euclidean distance to 6th order streams and above
# lf_exp315 Exposure to the NW (low = exposed (drier forests); high = sheltered (moister forests)).
# lf_logre10 Cold air drainage
# lf_tpi0360 Topographic position index using neighbourhood of 360m radius
# lf_tpi2000 Topographic position index using neighbourhood of 2000m radius
# sp_cly000_100prop Clay content proportionally combined depths from 0 to 100 cm
# sp_slt000_100prop Silt content proportionally combined depths from 0 to 100 cm
# sp_snd000_100prop Sand content proportionally combined depths from 0 to 100 cm

# list covariates
Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
                    "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
                    "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

# finalise data for model fitting
Data_Final <- TreeData %>% filter(include_site == 1) %>% dplyr::select(all_of(Species$PATNLabel), all_of(Covariates))

# fit regression trees

# test using GitKraken 

