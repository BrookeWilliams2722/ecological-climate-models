# using boosted regression trees to predict
# distribution probabilities of 15 eucalyptus species across NSW
# following Elith et al. 2008, A working guide to boosted regression trees 

# load required libraries
library(tidyverse)
library(gbm)
library(dismo)
library(reshape2)
library(raster)
library(rgdal)

# set working directory

setwd("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models")

# load Elith et al modified functions from 'gbm' original functions   

source("input/brt.functions.R")

################################################################################
# step 1 - get inputs ready for model fitting 
################################################################################

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

# copy data as.data.frame
data_trees <- as.data.frame(Data_Final)

# turn all species names to integer 
data_trees[1:15] <- lapply(data_trees[1:15], as.integer)

#########################################################################
### analyse input data - check total number of observations per spp ----
no_trees <- Data_Final %>% 
  dplyr::select(all_of(Species$PATNLabel)) %>% 
  gather(key = "spp", value = "observation") %>% 
  subset(observation == 1) %>% 
  group_by(spp) %>% 
  summarise(total_obs = n())
###----

######################################################################################################
# step 2 - optimise the number of trees required in the gbm to obtain the minimum predictive deviance
#          test different tree complexity and learning rate parameters following Elith et al. 2008
######################################################################################################

# the list of species 

# 1 Corygumm	Corymbia gummifera
# 2 Eucabanc	Eucalyptus bancroftii
# 3 Eucabosi	Eucalyptus bosistoana
# 4 Eucadean	Eucalyptus deanei
# 5 Eucaeuge	Eucalyptus eugenioides
# 6 Eucagran	Eucalyptus grandis
# 7 Eucalong	Eucalyptus longifolia
# 8 Eucapani	Eucalyptus paniculata
# 9 Eucaparr	Eucalyptus parramattensis
# 10 Eucaprop	Eucalyptus propinqua
# 11 Eucaresi	Eucalyptus resinifera
# 12 Eucarobu	Eucalyptus robusta
# 13 Eucasali	Eucalyptus saligna
# 14 Eucatric	Eucalyptus tricarpa
# 15 Melaquin	Melaleuca quinquenervia

# tests for each spp ----
# testing for sp 1 
argssp1 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp1 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33, # the variables
                           gbm.y = 1,     # the species 
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp1 <- pmap(argssp1, gbm.step.presp1)

# testing for sp 2 
argssp2 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp2 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 2,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp2 <- pmap(argssp2, gbm.step.presp2)

# testing for sp 3 
argssp3 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp3 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 3,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp3 <- pmap(argssp3, gbm.step.presp3)


# testing for sp 4 
argssp4 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp4 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 4,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp4 <- pmap(argssp4, gbm.step.presp4)

# testing for sp 5 
argssp5 <- list(tree.complexity = c(5,5,4,4), learning.rate = c(0.05,0.01,0.05,0.01)) 
gbm.step.presp5 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 5,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp5 <- pmap(argssp5, gbm.step.presp5)

# testing for sp 6 
argssp6 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp6 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 6,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp6 <- pmap(argssp6, gbm.step.presp6)

# testing for sp 7 
argssp7 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp7 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 7,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp7 <- pmap(argssp7, gbm.step.presp7)

# testing for sp 8 
argssp8 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp8 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 8,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp8 <- pmap(argssp8, gbm.step.presp8)

# testing for sp 9 
argssp9 <- list(tree.complexity = c(5,5), learning.rate = c(0.01,0.005)) 
gbm.step.presp9 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 9,
                           family = "bernoulli",
                           bag.fraction = 0.75)
results_sp9 <- pmap(argssp9, gbm.step.presp9)

# testing for sp 10 
argssp10 <- list(tree.complexity = c(5,5,4,4), learning.rate = c(0.05,0.01,0.05,0.01)) 
gbm.step.presp10 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 10,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp10 <- pmap(argssp10, gbm.step.presp10)

# testing for sp 11 
argssp11 <- list(tree.complexity = c(5,5,4,4), learning.rate = c(0.05,0.01,0.05,0.01)) 
gbm.step.presp11 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 11,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp11 <- pmap(argssp11, gbm.step.presp11)

# testing for sp 12 
argssp12 <- list(tree.complexity = c(5,5), learning.rate = c(0.01,0.005)) 
gbm.step.presp12 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 12,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp12 <- pmap(argssp12, gbm.step.presp12)

# testing for sp 13 
argssp13 <- list(tree.complexity = c(5,5,4,4), learning.rate = c(0.05,0.01,0.05,0.01)) 
gbm.step.presp13 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 13,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp13 <- pmap(argssp13, gbm.step.presp13)

# testing for sp 14 
argssp14 <- list(tree.complexity = c(5,5), learning.rate = c(0.01,0.005)) 
gbm.step.presp14 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 14,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp14 <- pmap(argssp14, gbm.step.presp14)

# testing for sp 15 
argssp15 <- list(tree.complexity = c(5,5,4,4), learning.rate = c(0.05,0.01,0.05,0.01)) 
gbm.step.presp15 <- partial(gbm.step,
                            data = data_trees,
                            gbm.x = 16:33,
                            gbm.y = 15,
                            family = "bernoulli",
                            bag.fraction = 0.75)
results_sp15 <- pmap(argssp15, gbm.step.presp15)

# ----

# create a list with all the results per spp

results_final <- list(results_sp1, results_sp2, results_sp3, results_sp4, results_sp5, results_sp6, results_sp7, results_sp8, results_sp9,
                      results_sp10, results_sp11, results_sp12, results_sp13, results_sp14, results_sp15)

# save all results as R object and load again (to free the global environment)

save(results_final, file = "output/gbm/results_final.RData")
load("output/gbm/results_final.RData")

############################################################################################################################
######## analyse results - extract $best.trees, $tree.complexity, $learning.rate, $deviance.mean, $deviance.se 
########                   from all tests to identify the learning rate and tree complexity that delivered the lowest deviance 

main_statistics_1_2 <- melt(lapply(lapply(lapply(results_final,"[[", 1),"[[",c("gbm.call")),"[",c("best.trees","tree.complexity","learning.rate"))) %>% 
  pivot_wider(names_from = L2, values_from = value) %>% 
  dplyr::select(-L1) %>% 
  bind_cols(melt(lapply(lapply(lapply(results_final,"[[", 1),"[[",c("cv.statistics")),"[",c("deviance.mean","deviance.se"))) %>% 
  pivot_wider(names_from = L2, values_from = value)) %>% 
  relocate(L1)

# other tests

 test <- results_final %>% flatten() %>% flatten() 
 test[names(test) %in% c("gbm.call","cv.statistics")]
 
 test.2 <- as.data.frame(test[["gbm.call"]][c("best.trees", "tree.complexity", "learning.rate")]) # this gets the data I want, but not for all the list 
 
 
 results_test <- data.frame()
 
 for (i in 1:length(test)) {
   results_test() <- test[i][["gbm.call"]][c("best.trees", "tree.complexity", "learning.rate")]
          }
   
 lapply(test, function(x) x[["gbm.call"]][c("best.trees", "tree.complexity", "learning.rate")])
 
 lapply(lapply(lapply(test,"[[", 1),"[[",c("gbm.call")),"[",c("best.trees","tree.complexity","learning.rate"))


# the not great way ----


# sp 1
test.results.sp1.gbm <- as.data.frame(results_final[[1]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[1]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
 # mutate(spp = 1)
 
test.results.sp1.cvstats <- as.data.frame(results_final[[1]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[1]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp1.gbm)

#sp 2

test.results.sp2.gbm <- as.data.frame(results_final[[2]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[2]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name" )])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = 2)
test.results.sp2.cvstats <- as.data.frame(results_final[[2]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[2]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[2]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp2.gbm)


#sp 3

test.results.sp3.gbm <- as.data.frame(results_final[[3]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[3]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = 3)
test.results.sp3.cvstats <- as.data.frame(results_final[[3]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[3]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[3]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp3.gbm)

# sp 4


test.results.sp4.gbm <- as.data.frame(results_final[[4]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[4]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = 4)
test.results.sp4.cvstats <- as.data.frame(results_final[[4]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[4]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[4]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp4.gbm)


# sp 5


test.results.sp5.gbm <- as.data.frame(results_final[[5]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[5]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = 5)
test.results.sp5.cvstats <- as.data.frame(results_final[[5]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[5]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[5]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp5.gbm)

# sp 6


test.results.sp6.gbm <- as.data.frame(results_final[[6]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[6]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = 6)
test.results.sp6.cvstats <- as.data.frame(results_final[[6]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[6]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[6]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp6.gbm)

# sp 7


test.results.sp7.gbm <- as.data.frame(results_final[[7]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[7]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[7]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[7]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 7)
test.results.sp7.cvstats <- as.data.frame(results_final[[7]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[7]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[7]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[7]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[7]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp7.gbm)

# sp 8


test.results.sp8.gbm <- as.data.frame(results_final[[8]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[8]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[8]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[8]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[8]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 8)
test.results.sp8.cvstats <- as.data.frame(results_final[[8]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[8]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[8]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[8]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[8]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp8.gbm)

# sp 9


test.results.sp9.gbm <- as.data.frame(results_final[[9]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[9]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[9]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[9]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[9]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 9)
test.results.sp9.cvstats <- as.data.frame(results_final[[9]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[9]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[9]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[9]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[9]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp9.gbm)

# sp 10


test.results.sp10.gbm <- as.data.frame(results_final[[10]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[10]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[10]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 10)
test.results.sp10.cvstats <- as.data.frame(results_final[[10]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[10]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[10]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp10.gbm)

# sp 11


test.results.sp11.gbm <- as.data.frame(results_final[[11]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[11]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[11]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[11]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[11]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 11)
test.results.sp11.cvstats <- as.data.frame(results_final[[11]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[11]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[11]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[11]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[11]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp11.gbm)

# sp 12


test.results.sp12.gbm <- as.data.frame(results_final[[12]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[12]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[12]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 12)
test.results.sp12.cvstats <- as.data.frame(results_final[[12]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[12]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[12]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp12.gbm)

# sp 13


test.results.sp13.gbm <- as.data.frame(results_final[[13]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[13]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[13]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[13]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[13]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 13)
test.results.sp13.cvstats <- as.data.frame(results_final[[13]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[13]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[13]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[13]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[13]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp13.gbm)

# sp 14


test.results.sp14.gbm <- as.data.frame(results_final[[14]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[14]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[14]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[[14]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[14]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 14)
test.results.sp14.cvstats <- as.data.frame(results_final[[14]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[14]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[14]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[[14]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[14]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp14.gbm)

# sp 15


test.results.sp15.gbm <- as.data.frame(results_final[[15]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[[15]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[[15]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = 15)
test.results.sp15.cvstats <- as.data.frame(results_final[[15]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[[15]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[[15]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp15.gbm)

main_statistics_allspp <- rbind(test.results.sp1.cvstats, test.results.sp2.cvstats, test.results.sp3.cvstats, test.results.sp4.cvstats,
                                test.results.sp5.cvstats, test.results.sp6.cvstats, test.results.sp7.cvstats, test.results.sp8.cvstats,
                                test.results.sp9.cvstats, test.results.sp10.cvstats, test.results.sp11.cvstats, test.results.sp12.cvstats,
                                test.results.sp13.cvstats, test.results.sp14.cvstats, test.results.sp15.cvstats)

# ----


# find the minimum trees, lr and tc based on minimum deviance.mean per spp

optimal <- main_statistics_allspp %>% 
  group_by(response.name) %>% 
  summarise(deviance.mean = min(deviance.mean)) 

optimal_results <- left_join(optimal, main_statistics_allspp, by = c("response.name", "deviance.mean"))

# save summary of the tree complexity and learning rate that result in the minimum deviance 

save(optimal_results, file = "output/gbm/optimal_results.RData")
write.csv(optimal_results, file = "output/gbm/optimal_results.csv")

# extract optimal models from results_final

eucabanc <- results_final[[1]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucabanc, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabanc.RData")
eucabosi <- results_final[[2]][[1]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucabosi, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabosi.RData")
eucadean <- results_final[[3]][[1]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucadean, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucadean.RData")
eucagran <- results_final[[4]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucagran, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucagran.RData")
eucalong <- results_final[[5]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucalong, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucalong.RData")
eucapani <- results_final[[6]][[3]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucapani, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucapani.RData")
corygumm <- results_final[[7]][[8]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(corygumm, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\corygumm.RData")
eucaprop <- results_final[[8]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucaprop, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaprop.RData")
eucaresi <- results_final[[9]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucaresi, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaresi.RData")
eucarobu <- results_final[[10]][[1]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucarobu, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucarobu.RData")
eucasali <- results_final[[11]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucasali, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucasali.RData")
eucatric <- results_final[[12]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucatric, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucatric.RData")
melaquin <- results_final[[13]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(melaquin, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\melaquin.RData")
eucaeuge <- results_final[[14]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucaeuge, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaeuge.RData")
eucaparr <- results_final[[15]][[2]]#[["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]
save(eucaparr, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaparr.RData")

# to obtain the importance of the different variables 
# example
summary(eucabanc) # this creates a graph
eucabanc[["contributions"]] # this prints the contributions of each variable


# assess change in predictive deviance when dropping variables  

eucabanc.simp <- gbm.simplify(eucabanc, n.drops = 10)
save(eucabanc.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabanc.simp.RData")
eucabosi.simp <- gbm.simplify(eucabosi, n.drops = 10)
save(eucabosi.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabosi.simp.RData")
eucadean.simp <- gbm.simplify(eucadean, n.drops = 10)
save(eucadean.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucadean.simp.RData")
eucagran.simp <- gbm.simplify(eucagran, n.drops = 10)
save(eucagran.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucagran.simp.RData")
eucalong.simp <- gbm.simplify(eucalong, n.drops = 10)
save(eucalong.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucalong.simp.RData")
eucapani.simp <- gbm.simplify(eucapani, n.drops = 10)
save(eucapani.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucapani.simp.RData")
corygumm.simp <- gbm.simplify(corygumm, n.drops = 10)
save(corygumm.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\corygumm.simp.RData")
eucaprop.simp <- gbm.simplify(eucaprop, n.drops = 10)
save(eucaprop.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaprop.simp.RData")
eucaresi.simp <- gbm.simplify(eucaresi, n.drops = 10)
save(eucaresi.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaresi.simp.RData")
eucarobu.simp <- gbm.simplify(eucarobu, n.drops = 10)
save(eucarobu.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucarobu.simp.RData")
eucasali.simp <- gbm.simplify(eucasali, n.drops = 10)
save(eucasali.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucasali.simp.RData")
eucatric.simp <- gbm.simplify(eucatric, n.drops = 10)
save(eucatric.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucatric.simp.RData")
melaquin.simp <- gbm.simplify(melaquin, n.drops = 10)
save(melaquin.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\melaquin.simp.RData")
eucaeuge.simp <- gbm.simplify(eucaeuge, n.drops = 10)
save(eucaeuge.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaeuge.simp.RData")
eucaparr.simp <- gbm.simplify(eucaparr, n.drops = 10)
save(eucaparr.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaparr.simp.RData")

# re-run models with variables dropped 
# based on parameters in  "R:\KPRIVATE19-A2212\analysis\ecological_climate_models\output\gbm\optimal_parameters.xlsx.csv"


corygumm.simp.simp <- gbm.step(data_trees, gbm.x =
                                 corygumm.simp$pred.list[[1]], gbm.y = 1, tree.complexity = 5, learning.rate =
                                 0.05)
save(corygumm.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\corygumm.simp.simp.RData")
eucabanc.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucabanc.simp$pred.list[[4]], gbm.y = 2, tree.complexity = 4, learning.rate =
                                 0.005)
save(eucabanc.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabanc.simp.simp.RData")
eucabosi.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucabosi.simp$pred.list[[9]], gbm.y = 3, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucabosi.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucabosi.simp.simp.RData")
eucadean.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucadean.simp$pred.list[[1]], gbm.y = 4, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucadean.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucadean.simp.simp.RData")
eucagran.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucagran.simp$pred.list[[1]], gbm.y = 6, tree.complexity = 5, learning.rate =
                                 0.005)
save(eucagran.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucagran.simp.simp.RData")
eucalong.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucalong.simp$pred.list[[4]], gbm.y = 7, tree.complexity = 5, learning.rate =
                                 0.005)
save(eucalong.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucalong.simp.simp.RData")
eucapani.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucapani.simp$pred.list[[9]], gbm.y = 8, tree.complexity = 4, learning.rate =
                                 0.01)
save(eucapani.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucapani.simp.simp.RData")
eucaprop.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucaprop.simp$pred.list[[2]], gbm.y = 10, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucaprop.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaprop.simp.simp.RData")
eucaresi.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucaresi.simp$pred.list[[2]], gbm.y = 11, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucaresi.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaresi.simp.simp.RData")
eucaeuge.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucaeuge.simp$pred.list[[1]], gbm.y = 5, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucaeuge.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucaeuge.simp.simp.RData")
eucarobu.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucarobu.simp$pred.list[[1]], gbm.y = 12, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucarobu.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucarobu.simp.simp.RData")
eucasali.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucasali.simp$pred.list[[1]], gbm.y = 13, tree.complexity = 5, learning.rate =
                                 0.01)
save(eucasali.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucasali.simp.simp.RData")
eucatric.simp.simp <- gbm.step(data_trees, gbm.x =
                                 eucatric.simp$pred.list[[2]], gbm.y = 14, tree.complexity = 5, learning.rate =
                                 0.005)
save(eucatric.simp.simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\eucatric.simp.simp.RData")


# load simplified models

load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\corygumm.simp.simp.RData")                             
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucabanc.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucabosi.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucadean.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucaeuge.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucagran.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucalong.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucapani.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucaparr.RData") # no variables to drop for eucaparr 
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucaprop.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucaresi.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucarobu.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucasali.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\eucatric.simp.simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\simplified_models\\melaquin.RData") # no variables to drop for melaquin 

# make predictions using 'predict' function form raster package

# save all simplified models as a list
models_simp <- list(corygumm.simp.simp,eucabanc.simp.simp,eucabosi.simp.simp,eucadean.simp.simp,eucaeuge.simp.simp,
                    eucagran.simp.simp,eucalong.simp.simp,eucapani.simp.simp,eucaparr,eucaprop.simp.simp,eucaresi.simp.simp,
                    eucarobu.simp.simp,eucasali.simp.simp,eucatric.simp.simp,melaquin)
save(models_simp, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\models_simp.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\models_simp.RData")


# to obtain relative importance of variables

# to obtain the relative contribution of the different variables 

models_simp[[1]][["contributions"]]
test <- models_simp[[10]][["contributions"]]
sum(test$rel.inf)


names(models_simp) <- species_names

contributions_models_simp <- melt(lapply(models_simp, '[', "contributions")) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  pivot_longer(cols = cw_precipwp:dl_strmdstge2, names_to = "variables", values_to = "relative_importance")


top_contributions <- contributions_models_simp %>% 
  dplyr::select(-c(variable, L2)) %>% 
  group_by(L1) %>% 
  slice_max(relative_importance, n=3)

frequencies <- top_contributions %>% 
  group_by(variables) %>% 
  summarise(frequencies = n())


all_frequencies <- contributions_models_simp %>% 
  group_by(variables) %>% 
  subset(!is.na(relative_importance)) %>% 
  summarise(frequencies = n())


# run the predictions for the baseline scenario

# read the tif files (variables) and save them as a raster stack for the baseline scenario

drx <- "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\input\\covariates_18_baseline"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

baseline_list_18 <- stack() 

for(i in 1:length(files.baseline)) { 
  baseline_list_18 <- stack(files.baseline[i], baseline_list_18)
}

save(baseline_list_18, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\baseline_list_18.RData")
load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\baseline_list_18.RData")

# add original data to the simplified model # check this in the code  ----
# angaus.tc5.lr005.simp <- gbm.step(Anguilla_train, gbm.x=angaus.simp$pred.list[[1]], gbm.y=2, tree.complexity=5, learning.rate=0.005)

corygumm.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees                          
eucabanc.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucabosi.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucadean.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucaeuge.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucagran.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucalong.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucapani.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucaparr
eucaprop.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucaresi.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucarobu.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucasali.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
eucatric.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
melaquin
# ----

# run all the predictions for the baseline scenario and save them as a raster stack 

species_names <- c("corygumm","eucabanc","eucabosi","eucadean","eucaeuge","eucagran","eucalong","eucapani","eucaparr","eucaprop","eucaresi",
                      "eucarobu","eucasali","eucatric","melaquin")

baseline_pred <- lapply(models_simp, function(x){predict(baseline_list_18,x,n.trees=x$gbm.call$best.trees, type="response")})

writeRaster(baseline_pred, filename="R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\predictions\\baseline.tif", bylayer=TRUE, suffix = species_names)


# run predictions for the other scenarios

load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\models_simp.RData")

species_names <- c("corygumm","eucabanc","eucabosi","eucadean","eucaeuge","eucagran","eucalong","eucapani","eucaparr","eucaprop","eucaresi",
                   "eucarobu","eucasali","eucatric","melaquin")

drx <- "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\input\\ECHAM_R1_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)]

Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
                "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
                "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R1_2039 <- stack() 

for(i in 1:length(files.baseline_18)) { 
  ECHAM_R1_2039 <- stack(files.baseline_18[i], ECHAM_R1_2039)
}

ECHAM_R1_2039_pred <- lapply(models_simp, function(x){predict(ECHAM_R1_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R1_2039_pred, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\predictions\\ECHAM_R1_2039\\ECHAM_R1_2039_pred.RData")
ECHAM_R1_2039_pred <- stack(ECHAM_R1_2039_pred)
writeRaster(ECHAM_R1_2039_pred, filename="R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\predictions\\ECHAM_R1_2039\\ECHAM_R1_2039.tif", bylayer=TRUE, suffix = species_names)


plot(ECHAM_R1_2039_pred[[1]])


CCCMA_R2_6079_pred <- stack(CCCMA_R2_6079_pred)
writeRaster(CCCMA_R2_6079_pred, filename="R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\predictions\\CCCMA_R2_6079\\CCCMA_R2_6079.tif", bylayer=TRUE, suffix = species_names, 
            overwrite = TRUE)



load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\predictions\\CCCMA_R2_6079\\CCCMA_R2_6079_pred.RData")

