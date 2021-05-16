# load required libraries
library(tidyverse)
library(gbm)
library(dismo)
library(reshape2)

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

###check total number of observations per spp ----
no_trees <- Data_Final %>% 
  dplyr::select(all_of(Species$PATNLabel)) %>% 
  gather(key = "spp", value = "observation") %>% 
  subset(observation == 1) %>% 
  group_by(spp) %>% 
  summarise(total_obs = n())
###----

# copy data as.data.frame
data_trees <- as.data.frame(Data_Final)

# turn all species names to integer 
data_trees[1:15] <- lapply(data_trees[1:15], as.integer)

# fit regression trees per spp ----

# testing for sp 1 
argssp1 <- list(tree.complexity = c(5,5,4,4,3,3,2,2,1,1), learning.rate = c(0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005,0.01,0.005)) 
gbm.step.presp1 <- partial(gbm.step,
                           data = data_trees,
                           gbm.x = 16:33,
                           gbm.y = 1,
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

save(results_final, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\results_final.RData")

load("R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\results_final.RData")

# extract $best.trees, $tree.complexity, $learning.rate, $deviance.mean, $deviance.se from all runs

main_statistics_1_2 <- melt(lapply(lapply(lapply(results_final,"[[", 1),"[[",c("gbm.call")),"[",c("best.trees","tree.complexity","learning.rate"))) %>% 
  pivot_wider(names_from = L2, values_from = value) %>% 
  dplyr::select(-L1) %>% 
  bind_cols(melt(lapply(lapply(lapply(results_final,"[[", 1),"[[",c("cv.statistics")),"[",c("deviance.mean","deviance.se"))) %>% 
  pivot_wider(names_from = L2, values_from = value)) %>% 
  relocate(L1)

# the not great way ----


# sp 1
test.results.sp1.gbm <- as.data.frame(results_final[["sp1"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp1"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
 # mutate(spp = "sp1")
 
test.results.sp1.cvstats <- as.data.frame(results_final[["sp1"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp1"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp1.gbm)

#sp 2

test.results.sp2.gbm <- as.data.frame(results_final[["sp2"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name" )])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = "sp2")
test.results.sp2.cvstats <- as.data.frame(results_final[["sp2"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp2"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp2"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp2.gbm)


#sp 3

test.results.sp3.gbm <- as.data.frame(results_final[["sp3"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = "sp3")
test.results.sp3.cvstats <- as.data.frame(results_final[["sp3"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp3"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp3"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp3.gbm)

# sp 4


test.results.sp4.gbm <- as.data.frame(results_final[["sp4"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = "sp4")
test.results.sp4.cvstats <- as.data.frame(results_final[["sp4"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp4"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp4"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp4.gbm)


# sp 5


test.results.sp5.gbm <- as.data.frame(results_final[["sp5"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = "sp5")
test.results.sp5.cvstats <- as.data.frame(results_final[["sp5"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp5"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp5"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp5.gbm)

# sp 6


test.results.sp6.gbm <- as.data.frame(results_final[["sp6"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # mutate(spp = "sp6")
test.results.sp6.cvstats <- as.data.frame(results_final[["sp6"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp6"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp6"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp6.gbm)

# sp 7


test.results.sp7.gbm <- as.data.frame(results_final[["sp7"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp7"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp7"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp7")
test.results.sp7.cvstats <- as.data.frame(results_final[["sp7"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp7"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp7"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp7"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp7"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp7.gbm)

# sp 8


test.results.sp8.gbm <- as.data.frame(results_final[["sp8"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp8"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp8"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp8"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp8")
test.results.sp8.cvstats <- as.data.frame(results_final[["sp8"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp8"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp8"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp8"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp8"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp8.gbm)

# sp 9


test.results.sp9.gbm <- as.data.frame(results_final[["sp9"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp9"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp9"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp9"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp9")
test.results.sp9.cvstats <- as.data.frame(results_final[["sp9"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp9"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp9"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp9"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp9"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp9.gbm)

# sp 10


test.results.sp10.gbm <- as.data.frame(results_final[["sp10"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp10"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp10")
test.results.sp10.cvstats <- as.data.frame(results_final[["sp10"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp10"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp10"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp10.gbm)

# sp 11


test.results.sp11.gbm <- as.data.frame(results_final[["sp11"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp11"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp11"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp11"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp11")
test.results.sp11.cvstats <- as.data.frame(results_final[["sp11"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp11"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp11"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp11"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp11"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp11.gbm)

# sp 12


test.results.sp12.gbm <- as.data.frame(results_final[["sp12"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp12"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp12")
test.results.sp12.cvstats <- as.data.frame(results_final[["sp12"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp12"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp12"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp12.gbm)

# sp 13


test.results.sp13.gbm <- as.data.frame(results_final[["sp13"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp13"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp13"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp13"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp13")
test.results.sp13.cvstats <- as.data.frame(results_final[["sp13"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp13"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp13"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp13"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp13"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp13.gbm)

# sp 14


test.results.sp14.gbm <- as.data.frame(results_final[["sp14"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp14"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp14"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  bind_rows(as.data.frame(results_final[["sp14"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp14")
test.results.sp14.cvstats <- as.data.frame(results_final[["sp14"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp14"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp14"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_rows(as.data.frame(results_final[["sp14"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp14"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp14.gbm)

# sp 15


test.results.sp15.gbm <- as.data.frame(results_final[["sp15"]][[1]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")]) %>% 
  bind_rows(as.data.frame(results_final[["sp15"]][[2]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) #%>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[3]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[4]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[5]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[6]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[7]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[8]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[9]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[10]][["gbm.call"]][c("best.trees","tree.complexity","learning.rate", "response.name")])) %>% 
  # mutate(spp = "sp15")
test.results.sp15.cvstats <- as.data.frame(results_final[["sp15"]][[1]][["cv.statistics"]][c("deviance.mean","deviance.se")]) %>%
  bind_rows(as.data.frame(results_final[["sp15"]][[2]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[3]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[4]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[5]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[6]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[7]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[8]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[9]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  # bind_rows(as.data.frame(results_final[["sp15"]][[10]][["cv.statistics"]][c("deviance.mean","deviance.se")])) %>% 
  bind_cols(test.results.sp15.gbm)

# ----

main_statistics_allspp <- rbind(test.results.sp1.cvstats, test.results.sp2.cvstats, test.results.sp3.cvstats, test.results.sp4.cvstats,
                                test.results.sp5.cvstats, test.results.sp6.cvstats, test.results.sp7.cvstats, test.results.sp8.cvstats,
                                test.results.sp9.cvstats, test.results.sp10.cvstats, test.results.sp11.cvstats, test.results.sp12.cvstats,
                                test.results.sp13.cvstats, test.results.sp14.cvstats, test.results.sp15.cvstats)

# find the minimum trees, lr and tc based on minimum deviance.mean per spp

optimal <- main_statistics_allspp %>% 
  group_by(response.name) %>% 
  summarise(deviance.mean = min(deviance.mean)) 

optimal_results <- left_join(optimal, main_statistics_allspp, by = c("response.name", "deviance.mean"))

# save results 

save(optimal_results, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\optimal_results.RData")

write.csv(optimal_results, file = "R:\\KPRIVATE19-A2212\\analysis\\ecological_climate_models\\output\\gbm\\optimal_results.csv")


