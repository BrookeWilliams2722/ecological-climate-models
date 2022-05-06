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

# set working directory to the main code directory

# setwd(" ")

# load Elith et al modified functions from 'gbm' original functions

source("brt.functions.R")

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
########                   from all tests to identify the learning rate and tree complexity that resulted in the lowest deviance

# flatten nested list to extract the different data

summary_data <- results_final %>% flatten() %>% flatten()
summary_data_gbm <- summary_data[grep("gbm.call", names(summary_data))] %>% flatten() %>% flatten()
summary_data_cv <- summary_data[grep("cv.statistics", names(summary_data))] %>% flatten() %>% flatten()

summary_data_best <- as.data.frame(summary_data_gbm[grep("best.trees", names(summary_data_gbm))])
summary_data_best <- data.frame(t(summary_data_best))

summary_data_tree <- as.data.frame(summary_data_gbm[grep("tree.complexity", names(summary_data_gbm))])
summary_data_tree <- data.frame(t(summary_data_tree))

summary_data_learning <- as.data.frame(summary_data_gbm[grep("learning.rate", names(summary_data_gbm))])
summary_data_learning <- data.frame(t(summary_data_learning))

summary_data_response <- as.data.frame(summary_data_gbm[grep("response.name", names(summary_data_gbm))])
summary_data_response <- data.frame(t(summary_data_response))

summary_data_deviance_mean <- as.data.frame(summary_data_cv[grep("deviance.mean", names(summary_data_cv))])
summary_data_deviance_mean <- data.frame(t(summary_data_deviance_mean))

summary_data_deviance_se <- as.data.frame(summary_data_cv[grep("deviance.se", names(summary_data_cv))])
summary_data_deviance_se<- data.frame(t(summary_data_deviance_se))

# join all data into a single dataframe

summary_complete <- cbind(summary_data_best, summary_data_tree, summary_data_learning, summary_data_response, summary_data_deviance_mean,
                          summary_data_deviance_se)

# find the minimum trees, learning rate and tree ccomplexity based on minimum deviance.mean per spp

optimal <- summary_complete %>%
  group_by(t.summary_data_response.) %>%
  summarise(t.summary_data_deviance_mean. = min(t.summary_data_deviance_mean.))

optimal_results <- left_join(optimal, summary_complete, by = c("t.summary_data_response.", "t.summary_data_deviance_mean."))

# save summary of the tree complexity and learning rate that result in the minimum deviance

save(optimal_results, file = "output/gbm/optimal_results.RData")
write.csv(optimal_results, file = "output/gbm/optimal_results.csv")

# extract optimal models from results_final based on summary above and save them as a separate element

eucabanc <- results_final[[1]][[2]] # first element in the list contains all the runs for x species, second element calls the model with the lowest deviance
save(eucabanc, file = "output/gbm/optimised_complete_models/eucabanc.RData")
eucabosi <- results_final[[2]][[1]]
save(eucabosi, file = "output/gbm/optimised_complete_models/eucabosi.RData")
eucadean <- results_final[[3]][[1]]
save(eucadean, file = "output/gbm/optimised_complete_models/eucadean.RData")
eucagran <- results_final[[4]][[2]]
save(eucagran, file = "output/gbm/optimised_complete_models/eucagran.RData")
eucalong <- results_final[[5]][[2]]
save(eucalong, file = "output/gbm/optimised_complete_models/eucalong.RData")
eucapani <- results_final[[6]][[3]]
save(eucapani, file = "output/gbm/optimised_complete_models/eucapani.RData")
corygumm <- results_final[[7]][[8]]
save(corygumm, file = "output/gbm/optimised_complete_models/corygumm.RData")
eucaprop <- results_final[[8]][[2]]
save(eucaprop, file = "output/gbm/optimised_complete_models/eucaprop.RData")
eucaresi <- results_final[[9]][[2]]
save(eucaresi, file = "output/gbm/optimised_complete_models/eucaresi.RData")
eucarobu <- results_final[[10]][[1]]
save(eucarobu, file = "output/gbm/optimised_complete_models/eucarobu.RData")
eucasali <- results_final[[11]][[2]]
save(eucasali, file = "output/gbm/optimised_complete_models/eucasali.RData")
eucatric <- results_final[[12]][[2]]
save(eucatric, file = "output/gbm/optimised_complete_models/eucatric.RData")
melaquin <- results_final[[13]][[2]]
save(melaquin, file = "output/gbm/optimised_complete_models/melaquin.RData")
eucaeuge <- results_final[[14]][[2]]
save(eucaeuge, file = "output/gbm/optimised_complete_models/eucaeuge.RData")
eucaparr <- results_final[[15]][[2]]
save(eucaparr, file = "output/gbm/optimised_complete_models/eucaparr.RData")

######################################################################################
########## analyse results -  to obtain the importance of the different variables use
##########                    example
summary(eucabanc) # this creates a graph
eucabanc[["contributions"]] # this prints the contributions of each variable

#######################################################################################
# step 3 - assess change in predictive deviance when dropping variables
#######################################################################################

eucabanc.simp <- gbm.simplify(eucabanc, n.drops = 10)
save(eucabanc.simp, file = "output/gbm/drop_variables_models/eucabanc.simp.RData")
eucabosi.simp <- gbm.simplify(eucabosi, n.drops = 10)
save(eucabosi.simp, file = "output/gbm/drop_variables_models/eucabosi.simp.RData")
eucadean.simp <- gbm.simplify(eucadean, n.drops = 10)
save(eucadean.simp, file = "output/gbm/drop_variables_models/eucadean.simp.RData")
eucagran.simp <- gbm.simplify(eucagran, n.drops = 10)
save(eucagran.simp, file = "output/gbm/drop_variables_models/eucagran.simp.RData")
eucalong.simp <- gbm.simplify(eucalong, n.drops = 10)
save(eucalong.simp, file = "output/gbm/drop_variables_models/eucalong.simp.RData")
eucapani.simp <- gbm.simplify(eucapani, n.drops = 10)
save(eucapani.simp, file = "output/gbm/drop_variables_models/eucapani.simp.RData")
corygumm.simp <- gbm.simplify(corygumm, n.drops = 10)
save(corygumm.simp, file = "output/gbm/drop_variables_models/corygumm.simp.RData")
eucaprop.simp <- gbm.simplify(eucaprop, n.drops = 10)
save(eucaprop.simp, file = "output/gbm/drop_variables_models/eucaprop.simp.RData")
eucaresi.simp <- gbm.simplify(eucaresi, n.drops = 10)
save(eucaresi.simp, file = "output/gbm/drop_variables_models/eucaresi.simp.RData")
eucarobu.simp <- gbm.simplify(eucarobu, n.drops = 10)
save(eucarobu.simp, file = "output/gbm/drop_variables_models/eucarobu.simp.RData")
eucasali.simp <- gbm.simplify(eucasali, n.drops = 10)
save(eucasali.simp, file = "output/gbm/drop_variables_models/eucasali.simp.RData")
eucatric.simp <- gbm.simplify(eucatric, n.drops = 10)
save(eucatric.simp, file = "output/gbm/drop_variables_models/eucatric.simp.RData")
melaquin.simp <- gbm.simplify(melaquin, n.drops = 10)
save(melaquin.simp, file = "output/gbm/drop_variables_models/melaquin.simp.RData")
eucaeuge.simp <- gbm.simplify(eucaeuge, n.drops = 10)
save(eucaeuge.simp, file = "output/gbm/drop_variables_models/eucaeuge.simp.RData")
eucaparr.simp <- gbm.simplify(eucaparr, n.drops = 10)
save(eucaparr.simp, file = "output/gbm/drop_variables_models/eucaparr.simp.RData")

#########################################################################################################
################## to see the results, i.e., change in predictive deviance when dropping variables, use:
eucaeuge.simp[["deviance.summary"]]

################## to see the list of variables used in each drop use
eucaparr.simp[["pred.list"]]

#################################################################################################################################
# step 4 - re-run models with the variables that did not decrease predictive deviance (dropped variables)
# based on summary of results in  "R:\KPRIVATE19-A2212\analysis\ecological_climate_models\output\gbm\optimal_parameters.xlsx.csv"
# no variables dropped for eucaparr and melaquin
# if the simplified model does not include the data frame with the original input data 'data_trees', then append it to the gbm,
#     e.g. corygumm.simp.simp[["gbm.call"]][["dataframe"]] <- data_trees
#################################################################################################################################

corygumm.simp.simp <- gbm.step(data_trees, gbm.x = corygumm.simp$pred.list[[1]], # the drop with the variables that did not decrease predictive deviance
                               gbm.y = 1, tree.complexity = 5, learning.rate = 0.05)
save(corygumm.simp.simp, file = "output/gbm/simplified_models/corygumm.simp.simp.RData")

eucabanc.simp.simp <- gbm.step(data_trees, gbm.x = eucabanc.simp$pred.list[[4]], gbm.y = 2, tree.complexity = 4, learning.rate = 0.005)
save(eucabanc.simp.simp, file = "output/gbm/simplified_models/eucabanc.simp.simp.RData")

eucabosi.simp.simp <- gbm.step(data_trees, gbm.x = eucabosi.simp$pred.list[[9]], gbm.y = 3, tree.complexity = 5, learning.rate = 0.01)
save(eucabosi.simp.simp, file = "output/gbm/simplified_models/eucabosi.simp.simp.RData")

eucadean.simp.simp <- gbm.step(data_trees, gbm.x = eucadean.simp$pred.list[[1]], gbm.y = 4, tree.complexity = 5, learning.rate = 0.01)
save(eucadean.simp.simp, file = "output/gbm/simplified_models/eucadean.simp.simp.RData")

eucagran.simp.simp <- gbm.step(data_trees, gbm.x = eucagran.simp$pred.list[[1]], gbm.y = 6, tree.complexity = 5, learning.rate = 0.005)
save(eucagran.simp.simp, file = "output/gbm/simplified_models/eucagran.simp.simp.RData")

eucalong.simp.simp <- gbm.step(data_trees, gbm.x = eucalong.simp$pred.list[[4]], gbm.y = 7, tree.complexity = 5, learning.rate = 0.005)
save(eucalong.simp.simp, file = "output/gbm/simplified_models/eucalong.simp.simp.RData")

eucapani.simp.simp <- gbm.step(data_trees, gbm.x = eucapani.simp$pred.list[[9]], gbm.y = 8, tree.complexity = 4, learning.rate = 0.01)
save(eucapani.simp.simp, file = "output/gbm/simplified_models/eucapani.simp.simp.RData")

eucaprop.simp.simp <- gbm.step(data_trees, gbm.x = eucaprop.simp$pred.list[[2]], gbm.y = 10, tree.complexity = 5, learning.rate = 0.01)
save(eucaprop.simp.simp, file = "output/gbm/simplified_models/eucaprop.simp.simp.RData")

eucaresi.simp.simp <- gbm.step(data_trees, gbm.x = eucaresi.simp$pred.list[[2]], gbm.y = 11, tree.complexity = 5, learning.rate = 0.01)
save(eucaresi.simp.simp, file = "output/gbm/simplified_models/eucaresi.simp.simp.RData")

eucaeuge.simp.simp <- gbm.step(data_trees, gbm.x = eucaeuge.simp$pred.list[[1]], gbm.y = 5, tree.complexity = 5, learning.rate = 0.01)
save(eucaeuge.simp.simp, file = "output/gbm/simplified_models/eucaeuge.simp.simp.RData")

eucarobu.simp.simp <- gbm.step(data_trees, gbm.x = eucarobu.simp$pred.list[[1]], gbm.y = 12, tree.complexity = 5, learning.rate = 0.01)
save(eucarobu.simp.simp, file = "output/gbm/simplified_models/eucarobu.simp.simp.RData")

eucasali.simp.simp <- gbm.step(data_trees, gbm.x = eucasali.simp$pred.list[[1]], gbm.y = 13, tree.complexity = 5, learning.rate = 0.01)
save(eucasali.simp.simp, file = "output/gbm/simplified_models/eucasali.simp.simp.RData")

eucatric.simp.simp <- gbm.step(data_trees, gbm.x = eucatric.simp$pred.list[[2]], gbm.y = 14, tree.complexity = 5, learning.rate = 0.005)
save(eucatric.simp.simp, file = "output/gbm/simplified_models/eucatric.simp.simp.RData")

#################################################################################################################################################
################## analyse results - to obtain summaries of the relative contribution of the different variables in the final, simplified, models

# save all simplified models as a list
models_simp <- list(corygumm.simp.simp,eucabanc.simp.simp,eucabosi.simp.simp,eucadean.simp.simp,eucaeuge.simp.simp,
                    eucagran.simp.simp,eucalong.simp.simp,eucapani.simp.simp,eucaparr,eucaprop.simp.simp,eucaresi.simp.simp,
                    eucarobu.simp.simp,eucasali.simp.simp,eucatric.simp.simp,melaquin)
save(models_simp, file = "output/gbm/models_simp.RData")

# add the names of the different spp to the list with the final, simplified, models
species_names <- c("corygumm","eucabanc","eucabosi","eucadean","eucaeuge","eucagran","eucalong","eucapani","eucaparr","eucaprop","eucaresi",
                   "eucarobu","eucasali","eucatric","melaquin")

names(models_simp) <- species_names

# list of the relative contribution of all variables for all spp

contributions_models_simp <- melt(lapply(models_simp, '[', "contributions")) %>%
  pivot_wider(names_from = var, values_from = value) %>%
  pivot_longer(cols = cw_precipwp:dl_strmdstge2, names_to = "variables", values_to = "relative_importance")

# filter to obtain the variables with the highest importance

top_contributions <- contributions_models_simp %>%
  dplyr::select(-c(variable, L2)) %>%
  group_by(L1) %>%
  slice_max(relative_importance, n=1)

# the number of models for which a top variable was the highest in importance
frequencies <- top_contributions %>%
  group_by(variables) %>%
  summarise(frequencies = n())

# the number of models for which a variable was used in the final, simplified, models
all_frequencies <- contributions_models_simp %>%
  group_by(variables) %>%
  subset(!is.na(relative_importance)) %>%
  summarise(frequencies = n())


##########################################################################
# step 5 - make predictions using 'predict' function form raster package
##########################################################################

######### baseline scenario ############

# read the tif files (variables) and save them as a raster stack

drx <- "input/covariates_18_baseline"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

baseline_list_18 <- stack()

for(i in 1:length(files.baseline)) {
  baseline_list_18 <- stack(files.baseline[i], baseline_list_18)
}

# run the predictions for all spp in the baseline scenario, save them as a list and as .tif

baseline_pred <- lapply(models_simp, function(x){predict(baseline_list_18,x,n.trees=x$gbm.call$best.trees, type="response")})
save(baseline_pred, file = "output/predictions/baseline/baseline_pred.RData")
baseline_pred <- stack(baseline_pred)
writeRaster(baseline_pred, filename="output/predictions/baseline/baseline.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R1_2039 ###############

drx <- "input/CCCMA_R1_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R1_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R1_2039 <- stack(files.baseline_18[i], CCCMA_R1_2039)
}

CCCMA_R1_2039_pred <- lapply(models_simp, function(x){predict(CCCMA_R1_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R1_2039_pred, file = "output/predictions/CCCMA_R1_2039/CCCMA_R1_2039_pred.RData")
CCCMA_R1_2039_pred <- stack(CCCMA_R1_2039_pred)
writeRaster(CCCMA_R1_2039_pred, filename="output/predictions/CCCMA_R1_2039/CCCMA_R1_2039.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R1_6079 ###############

drx <- "input/CCCMA_R1_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R1_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R1_6079 <- stack(files.baseline_18[i], CCCMA_R1_6079)
}

CCCMA_R1_6079_pred <- lapply(models_simp, function(x){predict(CCCMA_R1_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R1_6079_pred, file = "output/predictions/CCCMA_R1_6079/CCCMA_R1_6079_pred.RData")
CCCMA_R1_6079_pred <- stack(CCCMA_R1_6079_pred)
writeRaster(CCCMA_R1_6079_pred, filename="output/predictions/CCCMA_R1_6079/CCCMA_R1_6079.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R2_2039 ###############

drx <- "input/CCCMA_R2_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R2_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R2_2039 <- stack(files.baseline_18[i], CCCMA_R2_2039)
}

CCCMA_R2_2039_pred <- lapply(models_simp, function(x){predict(CCCMA_R2_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R2_2039_pred, file = "output/predictions/CCCMA_R2_2039/CCCMA_R2_2039_pred.RData")
CCCMA_R2_2039_pred <- stack(CCCMA_R2_2039_pred)
writeRaster(CCCMA_R2_2039_pred, filename="output/predictions/CCCMA_R2_2039/CCCMA_R2_2039.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R2_6079 ###############

drx <- "input/CCCMA_R2_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R2_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R2_6079 <- stack(files.baseline_18[i], CCCMA_R2_6079)
}

CCCMA_R2_6079_pred <- lapply(models_simp, function(x){predict(CCCMA_R2_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R2_6079_pred, file = "output/predictions/CCCMA_R2_6079/CCCMA_R2_6079_pred.RData")
CCCMA_R2_6079_pred <- stack(CCCMA_R2_6079_pred)
writeRaster(CCCMA_R2_6079_pred, filename="output/predictions/CCCMA_R2_6079/CCCMA_R2_6079.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R3_2039 ###############

drx <- "input/CCCMA_R3_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R3_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R3_2039 <- stack(files.baseline_18[i], CCCMA_R3_2039)
}

CCCMA_R3_2039_pred <- lapply(models_simp, function(x){predict(CCCMA_R3_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R3_2039_pred, file = "output/predictions/CCCMA_R3_2039/CCCMA_R3_2039_pred.RData")
CCCMA_R3_2039_pred <- stack(CCCMA_R3_2039_pred)
writeRaster(CCCMA_R3_2039_pred, filename="output/predictions/CCCMA_R3_2039/CCCMA_R3_2039.tif", bylayer=TRUE, suffix = species_names)

########### CCCMA_R3_6079 ###############

drx <- "input/CCCMA_R3_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CCCMA_R3_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CCCMA_R3_6079 <- stack(files.baseline_18[i], CCCMA_R3_6079)
}

CCCMA_R3_6079_pred <- lapply(models_simp, function(x){predict(CCCMA_R3_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CCCMA_R3_6079_pred, file = "output/predictions/CCCMA_R3_6079/CCCMA_R3_6079_pred.RData")
CCCMA_R3_6079_pred <- stack(CCCMA_R3_6079_pred)
writeRaster(CCCMA_R3_6079_pred, filename="output/predictions/CCCMA_R3_6079/CCCMA_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R1_6079 ###############

drx <- "input/CSIRO_R1_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R1_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R1_6079 <- stack(files.baseline_18[i], CSIRO_R1_6079)
}

CSIRO_R1_6079_pred <- lapply(models_simp, function(x){predict(CSIRO_R1_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R1_6079_pred, file = "output/predictions/CSIRO_R1_6079/CSIRO_R1_6079_pred.RData")
CSIRO_R1_6079_pred <- stack(CSIRO_R1_6079_pred)
writeRaster(CSIRO_R1_6079_pred, filename="output/predictions/CSIRO_R1_6079/CSIRO_R1_2039.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R1_2039 ###############

drx <- "input/CSIRO_R1_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R1_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R1_2039 <- stack(files.baseline_18[i], CSIRO_R1_2039)
}

CSIRO_R1_2039_pred <- lapply(models_simp, function(x){predict(CSIRO_R1_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R1_2039_pred, file = "output/predictions/CSIRO_R1_2039/CSIRO_R1_2039_pred.RData")
CSIRO_R1_2039_pred <- stack(CSIRO_R1_2039_pred)
writeRaster(CSIRO_R1_2039_pred, filename="output/predictions/CSIRO_R1_2039/CSIRO_R1_2039.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R2_2039 ###############

drx <- "input/CSIRO_R2_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R2_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R2_2039 <- stack(files.baseline_18[i], CSIRO_R2_2039)
}

CSIRO_R2_2039_pred <- lapply(models_simp, function(x){predict(CSIRO_R2_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R2_2039_pred, file = "output/predictions/CSIRO_R2_2039/CSIRO_R2_2039_pred.RData")
CSIRO_R2_2039_pred <- stack(CSIRO_R2_2039_pred)
writeRaster(CSIRO_R2_2039_pred, filename="output/predictions/CSIRO_R2_2039/CSIRO_R2_2039.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R2_6079 ###############

drx <- "input/CSIRO_R2_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R2_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R2_6079 <- stack(files.baseline_18[i], CSIRO_R2_6079)
}

CSIRO_R2_6079_pred <- lapply(models_simp, function(x){predict(CSIRO_R2_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R2_6079_pred, file = "output/predictions/CSIRO_R2_6079/CSIRO_R2_6079_pred.RData")
CSIRO_R2_6079_pred <- stack(CSIRO_R2_6079_pred)
writeRaster(CSIRO_R2_6079_pred, filename="output/predictions/CSIRO_R2_6079/CSIRO_R2_6079.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R3_2039 ###############

drx <- "input/CSIRO_R3_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R3_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R3_2039 <- stack(files.baseline_18[i], CSIRO_R3_2039)
}

CSIRO_R3_2039_pred <- lapply(models_simp, function(x){predict(CSIRO_R3_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R3_2039_pred, file = "output/predictions/CSIRO_R3_2039/CSIRO_R3_2039_pred.RData")
CSIRO_R3_2039_pred <- stack(CSIRO_R3_2039_pred)
writeRaster(CSIRO_R3_2039_pred, filename="output/predictions/CSIRO_R3_2039/CSIRO_R3_2039.tif", bylayer=TRUE, suffix = species_names)

########### CSIRO_R3_6079 ###############

drx <- "input/CSIRO_R3_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

CSIRO_R3_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  CSIRO_R3_6079 <- stack(files.baseline_18[i], CSIRO_R3_6079)
}

CSIRO_R3_6079_pred <- lapply(models_simp, function(x){predict(CSIRO_R3_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(CSIRO_R3_6079_pred, file = "output/predictions/CSIRO_R3_6079/CSIRO_R3_6079_pred.RData")
CSIRO_R3_6079_pred <- stack(CSIRO_R3_6079_pred)
writeRaster(CSIRO_R3_6079_pred, filename="output/predictions/CSIRO_R3_6079/CSIRO_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R1_2039 ###############

drx <- "input/ECHAM_R1_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R1_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R1_2039 <- stack(files.baseline_18[i], ECHAM_R1_2039)
}

ECHAM_R1_2039_pred <- lapply(models_simp, function(x){predict(ECHAM_R1_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R1_2039_pred, file = "output/predictions/ECHAM_R1_2039/ECHAM_R1_2039_pred.RData")
ECHAM_R1_2039_pred <- stack(ECHAM_R1_2039_pred)
writeRaster(ECHAM_R1_2039_pred, filename="output/predictions/ECHAM_R1_2039/ECHAM_R1_2039.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R1_6079 ###############

drx <- "input/ECHAM_R1_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R1_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R1_6079 <- stack(files.baseline_18[i], ECHAM_R1_6079)
}

ECHAM_R1_6079_pred <- lapply(models_simp, function(x){predict(ECHAM_R1_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R1_6079_pred, file = "output/predictions/ECHAM_R1_6079/ECHAM_R1_6079_pred.RData")
ECHAM_R1_6079_pred <- stack(ECHAM_R1_6079_pred)
writeRaster(ECHAM_R1_6079_pred, filename="output/predictions/ECHAM_R1_6079/ECHAM_R1_6079.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R2_2039 ###############

drx <- "input/ECHAM_R2_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R2_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R2_2039 <- stack(files.baseline_18[i], ECHAM_R2_2039)
}

ECHAM_R2_2039_pred <- lapply(models_simp, function(x){predict(ECHAM_R2_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R2_2039_pred, file = "output/predictions/ECHAM_R2_2039/ECHAM_R2_2039_pred.RData")
ECHAM_R2_2039_pred <- stack(ECHAM_R2_2039_pred)
writeRaster(ECHAM_R2_2039_pred, filename="output/predictions/ECHAM_R2_2039/ECHAM_R2_2039.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R2_6079 ###############

drx <- "input/ECHAM_R2_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R2_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R2_6079 <- stack(files.baseline_18[i], ECHAM_R2_6079)
}

ECHAM_R2_6079_pred <- lapply(models_simp, function(x){predict(ECHAM_R2_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R2_6079_pred, file = "output/predictions/ECHAM_R2_6079/ECHAM_R2_6079_pred.RData")
ECHAM_R2_6079_pred <- stack(ECHAM_R2_6079_pred)
writeRaster(ECHAM_R2_6079_pred, filename="output/predictions/ECHAM_R2_6079/ECHAM_R2_6079.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R3_2039 ###############

drx <- "input/ECHAM_R3_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R3_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R3_2039 <- stack(files.baseline_18[i], ECHAM_R3_2039)
}

ECHAM_R3_2039_pred <- lapply(models_simp, function(x){predict(ECHAM_R3_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R3_2039_pred, file = "output/predictions/ECHAM_R3_2039/ECHAM_R3_2039_pred.RData")
ECHAM_R3_2039_pred <- stack(ECHAM_R3_2039_pred)
writeRaster(ECHAM_R3_2039_pred, filename="output/predictions/ECHAM_R3_2039/ECHAM_R3_2039.tif", bylayer=TRUE, suffix = species_names)

########### ECHAM_R3_6079 ###############

drx <- "input/ECHAM_R3_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

ECHAM_R3_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  ECHAM_R3_6079 <- stack(files.baseline_18[i], ECHAM_R3_6079)
}

ECHAM_R3_6079_pred <- lapply(models_simp, function(x){predict(ECHAM_R3_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(ECHAM_R3_6079_pred, file = "output/predictions/ECHAM_R3_6079/ECHAM_R3_6079_pred.RData")
ECHAM_R3_6079_pred <- stack(ECHAM_R3_6079_pred)
writeRaster(ECHAM_R3_6079_pred, filename="output/predictions/ECHAM_R3_6079/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R1_2039 ###############

drx <- "input/MIROC_R1_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R1_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R1_2039 <- stack(files.baseline_18[i], MIROC_R1_2039)
}

MIROC_R1_2039_pred <- lapply(models_simp, function(x){predict(MIROC_R1_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R1_2039_pred, file = "output/predictions/MIROC_R1_2039/MIROC_R1_2039_pred.RData")
MIROC_R1_2039_pred <- stack(MIROC_R1_2039_pred)
writeRaster(MIROC_R1_2039_pred, filename="output/predictions/MIROC_R1_2039/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R1_6079 ###############

drx <- "input/MIROC_R1_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R1_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R1_6079 <- stack(files.baseline_18[i], MIROC_R1_6079)
}

MIROC_R1_6079_pred <- lapply(models_simp, function(x){predict(MIROC_R1_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R1_6079_pred, file = "output/predictions/MIROC_R1_6079/MIROC_R1_6079_pred.RData")
MIROC_R1_6079_pred <- stack(MIROC_R1_6079_pred)
writeRaster(MIROC_R1_6079_pred, filename="output/predictions/MIROC_R1_6079/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R2_2039 ###############

drx <- "input/MIROC_R2_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R2_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R2_2039 <- stack(files.baseline_18[i], MIROC_R2_2039)
}

MIROC_R2_2039_pred <- lapply(models_simp, function(x){predict(MIROC_R2_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R2_2039_pred, file = "output/predictions/MIROC_R2_2039/MIROC_R2_2039_pred.RData")
MIROC_R2_2039_pred <- stack(MIROC_R2_2039_pred)
writeRaster(MIROC_R2_2039_pred, filename="output/predictions/MIROC_R2_2039/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R2_6079 ###############

drx <- "input/MIROC_R2_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R2_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R2_6079 <- stack(files.baseline_18[i], MIROC_R2_6079)
}

MIROC_R2_6079_pred <- lapply(models_simp, function(x){predict(MIROC_R2_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R2_6079_pred, file = "output/predictions/MIROC_R2_6079/MIROC_R2_6079_pred.RData")
MIROC_R2_6079_pred <- stack(MIROC_R2_6079_pred)
writeRaster(MIROC_R2_6079_pred, filename="output/predictions/MIROC_R2_6079/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R3_2039 ###############

drx <- "input/MIROC_R3_2039"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R3_2039 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R3_2039 <- stack(files.baseline_18[i], MIROC_R3_2039)
}

MIROC_R3_2039_pred <- lapply(models_simp, function(x){predict(MIROC_R3_2039,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R3_2039_pred, file = "output/predictions/MIROC_R3_2039/MIROC_R3_2039_pred.RData")
MIROC_R3_2039_pred <- stack(MIROC_R3_2039_pred)
writeRaster(MIROC_R3_2039_pred, filename="output/predictions/MIROC_R3_2039/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)

########### MIROC_R3_6079 ###############

drx <- "input/MIROC_R3_6079"
files.baseline <- list.files(path = drx, pattern = "*.tif$", full.names = TRUE)

files.baseline_18 <- files.baseline[c(2:7,21,19,20,9,11,14,15,17,18,22,25,26)] # select only the 18 variables used in these models

# Covariates <- c("ce_radseas", "ct_tempann", "ct_tempiso", "ct_tempmtcp", "cw_precipdp", "cw_precipwp", "sp_awc000_100", "so_ph000_100", "so_soc000_100",
#                 "dl_strmdstge2", "dl_strmdstge6", "lf_exp315", "lf_logre10", "lf_tpi0360", "lf_tpi2000",
#                 "sp_cly000_100prop", "sp_slt000_100prop", "sp_snd000_100prop")

MIROC_R3_6079 <- stack()

for(i in 1:length(files.baseline_18)) {
  MIROC_R3_6079 <- stack(files.baseline_18[i], MIROC_R3_6079)
}

MIROC_R3_6079_pred <- lapply(models_simp, function(x){predict(MIROC_R3_6079,x,n.trees=x$gbm.call$best.trees, type="response")})
save(MIROC_R3_6079_pred, file = "output/predictions/MIROC_R3_6079/MIROC_R3_6079_pred.RData")
MIROC_R3_6079_pred <- stack(MIROC_R3_6079_pred)
writeRaster(MIROC_R3_6079_pred, filename="output/predictions/MIROC_R3_6079/ECHAM_R3_6079.tif", bylayer=TRUE, suffix = species_names)
