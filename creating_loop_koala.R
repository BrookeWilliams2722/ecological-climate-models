library(tidyverse)
library(gbm)
library(dismo)

#trying to create a loop

model.data <- read.csv("Z:/Documents/Jaramar/Otros proyectos academicos_ultima_version/RA_Koala/info_Elith_supp_mat/jane_1390_sm_appendixs3/data/model.data.csv")

model.data$Method <- as.factor(model.data$Method) # this is my incorporation 

system.time(
angaus.tc5.lr1 <- gbm.step(data=model.data, 
                            gbm.x = 3:14,
                            gbm.y = 2,
                            family = "bernoulli",
                            tree.complexity = 5,
                            learning.rate = 0.1,
                            bag.fraction = 0.5))

# Loop using pmap 
args <- list(tree.complexity = c(5,4,3), learning.rate = c(0.1,0.05,0.01))
gbm.step.pre <- partial(gbm.step,
                        data = model.data,
                        gbm.x = 3:14,
                        gbm.y = 2,
                        family = "bernoulli",
                        bag.fraction = 0.5)

results <- pmap(args, gbm.step.pre)

# Retrieve results from each model & spp 
results[[3]][["cv.statistics"]] [["gbm.call"]]

# Create a list of lists of the data you have (that you didnt create through pmap), per spp 
angus_list <- list(angaus.tc3.lr01, angaus.tc4.lr05)

# create vector with names
test_names <- c("tc3.lr01", "tc4.lr01")

# Change the names of the 'result list'
names(angus_list) <- test_names

# Save the list as an R object 

save(angus_list, file = "Z:\\Documents\\Jaramar\\Otros proyectos academicos_ultima_version\\RA_Koala\\ecological_climate_models_local\\ecological-climate-models\\output\\gbm\\angus_list.RData")

save.image(file = "Z:\\Documents\\Jaramar\\Otros proyectos academicos_ultima_version\\RA_Koala\\ecological_climate_models_local\\ecological-climate-models\\output\\gbm\\my_work_space.RData") # save the whole workspace 


# Create a table with the different no of trees, mean deviance and se for each spp

# call for the elements in a list inside a list
angus_list[[2]][[34]][[1]] 

# another way of calling the elements of a list inside a list
angus_list[[1]][[34]][c("deviance.mean", "deviance.se")]

# using lapply to apply the funciton to the entire list of lists and turn it into a data.frame
statistics <- as.data.frame(sapply(angus_list, function (x) x[["cv.statistics"]][c("deviance.mean", "deviance.se")]))

statistics.t <- as.data.frame(t(statistics)) %>% 
  mutate(test = test_names)


corygumm[[2]][["gbm.call"]]



# ------------------------------------------------


# the code I've been using for the other runs

trees.tc5.lr01.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 5,
                                    learning.rate = 0.01,
                                    bag.fraction = 0.75)

trees.tc5.lr05.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 5,
                                    learning.rate = 0.05,
                                    bag.fraction = 0.75)

trees.tc5.lr1.corygumm <- gbm.step(data=data_trees, 
                                   gbm.x = 16:33,
                                   gbm.y = 1,
                                   family = "bernoulli",
                                   tree.complexity = 5,
                                   learning.rate = 0.1,
                                   bag.fraction = 0.75)

# testing tc for lr 0.1
trees.tc4.lr1.corygumm <- gbm.step(data=data_trees, 
                                   gbm.x = 16:33,
                                   gbm.y = 1,
                                   family = "bernoulli",
                                   tree.complexity = 4,
                                   learning.rate = 0.1,
                                   bag.fraction = 0.75)

trees.tc3.lr1.corygumm <- gbm.step(data=data_trees, 
                                   gbm.x = 16:33,
                                   gbm.y = 1,
                                   family = "bernoulli",
                                   tree.complexity = 3,
                                   learning.rate = 0.1,
                                   bag.fraction = 0.75)

trees.tc2.lr1.corygumm <- gbm.step(data=data_trees, 
                                   gbm.x = 16:33,
                                   gbm.y = 1,
                                   family = "bernoulli",
                                   tree.complexity = 2,
                                   learning.rate = 0.1,
                                   bag.fraction = 0.75)

trees.tc1.lr1.corygumm <- gbm.step(data=data_trees, 
                                   gbm.x = 16:33,
                                   gbm.y = 1,
                                   family = "bernoulli",
                                   tree.complexity = 1,
                                   learning.rate = 0.1,
                                   bag.fraction = 0.75)

# testing tc for lr 0.05


trees.tc4.lr05.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 4,
                                    learning.rate = 0.05,
                                    bag.fraction = 0.75)

trees.tc3.lr05.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 3,
                                    learning.rate = 0.05,
                                    bag.fraction = 0.75)

trees.tc2.lr05.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 2,
                                    learning.rate = 0.05,
                                    bag.fraction = 0.75)

trees.tc1.lr05.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 1,
                                    learning.rate = 0.05,
                                    bag.fraction = 0.75)

# testing tc for lr 0.01

trees.tc1.lr01.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 1,
                                    learning.rate = 0.01,
                                    bag.fraction = 0.75)

trees.tc2.lr01.corygumm <- gbm.step(data=data_trees, 
                                    gbm.x = 16:33,
                                    gbm.y = 1,
                                    family = "bernoulli",
                                    tree.complexity = 2,
                                    learning.rate = 0.01,
                                    bag.fraction = 0.75)






# #Loop through all N values in vector
# for (i in 1:length(N_values)) {
#   
#   N = N_values[i]
#   
#   #Loop through all StDev values in vector for each 
#   #iteration of all N values
#   for (j in 1:length(StDv_values) {
#     
#     StDv = StDv_values[j]
#     
#     MyModel <- insert your model here... etc...
#     
#     
#     Output <- data.frame(etc...
#     )
#     
#     CompiledDF <- rbind(CompiledDF, Output)
#     
#     
#   }
# }
# 
# #Loop through all N values in vector ----
# for (i in 1:length(tree_complexity_val)) {
#   
#   tree.complexity = tree_complexity_val[i]
#   
#   #Loop through all StDev values in vector for each 
#   #iteration of all N values
#   for (j in 1:length(learning_rate_val) {
#     
#     learning.rate = learning_rate_val[j]
#     
#     MyModel <- insert your model here... etc...
#     
#     
#     Output <- data.frame(etc...
#     )
#     
#     CompiledDF <- rbind(CompiledDF, Output)
#     
#     
#   }
# }
# 
# 
# 
# # #######
# hsb2 <- read.csv("https://stats.idre.ucla.edu/stat/data/hsb2.csv")
# names(hsb2)
# 
# varlist <- names(hsb2)[8:11]
# 
# models <- lapply(varlist, function(x) {
#   lm(substitute(read ~ i, list(i = as.name(x))), data = hsb2)
# })
# 
# ## look at the first element of the list, model 1
# models[[1]]


