# Written by Felicity Charles
# Date:1/08/2023

##### Fire frequency analysis ----
# This script tests for correlations in the data, performs model selection, investigates spatial autocorrelation and produces training and testing data sets. 

# 1. Load required packages ----
library(MASS)
library(blockCV)
library(automap)
library(sf)
library(gstat)
library(dismo)


# 2. Read in the data and prepare for analysis steps ----
# 2.1 Point location data
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_resampled.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Rand_fire <- Rand_fire[, c(2,4, 3, 5:17)]


Background_data <- read.csv('./00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled.csv', header = T)
head(Background_data); dim(Background_data)
Background_data <- Background_data[, c(2:17)]

head(Background_data)
unique(is.na(Background_data))
# 2.3 Combine presence and background points into one dataframe
Pres_back <- rbind(Rand_fire, Background_data)
head(Pres_back); tail(Pres_back); dim(Pres_back)
unique(is.na(Pres_back))

# 3. Test for correlations ----
# 3.1 Produce a correlogram ----
cor1 <- ggstatsplot::ggcorrmat(Pres_back,
                      type = "non-parametric", # Assuming that we are looking at non-parametric data here. The data is not normally distributed
                      label = T,
                      cor.vars = c("QPWS_rand_firefreq", "Sentinel_rand_firefreq", "TWI", "tempseason", "precipseason", "diurnal_temp", "solar_radiation", "FPC", "soil_clay", "slope", "aspect", "topo_position", "elevation"),
                      size = 2)
cor1
# Insignificant correlations are shown by those with a cross through the box. No correlations appear to have a Spearman rho greater than 0.8, which is our cut-off value.

# Syphard et al. 2008 International Journal of Wildland Fire it was suggested to perform a bonferroni adjustment due to the large number of tests. Let's see what this returns
cor2 <- ggstatsplot::ggcorrmat(Rand_fire,
                              type = "non-parametric", # Assuming that we are looking at non-parametric data here. The data is not normally distributed
                              label = T,
                              cor.vars = c("QPWS_rand_firefreq", "Sentinel_rand_firefreq", "TWI", "tempseason", "precipseason", "diurnal_temp", "solar_radiation", "FPC", "soil_clay", "slope", "aspect", "topo_position", "elevation"),
                              p.adjust.methods = "bonferroni", 
                              size = 2)
cor2
# The output correlogram is not changed. Thus, correlations would not suggest that any variables need be removed due to high correlations. This would only be true if we lowered the threshold.


# 4. Determine which variables to include in the modelling ----
# Lets use a basic linear regression model to perform stepwise variable elimination. We will also investigate the results of stepwise elimination if it was fit on a boosted regression tree model using code produced by Jane Elith and John Leathwick for Elith, Leathwick & Hastie (2008) Journal of Animal Ecology. 
# 4.1 Stepwise elimination ----
# Following other papers on the topic we want to use AIC backwards stepwise elimination

full.model <- lm(QPWS_rand_firefreq ~ Sentinel_rand_firefreq + TWI + tempseason + precipseason + diurnal_temp + solar_radiation + FPC + soil_clay + slope + aspect + topo_position + elevation, data = Pres_back)

step.model <- stepAIC(full.model, direction = "backward")
summary(step.model)
# This suggests that some variabels may be dropped lm(formula = QPWS_rand_firefreq ~ Sentinel_rand_firefreq + tempseason + precipseason + diurnal_temp + solar_radiation + FPC + soil_clay +  slope + topo_position, data = Pres_back)


# Remove the dropped predictors from the dataset
head(Pres_back)
Pres_back <- Pres_back[,c(1:5,7:15)]
head(Pres_back); dim(Pres_back)


# Test for spatial autocorrelation ----
# 5.1 Use variograms to determine the extent of spatial autocorrelation ----
# Transform the data 
Pres_back_sf <- st_as_sf(Pres_back, coords = c("Lon", "Lat"), crs = 'EPSG:3577')
class(Pres_back_sf)

# Use the blockCV package to estimate extent of spatial autocorrelation
sac <- cv_spatial_autocor(x = Pres_back_sf, column = 'QPWS_rand_firefreq')
plot(sac$variograms[[1]])
# According to this if we were to fit our own empirical variogram the parameters should be nugget = 0.88, sill = 0.88, range = 35843, model = Ste


# Experimental variogram
vario1 <- variogram(QPWS_rand_firefreq ~ Sentinel_rand_firefreq + tempseason + precipseason + diurnal_temp + solar_radiation + FPC + soil_clay + slope + topo_position, data = Pres_back_sf)
plot(vario1)
summary(vario1)

# Fit the empirical variogram using the parameters suggested from the blockCV variogram.
vario.fit <- fit.variogram(vario1,
                           model = vgm(psill = 0.88,
                                       model = "Ste",
                                       range = 35843,
                                       nugget = 0.88))

vario.fit # Look at the result
# Parameter estimates can be adjusted further
plot(vario1, vario.fit)



# Make adjustments to the empirical variogram
vario.fit1 <- fit.variogram(vario1, 
                            model = vgm(psill = 0.00,
                                        model = "Ste",
                                        range = 35843,
                                        nugget = 0.954986))
vario.fit1
plot(vario1, vario.fit1) # Change is minimal but now we know what the block size should be for spatial blocking of the data



# 3.3 Spatially block the data ----
# Random spatial blocking
# While there are other methods for spatial blocking, the large number of points makes it computationally expensive and even on the remote desktop, fails to run. So in this case we will just continue with random spatial blocking. 

sb_folds <- cv_spatial(x = Pres_back_sf,
                       column = "QPWS_rand_firefreq", # The response column
                       k = 5L, # number of folds
                       size = 35843, # size of the blocks
                       selection = "random", # random blocks-to-fold
                       seed = 503, # Set a random seed for reproducibility
                       iteration = 50L) #  find evenly dispersed folds over 50 attempts
sb_folds$records # Splitting the data based on the fire frequency value, some folds will have no points in a particular frequency due to their 'rarity' across the landscape. The train_x relates to the fire frequency with train_0 = fire freq of 0. As fire frequency increases, their abundance in the landscape decreases.  

cv_plot(cv = sb_folds,
        x = Pres_back_sf)
# We can see that the data is distributed between training and testing folds randomly across the region of interest.

# Check environmental similarity between the training and testing folds
# This gives information on whether there is possible extrapolation in the testing folds by representing how similar a point in a testing fold is to a training fold. The negative values are the sites where at least one variable has a value that is outside the range of environments over the reference set (training folds), indicating novel environments.
predictors <- rast(c(tempseason, precipseason, diurnal_temp, solar_radiation, FPC, soil_clay, slope, aspect, topo_position))

cv_similarity(cv = sb_folds,
              x = Pres_back_sf,
              r = predictors)
# Look at https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgeb.13639&file=geb13639-sup-0001-AppendixS1.pdf section 9 to understand better

# These folds all have quite similar environments in the training and testing dataset. 


# 3.4 Extract the fold indices for training and testing data ----
folds <- sb_folds$folds_list
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
}

# Split these data into training, testing and validation data.
# We want 70% training, 15% validation, and 15% testing
length(folds[[k]][[1]]) # Length of the training set
length(folds[[k]][[2]]) # Length of the testing set

# This is an approximate 80 : 20 spit 
# We want to take 5% of the testing set and 10% of the training set to create the new validation set. This way we will have 70% training data, 15% testing data, and 15% validation data.

test <- Rand_fire[testSet,]
train <- Rand_fire[trainSet,]





# 5. Tune a boosted regression tree model ----
# Begin by training a BRT model with default parameters to determine if default parameter settings work. 
  # NOTE: parameters that are not set to the defaults are tree.complexity (because we have enough data points for a more complex tree) 

# For tuning the hyperparameters it is suggested to create 10 training and testing subsets so lets run the spatial blocking agian for this step

sb_tuning_folds <- cv_spatial(x = Pres_back_sf,
                              column = "QPWS_rand_firefreq", # The response column
                              k = 10L, # number of folds
                              size = 35843, # size of the blocks
                              selection = "random", # random blocks-to-fold
                              seed = 503, # Set a random seed for reproducibility
                              iteration = 50L,
                              plot = F)


tune_folds <- sb_tuning_folds$folds_list
for(k in seq_len(length(tune_folds))){
  ttrainSet <- unlist(tune_folds[[k]][1]) # Training set indices are the first element
}

train_tune <- Pres_back[ttrainSet,]
#test_tune <- Pres_back[ttestSet,]
head(train_tune)

# Clean up the environment
rm(Rand_fire, Background_data, Pres_back, Pres_back_sf, sb_tuning_folds, tune_folds, k, ttrainSet)
gc()
save.image('./02_Workspaces/004_predictive_modelling_prior_hypertuning_test.R')
#load('./02_Workspaces/004_predictive_modelling_prior_hypertuning_test.R')

# Test learning rate adjustments in increments of 0.5 from 0.1 to 0.0001 to find the optimum ----

mod.tc5.lr1 <- gbm.step(data = train_tune, # Use the rows of the data that have already been allocated for training the model.
                         gbm.x = c(3, 6:13), # The predictor variables
                         gbm.y = 2, # The response variable
                         family = "poisson",
                         tree.complexity = 5,
                         learning.rate = 0.1,
                         bag.fraction = 0.75)
# fitting final gbm model with a fixed number of 4100 trees for QPWS_rand_firefreq

#mean total deviance = 0.765 
#mean residual deviance = 0.064 

#estimated cv deviance = 0.308 ; se = 0.004 

#training data correlation = 0.967 
#cv correlation =  0.774 ; se = 0.004 

#elapsed time -  0.06 minutes 


mod.tc5.lr05 <- gbm.step(data = train_tune, # Use the rows of the data that have already been allocated for training the model.
                         gbm.x = c(3, 6:13), # The predictor variables
                         gbm.y = 2, # The response variable
                         family = "poisson",
                         tree.complexity = 5,
                         learning.rate = 0.05,
                         bag.fraction = 0.75)
# fitting final gbm model with a fixed number of 7150 trees for QPWS_rand_firefreq

#mean total deviance = 0.765 
#mean residual deviance = 0.075 

#estimated cv deviance = 0.309 ; se = 0.007 

#training data correlation = 0.962 
#cv correlation =  0.774 ; se = 0.007 

#elapsed time -  0.11 minutes 


mod.tc5.lr01 <- gbm.step(data = train_tune, # Use the rows of the data that have already been allocated for training the model.
                         gbm.x = c(3, 6:13), # The predictor variables
                         gbm.y = 2, # The response variable
                         family = "poisson",
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.75)
# maximum tree limit reached - results may not be optimal - refit with faster learning rate or increase maximum number of trees 

# mean total deviance = 0.765 
#mean residual deviance = 0.204 

#estimated cv deviance = 0.335 ; se = 0.007 

#training data correlation = 0.878 
#cv correlation =  0.756 ; se = 0.008 

#elapsed time -  0.33 minutes 

mod.tc5.lr005 <- gbm.step(data = train_tune,
                          gbm.x = c(3, 6:13), # The predictor variables
                          gbm.y = 2,
                          family = 'poisson',
                          tree.complexity = 5,
                          learning.rate = 0.005,
                          bag.fraction = 0.75)
# maximum tree limit reached - results may not be optimal. Higher learning rate may be better

#mean total deviance = 0.765 
#mean residual deviance = 0.284 

#estimated cv deviance = 0.371 ; se = 0.005 

#training data correlation = 0.816 
#cv correlation =  0.728 ; se = 0.003 

#elapsed time -  0.32 minutes 


mod.tc5.lr0001 <- gbm.step(data = train_tune,
                          gbm.x = c(3, 6:13), # The predictor variables
                          gbm.y = 2,
                          family = 'poisson',
                          tree.complexity = 5,
                          learning.rate = 0.0001,
                          bag.fraction = 0.75)
# maximum tree limit reached - results may not be optimal - refit with faster learning rate or increase maximum number of trees 

#mean total deviance = 0.765 
#mean residual deviance = 0.652 

#estimated cv deviance = 0.658 ; se = 0.007 

#training data correlation = 0.465 
#cv correlation =  0.448 ; se = 0.01 

#elapsed time -  0.49 minutes 
#save.image('./02_Workspaces/004_predictive_modelling.R')


# It would appear that decreasing the learning rate with the current settings increases the mean residual deviance and estimated cv deviance. A faster learning rate may be better. Will compare two of the best when considering the tree complexity also however. 



# Test tree complexity hyperparameter setting ----
# Note interactions will only be produced if the data can support this - let's try from 1 - 9 as we have 9 predictor variables with lr of 0.1 and 0.05. We can test combinations of learning rate and tree complexity to find the optimum

mod.tc1.lr1 <- gbm.step(data = train_tune, 
                        gbm.x = c(3, 6:13),
                        gbm.y = 2,
                        family = 'poisson',
                        tree.complexity = 1,
                        learning.rate = 0.1,
                        bag.fraction = 0.75)
# maximum tree limit reached - results may not be optimal - refit with faster learning rate or increase maximum number of trees 
#mean total deviance = 0.765 
#mean residual deviance = 0.403 

#estimated cv deviance = 0.484 ; se = 0.009 

#training data correlation = 0.699 
#cv correlation =  0.607 ; se = 0.009 

#elapsed time -  0.08 minutes 

mod.tc2.lr1 <- gbm.step(data = train_tune, 
                        gbm.x = c(3, 6:13),
                        gbm.y = 2,
                        family = 'poisson',
                        tree.complexity = 2,
                        learning.rate = 0.1,
                        bag.fraction = 0.75) 


## ONLY RUN TO HERE

mod.tc3.lr1 <- gbm.step(data = train_tune, 
                        gbm.x = c(3, 6:13),
                        gbm.y = 2,
                        family = 'poisson',
                        tree.complexity = 3,
                        learning.rate = 0.1,
                        bag.fraction = 0.75)

mod.tc4.lr1 <- gbm.step(data = train_tune, 
                        gbm.x = c(3, 6:13),
                        gbm.y = 2,
                        family = 'poisson',
                        tree.complexity = 4,
                        learning.rate = 0.1,
                        bag.fraction = 0.75)

mod.tc5.lr1 <- gbm.step(data = train_tune, 
                        gbm.x = c(3, 6:13),
                        gbm.y = 2,
                        family = 'poisson',
                        tree.complexity = 5,
                        learning.rate = 0.1,
                        bag.fraction = 0.75)


load('./02_Workspaces/004_predictive_modelling.R')


# Lets try a parallelisation solution to run hyperparameter tuning
# Following a solution from https://stackoverflow.com/questions/29873577/r-dismogbm-step-parameter-selection-function-in-parallel


library(doParallel)
library(foreach)
library(data.table)
# Create grid of hyperparameter and variable combinations
# Following suggested ranges in Elith, Leathwick and Hastie 2008 Journal of Animal Ecology paper.

hyper_grid <- expand.grid(learning.rate = c(0.0001,
                                            0.005,
                                            0.001,
                                            0.05,
                                            0.01,
                                            0.5,
                                            0.1),
                          tree.complexity = seq(from = 1, to = 5, by = 1),
                          bag.fraction = seq(0.45, 0.9, 0.05))

# Set up cluster to run in parallel
ncores <- detectCores()
cl <- makeCluster(ncores-1, outfile = "FullGBMGridSearchListening.txt") 
registerDoParallel(cl)

# Run grid search
print(Sys.time()) # #3:19pm 9/08/2024 start time

system.time(hyper_grid_res <- foreach (i = 1:nrow(hyper_grid), .packages = c("gbm", 'data.table'), .combine = rbind)
            
            %dopar% {
              
              TeachingDemos::char2seed("reproducibility", set = T) # set the seed for reproducibility
              
              # Train models
              gbm.tune <- dismo::gbm.step(data = train_tune,
                                   gbm.x = c(3, 6:13),
                                   gbm.y = 2,
                                   family = 'poisson',
                                   tree.complexity = hyper_grid$tree.complexity[i],
                                   learning.rate = hyper_grid$learning.rate[i],
                                   bag.fraction = hyper_grid$bag.fraction[i],
                                   plot.main = F)
              
              # Extract the info we need for model ranking
              tree.complexity <- hyper_grid$tree.complexity[i]
              learning.rate <- hyper_grid$learning.rate[i]
              bag.fraction <- hyper_grid$bag.fraction[i]
              n.trees <- gbm.tune$n.trees
              CV_correlation <- gbm.tune$cv.statistics$correlation.mean
              CV_deviance_explained <- (((gbm.tune$self.statistics$mean.null- gbm.tune$cv.statistics$deviance.mean)/gbm.tune$self.statistics$mean.null)*100)
              
              print(i) # Keep track of progress in listener
              
              # Combine desired outputs in a data.table
              data.table(tree.complexity = tree.complexity,
                         learning.rate = learning.rate,
                         bag.fraction = bag.fraction,
                         n.trees = n.trees,
                         CV_correlation = CV_correlation, 
                         CV_deviance_explained = CV_deviance_explained)
            }
            )




stopCluster(cl)

print(Sys.time())
# Order by deviance explained > number of predictors > n.trees > tree.complexity and place into a new object 

Full_BRT_parallel_grid.search.results <- hyper_grid_res %>% dplyr::arrange(desc(CV_deviance_explained), desc(n.trees), tree.complexity, learning.rate, bag.fraction)
