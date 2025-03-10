# Written by Felicity Charles
# Date:1/08/2023

##### Fire frequency analysis ----
# This script tests for correlations in the data, performs model selection, investigates spatial autocorrelation and produces training and testing data sets. 

# R version 4.3.1

# 1. Load required packages ----
library(MASS) # MASS_7.3-60
library(blockCV) # blockCV_3.1-4 
library(automap) #  automap_1.1-9
library(sf) # sf_1.0-14
library(gstat) # gstat_2.1-1
library(dismo) # dismo_1.3-14
library(terra) # terra_1.7-78 
library(ModelMetrics) # ModelMetrics_1.2.2.2
library(precrec) # precrec_0.14.4 
library(ggplot2) # ggplot2_3.5.1
library(Metrics) # Metrics_0.1.4 
library(mgcv) # mgcv_1.9-1
library(tidyterra) # tidyterra_0.6.1 
library(ggspatial) # ggspatial_1.1.9
library(caret) # caret_6.0-94
library(gbm) # gbm_2.1.9
library(doParallel) # doParallel_1.0.17


# Other attached packages not called directly
# iterators_1.0.14    
# foreach_1.5.2       
# lattice_0.21-8      
# nlme_3.1-162             
# raster_3.6-23        
# sp_2.0-0            



# 2. Read in the data and prepare for analysis steps ----
# 2.1 Point location data
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_resampled.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Rand_fire <- Rand_fire[, c(3:16)]


Background_data <- read.csv('./00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled.csv', header = T)
head(Background_data); dim(Background_data)
Background_data <- Background_data[, c(3:16)]

head(Background_data)
unique(is.na(Background_data))

# 2.3 Combine presence and background points into one dataframe
Pres_back <- rbind(Rand_fire, Background_data)
head(Pres_back); tail(Pres_back); dim(Pres_back)
unique(is.na(Pres_back))
unique(Pres_back$Sentinel_rand_firefreq)
str(Pres_back)

# 2.2 Environmental data
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')
Sentinel_ff <- round(Sentinel_ff)
environmental_preds <- rast('./00_Data/SDM_data/predictors.tif')
names(environmental_preds) <- c("QPWS_ff", "TWI", "Temp_season", "Precip_season", "Diurnal_temp", "FPC", "Soil_clay", "Slope", "Aspect", "TPI", "Elevation")


# 3. Test for correlations ----
# 3.1 Produce a correlogram ----
cor1 <- ggstatsplot::ggcorrmat(Pres_back,
                      type = "non-parametric", # Assuming that we are looking at non-parametric data here. The data is not normally distributed
                      label = T,
                      cor.vars = c("QPWS_rand_firefreq", "Sentinel_rand_firefreq", "TWI", "tempseason", "precipseason", "diurnal_temp", "FPC", "soil_clay", "slope", "aspect", "topo_position", "elevation"),
                      size = 2)
cor1
# Insignificant correlations are shown by those with a cross through the box. No correlations appear to have a Spearman rho greater than 0.8, which is our cut-off value.




# 4. Determine which variables to include in the modelling ----
# Lets use a basic linear regression model to perform stepwise variable elimination. 

# 4.1 Stepwise elimination ----
# Following other papers on the topic we want to use AIC backwards stepwise elimination

full.model <- lm(Sentinel_rand_firefreq ~ QPWS_rand_firefreq + TWI + tempseason + precipseason + diurnal_temp + FPC + soil_clay + slope + aspect + topo_position + elevation, data = Pres_back)

step.model <- stepAIC(full.model, direction = "backward")
summary(step.model)
# This suggests that some variables may be dropped lm(formula = Sentinel_rand_firefreq ~ QPWS_rand_firefreq + TWI + tempseason + precipseason + diurnal_temp + FPC + slope + aspect + topo_position + elevation, data = Pres_back)

head(Pres_back);dim(Pres_back)
Pres_back <- Pres_back[,c(1:9, 11:14)]
head(Pres_back); dim(Pres_back)
colnames(Pres_back) <- c("Sentinel_ff", "QPWS_ff", 'Lon', 'Lat', "TWI", "Temp_season", "Precip_season", "Diurnal_temp", "FPC", "Slope", "Aspect", "TPI", "Elevation")
head(Pres_back)


head(environmental_preds)
environmental_preds <- subset(environmental_preds, c(1:6, 8:11))
head(environmental_preds)


# Test for spatial autocorrelation ----
# 5.1 Use variograms to determine the extent of spatial autocorrelation ----
# Transform the data 
Pres_back_sf <- st_as_sf(Pres_back, coords = c("Lon", "Lat"), crs = 'EPSG:3577')
class(Pres_back_sf)

# Use the blockCV package to estimate extent of spatial autocorrelation
sac <- cv_spatial_autocor(x = Pres_back_sf, column = 'Sentinel_ff')
plot(sac$variograms[[1]])
# According to this if we were to fit our own empirical variogram the parameters should be nugget = 0.95, sill = 2.7, range = 10761, model = Ste


# Experimental variogram
vario1 <- variogram(Sentinel_ff ~ QPWS_ff + TWI + Temp_season + Precip_season + Diurnal_temp + FPC + Slope + Aspect + TPI + Elevation, data = Pres_back_sf)
plot(vario1)
summary(vario1)

# Fit the empirical variogram using the parameters suggested from the blockCV variogram.
vario.fit <- fit.variogram(vario1,
                           model = vgm(psill = 2.7,
                                       model = "Ste",
                                       range = 10761,
                                       nugget = 0.95))

vario.fit # Look at the result
# Parameter estimates can be adjusted further
plot(vario1, vario.fit)



# Make adjustments to the empirical variogram
vario.fit1 <- fit.variogram(vario1, 
                            model = vgm(psill = 1.245992,
                                        model = "Ste",
                                        range = 11926.02,
                                        nugget = 1.117106))
vario.fit1
plot(vario1, vario.fit1) # Change is minimal but now we know what the block size should be for spatial blocking of the data

vario.fit2 <- fit.variogram(vario1, 
                            model = vgm(psill = 1.246051,
                                        model = "Ste",
                                        range = 11924.58,
                                        nugget = 1.117028))
vario.fit2

vario.fit3 <- fit.variogram(vario1, 
                            model = vgm(psill = 1.246052,
                                        model = "Ste",
                                        range = 11924.55,
                                        nugget = 1.117027))
vario.fit3


# 3.3 Spatially block the data ----
# Random spatial blocking
# While there are other methods for spatial blocking, the large number of points makes it computationally expensive and even on the remote desktop, fails to run. So in this case we will just continue with random spatial blocking. 

sb_folds <- cv_spatial(x = Pres_back_sf,
                       column = "Sentinel_ff", # The response column
                       k = 5L, # number of folds
                       size = 11924, # size of the blocks
                       selection = "random", # random blocks-to-fold
                       seed = 503, # Set a random seed for reproducibility
                       iteration = 50L) #  find evenly dispersed folds over 50 attempts
sb_folds$records # Splitting the data based on the fire frequency value, some folds will have no points in a particular frequency due to their 'rarity' across the landscape. The train_x relates to the fire frequency with train_0 = fire freq of 0. As fire frequency increases, their abundance in the landscape decreases.  

cv_plot(cv = sb_folds,
        x = Pres_back_sf)
# We can see that the data is distributed between training and testing folds randomly across the region of interest.

# Check environmental similarity between the training and testing folds
# This gives information on whether there is possible extrapolation in the testing folds by representing how similar a point in a testing fold is to a training fold. The negative values are the sites where at least one variable has a value that is outside the range of environments over the reference set (training folds), indicating novel environments.


cv_similarity(cv = sb_folds,
              x = Pres_back_sf,
              r = environmental_preds)
# Look at https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgeb.13639&file=geb13639-sup-0001-AppendixS1.pdf section 9 to understand better

# These folds all have quite similar environments in the training and testing dataset. 


# 3.4 Extract the fold indices for training and testing data ----
folds <- sb_folds$folds_list
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
}




# Split these data into training, testing and validation data.
# We want 80% training, 20% testing
length(folds[[k]][[1]]) # Length of the training set
length(folds[[k]][[2]]) # Length of the testing set

head(Pres_back)
testing <- Pres_back[testSet, c(1,2, 5:ncol(Pres_back))]
training <- Pres_back[trainSet,c(1,2, 5:ncol(Pres_back))]


colnames(training) <- c("Sentinel_ff", "QPWS_ff", "TWI", "Temp_season", "Precip_season", "Diurnal_temp", "FPC", "Slope", "Aspect", "TPI", "Elevation")
head(training)


colnames(testing) <- c("Sentinel_ff", "QPWS_ff","TWI", "Temp_season", "Precip_season", "Diurnal_temp", "FPC", "Slope", "Aspect", "TPI", "Elevation")
testing$pres <- ifelse(testing$QPWS_ff == 0, 0, 1)
head(testing)
tail(testing)


save.image('./02_Workspaces/004_predictive_modelling_pre_hypertune.RData')




# 5. Tune a boosted regression tree model ----
# Use caret to tune a boosted regression tree model. Caret only has a gbm package boosted regression tree algorithm, however, this was the boosted regression tree which the dismo::gbm.step was based.
# The following steps to tune hyperparameters was run on the HPC
library(caret)
library(gbm)
library(doParallel)


## HPC Parallel processing
#cl <- makeCluster(32) # Requested 16 cores so have 32 threads
#registerDoParallel(cl)


# 5.1 Create grid of values for training the data ----
fitcontrol <- trainControl(method = "cv", number = 10, search = "grid", allowParallel = T)



# Note here that tree complexity is interaction.depth and shrinkage is the learning rate.
# Train a BRT model while tuning parameters. START: 2024-09-23 15:55:12 AEST; END 2024-09-24 21:04:05 AEST
gbmGrid <- expand.grid(n.trees = seq(from = 500, to = 10000, by = 100),
                       interaction.depth = seq(from = 1, to = 8, by = 1),
                       shrinkage = c(0.0001,
                                     0.005,
                                     0.001,
                                     0.05,
                                     0.01,
                                     0.5,
                                     0.1),
                       n.minobsinnode = c(50, 100, 200))



# 5.2 Tune hyperparemeters  ----
set.seed(243)
gbm_tune <- train(Sentinel_ff ~., 
            data = training, 
            method = "gbm", 
            trControl = fitcontrol, 
            tuneGrid = gbmGrid)
gbm_tune
stopCluster(cl)
# save.image('./Documents/004_predictive_model_hyperparameter_tuning.RData')

# RMSE was used to select the optimal model using the smallest value. The final values used for the model were n.trees = 10000, interaction.depth = 8, shrinkage = 0.1 and n.minobsinnode = 50.

#load('./02_Workspaces/004_predictive_model_hyperparameter_tuning.RData')

# Plot the resampling profile
p1 <- plot(gbm_tune, metric = 'RMSE')


# FINAL MODEL HYPERPARAMETER SETTINGS
# Without considering weighting of psuedoabsences/background points, the final model is n.trees = 10000, interaction.depth = 8, shrinkage = 0.1 and n.minobsinnode = 50.

# Lets now confirm the optimal number of trees using the dismo::gbm.step function as this is the function we shall be using for modelling as it was designed to work with this sort of presence/background data



### 6. Train a BRT model with optimised parameter settings ----

# As part of this, as we are using presence/background points we want to weight the points such that background points have a lower weighting then presences
# The method that has been accepted for presence only data would be an infinitely weighted logistic regression. The Valavi paper does not run IWLR BRT, so referring to the GAM and GLM implementations. The BRT method uses a method that they note is naive. We will compare both these methods. 

#load('./02_Workspaces/004_predictive_model_hyperparameter_tuning.RData')
load('./02_Workspaces/004_predictive_modelling_pre_hypertune.RData')


environmental_preds <- rast('./00_Data/SDM_data/predictors.tif')
names(environmental_preds) <- c("QPWS_ff", "TWI", "Temp_season", "Precip_season", "Diurnal_temp", "FPC", "Soil_clay", "Slope", "Aspect", "TPI", "Elevation")
head(environmental_preds)
environmental_preds <- subset(environmental_preds, c(1:6, 8:11))
head(environmental_preds)




# 6.1 Check data structure and add column for binary coding of presences and absences
head(Pres_back)
Pres_back <- Pres_back[, c(1:2, 5:ncol(Pres_back))]
Pres_back$pres <- ifelse(Pres_back$QPWS_ff == 0, 0, 1)
head(Pres_back);dim(Pres_back)


### 6.2 Model 1 - without case weights

for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  # Model with no weights 
  set.seed(480)
  fire_tc8lr.1<- gbm.step(Pres_back[trainSet,],
                          gbm.x = 2:10,
                          gbm.y = 1,
                          family = "poisson",
                          tree.complexity = 8,
                          learning.rate = 0.1)
}




# Write parameter file from model
param_file_m1 <- paste('./04_Results/Model_evaluation_statistics/Unweighted_model_lr_1.txt', sep = "")
write("Unweighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evluation metrics.", file = param_file_m1, sep = "")
write(paste("Predictors = ", fire_tc8lr.1$gbm.call$predictor.names, sep = ""), file = param_file_m1, append = T)
write(paste("Response = ", fire_tc8lr.1$gbm.call$response.name, sep = ""), file = param_file_m1, append = T)
write(paste("Family = ", fire_tc8lr.1$gbm.call$family, sep = ""), file = param_file_m1, append = T)
write(paste("Tree complexity = ", fire_tc8lr.1$gbm.call$tree.complexity, sep = ""), file = param_file_m1, append = T)
write(paste("Learning rate = ", fire_tc8lr.1$gbm.call$learning.rate, sep = ""), file = param_file_m1, append = T)
write(paste("CV folds = ", fire_tc8lr.1$gbm.call$cv.folds, sep = ""), file = param_file_m1, append = T)
write(paste("Best number of trees = ", fire_tc8lr.1$gbm.call$best.trees, sep = ""), file = param_file_m1, append = T)
write(paste("Time taken = ", fire_tc8lr.1$gbm.call$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m1, append = T)
write(paste("Mean total deviance = ", round(fire_tc8lr.1$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Mean residual deviance = ", round(fire_tc8lr.1$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Estimated cv deviance = ", round(fire_tc8lr.1$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("training data correlation = ", round(fire_tc8lr.1$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("CV correlation = ", round(fire_tc8lr.1$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m1, append = T)



# Mean total deviance fire_tc8lr.1$self.statistics$mean.null
# Mean residual deviance fire_tc8lr.1$self.statistics$mean.resid
# Estimated cv deviance accessed by fire_tc8lr.1$cv.statistics$deviance.mean and its standard error fire_tc8lr.1$cv.statistics$deviance.se
# Training data correlation fire_tc8lr.1$self.statistics$correlation
# cv corrleation fire_tc8lr.1$cv.statistics$correlation.mean and standard error fire_tc8lr.1$cv.statistics$correlation.se
# Training data AUC score fire_tc8lr.1$self.statistics$discrimination
# cv AUC score fire_tc8lr.1$cv.statistics$discrimination.mean and se fire_tc8lr.1$cv.statistics$discrimination.se


save.image('./02_Workspaces/004_predictive_modelling_post_hypertune.RData')





### 6.3 Model 2 - with case weights like in Valavi paper
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 12])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 12])["0"]) # Number of absences
  
  # Model with down-weighted background points 
  set.seed(480)
  fire_tc8lr.1_wt <- gbm.step(Pres_back[trainSet,],
                              gbm.x = 2:10,
                              gbm.y = 1,
                              family = "poisson",
                              tree.complexity = 8,
                              learning.rate = 0.1,
                              site.weights = ifelse(Pres_back[trainSet, 12] == 1, 1, prNum/bgNum))
  summary(fire_tc8lr.1_wt)
  

}


# Save model metrics
param_file_m2 <- paste('./04_Results/Model_evaluation_statistics/Down-weighted_model_lr_1.txt', sep = "")
write("Down-weighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evluation metrics.", file = param_file_m2, sep = "")
write(paste("Predictors = ", fire_tc8lr.1_wt$gbm.call$predictor.names, sep = ""), file = param_file_m2, append = T)
write(paste("Response = ", fire_tc8lr.1_wt$gbm.call$response.name, sep = ""), file = param_file_m2, append = T)
write(paste("Family = ", fire_tc8lr.1_wt$gbm.call$family, sep = ""), file = param_file_m2, append = T)
write(paste("Tree complexity = ", fire_tc8lr.1_wt$gbm.call$tree.complexity, sep = ""), file = param_file_m2, append = T)
write(paste("Learning rate = ", fire_tc8lr.1_wt$gbm.call$learning.rate, sep = ""), file = param_file_m2, append = T)
write(paste("CV folds = ", fire_tc8lr.1_wt$gbm.call$cv.folds, sep = ""), file = param_file_m2, append = T)
write(paste("Best number of trees = ", fire_tc8lr.1_wt$gbm.call$best.trees, sep = ""), file = param_file_m2, append = T)
write(paste("Time taken = ", fire_tc8lr.1$gbm.call_wt$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m2, append = T)
write(paste("Mean total deviance = ", round(fire_tc8lr.1_wt$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Mean residual deviance = ", round(fire_tc8lr.1_wt$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Estimated cv deviance = ", round(fire_tc8lr.1_wt$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1_wt$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("training data correlation = ", round(fire_tc8lr.1_wt$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("CV correlation = ", round(fire_tc8lr.1_wt$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1_wt$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m2, append = T)

# 6.4 Model 3 - with IWLR case weights - up weighting

for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  

  # Model with IWLR weights 
  set.seed(480)
  fire_tc8lr.1_IWLR <- gbm.step(Pres_back[trainSet,],
                              gbm.x = 2:10,
                              gbm.y = 1,
                              family = "poisson",
                              tree.complexity = 8,
                              learning.rate = 0.1,
                              site.weights = (10^6)^(1-Pres_back[trainSet, 12]))
  summary(fire_tc8lr.1_IWLR)
  
}


param_file_m3 <- paste('./04_Results/Model_evaluation_statistics/IWLR_model_lr_1.txt', sep = "")
write("ILWR weighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evluation metrics.", file = param_file_m3, sep = "")
write(paste("Predictors = ", fire_tc8lr.1_IWLR$gbm.call$predictor.names, sep = ""), file = param_file_m3, append = T)
write(paste("Response = ", fire_tc8lr.1_IWLR$gbm.call$response.name, sep = ""), file = param_file_m3, append = T)
write(paste("Family = ", fire_tc8lr.1_IWLR$gbm.call$family, sep = ""), file = param_file_m3, append = T)
write(paste("Tree complexity = ", fire_tc8lr.1_IWLR$gbm.call$tree.complexity, sep = ""), file = param_file_m3, append = T)
write(paste("Learning rate = ", fire_tc8lr.1_IWLR$gbm.call$learning.rate, sep = ""), file = param_file_m3, append = T)
write(paste("CV folds = ", fire_tc8lr.1_IWLR$gbm.call$cv.folds, sep = ""), file = param_file_m3, append = T)
write(paste("Best number of trees = ", fire_tc8lr.1_IWLR$gbm.call$best.trees, sep = ""), file = param_file_m3, append = T)
write(paste("Time taken = ", fire_tc8lr.1_IWLR$gbm.call$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m3, append = T)
write(paste("Mean total deviance = ", round(fire_tc8lr.1_IWLR$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Mean residual deviance = ", round(fire_tc8lr.1_IWLR$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Estimated cv deviance = ", round(fire_tc8lr.1_IWLR$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1_IWLR$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("training data correlation = ", round(fire_tc8lr.1_IWLR$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("CV correlation = ", round(fire_tc8lr.1_IWLR$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(fire_tc8lr.1_IWLR$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m3, append = T)

save.image('./02_Workspaces/004_predictive_modelling_post_hypertune.RData')



# 7. Model evaluations ----

## 7.1 Model fitted values ----
# Value above the plot indicates the weighted mean of the fitted values in relation to each non-factor predictor
dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m1 <- gbm.plot.fits(fire_tc8lr.1) 


dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m2 <- gbm.plot.fits(fire_tc8lr.1_wt)


dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m3 <- gbm.plot.fits(fire_tc8lr.1_IWLR)

# All models have similar fitted values


# 7.2 Partial dependence plots ----
# Plot the partial dependence of the response on the predictors
dev.new(width = 15, height = 7, noRStudioGD = T)
p_fun_m1 <- gbm.plot(fire_tc8lr.1, plot.layout = c(2,5))

dev.new(width = 15, height = 7, noRStudioGD = T)
p_fun_m2 <- gbm.plot(fire_tc8lr.1_wt, plot.layout = c(2,5))

dev.new(width = 15, height = 7, noRStudioGD = T)
pr_fun_m3 <- gbm.plot(fire_tc8lr.1_IWLR, plot.layout = c(2,5))

# The IWLR model had very minimal influence of aspect and TWI, with no influence of QPWS fire frequency. Temperature seasonality was the most influential predictor variable in all models, followed by diurnal temperature. All models ranked TWI and aspect among the lowest influential variables.



# 8. Produce spatial maps of predictions ----

# Make sure training data column names match the predictors
head(training)
head(environmental_preds)


# 8.2 Predict with terra ----

# 8.2.1 Model 1 - no weighting

### NOTE: Production of predictions for the models does not depend on the weights used during training. During the training phase the model has learned a function to give more importance to certain data points, therefore, we no longer need the weights for the prediction phase.
# Further to this, predictions are not returned as integer, this is because while Poisson distributions are used to model count data, it is a continuous distribution so the mean of the distribution (lambda) can be non-integer. Suggestions to handle this is rounding the data such as using round(), floor(), or ceiling(). 
# https://stackoverflow.com/questions/62912582/why-are-the-predictions-from-poisson-lasso-regression-model-in-glmnet-not-intege
  # Looking at the documentation for round(), which will round down for values between .1-.50 and round up for values between .51-.99

unweighted_pred <- terra::predict(object = environmental_preds,
                                   model = fire_tc8lr.1,
                                   type = "response",
                                   n.trees = fire_tc8lr.1$gbm.call$best.trees,
                                filename = './04_Results/Prediction_rasters/Unweighted_pred.tif', overwrite = T)

plot(unweighted_pred) # Check how this looks
unweighted_pred
# The value seems too high
round(unique(unweighted_pred$lyr1))

hist(unweighted_pred[unweighted_pred<10])
hist(unweighted_pred[unweighted_pred >1 & unweighted_pred<10])
table(round(unweighted_pred[unweighted_pred >1 & unweighted_pred <10]))



# 8.2.2 Model 2 - background points down-weighted
down_wt_pred <- terra::predict(object = environmental_preds,
                             model = fire_tc8lr.1_wt,
                             type = "response",
                             n.trees = fire_tc8lr.1_wt$gbm.call$best.trees,
                             filename = './04_Results/Prediction_rasters/Downweighted_pred.tif', overwrite = T)


plot(down_wt_pred)
down_wt_pred


hist(down_wt_pred[down_wt_pred >1 & down_wt_pred <10])
table(round(down_wt_pred[down_wt_pred >1 & down_wt_pred <10]))


# 8.2.3 Model 3 - IWLR
# Before running the prediction step we need to provide a weights argument and predict on the response scale to get the expected outcome, otherwise predictions returned do not follow the expected scale.
IWLR_pred <- terra::predict(object = environmental_preds,
                            model = fire_tc8lr.1_IWLR,
                            type = 'response',
                            n.trees = fire_tc8lr.1_IWLR$gbm.call$best.trees,
                            filename = './04_Results/Prediction_rasters/IWLR_pred.tif', overwrite = T)

plot(IWLR_pred)
IWLR_pred

QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
plot(QPWS_ff)

# Much lower predictions of fire frequency which does seem more similar to that of QPWS dat. However, this model only tended to build 50 trees, so it may not be performing comparatively well to the other BRT models which tend to overpredict.



save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')


# 9. Spatial predictive model performance evaluation ----
# 9.1 AUC and Precision-recall plots ----

# Note that ROC and PR curves work on the basis of probabilities (0,1)
# Need to extract predictions for the same coordinates as the testing data
Pres_back_crds <- rbind(Rand_fire, Background_data)
test_dat <- Pres_back_crds[testSet, c(1:9, 11:14)]
head(test_dat)
test_dat_crds <- test_dat[, 3:4]


preds_unweighted <- extract(unweighted_pred, test_dat_crds)
preds_unweighted <- preds_unweighted[,2]
preds_downwt <- extract(down_wt_pred, test_dat_crds)
preds_downwt <- preds_downwt[,2]
preds_IWLR <- extract(IWLR_pred, test_dat_crds)
preds_IWLR <- preds_IWLR[, 2]


# 9.1.1 Unweighted model 
sm1 <- mmdata(preds_unweighted, labels = ifelse(testing[,1] !=0, 1, 0))
sm1_eval <- evalmod(sm1, mode = 'rocprc')
sm1_eval
sm1_eval_basic <- evalmod(sm1, mode = 'basic')
sm1_eval_basic
# This returns numerous model evaluation measures which we may be interested in such as the error rate, accuracy, and precision. F1 score represents the harmonic mean of the precision and recall, provides a balanced measure of model performance. Matthews correlation coefficient is a reliable statistical rate which produces a high score only if the prediction obtained good results in all four confusion matrix categories (true pos, false neg, true neg, false pos)

# We can add these to the model evaluation file. We cannot extract the summary information provided by evalmod so have to copy and paste the values manually
# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m1, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm1_eval, "aucs")[1,4], digits = 3)), file = param_file_m1, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm1_eval, "aucs")[2,4], digits = 3)), file = param_file_m1, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m1, append = T)
write(paste("Classification error rate = ", round(attr(sm1_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Accuracy = ", round(attr(sm1_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Precision = ", round(attr(sm1_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm1_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm1_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m1, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall  = ", round(attr(sm1_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm1_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m1, append = T)


# 9.1.2 Down-weighted model
sm2 <- mmdata(preds_downwt, labels = ifelse(testing[,1] != 0, 1, 0))
sm2_eval <- evalmod(sm2)
sm2_eval
sm2_eval_basic <- evalmod(sm2, mode = 'basic')
sm2_eval_basic

# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m2, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm2_eval, "aucs")[1,4], digits = 3)), file = param_file_m2, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm2_eval, "aucs")[2,4], digits = 3)), file = param_file_m2, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m2, append = T)
write(paste("Classification error rate =", round(attr(sm2_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Accuracy = ", round(attr(sm2_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Precision = ", round(attr(sm2_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm2_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm2_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m2, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(sm2_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm2_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m2, append = T)



# 9.1.3 IWLR-weighted model
sm3 <- mmdata(preds_IWLR, labels = ifelse(testing[,1] != 0, 1, 0))
sm3_eval <- evalmod(sm3)
sm3_eval
sm3_eval_basic <- evalmod(sm3, mode = "basic")
sm3_eval_basic
# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m3, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm3_eval, "aucs")[1,4], digits = 3)), file = param_file_m3, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm3_eval, "aucs")[2,4], digits = 3)), file = param_file_m3, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m3, append = T)
write(paste("Classification error rate = ", round(attr(sm2_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Accuracy = ", round(attr(sm2_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Precision = ", round(attr(sm2_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm2_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm2_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m3, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(sm2_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm2_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m3, append = T)


# Based on these statistics the down-weighted model is likely performing better than the other two models.




# 9.2 Model prediction deviance from testing data and Pearsons correlation coefficient ----
# 9.2.1 For the unweighted model
# Calculate some further model metrics, note we calculate deviance from QPWS data as if it were calculated for Sentinel it returns an infinite number, we are also more interetsed in how well our model predicts for QPWS than Sentinel.
preds_unweighted[is.na(preds_unweighted)] <- 0

m1_mse <- mse(Pres_back[testSet, 1], preds_unweighted)
m1_r2 <- (1-m1_mse)/var(Pres_back[trainSet, 1])
m1_dev <- calc.deviance(testing$QPWS_ff, preds_unweighted, family = "poisson", calc.mean = T)

write(paste("Mean squared error = ", round(m1_mse, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("R-squared = ", round(m1_r2, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Deviance of observed and predicted values for QPWS data = ", round(m1_dev,digit = 3), sep = ""), file = param_file_m1, append = T)

# 9.2.2 For the down-weighted model
preds_downwt[is.na(preds_downwt)] <- 0

m2_mse <- mse(Pres_back[testSet, 1], preds_downwt)
m2_r2 <- (1-m2_mse)/var(Pres_back[trainSet, 1])
m2_dev <- calc.deviance(testing$QPWS_ff, preds_downwt, family = "poisson", calc.mean = T)


write(paste("Mean squared error = ", round(m2_mse, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("R-squared = ", round(m2_r2, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Deviance of observed and predicted values = ", round(m2_dev,digit = 3), sep = ""), file = param_file_m2, append = T)

# 9.2.3 IWLR weighting
preds_IWLR[is.na(preds_IWLR)] <- 0
m3_mse <- mse(Pres_back[testSet, 1], preds_IWLR)
m3_r2 <- (1-m3_mse)/var(Pres_back[trainSet, 1])
m3_dev <- calc.deviance(testing$QPWS_ff, preds_IWLR, family = "poisson", calc.mean = T)

write(paste("Mean squared error = ", round(m3_mse, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("R-squared = ", round(m3_r2, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Deviance of observed and predicted values = ", round(m3_dev,digit = 3), sep = ""), file = param_file_m3, append = T)


# Pearsons correlation coefficient (r2) is highest for the down-weighted model, followed by the unweighted model. Deviance of the observed and predicted values is also lowest for the down-weighted model, followed by the unweighted model.


# 9.3 Compare all models
all_scores <- join_scores(preds_unweighted, preds_downwt, preds_IWLR)
all_preds <- mmdata(all_scores, labels = ifelse(testing[,1] != 0, 1, 0), modnames = c("No weighting", "Down-weighted", "IWLR"))
preds_curves <- evalmod(all_preds)
preds_curves
preds_curves_basic <- evalmod(all_preds, mode = "basic")
preds_curves_basic
# We can see that all models perform quite similarly. Performance of the unweighted and down-weighted BRT is moderate as they fall within the range of 0.75 to 0.8., performance of the IWLR BRT is poor. There are only very small differences between the models performance based on AUC. The down-weighted model has marginally higher values, therefore is likely the best model for a BRT implementation.


#dev.new(height = 7, width = 5, dpi = 80)
par(mfrow = c(2,1), oma = c(0,0,0,0))
plot(preds_curves) # P is the number of positives and N is the number of negatives
# #y-axis sensitivity is the true positive rate, x-axis 1- specificity is the false positive rate. The curve is defined by how well it can identify areas with no fire compared to those with fire. 
# All models perform quite similarly. Precision-recall of the IWLR model is slightly lower to begin with but there were far fewer trees in this model.

save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')




# 10. Compare BRT model to GAM and GLM ----
# While we could run all the types of models that we have run previously for GAM and GLM, we will only run a down-weighted model as this is the best model that we have for the boosted regression tree implementation

#load('./02_Workspaces/004_predictive_modelling_predictions.RData')


# 10.1.1 GAM with down-weighting
# Here we are going to use the mgcv package as it has the bam function for producing generalised additive models with large datasets. Computation time is much faster using mgcv::bam() than mgcv::gam() or gamm4::gamm4().
# Using tensor product smooths as this type of smoothing is useful for when variables are not on the same scale, default basis function are cubic regression splines which are best for large datasets, only specifying  basis functions for the climatic data as it is cyclical, cc is also a type of cubic regression spline.

for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 13])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 13])["0"]) # Number of backgrounds
  
  
  # Model with down-weighted background points 
  set.seed(480)
  fire_gam <- bam(Sentinel_ff ~ te(QPWS_ff, k = 6) + te(TWI) + te(Temp_season, bs = 'cc', k = 6) + te(Precip_season, bs = 'cc', k = 6) + te(Diurnal_temp, bs = 'cc', k = 6) + te(FPC, k = 9) + + te(Slope) + te(Aspect, k = 8) + te(TPI) + te(Elevation, k = 12),
                  data = Pres_back[trainSet,],
                  family = "poisson",
                  weights = ifelse(Pres_back[trainSet, 12] == 1, 1, prNum/bgNum))
}


# Check how the gam looks
gam.check(fire_gam)
plot(fire_gam)
plot.gam(fire_gam, residuals = T)
gam.check(fire_gam)
summary(fire_gam)


# 10.2.1 GLM with down-weighting
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 13])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 13])["0"]) # Number of backgrounds
  
  
  # Model with down-weighted background points 
  set.seed(480)
  fire_glm <- glm(Sentinel_ff ~ QPWS_ff + TWI + Temp_season + Precip_season + Diurnal_temp + FPC + Slope + Aspect + TPI + Elevation,
                  data = Pres_back[trainSet,],
                  family = "poisson",
                  weights = ifelse(Pres_back[trainSet, 12] == 1, 1, prNum/bgNum))
}
summary(fire_glm)
summary.glm(fire_glm)



# 10.3 How does GLM and GAM compare to the BRT models? ----
# 10.3.1 Produce spatial predictions for GAM and GLM
# If we were to decide that a GAM implementation is the best we need to produce these predictions on a map, this will also assist with comparisons to other models spatial predictions
gam_pred <- terra::predict(object = environmental_preds,
                           type = 'response',
                           model = fire_gam,
                           na.rm = F,
                           filename = './04_Results/Prediction_rasters/GAM_pred.tif', overwrite = T)
plot(gam_pred)



glm_pred <- terra::predict(object = environmental_preds,
                           model = fire_glm,
                           type = 'response',
                           na.rm = F,
                           filename = './04_Results/Prediction_rasters/GLM_pred.tif', overwrite = T)
plot(glm_pred)


# Extract predictions for points in same manner as BRT predictions
preds_gam <- extract(gam_pred, test_dat_crds)
preds_gam <- preds_gam[,2]
preds_glm <- extract(glm_pred, test_dat_crds)
preds_glm <- preds_glm[,2]



# Evaluate model performance 
gam_mm <- mmdata(preds_gam, labels = ifelse(testing[,1] !=0, 1, 0))
gam_eval <- evalmod(gam_mm, mode = 'rocprc')
gam_eval
gam_eval_basic <- evalmod(gam_mm, mode = 'basic')
gam_eval_basic

# Write model metrics and evaluation statistics to file
param_file_m4 <- paste('./04_Results/Model_evaluation_statistics/GAM.txt', sep = "")
write("Generalised additve model for predicting fire frequency in South east Queensland. The following information provides details on model parameters and evluation metrics.", file = param_file_m4, sep = "")
write(paste("Model = ", fire_gam$call, sep = ""), file = param_file_m4, append = T)
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m4, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(gam_eval, "aucs")[1,4], digits = 3)), file = param_file_m4, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(gam_eval, "aucs")[2,4], digits = 3)), file = param_file_m4, append = T)
write(paste("Basic performance evaluation measures averages", file = param_file_m4, append = T))
write(paste("Classification error rate = ", round(attr(gam_eval_basic, "eval_summary")[4, 7], digits = 3)), file = param_file_m4, append = T)
write(paste("Accuracy = ", round(attr(gam_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m4, append = T)
write(paste("Precision = ", round(attr(gam_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m4, append = T)
write(paste("Specificity (TNR) = ", round(attr(gam_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m4, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(gam_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m4, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(gam_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m4, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(gam_eval_basic, "eval_summary")[9, 7], digits = 3)), file = param_file_m4, append = T)


# Calculating mean-squared error initially returns NULL for gam, so compare the predictions for the GAM to a BRT model
unique(is.na(preds_gam))
unique(is.na(preds_downwt))

# Need to replace NAs with 0 for gam
preds_gam[is.na(preds_gam)] <- 0

gamm_mse <- mse(Pres_back[testSet, 1], preds_gam)
gamm_r2 <- (1-gamm_mse)/var(Pres_back[trainSet, 1])
gamm_dev <- calc.deviance(testing$QPWS_ff, preds_gam, family = "gaussian", calc.mean = T)
gamm_dev # Returns infinite unless we change the family to guassian from poisson


write(paste("Mean squared error = ", round(gamm_mse, digit = 3), sep = ""), file = param_file_m4, append = T)
write(paste("R-squared = ", round(gamm_r2, digit = 3), sep = ""), file = param_file_m4, append = T)
write(paste("Deviance of observed and predicted values = ", round(gamm_dev,digit = 3), sep = ""), file = param_file_m4, append = T)


glm_mm <- mmdata(preds_glm, labels = ifelse(testing[,1] !=0, 1, 0))
glm_eval <- evalmod(glm_mm, mode = 'rocprc')
glm_eval
glm_eval_basic <- evalmod(glm_mm, mode = 'basic')
glm_eval_basic


param_file_m5 <- paste('./04_Results/Model_evaluation_statistics/GLM.txt', sep = "")
write("Generalised linear model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evluation metrics.", file = param_file_m5, sep = "")
write(paste("Model = ", fire_glm$call, sep = ""), file = param_file_m5, append = T)
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m5, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(glm_eval, "aucs")[1,4], digits = 3)), file = param_file_m5, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(glm_eval, "aucs")[1,4], digits = 3)), file = param_file_m5, append = T)
write(paste("Basic performance evaluation measures averages", file = param_file_m5, append = T))
write(paste("Classification error rate = ", round(attr(glm_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Accuracy = ", round(attr(glm_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Precision = ", round(attr(glm_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Specificity (TNR) = ", round(attr(glm_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(glm_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m5, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(glm_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(glm_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m5, append = T)


unique(is.na(preds_glm))
preds_glm[is.na(preds_glm)] <- 0

lm_mse <- mse(Pres_back[testSet, 1], preds_glm)
lm_r2 <- (1-lm_mse)/var(Pres_back[trainSet, 1])
lm_dev <- calc.deviance(testing$QPWS_ff, preds_glm, family = 'gaussian', calc.mean = T) # Would usually specify family as poisson but this will return a value of infinite so instead specifying guassian

write(paste("Mean squared error = ", round(lm_mse, digit = 3), sep = ""), file = param_file_m5, append = T)
write(paste("R-squared = ", round(lm_r2, digit = 3), sep = ""), file = param_file_m5, append = T)
write(paste("Deviance of observed and predicted values = ", round(lm_dev,digit = 3), sep = ""), file = param_file_m5, append = T)

save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')





# Compare models
all_models <- join_scores(preds_unweighted, preds_downwt, preds_IWLR, preds_gam, preds_glm)
all_models_mm <- mmdata(all_models, labels = ifelse(testing[,1] !=0, 1, 0), modnames = c("Unweighted", "Down-weighted", "IWLR", "GAM", "GLM"))
all_models_eval <- evalmod(all_models_mm)
all_models_eval

# The down-weighted BRT model does perform marginally better than a down-weighted GAM. A GLM performs the worst.


# Lets compare some further model metrics
gam_eval_basic
sm2_eval_basic
# These statistics also show that the down-weighted BRT performs marginally better than a GAM

#dev.new(height = 7, width = 5, dpi = 80)
dev.new(height = 10, width = 20)
par(mfrow = c(2,1), oma = c(0,0,0,0))
plot(all_models_eval)

# The GAM ROC is similar to that of down-weighted BRT but its precision-recall curve is lower than down-weighted and unweighted BRT models. GLM is the worst performing model. 


save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')



# 11. Plot the models ----
unweighted_pred <- rast('./04_Results/Prediction_rasters/Unweighted_pred.tif')
down_wt_pred <- rast('./04_Results/Prediction_rasters/Downweighted_pred.tif')
IWLR_pred <- rast('./04_Results/Prediction_rasters/IWLR_pred.tif')
gam_pred <- rast('./04_Results/Prediction_rasters/GAM_pred.tif')
glm_pred <- rast('./04_Results/Prediction_rasters/GLM_pred.tif')
protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(down_wt_pred)
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif') 

#Bulimbah <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-bulimbah-nature-refuge.kmz", destfile = './00_Data/Spatial data/Bulimbah_nature_refuge.kmz', mode = "wb", cacheOK = F)
#unzip(zipfile = './00_Data/Spatial data/Bulimbah_nature_refuge.kmz', exdir = './00_Data/Spatial data/')

#Gillies <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-gillies-ridge-nature-refuge.kmz", destfile = './00_Data/Spatial data/Gillies_nature_refuge.kmz', mode = "wb", cacheOK = F)
#unzip(zipfile = './00_Data/Spatial data/Gillies_nature_refuge.kmz', exdir = './00_Data/Spatial data/')

#Bartopia <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-bartopia-nature-refuge.kmz", destfile = './00_Data/Spatial data/Bartopia.kmz', mode = "wb", cacheOK = F)
#unzip(zipfile = './00_Data/Spatial data/Bartopia.kmz', exdir = './00_Data/Spatial data/')

#unzip(zipfile = './00_Data/Spatial data/Entire_OHV_perimeter.kmz', exdir='./00_Data/Spatial data/') # Note: this needs to be renamed

# Read in nature refuge files
HV <- vect('./00_Data/Spatial data/HV.kml') %>% 
  project('EPSG:3577')
Gillies <- vect('./00_Data/Spatial data/nature-refuge-gillies-ridge-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bulimbah <- vect('./00_Data/Spatial data/nature-refuge-bulimbah-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bartopia <- vect('./00_Data/Spatial data/nature-refuge-bartopia-nature-refuge.kml') %>% 
  project('EPSG:3577')

#Prior to plotting we need to mask the areas off the coast of Australia, these areas have been given 0s but we would want these to be NA values anyways.
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp') %>% 
  project('EPSG:3577')

unweighted_pred <- mask(unweighted_pred, Aus)
down_wt_pred <- mask(down_wt_pred, Aus)
IWLR_pred <- mask(IWLR_pred, Aus)

# Masking may not be necessary for the GAM and GLM but do this anyway
gam_pred <-  mask(gam_pred, Aus)
glm_pred <- mask(glm_pred, Aus)


SEQ <- crop(Aus, down_wt_pred)
plot(SEQ)

# Mask Sentinel by QPWS 
Sent_m <- mask(Sentinel_ff, QPWS_ff, inverse = T)
plet(Sent_m)

# When plotting, rather than replacing 0s with NA we need to identify method to plot 0s as white

unweighted <- ggplot() + 
  geom_spatraster(data = unweighted_pred) +
  theme_minimal()+
  scale_fill_gradient(low = "#FFF5F0", high = "darkred", limits = c(1,30), na.value = "white")



downweighted <- ggplot() + 
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = down_wt_pred) +
  theme_minimal()+
  scale_fill_gradient(low = "#FFF5F0", high = "darkred", limits = c(1,40), breaks = c(5,10,15,20,25, 30, 35, 40), na.value = "transparent")+
  geom_spatvector(data = protected_land, fill = 'transparent') +
  geom_spatvector(data = HV, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bulimbah, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bartopia, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Gillies, fill = 'transparent', col = 'gray40')+
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks')+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(face = 'bold'))


IWLR <- ggplot() + 
  geom_spatraster(data = IWLR_pred) +
  theme_minimal()+
  scale_fill_gradient(low = "#FFF5F0", high = "darkred", limits = c(1,30), na.value = "white")


QPWS <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = environmental_preds$QPWS_ff) +
  theme_minimal()+
  scale_fill_gradient(low = "#FFF5F0", high = "darkred", limits = c(1,40), breaks = c(5,10,15,20,25, 30, 35, 40), na.value = "transparent")+
  geom_spatvector(data = protected_land, fill = 'transparent') +
  geom_spatvector(data = HV, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bulimbah, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bartopia, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Gillies, fill = 'transparent', col = 'gray40')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks')+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(face = 'bold'))




GAM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = gam_pred) +
  theme_minimal() +
  scale_fill_gradient(low = '#FFF5F0', high = 'darkred', limits = c(1,30), breaks = c(5,10,15,20,25, 30), na.value = 'transparent')+
  geom_spatvector(data = protected_land, fill = 'transparent') +
  geom_spatvector(data = HV, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bulimbah, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bartopia, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Gillies, fill = 'transparent', col = 'gray40')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks')+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(face = 'bold'))




GLM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = glm_pred) +
  theme_minimal() +
  scale_fill_gradient(low = '#FFF5F0', high = 'darkred', limits = c(1,40), breaks = c(5,10,15,20,25, 30, 35, 40), na.value = 'transparent')+
  geom_spatvector(data = protected_land, fill = 'transparent') +
  geom_spatvector(data = HV, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bulimbah, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bartopia, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Gillies, fill = 'transparent', col = 'gray40')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks')+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(face = 'bold'))



Sent <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = environmental_preds$QPWS_ff) +
  geom_spatraster(data = Sent_m) +
  theme_minimal() +
  scale_fill_gradient(low = '#FFF5F0', high = 'darkred', limits = c(1,30), breaks = c(5,10,15,20,25, 30), na.value = 'transparent')+
  geom_spatvector(data = protected_land, fill = 'transparent') +
  geom_spatvector(data = HV, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bulimbah, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Bartopia, fill = 'transparent', col = 'gray40')+
  geom_spatvector(data = Gillies, fill = 'transparent', col = 'gray40')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks')+
  theme(legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(face = 'bold'))



save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')
