# Written by Felicity Charles
# Date:1/08/2023
# Updated 25/05/2026

##### Fire frequency analysis ----
# This script tests for correlations in the data, performs model selection, investigates spatial autocorrelation and produces training and testing data sets. 

# R version 4.5.1

# 1. Load required packages ----
# To install blockCV with latest version of R run the below code 
#remotes::install_github("rvalavi/blockCV", dependencies = TRUE)

library(MASS) # MASS_7.3-60 # Updated to 7.3-65
library(blockCV) # blockCV_3.1-4 # Updated to 3.2-0
library(automap) #  automap_1.1-9 # Updated to 1.1-20
library(sf) # sf_1.0-14 # Updated to 1.0-21
library(gstat) # gstat_2.1-1 # Updated to 2.1-4
library(dismo) # dismo_1.3-14 # Updated to 1.3-16
library(terra) # terra_1.7-78 # Updated to 1.8-60
library(ModelMetrics) # ModelMetrics_1.2.2.2
library(precrec) # precrec_0.14.4  # Updated to 0.14.5
library(ggplot2) # ggplot2_3.5.1 # Updated to 3.5.2
library(Metrics) # Metrics_0.1.4 
library(mgcv) # mgcv_1.9-1 # Updated to 1.9-3
library(tidyterra) # tidyterra_0.6.1 # Updated to 0.7.2
library(ggspatial) # ggspatial_1.1.9
library(caret) # caret_6.0-94 # Updated to 7.0-1
library(gbm) # gbm_2.1.9 # Updated to 2.2.2
library(doParallel) # doParallel_1.0.17
library(glmm.hp) # 0.1-7 # Updated to 0.1-8
library(gam.hp) # 0.0-3
library(ggstatsplot) # 0.13.3
# Other attached packages not called directly
# iterators_1.0.14    
# foreach_1.5.2       
# lattice_0.21-8 # Updated to 0.22-7  
# nlme_3.1-162      # Updated to 3.1-168        
# raster_3.6-23 # Updated to 3.6-32       
# sp_2.0-0 # Updated to 2.2-0
# vegan_2.6-4 # Updated to 2.7-1
# permute_0.9-7 # Updated to 0.9-8
# MuMIn_1.47.5 # Updated to 1.48.11



# 2. Read in the data and prepare for analysis steps ----
# 2.1 Point location data
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_IBRA_resampled.csv', header = T)
head(Rand_fire); dim(Rand_fire)


Background_data <- read.csv('./00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled_IBRA.csv', header = T)
head(Background_data)
Background_data$QPWS_rand_firefreq <- round(Background_data$QPWS_rand_firefreq)
head(Background_data); dim(Background_data)


# 2.3 Combine presence and background points into one dataframe
Pres_back <- rbind(Rand_fire, Background_data)
head(Pres_back); tail(Pres_back); dim(Pres_back)
colnames(Pres_back) <- c("Sentinel_ff", "QPWS_ff", "Lon", "Lat", "TWI", "Tempseason", "Precipseason", 'FPC', 'Soil_clay', 'Slope', 'Aspect', 'TPI', 'Elevation', 'BVG')
unique(is.na(Pres_back))
unique(Pres_back$Sentinel_rand_firefreq)
str(Pres_back)

# 2.2 Environmental data
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif')
Sentinel_ff <- round(Sentinel_ff)
environmental_preds <- rast('./00_Data/SDM_data/predictors_SEQ_IBRA.tif')
names(environmental_preds)
names(environmental_preds)[12] <- "BVG"
names(environmental_preds) <- c("QPWS_ff", "TWI", "Tempseason", "Precipseason", 'FPC', 'Soil_clay', 'Slope', 'Aspect', 'TPI', 'Elevation', 'BVG')



# 3. Test for correlations ----
# 3.1 Produce a correlogram ----
cor1 <- ggstatsplot::ggcorrmat(Pres_back,
                               type = "non-parametric", # Assuming that we are looking at non-parametric data here. The data is not normally distributed
                               label = T,
                               cor.vars = c("QPWS_ff", "Sentinel_ff", "TWI", "Tempseason", "Precipseason", "FPC", "Soil_clay", "Slope", "Aspect", "TPI", "Elevation", "BVG"),
                               size = 2)
cor1
# Insignificant correlations are shown by those with a cross through the box. No correlations appear to have a Spearman rho greater than 0.8, which is our cut-off value.




# 4. Determine which variables to include in the modelling ----
# Lets use a basic linear regression model to perform stepwise variable elimination. 

# 4.1 Stepwise elimination ----
# Following other papers on the topic we want to use AIC backwards stepwise elimination

full.model <- lm(Sentinel_ff ~ QPWS_ff + TWI + Tempseason + Precipseason + FPC + Soil_clay + Slope + Aspect + TPI + Elevation + BVG, data = Pres_back)

step.model <- stepAIC(full.model, direction = "backward")
summary(step.model)
# This suggests that no variables need to be dropped


# 5.1 Test for spatial autocorrelation ----
# Use variograms to determine the extent of spatial autocorrelation ----
# Transform the data 
Pres_back_sf <- st_as_sf(Pres_back, coords = c("Lon", "Lat"), crs = 'EPSG:3577')
class(Pres_back_sf)

# Use the blockCV package to estimate extent of spatial autocorrelation
sac <- cv_spatial_autocor(x = Pres_back_sf, column = 'Sentinel_ff')
plot(sac$variograms[[1]])
# According to this if we were to fit our own empirical variogram the parameters should be nugget = 0.34, sill = 2.3, range = 8216, model = Ste


# Experimental variogram
vario1 <- variogram(Sentinel_ff ~ QPWS_ff + TWI + Tempseason + Precipseason + FPC + Slope + Aspect + TPI + Elevation + BVG, data = Pres_back_sf)
plot(vario1)
summary(vario1)

# Fit the empirical variogram using the parameters suggested from the blockCV variogram.
vario.fit <- fit.variogram(vario1,
                           model = vgm(psill = 2.3,
                                       model = "Ste",
                                       range = 8216,
                                       nugget = 0.34))

vario.fit # Look at the result
# Parameter estimates can be adjusted further
plot(vario1, vario.fit)



# Make adjustments to the empirical variogram
vario.fit1 <- fit.variogram(vario1, 
                            model = vgm(psill = 0.4160198,
                                        model = "Ste",
                                        range = 29108.62,
                                        nugget = 1.2800346))
vario.fit1
plot(vario1, vario.fit1) # Change is minimal 


vario.fit2 <- fit.variogram(vario1, 
                            model = vgm(psill = 0.4160198,
                                        model = "Ste",
                                        range = 29108.61,
                                        nugget = 1.2800346))
vario.fit2
plot(vario1, vario.fit2) # No change. Now we know what the block size should be for spatial blocking of the data

# 3.3 Spatially block the data ----
# Random spatial blocking
# While there are other methods for spatial blocking, the large number of points makes it computationally expensive and even on the remote desktop, fails to run. So in this case we will just continue with random spatial blocking. 

sb_folds <- cv_spatial(x = Pres_back_sf,
                       column = "Sentinel_ff", # The response column
                       k = 5L, # number of folds
                       size = 29109, # size of the blocks
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
head(training)
testing$pres <- ifelse(testing$QPWS_ff == 0, 0, 1)
head(testing);tail(testing)
unique(testing$pres)

save.image('./02_Workspaces/004_predictive_modelling_pre_hypertune_SEQ_IBRA.RData')




# 5. Tune a boosted regression tree model ----
# NOTE: Tuning the boosted regression tree model was not re-run for the SEQ IBRA re-analysis, we are using the same model settings as were originally identified.
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

# NOTE: when running the original analysis we loaded the previously created workspaces at this stage, we won't do so for the re-analysis as we are not re-running the hyperparameter tuning.
#load('./02_Workspaces/004_predictive_model_hyperparameter_tuning.RData')
#load('./02_Workspaces/004_predictive_modelling_pre_hypertune.RData')

# Re-load environmental predictors
#environmental_preds <- rast('./00_Data/SDM_data/predictors_SEQ_IBRA.tif.tif')
#names(environmental_preds); head(environmental_preds)
#names(environmental_preds) <- c("QPWS_ff", "TWI", "Tempseason", "Precipseason", "FPC", "Soil_clay", "Slope", "Aspect", "TPI", "Elevation", "BVG")
#names(environmental_preds)

# 6.1 Check data structure and add column for binary coding of presences and absences
head(Pres_back)
Pres_back <- Pres_back[, c(1:2, 5:ncol(Pres_back))] # Remove the coordinate columns
Pres_back$pres <- ifelse(Pres_back$QPWS_ff == 0, 0, 1)
head(Pres_back);dim(Pres_back)


### 6.2 Model 1 - without case weights
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  # Model with no weights 
  set.seed(480)
  unwt_fpc <- gbm.step(Pres_back[trainSet,],
                       gbm.x = c(2:6, 8:13), # Exclude the binary presence/absence column
                       gbm.y = 1,
                       family = "poisson",
                       tree.complexity = 8,
                       learning.rate = 0.1)
}




# Write parameter file from model
param_file_m1 <- paste('./04_Results/Model_evaluation_statistics/SEQ_IBRA_Unweighted_model_FPC.txt', sep = "")
write("Unweighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evaluation metrics.", file = param_file_m1, sep = "")
write(paste("Predictors = ", unwt_fpc$gbm.call$predictor.names, sep = ""), file = param_file_m1, append = T)
write(paste("Response = ", unwt_fpc$gbm.call$response.name, sep = ""), file = param_file_m1, append = T)
write(paste("Family = ", unwt_fpc$gbm.call$family, sep = ""), file = param_file_m1, append = T)
write(paste("Tree complexity = ", unwt_fpc$gbm.call$tree.complexity, sep = ""), file = param_file_m1, append = T)
write(paste("Learning rate = ", unwt_fpc$gbm.call$learning.rate, sep = ""), file = param_file_m1, append = T)
write(paste("CV folds = ", unwt_fpc$gbm.call$cv.folds, sep = ""), file = param_file_m1, append = T)
write(paste("Best number of trees = ", unwt_fpc$gbm.call$best.trees, sep = ""), file = param_file_m1, append = T)
write(paste("Time taken = ", unwt_fpc$gbm.call$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m1, append = T)
write(paste("Mean total deviance = ", round(unwt_fpc$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Mean residual deviance = ", round(unwt_fpc$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Estimated cv deviance = ", round(unwt_fpc$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(unwt_fpc$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("training data correlation = ", round(unwt_fpc$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("CV correlation = ", round(unwt_fpc$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(unwt_fpc$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m1, append = T)




# Mean total deviance fire_tc8lr.1$self.statistics$mean.null
# Mean residual deviance fire_tc8lr.1$self.statistics$mean.resid
# Estimated cv deviance accessed by fire_tc8lr.1$cv.statistics$deviance.mean and its standard error fire_tc8lr.1$cv.statistics$deviance.se
# Training data correlation fire_tc8lr.1$self.statistics$correlation
# cv corrleation fire_tc8lr.1$cv.statistics$correlation.mean and standard error fire_tc8lr.1$cv.statistics$correlation.se
# Training data AUC score fire_tc8lr.1$self.statistics$discrimination
# cv AUC score fire_tc8lr.1$cv.statistics$discrimination.mean and se fire_tc8lr.1$cv.statistics$discrimination.se





### 6.3 Model 2 - with case weights like in Valavi paper
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 14])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 14])["0"]) # Number of absences
  
  # Model with down-weighted background points 
  set.seed(480)
  dwt_fpc <- gbm.step(Pres_back[trainSet,],
                      gbm.x = c(2:6, 8:13),
                      gbm.y = 1,
                      family = "poisson",
                      tree.complexity = 8,
                      learning.rate = 0.1,
                      site.weights = ifelse(Pres_back[trainSet, 14] == 1, 1, prNum/bgNum))
}


# Save model metrics
param_file_m2 <- paste('./04_Results/Model_evaluation_statistics/SEQ_IBRA_Down-weighted_model_FPC.txt', sep = "")
write("Down-weighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evaluation metrics.", file = param_file_m2, sep = "")
write(paste("Predictors = ", dwt_fpc$gbm.call$predictor.names, sep = ""), file = param_file_m2, append = T)
write(paste("Response = ", dwt_fpc$gbm.call$response.name, sep = ""), file = param_file_m2, append = T)
write(paste("Family = ", dwt_fpc$gbm.call$family, sep = ""), file = param_file_m2, append = T)
write(paste("Tree complexity = ", dwt_fpc$gbm.call$tree.complexity, sep = ""), file = param_file_m2, append = T)
write(paste("Learning rate = ", dwt_fpc$gbm.call$learning.rate, sep = ""), file = param_file_m2, append = T)
write(paste("CV folds = ", dwt_fpc$gbm.call$cv.folds, sep = ""), file = param_file_m2, append = T)
write(paste("Best number of trees = ", dwt_fpc$gbm.call$best.trees, sep = ""), file = param_file_m2, append = T)
write(paste("Time taken = ", dwt_fpc$gbm.call_wt$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m2, append = T)
write(paste("Mean total deviance = ", round(dwt_fpc$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Mean residual deviance = ", round(dwt_fpc$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Estimated cv deviance = ", round(dwt_fpc$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(dwt_fpc$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("training data correlation = ", round(dwt_fpc$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("CV correlation = ", round(dwt_fpc$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(dwt_fpc$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m2, append = T)





# 6.4 Model 3 - with IWLR case weights - up weighting
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  # Model with IWLR weights 
  set.seed(480)
  iwlr_fpc <- gbm.step(Pres_back[trainSet,],
                       gbm.x = c(2:6, 8:13),
                       gbm.y = 1,
                       family = "poisson",
                       tree.complexity = 8,
                       learning.rate = 0.1,
                       site.weights = (10^6)^(1-Pres_back[trainSet, 14]))
}


param_file_m3 <- paste('./04_Results/Model_evaluation_statistics/SEQ_IBRA_IWLR_model_FPC.txt', sep = "")
write("ILWR weighted model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evaluation metrics.", file = param_file_m3, sep = "")
write(paste("Predictors = ", iwlr_fpc$gbm.call$predictor.names, sep = ""), file = param_file_m3, append = T)
write(paste("Response = ", iwlr_fpc$gbm.call$response.name, sep = ""), file = param_file_m3, append = T)
write(paste("Family = ", iwlr_fpc$gbm.call$family, sep = ""), file = param_file_m3, append = T)
write(paste("Tree complexity = ", iwlr_fpc$gbm.call$tree.complexity, sep = ""), file = param_file_m3, append = T)
write(paste("Learning rate = ", iwlr_fpc$gbm.call$learning.rate, sep = ""), file = param_file_m3, append = T)
write(paste("CV folds = ", iwlr_fpc$gbm.call$cv.folds, sep = ""), file = param_file_m3, append = T)
write(paste("Best number of trees = ", iwlr_fpc$gbm.call$best.trees, sep = ""), file = param_file_m3, append = T)
write(paste("Time taken = ", iwlr_fpc$gbm.call$elapsed.time.minutes, " minutes", sep = ""), file = param_file_m3, append = T)
write(paste("Mean total deviance = ", round(iwlr_fpc$self.statistics$mean.null, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Mean residual deviance = ", round(iwlr_fpc$self.statistics$mean.resid, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Estimated cv deviance = ", round(iwlr_fpc$cv.statistics$deviance.mean, digit = 3), " and standard error = ", round(iwlr_fpc$cv.statistics$deviance.se, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("training data correlation = ", round(iwlr_fpc$self.statistics$correlation, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("CV correlation = ", round(iwlr_fpc$cv.statistics$correlation.mean, digit = 3), " and standard error = ", round(iwlr_fpc$cv.statistics$correlation.se, digit = 3), sep = ""), file = param_file_m3, append = T)


# 7. Model evaluations ----

## 7.1 Model fitted values ----
# Value above the plot indicates the weighted mean of the fitted values in relation to each non-factor predictor
dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m1_fpc <- gbm.plot.fits(unwt_fpc)

dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m2_fpc <- gbm.plot.fits(dwt_fpc)

dev.new(width = 7, height = 10, noRStudioGD = T)
p_fits_m3_fpc <- gbm.plot.fits(iwlr_fpc)

# All models have similar fitted values


# 7.2 Marginal effects plots ----
# Plot the partial dependence of the response on the predictors
dev.new(width = 15, height = 7, noRStudioGD = T)
pr_fun_m1_fpc <- gbm.plot(unwt_fpc, plot.layout = c(3,5))

dev.new(width = 15, height = 7, noRStudioGD = T)
pr_fun_m2_fpc <- gbm.plot(dwt_fpc, plot.layout = c(3,5))


dev.new(width = 15, height = 7, noRStudioGD = T)
pr_fun_m3_fpc <- gbm.plot(iwlr_fpc, plot.layout = c(3,5))


# 7.3 Relative influence data frames
unwt_rel_inf_fpc <- as.data.frame(unwt_fpc$contributions[2])
unwt_rel_inf_fpc$Variable <- as.factor(row.names(unwt_rel_inf_fpc))


dnwt_rel_inf_fpc <- as.data.frame(dwt_fpc$contributions[2])
dnwt_rel_inf_fpc$Variable <- as.factor(row.names(dnwt_rel_inf_fpc))


iwlr_rel_inf_fpc <- as.data.frame(iwlr_fpc$contributions[2])
iwlr_rel_inf_fpc$Variable <- as.factor(row.names(iwlr_rel_inf_fpc))


# The IWLR model had very minimal influence of aspect and TWI, with no influence of QPWS fire frequency. Temperature seasonality was the most influential predictor variable in all models, followed by diurnal temperature. All models ranked TWI and aspect among the lowest influential variables.



# 8. Produce spatial maps of predictions ----

# Make sure training data column names match the predictors
head(training)
environmental_preds

env_fpc <- c(environmental_preds$QPWS_ff, environmental_preds$TWI, environmental_preds$Tempseason, environmental_preds$Precipseason, environmental_preds$FPC, environmental_preds$Soil_clay, environmental_preds$Slope, environmental_preds$Aspect, environmental_preds$TPI, environmental_preds$Elevation, environmental_preds$BVG)


# 8.2 Predict with terra ----

# 8.2.1 Model 1 - no weighting

### NOTE: Production of predictions for the models does not depend on the weights used during training. During the training phase the model has learned a function to give more importance to certain data points, therefore, we no longer need the weights for the prediction phase.
# Further to this, predictions are not returned as integer, this is because while Poisson distributions are used to model count data, it is a continuous distribution so the mean of the distribution (lambda) can be non-integer. Suggestions to handle this is rounding the data such as using round(), floor(), or ceiling(). 
# https://stackoverflow.com/questions/62912582/why-are-the-predictions-from-poisson-lasso-regression-model-in-glmnet-not-intege
# Looking at the documentation for round(), which will round down for values between .1-.50 and round up for values between .51-.99

unweighted_pred <- terra::predict(object = environmental_preds,
                                  model = unwt_fpc,
                                  type = "response",
                                  n.trees = unwt_fpc$gbm.call$best.trees,
                                  filename = './04_Results/Prediction_rasters/SEQ_IBRA_Unweighted_pred_FPC.tif', overwrite = T)

plot(unweighted_pred) # Check how this looks
unweighted_pred
# The value seems too high
round(unique(unweighted_pred$lyr1))

hist(unweighted_pred[unweighted_pred<10])
hist(unweighted_pred[unweighted_pred >1 & unweighted_pred<10])
table(round(unweighted_pred[unweighted_pred >1 & unweighted_pred <10]))
hist(unweighted_pred[unweighted_pred >1])




# 8.2.2 Model 2 - background points down-weighted
downweighted_pred_fpc <- terra::predict(object = env_fpc,
                                        model = dwt_fpc,
                                        type = "response",
                                        n.trees = dwt_fpc$gbm.call$best.trees,
                                        filename = './04_Results/Prediction_rasters/SEQ_IBRA_Downweighted_FPC_pred.tif', overwrite = T)


# 8.2.3 Model 3 - IWLR
# Before running the prediction step we need to provide a weights argument and predict on the response scale to get the expected outcome, otherwise predictions returned do not follow the expected scale.
IWLR_pred_fpc <- terra::predict(object = environmental_preds,
                                model = iwlr_fpc,
                                type = 'response',
                                n.trees = iwlr_fpc$gbm.call$best.trees,
                                filename = './04_Results/Prediction_rasters/SEQ_IBRA_IWLR_FPC_pred.tif', overwrite = T)



save.image('./02_Workspaces/004_predictive_modelling_predictions_SEQ_IBRA.RData')


# 9. Spatial predictive model performance evaluation ----
# 9.1 AUC and Precision-recall plots ----

# Note that ROC and PR curves work on the basis of probabilities (0,1)
# Need to extract predictions for the same coordinates as the testing data
Pres_back_crds <- rbind(Rand_fire, Background_data)
test_dat <- Pres_back_crds[testSet, c(1:9, 11:14)]
head(test_dat)
test_dat_crds <- test_dat[, 3:4]


preds_unweighted_fpc <- extract(unweighted_pred, test_dat_crds)
preds_unweighted_fpc <- preds_unweighted_fpc[,2]

preds_downwt_fpc <- extract(downweighted_pred_fpc, test_dat_crds)
preds_downwt_fpc <- preds_downwt_fpc[,2]

preds_iwlr_fpc <- extract(IWLR_pred_fpc, test_dat_crds)
preds_iwlr_fpc <- preds_iwlr_fpc[,2]

# 9.1.1 Unweighted FPC model 
sm1_fpc <- mmdata(preds_unweighted_fpc, labels = ifelse(testing[,1] !=0, 1, 0))
sm1_fpc_eval <- evalmod(sm1_fpc, mode = 'rocprc')
sm1_fpc_eval
sm1_fpc_eval_basic <- evalmod(sm1_fpc, mode = 'basic')
sm1_fpc_eval_basic
# This returns numerous model evaluation measures which we may be interested in such as the error rate, accuracy, and precision. F1 score represents the harmonic mean of the precision and recall, provides a balanced measure of model performance. Matthews correlation coefficient is a reliable statistical rate which produces a high score only if the prediction obtained good results in all four confusion matrix categories (true pos, false neg, true neg, false pos)

# We can add these to the model evaluation file. We cannot extract the summary information provided by evalmod so have to copy and paste the values manually
# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m1, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm1_fpc_eval, "aucs")[1,4], digits = 3)), file = param_file_m1, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm1_fpc_eval, "aucs")[2,4], digits = 3)), file = param_file_m1, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m1, append = T)
write(paste("Classification error rate = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Accuracy = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Precision = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m1, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall  = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m1, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm1_fpc_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m1, append = T)




# 9.1.3 Down-weighted FPC model
sm2_fpc <- mmdata(preds_downwt_fpc, labels = ifelse(testing[,1] !=0, 1, 0))
sm2_fpc_eval <- evalmod(sm2_fpc, mode = 'rocprc')
sm2_fpc_eval
sm2_fpc_eval_basic <- evalmod(sm2_fpc, mode = 'basic')
sm2_fpc_eval_basic

# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m2, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm2_fpc_eval, "aucs")[1,4], digits = 3)), file = param_file_m2, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm2_fpc_eval, "aucs")[2,4], digits = 3)), file = param_file_m2, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m2, append = T)
write(paste("Classification error rate =", round(attr(sm2_fpc_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Accuracy = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Precision = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m2, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m2, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm2_fpc_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m2, append = T)



# 9.1.5 IWLR-weighted FPC model
sm3_fpc <- mmdata(preds_iwlr_fpc, labels = ifelse(testing[,1] !=0, 1, 0))
sm3_fpc_eval <- evalmod(sm3_fpc, mode = 'rocprc')
sm3_fpc_eval
sm3_fpc_eval_basic <- evalmod(sm3_fpc, mode = 'basic')
sm3_fpc_eval_basic
# Add measures to file
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m3, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(sm3_fpc_eval, "aucs")[1,4], digits = 3)), file = param_file_m3, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(sm3_fpc_eval, "aucs")[2,4], digits = 3)), file = param_file_m3, append = T)
write(paste("Basic performance evaluation measures averages"), file = param_file_m3, append = T)
write(paste("Classification error rate = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Accuracy = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Precision = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Specificity (TNR) = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m3, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m3, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(sm3_fpc_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m3, append = T)



# Compare performances
sm1_fpc_eval # Performs best of unweighted

sm2_fpc_eval # Performs similarly to model with NDVI and FPC

sm3_fpc_eval # Performs best of IWLR



sm1_fpc_eval; sm2_fpc_eval; sm3_fpc_eval # Downweighted performs best but only marginally more so than unweighted model. 




# 9.2 Model prediction deviance from testing data and Pearsons correlation coefficient ----
# 9.2.1 For the unweighted model
# Calculate some further model metrics, note we calculate deviance from QPWS data as if it were calculated for Sentinel it returns an infinite number, we are also more interetsed in how well our model predicts for QPWS than Sentinel.
preds_unweighted_fpc[is.na(preds_unweighted_fpc)] <- 0

m1_fpc_mse <- mse(Pres_back[testSet, 1], preds_unweighted_fpc)
m1_fpc_r2 <- (1-m1_fpc_mse)/var(Pres_back[trainSet, 1])
m1_fpc_dev <- calc.deviance(testing$QPWS_ff, preds_unweighted_fpc, family = "poisson", calc.mean = T)
write(paste("Mean squared error = ", round(m1_fpc_mse, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("R-squared = ", round(m1_fpc_r2, digit = 3), sep = ""), file = param_file_m1, append = T)
write(paste("Deviance of observed and predicted values for QPWS data = ", round(m1_fpc_dev,digit = 3), sep = ""), file = param_file_m1, append = T)



# 9.2.2 For the down-weighted model
preds_downwt_fpc[is.na(preds_downwt_fpc)] <- 0

m2_fpc_mse <- mse(Pres_back[testSet, 1], preds_downwt_fpc)
m2_fpc_r2 <- (1-m2_fpc_mse)/var(Pres_back[trainSet, 1])
m2_fpc_dev <- calc.deviance(testing$QPWS_ff, preds_downwt_fpc, family = "poisson", calc.mean = T)
write(paste("Mean squared error = ", round(m2_fpc_mse, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("R-squared = ", round(m2_fpc_r2, digit = 3), sep = ""), file = param_file_m2, append = T)
write(paste("Deviance of observed and predicted values = ", round(m2_fpc_dev,digit = 3), sep = ""), file = param_file_m2, append = T)




# 9.2.3 IWLR weighting
preds_iwlr_fpc[is.na(preds_iwlr_fpc)] <- 0


m3_fpc_mse <- mse(Pres_back[testSet, 1], preds_iwlr_fpc)
m3_fpc_r2 <- (1-m3_fpc_mse)/var(Pres_back[trainSet, 1])
m3_fpc_dev <- calc.deviance(testing$QPWS_ff, preds_iwlr_fpc, family = "poisson", calc.mean = T)
write(paste("Mean squared error = ", round(m3_fpc_mse, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("R-squared = ", round(m3_fpc_r2, digit = 3), sep = ""), file = param_file_m3, append = T)
write(paste("Deviance of observed and predicted values = ", round(m3_fpc_dev,digit = 3), sep = ""), file = param_file_m3, append = T)


# Pearsons correlation coefficient (r2) is highest for the down-weighted model, followed by the unweighted model.  Deviance of the observed and predicted values is also lowest for the down-weighted model, followed by the unweighted model for models with all enviro vars.


# 9.3 Compare all models
# P is the number of positives and N is the number of negatives
# #y-axis sensitivity is the true positive rate, x-axis 1- specificity is the false positive rate. The curve is defined by how well it can identify areas with no fire compared to those with fire. 
# All models perform quite similarly. Precision-recall of the IWLR model is slightly lower to begin with but there were far fewer trees in this model.

# FPC models
all_scores_fpc <- join_scores(preds_unweighted_fpc, preds_downwt_fpc, preds_iwlr_fpc)
all_preds_fpc <- mmdata(all_scores_fpc, labels = ifelse(testing[,1] != 0, 1, 0), modnames = c("No weighting", "Down-weighted", "IWLR"))
preds_curves_fpc <- evalmod(all_preds_fpc)
preds_curves_fpc
preds_curves_basic_fpc <- evalmod(all_preds_fpc, mode = "basic")
preds_curves_basic_fpc
# We can see that all models perform quite similarly. Performance of the unweighted and down-weighted BRT is moderate as they fall within the range of 0.75 to 0.8., performance of the IWLR BRT is poor. There are only very small differences between the models performance based on AUC. The down-weighted model has marginally higher values, therefore is likely the best model for a BRT implementation.


#dev.new(height = 7, width = 5, dpi = 80)
par(mfrow = c(2,1), oma = c(0,0,0,0))
plot(preds_curves_fpc) # P is the number of positives and N is the number of negatives



# All models
allsc <- join_scores(preds_unweighted_fpc, preds_downwt_fpc, preds_iwlr_fpc)
all_pr <- mmdata(allsc, labels = ifelse(testing[,1] != 0, 1, 0), modnames = c("Non", "Down", "IWLR"))
pr_curves <- evalmod(all_pr)
pr_curves
pr_curves_basic <- evalmod(all_pr, mode = 'basic')
pr_curves_basic
plot(pr_curves)

save.image('./02_Workspaces/004_predictive_modelling_predictions_SEQ_IBRA.RData')




# 10. Compare BRT model to GAM and GLM ----
# While we could run all the types of models that we have run previously for GAM and GLM, we will only run a down-weighted model as this is the best model that we have for the boosted regression tree implementation

#load('./02_Workspaces/004_predictive_modelling_predictions_SEQ_IBRA.RData')

# 10.1.1 GAM with down-weighting
# Here we are going to use the mgcv package as it has the bam function for producing generalised additive models with large datasets. Computation time is much faster using mgcv::bam() than mgcv::gam() or gamm4::gamm4().
# Using tensor product smooths as this type of smoothing is useful for when variables are not on the same scale, default basis function are cubic regression splines which are best for large datasets as they have modest sized sets of knots spread evenly through the covariate values. So only specifying  basis functions for the climatic data as it is cyclical, cc is also a type of cubic regression spline.

for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 14])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 14])["0"]) # Number of backgrounds
  
  
  # Model with down-weighted background points 
  set.seed(480)
  fpc_gam <- bam(Sentinel_ff ~ te(QPWS_ff) + te(TWI, k = 4) + te(Tempseason, bs = "cc", k = 6) + te(Precipseason, bs = "cc", k = 6) + te(FPC, k = 6) + te(Soil_clay) + te(Slope) + te(Aspect, k = 8) + te(TPI) + te(Elevation, k = 8) + te(BVG, bs = 're'),
                 data = Pres_back[trainSet,],
                 family = "poisson",
                 weights = ifelse(Pres_back[trainSet, 14] == 1, 1, prNum/bgNum))
} 


# Check how the gam looks
par(mfrow = c(2,2)); gam.check(fpc_gam)
par(mfrow = c(3,4)); plot(fpc_gam)
plot.gam(fpc_gam, residuals = T)
summary(fpc_gam)

fpc_gam_inf <- gam.hp(fpc_gam)

fpc_gam_rel_inf<- as.data.frame(fpc_gam_inf[,4])
colnames(fpc_gam_rel_inf) <-  "I.perc"
fpc_gam_rel_inf$Variable <- as.factor(rownames(fpc_gam_inf))
fpc_gam_rel_inf
str(fpc_gam_rel_inf)




# 10.2 GLM model
for(k in seq_len(length(folds))){
  trainSet <- unlist(folds[[k]][1]) # Training set indices are the first element
  testSet <- unlist(folds[[k]][2]) # Testing set indices are the second element
  
  
  prNum <- as.numeric(table(Pres_back[trainSet, 14])["1"]) # Number of presences
  bgNum <- as.numeric(table(Pres_back[trainSet, 14])["0"]) # Number of backgrounds
  
  
  # Model with down-weighted background points 
  set.seed(480)
  fpc_glm <- glm(Sentinel_ff ~ QPWS_ff + TWI + Tempseason + Precipseason + FPC + Soil_clay + Slope + Aspect + TPI + Elevation + BVG,
                 data = Pres_back[trainSet,],
                 family = "poisson",
                 weights = ifelse(Pres_back[trainSet, 14] == 1, 1, prNum/bgNum))
}
summary(fpc_glm)
summary.glm(fpc_glm)


# Calculate relative influence of variables for the glm
fpc_glm_inf <- glmm.hp(fpc_glm)

# We want the values for I.perc(%) as this is comparable to the gbm.name$contributions table. 
fpc_glm_rel_inf <- as.data.frame(fpc_glm_inf$delta[,4])
colnames(fpc_glm_rel_inf) <-  "I.perc"
fpc_glm_rel_inf$Variable <- as.factor(rownames(fpc_glm_rel_inf))
str(fpc_glm_rel_inf)




# 10.3 How does GLM and GAM compare to the BRT models? ----
# 10.3.1 Produce spatial predictions for GAM and GLM
# If we were to decide that a GAM implementation is the best we need to produce these predictions on a map, this will also assist with comparisons to other models spatial predictions
gam_fpc_pred <- terra::predict(object = env_fpc,
                               type = 'response',
                               model = fpc_gam,
                               na.rm = F,
                               filename = './04_Results/Prediction_rasters/SEQ_IBRA_GAM_FPC_pred.tif', overwrite = T)
plot(gam_fpc_pred) 




glm_fpc_pred <- terra::predict(object = env_fpc,
                               model = fpc_glm,
                               type = 'response',
                               na.rm = F,
                               filename = './04_Results/Prediction_rasters/SEQ_IBRA_GLM_FPC_pred.tif', overwrite = T)
plot(glm_fpc_pred)


# Extract predictions for points in same manner as BRT predictions
preds_fpc_gam <- extract(gam_fpc_pred, test_dat_crds)
preds_fpc_gam <- preds_fpc_gam[,2]

preds_fpc_glm <- extract(glm_fpc_pred, test_dat_crds)
preds_fpc_glm <- preds_fpc_glm[,2]

# Evaluate model performance 
# GAM FPC model
gam_mm <- mmdata(preds_fpc_gam, labels = ifelse(testing[,1] !=0, 1, 0))
gam_eval <- evalmod(gam_mm, mode = 'rocprc')
gam_eval
gam_eval_basic <- evalmod(gam_mm, mode = 'basic')
gam_eval_basic

# Write model metrics and evaluation statistics to file
param_file_m4 <- paste('./04_Results/Model_evaluation_statistics/SEQ_IBRA_GAM_FPC.txt', sep = "")
write("Generalised additve model for predicting fire frequency in South east Queensland. The following information provides details on model parameters and evluation metrics.", file = param_file_m4, sep = "")
write(paste("Model = ", fpc_gam$call, sep = ""), file = param_file_m4, append = T)
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
unique(is.na(preds_fpc_gam))
unique(is.na(preds_downwt_fpc))

# Need to replace NAs with 0 for gam
preds_fpc_gam[is.na(preds_fpc_gam)] <- 0

gamm_mse <- mse(Pres_back[testSet, 1], preds_fpc_gam)
gamm_r2 <- (1-gamm_mse)/var(Pres_back[trainSet, 1])
gamm_dev <- calc.deviance(testing$QPWS_ff, preds_fpc_gam, family = "gaussian", calc.mean = T)
gamm_dev # Returns infinite unless we change the family to guassian from poisson


write(paste("Mean squared error = ", round(gamm_mse, digit = 3), sep = ""), file = param_file_m4, append = T)
write(paste("R-squared = ", round(gamm_r2, digit = 3), sep = ""), file = param_file_m4, append = T)
write(paste("Deviance of observed and predicted values = ", round(gamm_dev,digit = 3), sep = ""), file = param_file_m4, append = T)



# GLM FPC model
glm_mm <- mmdata(preds_fpc_glm, labels = ifelse(testing[,1] !=0, 1, 0))
glm_eval <- evalmod(glm_mm, mode = 'rocprc')
glm_eval
glm_eval_basic <- evalmod(glm_mm, mode = 'basic')
glm_eval_basic


param_file_m5 <- paste('./04_Results/Model_evaluation_statistics/SEQ_IBRA_GLM_FPC.txt', sep = "")
write("Generalised linear model for predicting fire frequency in South east Queensland using a Boosted regression tree. The following information provides details on model parameters and evluation metrics.", file = param_file_m5, sep = "")
write(paste("Model = ", fpc_glm$call, sep = ""), file = param_file_m5, append = T)
write(paste("The following model evaluation measures were calculated using precrec::evalmod(), by including mode = 'basic' this returns further measures beyond AUC ROC and precision-recall curves."), file = param_file_m5, append = T)
write(paste("Area Under the Reciever Operating Characteristic Curve (AUC ROC)  = ", round(attr(glm_eval, "aucs")[1,4], digits = 3)), file = param_file_m5, append = T)
write(paste("Precision-Recall curve (PRC) = ", round(attr(glm_eval, "aucs")[2,4], digits = 3)), file = param_file_m5, append = T)
write(paste("Basic performance evaluation measures averages", file = param_file_m5, append = T))
write(paste("Classification error rate = ", round(attr(glm_eval_basic, "eval_summary")[4,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Accuracy = ", round(attr(glm_eval_basic, "eval_summary")[5,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Precision = ", round(attr(glm_eval_basic, "eval_summary")[8,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Specificity (TNR) = ", round(attr(glm_eval_basic, "eval_summary")[6,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Sensitivity (TPR) = ", round(attr(glm_eval_basic, "eval_summary")[7,7], digits = 3)), file = param_file_m5, append = T)
write(paste("F-score, a balanced measure of model performance based on precision and recall = ", round(attr(glm_eval_basic, "eval_summary")[10,7], digits = 3)), file = param_file_m5, append = T)
write(paste("Matthews correlation coefficient = ", round(attr(glm_eval_basic, "eval_summary")[9,7], digits = 3)), file = param_file_m5, append = T)


unique(is.na(preds_fpc_glm))
preds_fpc_glm[is.na(preds_fpc_glm)] <- 0

lm_mse <- mse(Pres_back[testSet, 1], preds_fpc_glm)
lm_r2 <- (1-lm_mse)/var(Pres_back[trainSet, 1])
lm_dev <- calc.deviance(testing$QPWS_ff, preds_fpc_glm, family = 'gaussian', calc.mean = T) # Would usually specify family as poisson but this will return a value of infinite so instead specifying guassian

write(paste("Mean squared error = ", round(lm_mse, digit = 3), sep = ""), file = param_file_m5, append = T)
write(paste("R-squared = ", round(lm_r2, digit = 3), sep = ""), file = param_file_m5, append = T)
write(paste("Deviance of observed and predicted values = ", round(lm_dev,digit = 3), sep = ""), file = param_file_m5, append = T)




save.image('./02_Workspaces/004_predictive_modelling_predictions_SEQ_IBRA.RData')



# Compare models
all_models <- join_scores(preds_unweighted_fpc, preds_downwt_fpc, preds_iwlr_fpc, preds_fpc_gam, preds_fpc_glm)
all_models_mm <- mmdata(all_models, labels = ifelse(testing[,1] !=0, 1, 0), modnames = c("Unweighted", "Down-weighted", "IWLR", "GAM", "GLM"))
all_models_eval <- evalmod(all_models_mm)
all_models_eval
# The down-weighted BRT model does perform marginally better than an unweighted BRT or down-weighted GLM and GAM. The IWLR model performs the worst.


#dev.new(height = 7, width = 5, dpi = 80)
dev.new(height = 10, width = 20)
par(mfrow = c(2,1), oma = c(0,0,0,0))
plot(all_models_eval)

# The GAM ROC is similar to that of down-weighted BRT, the GLM also performs quite well. 




save.image('./02_Workspaces/004_predictive_modelling_predictions.RData')





# 11. Compare the relative influence of variables for all models ----
# GLM, GAM, down-weighted BRT, unweighted BRT, Infinite BRT
fpc_glm_rel_inf
str(fpc_glm_rel_inf)

quartz()
dev.new(width =15, height = 8, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(1, 2), mar = c(9,4,3.5,1))

barplot(fpc_glm_rel_inf$I.perc ~ fpc_glm_rel_inf$Variable, ylim = c(0,100), ylab = "", xlab = "", las = 2, names.arg = c("Aspect", "BVG", "Elevation", "FPC", "Precipitation \n seasonality", "Public land \n fire frequency", "Slope", "% soil clay", "Temperature \n seasonality", "TPI", "TWI"), yaxt = "n", cex.names = 1.2)
mtext(expression(bold("Environmental predictor")), side = 1, line = 8, cex = 1.5)
axis(side = 2, at = seq(0,100, 10), las = 1, line = -0.5, cex.axis = 1.2)
axis(side = 2, at = seq(0,100, 5), labels = F, line = -0.5)
mtext(expression(bold("Relative contribution (%)")), side = 2, line = 1.8, cex = 1.5)
mtext(expression(bold("(a) GLM")), line = 1, cex = 2)

barplot(fpc_gam_rel_inf$I.perc ~ fpc_gam_rel_inf$Variable, ylim = c(0,100), ylab = "", xlab = "", las = 2, names.arg = c("Aspect", "BVG", "Elevation", "FPC", "Precipitation \n seasonality", "Public land \n fire frequency", "Slope", "% soil clay", "Temperature \n seasonality", "TPI", "TWI"), yaxt = "n", cex.names = 1.2)
mtext(expression(bold("Environmental predictor")), side = 1, line = 8, cex = 1.5)
axis(side = 2, at = seq(0,100, 10), las = 1, line = -0.5, cex.axis = 1.2)
axis(side = 2, at = seq(0,100, 5), labels = F, line = -0.5)
mtext(expression(bold("Relative contribution (%)")), side = 2, line = 1.8, cex = 1.5)
mtext(expression(bold("(b) GAM")), line = 1, cex = 2)
quartz.save(file = './04_Results/Plots/Variable_contributions/Main_FPC.png', type = 'png', dpi = 80)

# BRT as subplots
quartz()
dev.new(width = 16, height = 6, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(1, 3), mar = c(11,4,3.5,1))

barplot(dnwt_rel_inf_fpc$rel.inf ~ dnwt_rel_inf_fpc$Variable, ylim = c(0,100), ylab = "", xlab = "", las = 2, names.arg = c("Aspect", "BVG", "Elevation", "FPC", "Precipitation \n seasonality", "Public land \n fire frequency", "Slope", "% soil clay", "Temperature \n seasonality", "TPI", "TWI"), yaxt = "n", cex.names = 1.5)
mtext(expression(bold("Environmental predictor")), side = 1, line = 10, cex = 1.5)
axis(side = 2, at = seq(0,100, 10), las = 1, line = -0.5, cex.axis = 1.5)
axis(side = 2, at = seq(0,100, 5), labels = F, line = -0.5)
mtext(expression(bold("Relative contribution (%)")), side = 2, line = 2, cex = 1.5)
mtext(expression(bold("(c) Down-weighted BRT")), line = 1, cex = 2)



barplot(unwt_rel_inf_fpc$rel.inf ~ unwt_rel_inf_fpc$Variable, ylim = c(0,100), ylab = "", xlab = "", las = 2, names.arg = c("Aspect", "BVG", "Elevation", "FPC", "Precipitation \n seasonality", "Public land \n fire frequency", "Slope", "% soil clay", "Temperature \n seasonality", "TPI", "TWI"), yaxt = "n", cex.names = 1.5)
mtext(expression(bold("Environmental predictor")), side = 1, line = 10, cex = 1.5)
axis(side = 2, at = seq(0,100, 10), las = 1, line = -0.5, cex.axis = 1.5)
axis(side = 2, at = seq(0,100, 5), labels = F, line = -0.5)
mtext(expression(bold("Relative contribution (%)")), side = 2, line = 2, cex = 1.5)
mtext(expression(bold("(d) Unweighted BRT")), line = 1, cex = 2)


barplot(iwlr_rel_inf_fpc$rel.inf ~ iwlr_rel_inf_fpc$Variable, ylim = c(0,100), ylab = "", xlab = "", las = 2, names.arg = c("Aspect", "BVG", "Elevation", "FPC", "Precipitation \n seasonality", "Public land \n fire frequency", "Slope", "% soil clay", "Temperature \n seasonality", "TPI", "TWI"), yaxt = "n", cex.names = 1.5)
mtext(expression(bold("Environmental predictor")), side = 1, line = 10, cex = 1.5)
axis(side = 2, at = seq(0,100, 10), las = 1, line = -0.5, cex.axis = 1.5)
axis(side = 2, at = seq(0,100, 5), labels = F, line = -0.5)
mtext(expression(bold("Relative contribution (%)")), side = 2, line = 2, cex = 1.5)
mtext(expression(bold("(e) Infinite BRT")), line = 1, cex = 2)
quartz.save(file = './04_Results/Plots/Variable_contributions/BRT_FPC.png', type = 'png', dpi = 80)



save.image('./02_Workspaces/004_predictive_modelling_predictions_SEQ_IBRA.RData')


