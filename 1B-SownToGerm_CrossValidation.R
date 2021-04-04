########################
### CROSS-VALIDATION ###
########################
# CODE FROM KENEN

#### LOAD MODEL & DATA ####
# NOTE: LOAD SownSurv MODEL STRUCTURE FROM MODEL SCRIPT 1-SownToGerm...
# load model output
load('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/03_Model_Output_Data/01_Seed_To_Germination_Model/GermModel_03Apr2021.RDS')
# full dataset
in_germ_data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/02_Processed_Data_For_Model_Input/01_Seed_To_Germination_Model/GerminationMatrix_FINAL_15DEC2020.csv', header=TRUE, stringsAsFactors=FALSE)
# subset dataset (for testing)
in_germ_data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixGermData_SourceSubset_updatedChambers.csv', header=TRUE, stringsAsFactors=FALSE)
  
#### LOAD DEPENDENCIES ####
library(dplyr)
library(ggplot2)
source('C:/Users/Maggie/Documents/WetlandEcology/SBWR/00_CrossValidationFunctions.R')

#### DATA PREPROCESSING ####
## save data as appropriate list for JAGS
germData <- list(germ=in_germ_data[,7:43],
                 species=as.integer(as.factor(germination_matrix$Species)),
                 source=as.integer(as.factor(germination_matrix$Source)),
                 seedmass=germination_matrix$SeedMass,
                 #sct=germination_matrix$SCT,
                 wp=germination_matrix$WP,
                 temp=germination_matrix$Temp,
                 cup=as.integer(as.factor(germination_matrix$CupNo)),
                 chamber=as.integer(as.factor(germination_matrix$Chamber)),
                 rep=as.integer(as.factor(germination_matrix$Trial)),
                 NSown=nrow(germination_matrix),
                 NSpecies = length(unique(germination_matrix$Species)),
                 NSources = length(unique(germination_matrix$Source)),
                 NCups = length(unique(germination_matrix$CupNo)),
                 NTimes = 37, #ncol(germData$germ)
                 NChambers = length(unique(germination_matrix$Chamber)))

## Apply representative 10-fold CV subset function to data
names(in_germ_data) # find columns with categorical data
in_germ_data$TenfoldID <- protected_CV_subset(CV_data=in_germ_data, colnum=c(1, 2, 3, 4, 5, 6, 52))


#### RUN 10-FOLD CROSS-VALIDATION MODELS ####
## Run CV parameter estimation function
# grab # of parameters from model output
num_pars <- ncol(out_model[[1]])
# get rownames from summary
summary <- summary(out_model)
par_names <- rownames(summary[[1]])
# run 10-fold CV 
CV_par_ests <- CV_parameter_estimates(CV_data='in_germ_data',
                                        number_parameters=num_pars,
                                        variable_names=par_names,
                                        data_statement=germData,
                                        model_structure=SownSurv,
                                        gelman_threshold_1=1.05,
                                        gelman_threshold_2=1.10,
                                        iters=10000)

## Copy tenfold ID to data list
germData$TenfoldID <- in_germ_data$TenfoldID 


#### CALCULATE PREDICTED VALUES ####
predicted_values <- as.data.frame(matrix(NA, nrow=nrow(germData$germ), ncol=germData$NTimes-1)) 

for(i in 1:nrow(germData$germ)){ 
  # Subset full parameter data frame to just those parameters relevant to the observation
  parameter_sub <- CV_par_ests[which(CV_par_ests$TenfoldIDExcluded==germData$TenfoldID[i]),]
  ## Calculate estimate at each timestep
  for (t in 1:(germData$NTimes-1)){
    # specify parameters to be used based on observation categories and timestep
    beta1_var <- paste0("beta1[", as.character(germData$species[i]), ",", as.character(t), "]")
    beta2_var <- paste0("beta2[", as.character(germData$source[i]), ",", as.character(t), "]")
    beta3_var <- paste0("beta3[", as.character(germData$wp[i]), ",", as.character(t), "]")
    beta4_var <- paste0("beta4[", as.character(germData$temp[i]), ",", as.character(t), "]")
    beta6_var <- paste0("beta6[", as.character(germData$cup[i]), "]")
    beta7_var <- paste0("beta7[", as.character(germData$chamber[i]), "]")
    beta8_var <- paste0("beta8[", as.character(germData$rep[i]), "]")
    
    # find categorical effects
    beta1 <- parameter_sub[1, colnames(parameter_sub)==beta1_var]
    beta2 <- parameter_sub[1, colnames(parameter_sub)==beta2_var]
    beta3 <- parameter_sub[1, colnames(parameter_sub)==beta3_var]
    beta4 <- parameter_sub[1, colnames(parameter_sub)==beta4_var]
    beta6 <- parameter_sub[1, colnames(parameter_sub)==beta6_var]
    beta7 <- parameter_sub[1, colnames(parameter_sub)==beta7_var]
    beta8 <- parameter_sub[1, colnames(parameter_sub)==beta8_var]
    
    # calculate fixed effects
    beta5 <- parameter_sub[1, colnames(parameter_sub)=="beta5"]*germData$seedmass[i]
    beta9 <- parameter_sub[1, colnames(parameter_sub)=="beta9"]*((t*2)-1)
    
    # calculate predicted value and add to dataframe
    predicted_values[i,t] <- beta1+beta2+beta3+beta4+beta5+beta6+beta7+beta8+beta9
  }#for t
}#for i


#### PLOT PREDICTED VS. OBSERVED RESULTS ####
## make combined pred-obs dataframe
pred_obs_df <- as.data.frame(matrix(NA,nrow=nrow(predicted_values), ncol=(germData$NTimes-1)*2))
# add labels
for(i in 1:(germData$NTimes-1)){
  names(pred_obs_df)[i] <- paste0("Predicted_", as.character(i))
  names(pred_obs_df)[i+germData$NTimes-1] <- paste0("Observed_", as.character(i))
}#for i
# add predicted data
pred_obs_df[,1:(germData$NTimes-1)] <- exp(predicted_values) / (1+exp(predicted_values))
# add observed data
pred_obs_df[,(germData$NTimes):ncol(pred_obs_df)] <- germData$germ[,1:germData$NTimes-1]

## convert pred v. obs data to long data format
pred_list <- c()
obs_list <- c()
time_list <- c()
for (t in 1:(germData$NTimes-1)){
  pred_list <- c(pred_list, pred_obs_df[,t])
  obs_list <- c(obs_list, pred_obs_df[,(t+germData$NTimes-1)])
  time_list <- c(time_list, rep(t, nrow(pred_obs_df)))
}#for t
long_pred_obs <- as.data.frame(matrix(NA, nrow=nrow(pred_obs_df)*(germData$NTimes-1), ncol=3))
names(long_pred_obs) <- c("Obs", "Pred", "Time")
long_pred_obs$Obs <- obs_list
long_pred_obs$Pred <- pred_list
long_pred_obs$Time <- as.factor(time_list)

## Visualize CV results
# skipping time steps 1 and 2 because none germinate
glm_pred_v_obs <- glm(Obs ~ Pred, data=long_pred_obs, family="binomial", na.action=na.omit)
summary(glm_pred_v_obs)
r <- cor(x=long_pred_obs$Pred,y=long_pred_obs$Obs)
ggplot(data=long_pred_obs, aes(x=Pred, y=Obs, col=Time))+
  #geom_abline(intercept=0, slope=1, color='black')+
  geom_point(alpha=0.)+
  ggtitle('Predicted vs. Observed')+
  #geom_smooth(method=lm,fill=NA)+
  labs(subtitle=paste0('r = ', round(r,3), "; R_sq = ", round(r^2,3)),
       xlab="Predicted",
       ylab="Observed")
  #, "; Slope = ", round(Kn_pred_v_actual$coefficients[2],3)))
plot(long_pred_obs$Pred, long_pred_obs$Obs, xlab="Predicted P(germ)", ylab="Probability of germination")
curve(predict(glm_pred_v_obs, data.frame(Pred=x), type="resp"), add=TRUE)

## Plot residuals
resids <- long_pred_obs$Obs - long_pred_obs$Pred
hist(resids, main="Histogram of Residuals (Obs-Pred)", freq=FALSE, col="grey")


#### Check "classification accuracy" ####
## Make confusion matrix
for(i in 1:nrow(long_pred_obs)){
  long_pred_obs$Pred[i] <- rbinom(1, 1, long_pred_obs$Pred[i])
}#i
# make classification df
check_class_acc <- as.data.frame(matrix(NA, nrow=2, ncol=2))
colnames(check_class_acc) <- c("Pred=0", "Pred=1")
rownames(check_class_acc) <- c("Obs=0", "Obs=1")
check_class_acc[2,2] <- length(which(long_pred_obs$Obs==1 & long_pred_obs$Pred==1))
check_class_acc[2,1] <- length(which(long_pred_obs$Obs==1 & long_pred_obs$Pred==0))
check_class_acc[1,1] <- length(which(long_pred_obs$Obs==0 & long_pred_obs$Pred==0))
check_class_acc[1,2] <- length(which(long_pred_obs$Obs==0 & long_pred_obs$Pred==1))
# convert to percentages
check_class_acc_prop <- round((check_class_acc / nrow(long_pred_obs))*100, 0)
# print confusion matrix as # obs and percents
check_class_acc_prop; check_class_acc

## Calculate transition times
# predict state based on transition probabilities
for (i in 1:nrow(pred_obs_df)){
  pred_obs_df[i,1] <- rbinom(1, 1, pred_obs_df[i,1])
  for (t in 2:(germData$NTimes-1)){
    pred_obs_df[i,t] <- pred_obs_df[i,t-1] * rbinom(1, 1, pred_obs_df[i,t])
  }#for t
}#for i
# calculate predicted and observed germinated timestep for each seed
pred_obs_df$PredGermTime <- rep(NA, nrow(pred_obs_df))
pred_obs_df$ObsGermTime <- rep(NA, nrow(pred_obs_df))
for (i in 1:nrow(pred_obs_df[,1:(germData$NTimes-1)])){
  pred_time <- which(pred_obs_df[i,1:(germData$NTimes-1)]==0)
  obs_time <- which(pred_obs_df[i,germData$NTimes:((germData$NTimes-1)*2)]==0)
  ifelse(length(pred_time)==0, 
         pred_obs_df$PredGermTime[i] <- 50, # give arbitrarily large value if never germinated
         pred_obs_df$PredGermTime[i] <- min(pred_time))
  ifelse(length(obs_time)==0, 
         pred_obs_df$ObsGermTime[i] <- 50, # give arbitrarily large value if never germinated
         pred_obs_df$ObsGermTime[i] <- min(obs_time))
}#i

## Check obs v pred transition times
lm_times <- lm(ObsGermTime~PredGermTime, data=pred_obs_df)
summary(lm_times)
plot(pred_obs_df$PredGermTime, pred_obs_df$ObsGermTime, 
     pch=20, col=alpha("blue", 0.02),
     main="Predicted vs. Observed Germination Times",
     sub=,
     xlab="Predicted",
     ylab="Observed")
abline(coef = c(0,1)) # comparison 1:1 line
abline(coef = lm_times$coefficients, col="blue") # regression line

## Check obs vs. pred without no germination ones
pred_obs_df_sub <- pred_obs_df[which(pred_obs_df$ObsGermTime!=50 & pred_obs_df$PredGermTime!=50),]
lm_times <- lm(ObsGermTime~PredGermTime, data=pred_obs_df_sub)
summary(lm_times)
plot(pred_obs_df_sub$PredGermTime, pred_obs_df_sub$ObsGermTime, 
     pch=20, col=alpha("blue", 0.02),
     main="Predicted vs. Observed Germination Times",
     sub=,
     xlab="Predicted",
     ylab="Observed")
abline(coef = c(0,1))
abline(coef = lm_times$coefficients, col="blue")
