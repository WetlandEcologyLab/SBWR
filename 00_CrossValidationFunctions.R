################################
## CROSS VALIDATION FUNCTIONS ##
################################
## Credited to Kenen Goodwin
library(rjags)

############################
## SUBSETTING DATA FOR CV ##
############################
## Create function for subsetting 10-fold CV data such that no levels of
## categorical predictors are lost. 
## i.e., enables representative sample is taken from categories
protected_CV_subset <- function(CV_data, colnum){
## CV_data: data to cross validate on.
## colnum: column numbers for which at least one observation is desired for each level.
## Returns vector of tenfold CV IDs.
  # Create empty storage vector.
  good <- c() 
  # Repeat until all subsets are good.
  while(length(good)==0 | 'Discard' %in% good){ 
    # Clear empty storage vector.
    good <- c() 
    # Label 10-fold CV subsets.
    CV_data$TenfoldID <- sample(rep(1:10,length=nrow(CV_data)), size=nrow(CV_data), replace=F)
    # Loop through all 10-fold CV subsets.
    for(i in 1:10){ 
      # Subset by 10-fold CV ID.
      sub <- subset(CV_data,TenfoldID!=i) 
      # Loop through specified columns.
      for(j in colnum){ 
        # Ensure that the same number of levels are present in the original and subsetted data.
        good <- c(good,ifelse(length(unique(sub[,j]))==length(unique(CV_data[,j])),'Good','Discard'))
      }#for j
    }#for i
  }#while
  
  return(CV_data$TenfoldID) # Return good subset label numbers.
}#function


###################################
## DEFINING PARAMETER ESTIMATION ##
###################################
# Define CV parameter estimation function.
CV_parameter_estimates <- function(CV_data, number_parameters, variable_names, data_statement, model_structure, gelman_threshold_1=1.05, gelman_threshold_2=1.10, iters=10000){
# CV_data: name as text of data frame to cross validate on. Must have a 'TenfoldID' field.
# number_parameters: number of parameters being estimated.
# variable_names: names of parameters being estimated. For beta4[1], beta4[2], etc., just put beta4.
# data_statement: JAGS data list code stored as text.
# model_structure: JAGS model structure stored as text, not as a text connection.
# gelman_threshold_1: first threshold for Gelman diagnostic
# gelman_Threshold_2: second higher threshold for Gelman diagnostic 
# iters: # of iterations to run through coda.samples
# Returns data frame of parameter point estimates for each tenfold CV ID.

  #### SET UP TO RUN MODELS ####
  # Create empty storage vector for Gelman Diagnostics beyond thresholds
  gelman_gt_1 <- c() 
  gelman_gt_2 <- c() 
  # set up empty storage vector for saving parameter estimates
  point_estimates_df_full <- as.data.frame(matrix(ncol=number_parameters,nrow=0))
  names(point_estimates_df_full) <- variable_names
  # Change data statement to work with subsets
  data_statement <- gsub(CV_data, 'sub', data_statement) 
  # Get data frame by name
  CV_data <- eval(parse(text=CV_data)) 
  
  # Evaulate model for each data subset
  for(i in 1:10){
    #### RUN MODEL ####
    # Subset data by 10-fold CV ID
    sub <- subset(CV_data, TenfoldID!=i) 
    # Get model data
    mod_data <- eval(parse(text=data_statement)) 
    mod_structure_eval <- eval(parse(text=model_structure))
    # Compile model
    model <- jags.model(mod_structure_eval, data=mod_data, n.chains=3) 
    # Run model
    update(model,n.iter=1000) 
    # Initiate model
    CV_Out <- coda.samples(model=model, variable.names=variable_names, n.iter=iters, thin=3) 
    
    #### EVALUATE GELMAN DIAGNOSTICS ####
    # Calc # parameters with Gelman Diagnostic greater than threshold 1
    gelman_gt_1 <- c(gelman_gt_1.05, sum(gelman.diag(CV_Out,multivariate=F)$psrf[,2] > gelman_threshold_1))
    # Calc # parameters with Gelman Diagnostic greater than threshold 2
    gelman_gt_2 <- c(gelman_gt_1.10, sum(gelman.diag(CV_Out,multivariate=F)$psrf[,2] > gelman_thresdhold_2))
    
    #### SAVE POINT ESTIMATES ####
    # Combine MCMC chains
    MCMC_Out <- rbind(CV_Out[[1]], CV_Out[[2]], CV_Out[[3]]) 
    point_estimates <- c() 
    # Loop through parameters and store point estimates
    out_summary <- summary(CV_Out)[[1]]
    out_summary_df <- as.data.frame(out_summary[,1], ncol=1)
    names(out_summary_df) <- paste('CV_', as.character(i), sep='')
    point_estimates <- out_summary_df#cbind(point_estimates, out_summary_df)
    # Turn point estimates into a dataframe
    point_estimates_df <- as.data.frame(t(point_estimates),nrow=1) 
    # Store point estimates
    point_estimates_df_full <- rbind(point_estimates_df_full, point_estimates_df)
  }#for i
  
  #### SAVE FULL OUTPUT ####
  # Give a 10-fold CV ID
  point_estimates_df_full$TenfoldIDExcluded <- 1:nrow(point_estimates_df_full) 
  # Rearrange columns within dataframe. Place the 10-fold CV ID in first column
  point_estimates_df_full <- point_estimates_df_full[,c(ncol(point_estimates_df_full), 1:(ncol(point_estimates_df_full)-1))]
  
  #### PRINT GELMAN STATISTICS ####
  print('##############################')
  print('Gelman Convergence Diagnostic:')
  for(i in 1:10){
    print('--------')
    print(paste0('Fold ',i,':'))
    print('--------')
    print(paste0('>', as.character(gelman_threshold_1),' = ',gelman_gt_1.05[i]))
    print(paste0('>', as.character(gelman_threshold_2),' = ',gelman_gt_1.05[i]))
  }#for i
  
  return(point_estimates_df_full) # Return point estimates as data frame.
  
}#function



