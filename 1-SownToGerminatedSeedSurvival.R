####################################################
### Transition from sown seed to germinated seed ###
####################################################

#### SET FILE PATHS ####
## path to input germination matrix
in_data_path <- read.csv("./00_ModelData/01_Raw_Checked_Data/01_Seed_To_Germination_Model/FullGerminationData_FINAL.csv")
## folder to save outputs to
output_folder <- './03_Model_Output_Data/01_Seed_To_Germination_Model'
## model name for output files
model_name <- 'GermModel_01.2022'

#### LOAD PACKAGES ####
#rm(list=ls()) # cleanup
library(rjags)
require(devtools)
library(jagsUI)
library(coda)
library(nimble)


#### DEFINE MODEL ####
SownSurv <- 
# use this line if running in nimble (remember to remove quotes at end bracket)
#nimbleCode({
# use this line instead if running in rjags
textConnection("model{ # has the same effect as writing all this in a separate text file
  ###########################
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter (prob of NOT transitioning from sown to germinated)
  #     based on each observation
  # germ = response matrix of binary germination observations
  #   rows=inds, cols=timesteps
  # NSown = # individual seeds sown (nrows in germination matrix)
  # NTimes = # timesteps (ncols in germination matrix)
  # beta1, 2, 3, 4, 7, 8 have multiple estimates for each species category
  # beta5, 6, 9 have single coefficients and use each observations hydrothermal time and seed mass value as a continuous parameter 
  
  for (i in 1:NSown){#i = individual seed
    for (t in 1:NTimes){#t = timestep

      ## Define linear model
      ##NOTE: beta9 defined by ((t-1)*2) because timesteps start at 0 (t-1)
      ## and proceed by intervals of 2
      sown_linear[i,t] <- beta1[species[i],t]+ beta2[source[i],t]+ beta3[wp[i],t] + beta4[temp[i],t] + beta5*seedmass[i] + beta6[cup[i]] + beta8[rep[i]]  + beta9*((t*2)-1) #+ beta7[chamber[i]]

    ## Transformation to logit scale for logistic regression
    prob_sown[i,t] <- exp(sown_linear[i,t]) / (1+exp(sown_linear[i,t]))

    }#t
  }#i
  
  ## Likelihood parameter depends on previous observation
  for (i in 1:NSown){
    for (t in 2:NTimes){
      germ[i,t] ~ dbern(prob_sown[i,t] * germ[i,t-1])
    }#i
  }#t

  ##############################
  ### PRIORS ON MAIN EFFECTS ###
  ##############################  

  ## species effect on germination probability - categorical - varies by time
  for (t in 1:NTimes){
    for (spp in 1:NSpecies){ 
      beta1[spp,t] ~ dnorm(kappa1[spp], pow(tau1[spp],2))
  }#spp

  ## sum-to-zero constraint not needed because no beta0- now beta1 species effect is taking up some of the intercept term
  #beta1[NSpecies,t] <- -1 * sum(beta1[1:(NSpecies-1),t]) # sum-to-zero constraint on species axis

  }#t

  # conditional prior on species means
  for (spp in 1:NSpecies){
    kappa1[spp] ~ dnorm(mu1, pow(sigma1,2)) # conditional prior
  }#spp

  # independent prior on species variances
  for (spp in 1:NSpecies){
    tau1[spp] ~ dt(0, pow(5,-2), 1)T(0,) # independent prior
  }#spp
  
  # uninformative hyperpriors on species mean
  mu1 ~ dnorm(0, 0.01)
  sigma1 ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## source effect on germination probability - hierarchical - varies by time
  for (t in 1:NTimes){
    for (src in 1:(NSources-1)){ 
      beta2[src,t] ~ dnorm(kappa2[src], pow(tau2[src],2)) 
    }#src
  beta2[NSources,t] <- -1 * sum(beta2[1:(NSources-1),t]) # sum-to-zero constraint
  }#t

  # conditional priors on source means
  for (src in 1:(NSources-1)){
    kappa2[src] ~ dnorm(mu2, pow(sigma2,2))
  }#src
    
  # independent priors on source variances
  for (src in 1:(NSources-1)){
    tau2[src] ~ dt(0, pow(2.5,-2), 1)T(0,)
  }#src

  # uninformative hyperpriors on source mean
  mu2 ~ dnorm(0, 0.01)
  sigma2 ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## wp effect on germination - hierarchical - varies by time
  for (t in 1:NTimes){
    beta3[1,t] <- 0
    beta3[2,t] ~ dnorm(kappa3, pow(tau3,2))
    #beta3[2,t] <- -1 * beta3[1,t] # sum-to-zero constraint
  }#t

  # independent priors on wp mean and variance
  kappa3 ~ dnorm(0, 0.01)
  tau3 ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## temp effect on germination - hierarchical - varies by time
  for (t in 1:NTimes){
    for (tmp in 1:2){
      beta4[tmp,t] ~ dnorm(kappa4[tmp], pow(tau4[tmp],2))
    }#tmp
    beta4[3,t] <- -1 * sum(beta4[1:2,t]) # sum-to-zero constraint
  }#t
  
  # independent priors on temp means
  for (tmp in 1:2){
    kappa4[tmp] ~ dnorm(0, 0.01)
  }#tmp

  # independent priors on temp variances
  for (tmp in 1:2){
    tau4[tmp] ~ dt(0, pow(2.5,-2), 1)T(0,)
  }

  ## hydrothermal time effect on germination - fixed linear - slightly constrained variance prior
  #beta3 ~ dnorm(0, 0.01)

  ## seed mass effect on germination - fixed linear - slightly constrained variance prior
  beta5 ~ dnorm(0, 0.01)

  ## seed coat thickness effect on germination - fixed linear - slightly constrained variance prior
  # NOTE: removed because estimate was centered on zero
  #beta6 ~ dnorm(0, 0.01)

  ## cup effect on germination - hierarchical
  for (c in 1:(NCups-1)){ 
    beta6[c] ~ dnorm(0, pow(nu[1],2)) 
  }#cup
  beta6[NCups] <- -1 * sum(beta6[1:(NCups-1)]) # sum-to-zero constraint

  ## chamber effect on germination - hierarchical
  #for (chm in 1:(NChambers-1)){ 
  #  beta7[chm] ~ dnorm(0, pow(nu[2],2)) 
  #}#chm
  #beta7[NChambers] <- -1 * sum(beta7[1:(NChambers-1)]) # sum-to-zero constraint

  ## rep effect on germination - hierarchical
  for (rep in 1:2){ 
    beta8[rep] ~ dnorm(0, pow(nu[3],2)) 
  }
  beta8[3] <- sum(beta8[1:2]) # sum-to-zero constraint

  ## days effect - fixed linear - slightly constrained
  beta9 ~ dnorm(0, 0.01)

  ##############################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIATION ###
  ##############################################
  nu[1] ~ dt(0, pow(2.5,-2), 1)T(0,)
  #nu[2] ~ dt(0, pow(2.5,-2), 1)T(0,)
  nu[3] ~ dt(0, pow(2.5,-2), 1)T(0,)

}")#end model block


#### LOAD DATA ####
## read in data- matrix format
germination_matrix <- read.csv(in_data_path, header=TRUE, stringsAsFactors=FALSE)
head(germination_matrix) # check


#### RUN MODEL - JAGS FORMAT ####
## save data as appropriate list
germData <- list(germ=germination_matrix[,7:43], 
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

## Initialize model
germModel <- jags.model(file=SownSurv, data=germData, n.chains=3, n.adapt=1000)

## Run model
update(germModel, n.iter=1000)
out_model <- coda.samples(model=germModel, variable.names = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "kappa1", "tau1", "kappa2", "tau2", "kappa3", "tau3", "kappa4", "tau4", "mu1", "mu2", "sigma1", "sigma2", "nu"), thin=10, n.iter=10000)


#### Save outputs ####
setwd(output_folder)

## Save output model
saveRDS(out_model, file=paste(model_name,".RDS", sep=""))

## Save plots to file
pdf(file=paste(model_name,".pdf", sep=""))
orig.mar <- par('mar') # save default plot margins to reset later
par(mar=c(1,1,1,1)) # change plot margins
plot(out_model) # plot outputs - chains and posteriors
dev.off() # close PDF
par(mar=orig.mar) # return plot margins to normal

## Save summary to CSV
summary <- summary(out_model)
write.csv(summary[[1]], paste(model_name, "_Summary.csv", sep=""))
write.csv(summary[[2]], paste(model_name, "_Quantiles.csv", sep=""))

#### Assess Convergence ####
## Check effective sample size- should be on the order of 1000s, >300 generally ok
effectiveSize(out_model) 

## Check Gelman-Rubin diagnostic plots
orig.mar <- par('mar') # save normal plot margins
par(mar=c(1,1,1,1)) # change plot margins
gelman.diag(out_model, multivariate=F) # should be close to 1
gelman.plot(out_model) # plot diagnostic plots
par(mar=orig.mar) # return plot margins to normal




#### RUN MODEL - NIMBLE FORMAT ####
{## GOOD WEBSITE FOR NIMBLE TUTORIAL: https://r-nimble.org/html_manual/cha-lightning-intro.html

## set up data
germData <- list(germ=germination_matrix[,5:10],
                 species=as.integer(as.factor(germination_matrix$Species)), 
                 source=as.integer(as.factor(germination_matrix$Source)),
                 #htt=germination_matrix$HTT,
                 temp=germination_matrix$Temp,
                 wp=germination_matrix$WP,
                 #rep=germination_matrix[,39],
                 cup=as.integer(as.factor(germination_matrix$CupNo)),
                 chamber=as.integer(as.factor(germination_matrix$Chamber)),
                 seedmass=germination_matrix$SeedMass,
                 sct=germination_matrix$SCT,
                 NSown=6185, #nrow(germination_matrix),
                 NSpecies=7, #length(unique(germination_matrix[,35])), 
                 NSources=38, #length(unique(germination_matrix[,36])),
                 NCups=216, #length(unique(germination_matrix[,40])),
                 NChambers=3, #length(unique(germination_matrix[,42])),
                 NTimes=5)

## create model
germModel <- nimbleModel(code = SownSurv,
                          name = "Survival to Germination",
                          constants = germData, # just like your usual jags.data
                          inits = 1000) # just like your usual inits function

## check out model components
germModel$getNodeNames()
germ$plotGraph() # directed acyclic graph

## show dependencies
germModel$getDependencies(c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta8", "beta9")) # all dependencies to betas
germModel$getDependencies(c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta8", "beta9"), determOnly=TRUE) # only deterministic dependencies

## check that the lifted nodes were initialized
germModel[[]]

## generate MCMC
germ_mcmc <- buildMCMC(germModel,
                          niter = 500,
                          nburnin = 250,
                          nchain = 3)

compiled_germ_model <- compileNimble(germModel)
compiled_germ_mcmc <- compileNimble(germ_mcmc)

samples <- runMCMC(compiled_germ_mcmc,
                   niter = 10000,#30000,
                   nburnin = 5000,#15000,
                   nchain = 3)

samples$chain1 <- mcmc(samples$chain1)
samples$chain2 <- mcmc(samples$chain2)
samples$chain3 <- mcmc(samples$chain3)

## view model outputs
samples_coda <- mcmc.list(samples)
sum_code <- summary(samples_coda)

summary(samples_coda)
gelman.diag(samples_coda)
plot(samples_coda)

}#group Nimble format


#test