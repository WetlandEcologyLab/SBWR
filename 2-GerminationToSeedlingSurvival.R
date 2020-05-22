#### LOAD PACKAGES ####
rm(list=ls()) # cleanup
library(rjags)
require(devtools)
library(jagsUI)
library(coda)
library(nimble)

##############################################################
### DEFINE SURVIVAL PROBABILITY HIERARCHICAL PROBABILITIES ###
##############################################################

#########################################################
### Transition from germinated seed to 7-day seedling ###
#########################################################
GermSurv <- textConnection("model{ # has the same effect as writing all this in a separate text file

  ########################### 
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
  for (i in 1:NGerm){
    for (t in 1:NTimes){
      germ_linear[i,t] <- beta1[species[i],t] + beta2[source[i],t] + beta3[wp[i],t]+ beta4[temp[i],t] + beta5*pc1[i] + beta6*pc2[i] + beta7*pc3[i] + beta8*pc4[i] + beta9[cup[i]] + beta10[chamber[i]] + beta11[rep[i]] + beta12*((t-1)*2)
    
      # transformation on logit scale for logistic regression
      prob_early[i,t] <- exp(germ_linear[i,t]) / (1+exp(germ_linear[i,t]))
    
    }#t
  }#i
  
  ## Likelihood parameter depends on previous observation
  for (i in 1:NGerm){ 
    for (t in 2:NTimes){
      early_survival[i,t] ~ dbern(prob_early[i,t] * early_survival[i,t-1]) 
    }#t
  }#i

  ##############################
  ### PRIORS ON MAIN EFFECTS ###
  ##############################
  ## species effect - categorical
  for (t in 1:NTimes){
    for (spp in 1:NSpecies){
      beta1[spp,t] ~ dnorm(kappa1[spp], tau1[spp]) 
    }#spp
  }#t

  # conditional prior on mean species effect
  for (spp in 1:NSpecies){
    kappa1[spp] ~ dnorm(mu1, sigma1)
  }#spp

  # independent prior on variance of species effect
  for (spp in 1:NSpecies){
    tau1[spp] ~ dt(0, pow(5,-2), 1)T(0,)
  }#spp

  # hyperpriors on overall species mean
  mu1 ~ dnorm(0, 0.01)
  sigma1 ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## source effect - categorical
  for (t in 1:NTimes){ 
    for (src in 1:(NSources-1)){
      beta2[src,t] ~ dnorm(kappa2[src], tau2[src])
    }#src
    beta2[NSources,t] <- -1 * sum(beta2[1:(NSources-1),t]) #sum-to-zero constraint at each timestep
  }#t

  # conditional prior on mean source effect
  for (src in 1:(NSources-1)){
    kappa2[src] ~ dnorm(mu2, sigma2)
  }#src

  # independent priors on variance of source effect
  for (src in 1:(NSources-1)){
    tau2[src] ~ dt(0, pow(2.5,-2), 1)T(0,)
  }#src
  
  # hyperpriors on overall mean source effect
  mu2 ~ dnorm(0, 0.01)
  sigma2 ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## water potential effect - categorical
  for (t in 1:NTimes){
    beta3[1,t] <- 0
    beta3[2,t] ~ dnorm(kappa3, tau3)
  }#t

  # hyperpriors on wp2 mean effect
  kappa3 ~ dnorm(0, 0.01)
  tau3 ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## temperature effect - categorical
  for (t in 1:NTimes){
    for (tmp in 1:2){ 
      beta4[tmp,t] ~ dnorm(kappa4[tmp], tau4[tmp]) 
    }#tmp
    beta4[3,t] <- -1 * sum(beta4[1:2,t]) #sum-to-zero constraint
  }#t

  # conditional prior on temp means
  for (tmp in 1:2){
    kappa4[tmp] ~ dnorm(mu3, sigma3)
  }#tmp
  
  # independent priors on temp variances
  for (tmp in 1:2){
    tau4[tmp] ~ dt(0, pow(2.5,-2), 1)T(0,)
  }#tmp

  # hyperpriors for overall temp mean
  mu3 ~ dnorm(0, 0.01)
  sigma3 ~ dt(0, pow(2.5,-1), 1)T(0,)

  ## imputation model for PC1
  #for (i in 1:NGerm){
  #  # linear model on PC1 value
  #  pc1_linear[i] <- delta0 + delta1[species[i]] + delta2[source[i]] + delta3[wp[i]] + delta4[temp[i]]
  #  # transformation
  #  prob_pc1[i] <- exp(pc1_linear[i]) / (1+exp(pc1_linear[i]))
  #}#i
  # PC1 likelihood
  #for (i in 1:NGerm){
  #  pc1[i] ~ dbern(prob_pc1[i])
  #}#i
  # priors on effects on PC1
  #delta0 ~ dnorm(0, 0.01)
  #for (spp in 1:(NSpecies-1)){
  #  delta1[spp] ~ dnorm(0, 0.01)
  #}#spp
  #delta1[NSpecies] <- -1 * sum(delta1[1:(NSpecies-1)]) #sum-to-zero constraint
  #for (src in 1:(NSources-1)){
  #  delta2[src] ~ dnorm(0, 0.01)
  #}#src
  #delta2[NSources] <- -1 * sum(delta2[1:(NSources-1)]) #sum-to-zero constraint
  #delta3[1] ~ dnorm(0, 0.01)
  #delta3[2] <- -1 * sum(delta3[1:2]) #sum-to-zero constraint
  #for (tmp in 1:2){
  #  delta4[tmp] ~ dnorm(0, 0.01)
  #}#tmp
  #delta4[3] <- -1 * sum(delta4[1:2]) #sum-to-zero constraint

  # prior on PC1 - continuous fixed effect
  beta5 ~ dnorm(0, 0.01)

  ## imputation model for PC2
  #for (i in 1:NGerm){
  #  # linear model on PC2 value
  #  pc2_linear[i] <- eta0 + eta1[species[i]] + eta2[source[i]] + eta3[wp[i]] + eta4[temp[i]]
  #  # transformation
  #  prob_pc2[i] <- exp(pc2_linear[i]) / (1+exp(pc2_linear[i]))
  #}#i
  # PC2 likelihood
  #for (i in 1:NGerm){
  #  pc2[i] ~ dbern(prob_pc2[i])
  #}#i
  # priors on effects on PC2
  #eta0 ~ dnorm(0, 0.01)
  #for (spp in 1:(NSpecies-1)){
  #  eta1[spp] ~ dnorm(0, 0.01)
  #}#spp
  #eta1[NSpecies] <- -1 * sum(eta1[1:(NSpecies-1)]) #sum-to-zero constraint
  #for (src in 1:(NSources-1)){
  #  eta2[src] ~ dnorm(0, 0.01)
  #}#src
  #eta2[NSources] <- -1 * sum(eta2[1:(NSources-1)]) #sum-to-zero constraint
  #eta3[1] ~ dnorm(0, 0.01)
  #eta3[2] <- -1 * sum(eta3[1]) #sum-to-zero constraint
  #for (tmp in 1:2){
  #  eta4[tmp] ~ dnorm(0, 0.01)
  #}#tmp
  #eta4[3] <- -1 * sum(eta4[1:2]) #sum-to-zero constraint

  ## prior on PC2 - continuous fixed effect
  beta6 ~ dnorm(0, 0.01)

  ## imputation model for PC3
  #for (i in 1:NGerm){
  #  # linear model on PC3 value
  #  pc3_linear[i] <- iota0 + iota1[species[i]] + iota2[source[i]] + iota3[wp[i]] + iota4[temp[i]]
  #  # transformation
  #  prob_pc3[i] <- exp(pc3_linear[i]) / (1+exp(pc3_linear[i]))
  #}#i
  # PC3 likelihood
  #for (i in 1:NGerm){
  #  pc3[i] ~ dbern(prob_pc3[i])
  #}#i
  # priors on effects on PC3
  #iota0 ~ dnorm(0, 0.01)
  #for (spp in 1:(NSpecies-1)){
  #  iota1[spp] ~ dnorm(0, 0.01)
  #}#spp
  #iota1[NSpecies] <- -1 * sum(iota1[1:(NSpecies-1)]) #sum-to-zero constraint
  #for (src in 1:(NSources-1)){
  #  iota2[src] ~ dnorm(0, 0.01)
  #}#src
  #iota2[NSources] <- -1 * sum(iota2[1:(NSources-1)]) #sum-to-zero constraint
  #iota3[1] ~ dnorm(0, 0.01)
  #iota3[2] <- -1 * sum(iota3[1]) #sum-to-zero constraint
  #for (tmp in 1:2){
  #  iota4[tmp] ~ dnorm(0, 0.01)
  #}#tmp
  #iota4[3] <- -1 * sum(iota4[1:2]) #sum-to-zero constraint

  ## prior on PC3 - continuous fixed effect
  beta7 ~ dnorm(0, 0.01)

  ## imputation model for PC4
  #for (i in 1:NGerm){
  #  # linear model on PC4 value
  #  pc4_linear[i] <- zeta0 + zeta1[species[i]] + zeta2[source[i]] + zeta3[wp[i]] + zeta4[temp[i]]
    # transformation
  #  prob_pc4[i] <- exp(pc4_linear[i]) / (1+exp(pc4_linear[i]))
  #}#i
  # PC4 likelihood
  #for (i in 1:NGerm){
  #  pc4[i] ~ dbern(prob_pc4[i])
  #}#i
  # priors on effects on PC4
  #zeta0 ~ dnorm(0, 0.01)
  #for (spp in 1:(NSpecies-1)){
  #  zeta1[spp] ~ dnorm(0, 0.01)
  #}#spp
  #zeta1[NSpecies] <- -1 * sum(zeta1[1:(NSpecies-1)]) #sum-to-zero constraint
  #for (src in 1:(NSources-1)){
  #  zeta2[src] ~ dnorm(0, 0.01)
  #}#src
  #zeta2[NSources] <- -1 * sum(zeta2[1:(NSources-1)]) #sum-to-zero constraint
  #zeta3[1] ~ dnorm(0, 0.01)
  #zeta3[2] <- -1 * sum(zeta3[1]) #sum-to-zero constraint
  #for (tmp in 1:2){
  #  zeta4[tmp] ~ dnorm(0, 0.01)
  #}#tmp
  #zeta4[3] <- -1 * sum(zeta4[1:2]) #sum-to-zero constraint

  ## prior on PC4 - continuous fixed effect
  beta8 ~ dnorm(0, 0.01)

  ## cup effect - hierarchical
  for (c in 1:(NCups-1)){ 
    beta9[c] ~ dnorm(0, nu[1]) 
  }#c
  beta9[NCups] <- -1 * sum(beta9[1:(NCups-1)]) #sum-to-zero constraint

  ## chamber effect - hierarchical
  for (chm in 1:(NChambers-1)){ 
    beta10[chm] ~ dnorm(0, nu[2]) 
  }#chm
  beta10[NChambers] <- -1 * sum(beta10[1:(NChambers-1)]) #sum-to-zero constraint

  ## rep effect - categorical
  for (rep in 1:2){ 
    beta11[rep] ~ dnorm(0, nu[3]) 
  }#rep
  beta11[3] <- -1 * sum(beta11[1:2]) #sum-to-zero constraint

  ## days effect - continuous fixed
  beta12 ~ dnorm(0, 0.01)

  ################################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIABILITY ###
  ################################################
  nu[1] ~ dt(0, pow(2.5,-2), 1)T(0,)
  nu[2] ~ dt(0, pow(2.5,-2), 1)T(0,)
  nu[3] ~ dt(0, pow(2.5,-2), 1)T(0,)

}")#end model


#########################################
### RUN PRINCIPAL COMPONENTS ANALYSES ###
#########################################
# from Emily, include the following: htt, synchrony(pop), mean germination rate(pop), time to germination, seed mass(pop), sct(pop), total mass, total length, shoot mass ratio, shoot length ratio, shoot elongation rate(pop), root elongation rate(pop)
# currently included: htt, seedmass, sct, germ time percentile, root length, shoot length, root SA, shoot, SA, total length, total SA, shoot elongation rate, root elongation rate, shoot length ratio, above ground fresh mass, below ground fresh mass, above ground dry mass, below ground dry mass, total dry mass, shoot mass ratio
# missing: synchrony, mean germination rate

# read in data
early_seedling_data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixEarlySeedlingData_21May2020.csv', header=TRUE, stringsAsFactors=FALSE)
# pull variables for pca
seedling_traits <- early_seedling_data[,c(50:52, 55:70)]
# remove observations with no data
seedling_traits <- seedling_traits[which(!is.na(seedling_traits$RootLength)),]

# run PCA
seedling_pca <- prcomp(seedling_traits)
summary(seedling_pca)
# using the first 4 terms from this PCA which explain 98.8% of the variation,
# using cutoff of > 1% variance explained for a single PC
# this cutoff is completely arbitrary though

# pull coefficients for principal components that describe most of variation
pc1_rotation <- as.matrix(seedling_pca$rotation[,1])
pc2_rotation <- as.matrix(seedling_pca$rotation[,2])
pc3_rotation <- as.matrix(seedling_pca$rotation[,3])
pc4_rotation <- as.matrix(seedling_pca$rotation[,4])

# calculate principal components from coefficients
pc1 <- as.matrix(seedling_traits) %*% pc1_rotation
pc2 <- as.matrix(seedling_traits) %*% pc2_rotation
pc3 <- as.matrix(seedling_traits) %*% pc3_rotation
pc4 <- as.matrix(seedling_traits) %*% pc4_rotation

# merge pc values back into data frame
early_seedling_data$PC1 <- NA
early_seedling_data$PC1[which(!is.na(early_seedling_data$RootLength))] <- pc1
early_seedling_data$PC2 <- NA
early_seedling_data$PC2[which(!is.na(early_seedling_data$RootLength))] <- pc2
early_seedling_data$PC3 <- NA
early_seedling_data$PC3[which(!is.na(early_seedling_data$RootLength))] <- pc3
early_seedling_data$PC4 <- NA
early_seedling_data$PC4[which(!is.na(early_seedling_data$RootLength))] <- pc4

#################
### RUN MODEL ###
#################
#### JAGS FORMAT ####
# save data as appropriate list
seedling7Data <- list(early_survival = early_seedling_data[,1:41], 
                     species = as.integer(as.factor(early_seedling_data$Species)),
                     source = as.integer(as.factor(early_seedling_data$Source)), 
                     wp = early_seedling_data$WP, 
                     temp = early_seedling_data$Temp,
                     cup = as.integer(as.factor(early_seedling_data$CupNo)), 
                     chamber = as.integer(as.factor(early_seedling_data$Chamber)),
                     rep = as.integer(early_seedling_data$Trial),
                     pc1 = as.numeric(early_seedling_data$PC1),
                     pc2 = as.numeric(early_seedling_data$PC2),
                     pc3 = as.numeric(early_seedling_data$PC3),
                     pc4 = as.numeric(early_seedling_data$PC4),
                     NGerm = nrow(early_seedling_data),
                     NSpecies = length(unique(germination_matrix$Species)), 
                     NSources = length(unique(germination_matrix$Source)),
                     NCups = length(unique(germination_matrix$CupNo)),
                     NTimes = 41,
                     NChambers = length(unique(germination_matrix$Chamber)))

## Initialize model
seedling7Model <- jags.model(file=GermSurv, data=seedling7Data, n.chains = 3)#adapt=500
## Run model
update(seedling7Model, n.iter=1000)
out_model <- coda.samples(model=germModel, variable.names = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10", "beta11", "beta12", "kappa1", "kappa2", "kappa3", "kapp4", "tau1", "tau2", "tau3", "tau4", "mu1", "mu2", "sigma1", "sigma2", "nu"), thin=1, n.iter = 1000)

#### Assess convergence - save to PDF (R plot console runs out of space) ####
## Check out convergence plots
# set up PDF file
#library(grid)
#pdf(file='C:/Users/Maggie/Documents/WetlandEcology/RevegModel/GermModel_4May2020_2.pdf')
# add text to top of PDF
#title <- "MODEL PARAMETRIZATION"
#a <- "Linear Model: beta0 + beta1[species[i],t]+beta2[source[i],t]+\nbeta3[wp[i],t]+beta4[temp[i],t]+beta5*seedmass[i]+beta6*sct[i]\n+beta7[cup[i]]+beta8[chamber[i]]"
#b <- "Transformation: prob_sown[i] <- exp(sown_linear[i]) / (1+exp(sown_linear[i])"
#c <- "Response likelihood: germ[i] ~ dbern(prob_sown[i])"
#d <- "CONTENTS"
#e <- "  - Convergence plots"
#content <- paste(title, "\n", "\n", a, "\n", b, "\n", c, "\n", "\n", d, "\n", e, "\n")
#grid.text(content, x = 0.05, hjust = 0, gp = gpar(fontsize = 11))
# plot trace plots to PDF
orig.mar <- par('mar') # save default plot margins to reset later
par(mar=c(1,1,1,1)) # change plot margins
plot(out_model, ask=TRUE) # plot chains
dev.off() # close PDF
par(mar=orig.mar) # return plot margins to normal

## check effective sample size- should be on the order of 1000s, >300 generally ok
effectiveSize(out_model) 

## check Gelman-Rubin diagnostic plots
orig.mar <- par('mar') # save normal plot margins
par(mar=c(1,1,1,1)) # change plot margins
gelman.diag(out_model, multivariate=F) # should be close to 1
gelman.plot(out_model) # plot diagnostic plots
par(mar=orig.mar) # return plot margins to normal

###########################################
### BACK-CALCULATE PRINCIPAL COMPONENTS ###
###########################################
rho1 <- seedling_pca[1,1] + seedling_pca[1,2] + seedling_pca[1,3] + seedling_pca[1,4]
rho2 <- seedling_pca[2,1] + seedling_pca[2,2] + seedling_pca$ rotation[2,3] + seedling_pca[2,4]
rho3 <- seedling_pca[3,1] + seedling_pca[3,2] + seedling_pca$ rotation[3,3] + seedling_pca[3,4]
rho4 <- seedling_pca[4,1] + seedling_pca[4,2] + seedling_pca[4,3] + seedling_pca[4,4]
rho5 <- seedling_pca[5,1] + seedling_pca[5,2] + seedling_pca[5,3] + seedling_pca[5,4]
rho6 <- seedling_pca[6,1] + seedling_pca[6,2] + seedling_pca[6,3] + seedling_pca[6,4]
rho7 <- seedling_pca[7,1] + seedling_pca[7,2] + seedling_pca[7,3] + seedling_pca[7,4]
rho8 <- seedling_pca[8,1] + seedling_pca[8,2] + seedling_pca[8,3] + seedling_pca[8,4]
rho9 <- seedling_pca[9,1] + seedling_pca[9,2] + seedling_pca[9,3] + seedling_pca[9,4]
rho10 <- seedling_pca[10,1] + seedling_pca[10,2] + seedling_pca[10,3] + seedling_pca[10,4]
rho11 <- seedling_pca[11,1] + seedling_pca[11,2] + seedling_pca[11,3] + seedling_pca[11,4]
rho12 <- seedling_pca[12,1] + seedling_pca[12,2] + seedling_pca[12,3] + seedling_pca[12,4]
rho13 <- seedling_pca[13,1] + seedling_pca[13,2] + seedling_pca[13,3] + seedling_pca[13,4]
rho14 <- seedling_pca[14,1] + seedling_pca[14,2] + seedling_pca[14,3] + seedling_pca[14,4]
rho15 <- seedling_pca[15,1] + seedling_pca[15,2] + seedling_pca[15,3] + seedling_pca[15,4]
rho16 <- seedling_pca[16,1] + seedling_pca[16,2] + seedling_pca[16,3] + seedling_pca[16,4]
rho17 <- seedling_pca[17,1] + seedling_pca[17,2] + seedling_pca[17,3] + seedling_pca[17,4]
rho18 <- seedling_pca[18,1] + seedling_pca[18,2] + seedling_pca[18,3] + seedling_pca[18,4]
rho19 <- seedling_pca[19,1] + seedling_pca[19,2] + seedling_pca[19,3] + seedling_pca[19,4]
rho20 <- seedling_pca[20,1] + seedling_pca[20,2] + seedling_pca[20,3] + seedling_pca[20,4]

#### View estimates ####
## summarize outputs
summary <- summary(out_model)
## compare single species-timestep means to modeled species-timestep effects...
chain1 <- out_model[[1]]
chain2 <- out_model[[2]]
chain3 <- out_model[[3]]

