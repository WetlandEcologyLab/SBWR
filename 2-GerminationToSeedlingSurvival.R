#### LOAD PACKAGES ####
rm(list=ls()) # cleanup
library(rjags)
require(devtools)
library(jagsUI)
library(coda)
library(nimble)

##############################
### CONVERT DATA TO MATRIX ###
##############################
data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/', header=TRUE, stringsAsFactors=FALSE)
# select only 7-day seedling data

# convert to time-to-event

# convert to matrix

# add covariates

#########################################
### RUN PRINCIPAL COMPONENTS ANALYSES ###
#########################################
# htt, synchrony(pop), mean germination rate(pop), time to germination, seed mass(pop), sct(pop), total mass, total length, shoot mass ratio, shoot length ratio, shoot elongation rate(pop), root elongation rate(pop)

# pull variables for pca

seedling_traits <- data[,c(11:37)]
# run PCA
seedling_pca <- prcomp(seedling_traits)
#summary(seedling_pca)
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
      germ_linear[i,t] <- beta1[species[i],t] + beta2[source[i],t] + beta3[wp[i],t]+ beta4[temp[i],t] + beta5*pc1[i,t] + beta6*pc2[i] + beta7*pc3[i] + beta8*pc4[i] + beta9[cup[i]] + beta10[chamber[i]] + beta11[rep[i]] + beta12*GermDays[i] + beta13*((t-1)*2)
    
      # transformation on logit scale for logistic regression
      prob_germ[i,t] <- exp(germ_linear[i,t]) / (1+exp(germ_linear[i,t]))
    
    }#t
  }#i
  
  ## Likelihood parameter depends on previous observation
  for (i in 1:NGerm){ 
    for (t in 2:NTimes){
      transition[i,t] ~ dbern(prob_germ[i,t] * transition[i,t-1]) 
    }#t
  }#i

  ##############################
  ### PRIORS ON MAIN EFFECTS ###
  ##############################
  ## species effect - categorical
  for (1 in 1:NTimes){
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
    kappa2 ~ dnorm(mu2, sigma2)
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
    beta3[2,t] <- dnorm(kappa3, tau3)
  }#t

  # hyperpriors on wp2 mean effect
  kappa3 ~ dnorm(0, 0.01)
  tau3 ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  ## temperature effect - categorical
  for (t in 1:NTimes){
    for (tmp in c(1,3)){ 
      beta4[tmp,t] ~ dnorm(kappa4[tmp], tau4[tmp]) 
    }#tmp
    beta4[2] <- -1 * sum(beta4[c(1,3),t]) #sum-to-zero constraint
  }#t

  # conditional prior on temp means
  for (tmp in c(1,3)){
    kappa4[tmp] ~ dnorm(mu3, tau3)
  }#tmp
  
  # independent priors on temp variances
  for (tmp in c(1,3)){
    tau4[tmp] ~ dt(0, pow(2.5, -2), 1)T(0,)
  }#tmp

  # hyperpriors for overall temp mean
  mu3 ~ dnorm(0, 0.01)
  tau3 ~ dgamma(0.01, 0.01)

  ## prior on PC1 - continuous fixed effect
  beta5 ~ dnorm(0, 0.01)

  ## prior on PC2 - continuous fixed effect
  beta6 ~ dnorm(0, 0.01)

  ## prior on PC3 - continuous fixed effect
  beta7 ~ dnorm(0, 0.01)

  ## prior on PC4 - continuous fixed effect
  beta8 ~ dnorm(0, 0.01)

  ## cup effect - hierarchical
  for (c in 1:(NCups-1)){ 
    beta9[cup] ~ dnorm(0, nu[1]) 
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
  }
  beta11[3] <- -1 * sum(beta11[1:2]) #sum-to-zero constraint

  ## germination days effect - continuous fixed
  beta12 ~ dnorm(0, 0.01)

  ## days to seedling effect 
  beta13 ~ dnorm(0, 0.01)

  ################################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIABILITY ###
  ################################################
  nu[1] ~ dgamma(0.01, 0.01)
  nu[2] ~ dgamma(0.01, 0.01)
  nu[3] ~ dgamma(0.01, 0.01)

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
                           
}")


#################
### RUN MODEL ###
#################
#### JAGS FORMAT ####
# save data as appropriate list
seedling7Data <- list(germ=germination_matrix[,5:10], 
                     species=as.integer(as.factor(germination_matrix$Species)),
                     source=as.integer(as.factor(germination_matrix$Source)), 
                     #seedmass=germination_matrix$SeedMass, 
                     #sct=germination_matrix$SeedMass, 
                     wp=germination_matrix$WP, 
                     temp=germination_matrix$Temp,
                     cup=as.integer(as.factor(germination_matrix$CupNo)), 
                     chamber=as.integer(as.factor(germination_matrix$Chamber)),
                     NSown=nrow(germination_matrix),
                     NSpecies = length(unique(germination_matrix$Species)), 
                     NSources = length(unique(germination_matrix$Source)),
                     NCups = length(unique(germination_matrix$CupNo)),
                     NTimes = 5,
                     NChambers = length(unique(germination_matrix$Chamber)))
                    # use if in long format:
                    #days=as.numeric(as.factor(germination_matrix$)))

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

#### View estimates ####
## summarize outputs
summary <- summary(out_model)
## compare single species-timestep means to modeled species-timestep effects...
chain1 <- out_model[[1]]
chain2 <- out_model[[2]]
chain3 <- out_model[[3]]

