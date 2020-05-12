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

####################################################
### Transition from sown seed to germinated seed ###
####################################################
SownSurv <- 
# use this line if running in nimble (remember to remove quotes at end bracket)
#nimbleCode({
# use this line instead if running in rjags
textConnection("model{ # has the same effect as writing all this in a separate text file
  ###########################
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
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
      sown_linear[i,t] <- beta1[species[i],t]+ beta2[source[i],t]+ beta3[wp[i],t] + beta4[temp[i],t] + beta5*seedmass[i] + beta6[cup[i]] + beta7[chamber[i]] + beta9*((t-1)*2)
      # + beta8[rep[i]]
    
    ## Transformation to logit scale
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
  for (t in 1:NTimes) {
    for (spp in 1:NSpecies) { 
      beta1[spp,t] ~ dnorm(kappa1[spp], tau1[spp]) 
    }#spp
  ## sum-to-zero constraint not needed because no beta0- now beta1 species effect is taking up some of the intercept term
  #beta1[NSpecies,t] <- -1 * sum(beta1[1:(NSpecies-1),t]) # sum-to-zero constraint on species axis
  }#t

  # conditional prior on species means
  for (spp in 1:NSpecies){
    kappa1[spp] ~ dnorm(mu1, eps1) # conditional prior
  }#spp

  # independent prior on species variances
  for (spp in 1:NSpecies){
    tau1[spp] ~ dgamma(0.01, 0.01) # independent prior
  }#spp
  
  # uninformative hyperpriors on species mean
  mu1 ~ dnorm(0, 0.01)
  eps1 ~ dgamma(0.01, 0.01)

  ## source effect on germination probability - hierarchical - varies by time
  for (t in 1:NTimes) {
    for (src in 1:(NSources-1)) { 
      beta2[src,t] ~ dnorm(kappa2[src], tau2[src]) 
    }#src
  beta2[NSources,t] <- -1 * sum(beta2[1:(NSources-1),t]) # sum-to-zero constraint
  }#t

  # conditional priors on source means
  for (src in 1:(NSources-1)){
    kappa2[src] ~ dnorm(mu2, eps2)
  }
    
  # independent priors on source variances
  for (src in 1:(NSources-1)){
    tau2[src] ~ dgamma(0.01, 0.01)
  }

  # uninformative hyperpriors on source mean
  mu2 ~ dnorm(0, 0.01)
  eps2 ~ dgamma(0.01, 0.01)

  ## wp effect on germination - hierarchical - varies by time
  for (t in 1:NTimes) {
    beta3[1,t] ~ dnorm(kappa3, tau3)
    beta3[2,t] <- -1 * beta3[1,t] # sum-to-zero constraint
  }#t

  # independent priors on wp mean and variance
  # NOTE: not conditional for the mean because there's only one level (i.e. WP 1) being drawn from the distribution
  kappa3 ~ dnorm(0, 0.01)
  tau3 ~ dgamma(0.01, 0.01)

  ## temp effect on germination - hierarchical - varies by time
  for (t in 1:NTimes) {
    for (tmp in 1:2){
      beta4[tmp,t] ~ dnorm(kappa4[tmp], tau4[tmp])
    }#tmp
  beta4[3,t] <- -1 * sum(beta4[1:2,t]) # sum-to-zero constraint
  }#t
  
  # independent priors on temp means
  # NOTE: Poor convergence with conditional hyperparameters kappa4 and tau4- not enough levels maybe?
  for (tmp in 1:2){
    kappa4[tmp] ~ dnorm(0, 0.01)
  }

  # independent priors on temp variances
  for (tmp in 1:2){
    tau4[tmp] ~ dgamma(0.01, 0.01)
  }

  ## hydrothermal time effect on germination - fixed linear - slightly constrained variance prior
  #beta3 ~ dnorm(0, 0.01)

  ## seed mass effect on germination - fixed linear - slightly constrained variance prior
  beta5 ~ dnorm(0, 0.01)

  ## seed coat thickness effect on germination - fixed linear - slightly constrained variance prior
  # NOTE: removed because estimate was centered on zero
  #beta6 ~ dnorm(0, 0.01)

  ## cup effect on germination - hierarchical
  for (c in 1:(NCups-1)) { 
    beta6[c] ~ dnorm(0, nu[1]) 
  }#cup
  beta6[NCups] <- -1 * sum(beta6[1:(NCups-1)]) # sum-to-zero

  ## chamber effect on germination - hierarchical
  for (chm in 1:(NChambers-1)) { 
    beta7[chm] ~ dnorm(0, nu[2]) 
  }#chm
  beta7[NChambers] <- -1 * sum(beta7[1:(NChambers-1)]) # sum-to-zero constraint

  ## rep effect on germination - hierarchical
  #for (rep in 1:2)   { 
  #  beta9[rep] ~ dnorm(0, nu[3]) 
  #}
  #beta9[3] <- sum(beta9[1:2]) # sum-to-zero constraint

  ## days effect - fixed linear - slightly constrained
  beta9 ~ dnorm(0, 0.01)

  ##############################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIATION ###
  ##############################################
  nu[1] ~ dgamma(0.01, 0.01)
  nu[2] ~ dgamma(0.01, 0.01)
  nu[3] ~ dgamma(0.01, 0.01)

}")#end model block

##################
### RUN MODELS ###
##################
# read in data- long format
#data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/DummyGerminationData_Trial3_23Apr20.csv',header=T,stringsAsFactors = FALSE)
# read in data- matrix format
germination_matrix <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixGermDataTrial3_28Apr_2020.csv', header=TRUE, stringsAsFactors=FALSE)

#### NIMBLE FORMAT ####
## GOOD WEBSITE FOR NIMBLE TUTORIAL: https://r-nimble.org/html_manual/cha-lightning-intro.html

# set up data
GermData <- list(germ=germination_matrix[,5:10],
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

# create model
germModel <- nimbleModel(code = SownSurv2,
                          name = "Survival to Germination",
                          constants = GermData, # just like your usual jags.data
                          inits = 1000) # just like your usual inits function

# check out model components
germModel$getNodeNames()
germ$plotGraph() # directed acyclic graph

# show dependencies
germModel$getDependencies(c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta8", "beta9")) # all dependencies to betas
germModel$getDependencies(c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta8", "beta9"), determOnly=TRUE) # only deterministic dependencies

# check that the lifted nodes were initialized
germModel[[]]

# generate MCMC
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

samples_coda <- mcmc.list(samples)
sum_code <- summary(samples_coda)

summary(samples_coda)
gelman.diag(samples_coda)
plot(samples_coda)



#### JAGS FORMAT ####
# save data as appropriate list
germData <- list(germ=germination_matrix[,5:10], 
                     species=as.integer(as.factor(germination_matrix$Species)),
                     source=as.integer(as.factor(germination_matrix$Source)), 
                     seedmass=germination_matrix$SeedMass, 
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
germModel <- jags.model(file=SownSurv, data=germData, n.chains = 3)#adapt=500
## Run model
update(germModel, n.iter=1000)
out_model26<- coda.samples(model=germModel, variable.names = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta9", "kappa1", "tau1", "kappa2", "tau2", "kappa3", "tau3", "mu1", "mu2", "eps1", "eps2", "nu"), thin=1, n.iter = 1000)

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
plot(out_model21, ask=TRUE) # plot chains
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


#### CHECKING COEFFICIENTS THAT ARE NOT CONVERGING ####
## Checking beta0 and beta9 against each other 
## NOTE: Subscripts for sums may need to be changed depending on model parameterization
# pull out coda samples from each chain
chain1 <- out_model16[[1]]
chain2 <- out_model16[[2]]
chain3 <- out_model16[[3]]
# take the sum of beta0 and beta9 for each chain
sum_chain1 <- chain1[,1] + chain1[,473]; traceplot(sum_chain1, col="black")
sum_chain2 <- chain2[,1] + chain2[,473]; traceplot(sum_chain2,col="red")
sum_chain3 <- chain3[,1] + chain3[,473]; plot(sum_chain3)
# return traceplots for each chain for the sum of beta0 and beta9
traceplot(sum_chain1, col="black")
traceplot(sum_chain2, col="red", add=TRUE)
traceplot(sum_chain3, col="green", add=TRUE)
# check the correlation betweeen beta0 and beta9
chain1_beta1 <- chain1[1:nrow(chain1), 1, drop=FALSE]
chain1_beta9 <- chain1[1:nrow(chain1), 473, drop=FALSE]
chain2_beta1 <- chain1[1:nrow(chain2), 1, drop=FALSE]
chain2_beta9 <- chain1[1:nrow(chain2), 473, drop=FALSE]
chain3_beta1 <- chain1[1:nrow(chain3), 1, drop=FALSE]
chain3_beta9 <- chain1[1:nrow(chain3), 473, drop=FALSE]
plot(chain1_beta1, chain1_beta9)
plot(chain2_beta1, chain2_beta9)
plot(chain3_beta1, chain3_beta9)
cor(chain1_beta1, chain1_beta9) # -0.879
cor(chain2_beta1, chain2_beta9) # -0.879
cor(chain3_beta1, chain3_beta9) # -0.879
# CONCLUSION: Beta0 and beta9 are highly correlated, indicating that the intercept and time coefficient are trying to account for the same variation. Kezia recommends dropping the intercept, especially since the time component has meaning.

## Checking seedmass (beta0) and seed coat thickness (beta6)
# NOTE: Subscripts may need to be changed depending on model parameterization
# check correlation of posteriors
chain1 <- out_model[[1]]
chain2 <- out_model[[2]]
chain3 <- out_model[[3]]
chain1_seedmass <- chain1[1:nrow(chain1),253, drop=FALSE]
chain1_sct <- chain1[1:nrow(chain1),252, drop=FALSE]
plot(chain1_seedmass, chain1_sct)
cor(chain1_sct,chain1_seedmass) # 0.113
# check correlation of original values
plot(GermData$seedmass, GermData$sct)
cor(GermData$seedmass, GermData$sct) # 0.073
## NO CORRELATION! If rerun and everything else  converges, try removing seed coat thickness. Currently the posterior for sct is also centered on 0 meaning that it may not be providing any information on germination.


####################
### PLOT OUTPUTS ###
####################
library(ggplot2)
library(vioplot)
summary <- summary(out_model26)
summary_df <- as.data.frame(summary$quantiles)
rownames(summary_df) <- rownames(summary$quantiles)
chains_merged <- rbind(out_model[[1]], out_model[[2]], out_model[[3]])

#### mu estimates ####
# overall species effect (mu1)
ggplot(mapping=aes(x="mu1",
                   ymin=summary_df[521,1],
                   lower=summary_df[521,2],
                   middle=summary_df[521,3],
                   upper=summary_df[521,4],
                   ymax=summary_df[521,5]))+
  geom_boxplot(stat="identity")+
  ggtitle("Overall species effect")+
  theme_classic()+
  xlab("Mu1")

#### Kappa1 estimates ####
# species means (kappa1[spp]) boxplot
ggplot(mapping=aes(x=levels(as.factor(germination_matrix$Species)),
                   ymin=summary_df[474:480,1], # 2.5% percentile
                   lower=summary_df[474:480,2], # 25% percentile
                   middle=summary_df[474:480,3], # median
                   upper=summary_df[474:480,4], # 75% percentile
                   ymax=summary_df[474:480,5]))+ # 97.5% percentile
  geom_boxplot(stat="identity")+
  ggtitle("Kappa1 estimates: Mean species effects")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Species")+
  geom_abline(aes(slope=0, intercept=0), col="red")
# species means (kappa1[spp]) violin plots
kappa1_1 <- chains_merged[,474]
kappa1_2 <- chains_merged[,475]
kappa1_3 <- chains_merged[,476]
kappa1_4 <- chains_merged[,477]
kappa1_5 <- chains_merged[,478]
kappa1_6 <- chains_merged[,479]
kappa1_7 <- chains_merged[,480]
vioplot(kappa1_1, kappa1_2, kappa1_3, kappa1_4, kappa1_5, kappa1_6, kappa1_7, col = "gray", names=levels(as.factor(germination_matrix$Species)))
mtext("Species", 1, cex = 1.4, line = 2.5)
mtext("Species mean", 2, cex = 1.4, line = 2.)
abline(h = 0, lty = 2)

#### kappa2 estimates ####
# source means (kappa2[src])
## NOTE: missing kappa2[NSources] due to sum-to-zero constraint
ggplot(mapping=aes(x=levels(as.factor(germination_matrix$Source))[1:germData$NSources-1],
                   ymin=summary_df[481:517, 1], # 2.5% percentile
                   lower=summary_df[481:517, 2], # 25% percentile
                   middle=summary_df[481:517, 3], # median
                   upper=summary_df[481:517, 4], # 75% percentile
                   ymax=summary_df[481:517, 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa2 estimates: Mean source effects")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### kappa3 estimates ####
# wp means
ggplot(mapping=aes(x="WP=1",
                   ymin=summary_df[518, 1], # 2.5% percentile
                   lower=summary_df[518, 2], # 25% percentile
                   middle=summary_df[518, 3], # median
                   upper=summary_df[518, 4], # 75% percentile
                   ymax=summary_df[518, 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa3 estimate: Mean WP effect for wp=1")+
    theme_classic()+
    #theme(axis.text.x=element_text(angle=90))+
    xlab("Water potential")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### kappa4 estimates ####
# temp means
ggplot(mapping=aes(x=c("28-10", "32-15"),
                   ymin=summary_df[519:520, 1], # 2.5% percentile
                   lower=summary_df[519:520, 2], # 25% percentile
                   middle=summary_df[519:520, 3], # median
                   upper=summary_df[519:520, 4], # 75% percentile
                   ymax=summary_df[519:520, 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa4 estimates: Mean temp effects")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### beta1 estimates: Species effects ####
# BOMA over time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=1 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: BOMA through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# DISP over time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=2 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: DISP through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# ELPA over time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=3 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: ELPA through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# JUBA through time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=4 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: JUBA through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# PHAU through time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=5 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: PHAU through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# SCAC through time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=6 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: SCAC through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

# SCAM through time
# set these based on data....
timesteps=germData$NTimes # number of timesteps
spp_level=7 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary$quantiles[index,1],
                   lower=summary$quantiles[index,2],
                   middle=summary$quantiles[index,3],
                   upper=summary$quantiles[index,4],
                   ymax=summary$quantiles[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: SCAM through time", sep=''))+
  xlab("Parameter name")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta2 estimates: Source effects ####
# beta2 estimates
source_level=1
timesteps=germData$NTimes
start=36+source_level
end=start+(38*timesteps)
index=seq(start, end-1, 38)
ggplot(mapping=aes(x=rownames(summary_df)[index],
                   ymin=summary_df[index,1],
                   lower=summary_df[index,2],
                   middle=summary_df[index,3],
                   upper=summary_df[index,4],
                   ymax=summary_df[index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta2 estimates: Source effects at t=0")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Source")

#### beta3 estimates: WP ####
ggplot(mapping=aes(x=rownames(summary_df)[226:235],
                   ymin=summary_df[226:235,1],
                   lower=summary_df[226:235,2],
                   middle=summary_df[226:235,3],
                   upper=summary_df[226:235,4],
                   ymax=summary_df[226:235,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta3 estimates: Water potential effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Water Potential")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta4 estimate: temp effect ####
ggplot(mapping=aes(x=rownames(summary_df)[236:250],
                   ymin=summary_df[236:250,1],
                   lower=summary_df[236:250,2],
                   middle=summary_df[236:250,3],
                   upper=summary_df[236:250,4],
                   ymax=summary_df[236:250,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta3 estimates: Water potential effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Temp effects")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta5 estimate: seed mass ####
ggplot(mapping=aes(x=rownames(summary_df)[251],
                   ymin=summary_df[251,1],
                   lower=summary_df[251,2],
                   middle=summary_df[251,3],
                   upper=summary_df[251,4],
                   ymax=summary_df[251,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta5 estimates: Seed mass effect")+
  xlab("Beta5")

#### beta6 estimates: cup effect #####
ggplot(mapping=aes(x=rownames(summary$quantiles)[252:467],
                   ymin=summary$quantiles[252:467,1],
                   lower=summary$quantiles[252:467,2],
                   middle=summary$quantiles[252:467,3],
                   upper=summary$quantiles[252:467,4],
                   ymax=summary$quantiles[252:467,5]))+ 
  geom_boxplot(stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  geom_abline(aes(slope=0, intercept=0), col="red")
  ggtitle("Beta6 estimates: Cup effect")+
  xlab("Beta6")

#### beta7 estimates: chamber effect #####
ggplot(mapping=aes(x=rownames(summary_df)[468:470],
                   ymin=summary_df[468:470,1],
                   lower=summary_df[468:470,2],
                   middle=summary_df[468:470,3],
                   upper=summary_df[468:470,4],
                   ymax=summary_df[468:470,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta7 estimates: Chamber effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Chamber")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta8 estimates: rep effect ####
ggplot(mapping=aes(x=rownames(summary$quantiles)[...],
                   ymin=summary$quantiles[:,1],
                   lower=summary$quantiles[:,2],
                   middle=summary$quantiles[:,3],
                   upper=summary$quantiles[:,4],
                   ymax=summary$quantiles[:,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta8 estimates: Chamber effects")+
  xlab("Chamber * rep number")+
  theme_classic()

#### beta9 estimates ####
ggplot(mapping=aes(x=rownames(summary$quantiles)[471],
                   ymin=summary$quantiles[471,1],
                   lower=summary$quantiles[471,2],
                   middle=summary$quantiles[471,3],
                   upper=summary$quantiles[471,4],
                   ymax=summary$quantiles[471,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta9 estimate: Day effect")+
  theme_classic()+
  xlab("Parameter")+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### Variance estimates #####
# nu1 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[286],
                   ymin=summary$quantiles[286,1],
                   lower=summary$quantiles[286,2],
                   middle=summary$quantiles[286,3],
                   upper=summary$quantiles[286,4],
                   ymax=summary$quantiles[286,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu1: germ variation across species estimates")+
  xlab("Nu1")

# nu2 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[287],
                   ymin=summary$quantiles[287,1],
                   lower=summary$quantiles[287,2],
                   middle=summary$quantiles[287,3],
                   upper=summary$quantiles[287,4],
                   ymax=summary$quantiles[287,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu2: variation across source estimates")+
  xlab("Nu2")
# nu4 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[288],
                   ymin=summary$quantiles[288,1],
                   lower=summary$quantiles[288,2],
                   middle=summary$quantiles[288,3],
                   upper=summary$quantiles[288,4],
                   ymax=summary$quantiles[288,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu4: variation across temp estimates")+
  xlab("Nu4")
# nu7 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[289],
                   ymin=summary$quantiles[289,1],
                   lower=summary$quantiles[289,2],
                   middle=summary$quantiles[289,3],
                   upper=summary$quantiles[289,4],
                   ymax=summary$quantiles[289,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu7: variation across cup estimates")+
  xlab("Nu7")
# nu8 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[290],
                   ymin=summary$quantiles[290,1],
                   lower=summary$quantiles[290,2],
                   middle=summary$quantiles[290,3],
                   upper=summary$quantiles[290,4],
                   ymax=summary$quantiles[290,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu8: variation across chamber estimates")+
  xlab("Nu8")


#############################
### SIMULATING POPULATION ###
#############################
# set input values based on factor levels
species <- 2 # options: BOMA, DISP, ELPA, JUBA, PHAU, SCAC, SCAM
source <- 10 # options restricted by species
source_name <- "DIST"
wp <- 1
temp <- 2
NSown <- 100
seedmass_csv <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedMass_FINAL_Processed.csv', header=TRUE)
names(seedmass_csv)[1] <- "Source"

# pull matrix of coefficients
chains_merged <- rbind(out_model[[1]], out_model[[2]], out_model[[3]])
sampled_rows <- sample(1:nrow(chains_merged), NSown, replace=TRUE)
coeff_matrix <- chains_merged[sampled_rows,]

# find parameter index for input parameters
## NOTE: Check these using rownames(summary_df)[beta1_cols]
beta1_cols <- seq(species, 35, 7)
beta2_cols <- seq(35+source, 225, 38)
beta3_cols <- seq(225+wp, 235, 2)
beta4_cols <- seq(235+temp, 250, 3)
seedmass <- mean(seedmass_csv$AvgSeedMass[which(seedmass_csv$Site==source_name)])
beta6_cols <- sample(

# simulate by timesteps
for (t in 1:germData$NTimes){
  
}

########################
### CROSS-VALIDATION ###
########################
# Combine all MCMC estimates.
MCMC<-rbind(out[[1]],out[[2]],out[[3]])

# Plotting parameter estimates for Stratum and Time.
# Get just the beta4 and beta5 estimates.
st_betas<-MCMC[,grepl('beta1',colnames(MCMC)) | grepl('beta2',colnames(MCMC)) | grepl('beta3',colnames(MCMC)) | grepl('beta4',colnames(MCMC)) | grepl('beta6', colnames(MCMC))]

# Label beta4s and beta5s as their respective strata and dates.
colnames(st_betas)<-c(as.character(unique(data$Species)),as.character(unique(data$Source)), c("WP1", "WP2"), as.character(unique(data$Temp)), c("Trial1", "Trial2", "Trial3"))

# Subsample the MCMC estimates for plotting.
library(dplyr)
st_reduced<-sample_n(as.data.frame(st_betas),500) # Subsample size of 500 each.

# Convert MCMC output into a data frame structure ggplot2 can understand.
st_data<-data.frame(matrix(NA,ncol=2,nrow=nrow(st_reduced)*ncol(st_reduced))) # Empty storage df.
colnames(st_data)<-c('Species','Estimate') # Name storage df columns.
a<-0 # Set counter to 0.
for(i in 1:ncol(st_reduced)){ # Loop through all estimates.
  for(j in 1:nrow(st_reduced)){
    a<-a+1 # Increase counter count.
    st_data$Species[a]<-as.character(colnames(st_reduced)[i]) # Store stratum/date.
    st_data$Estimate[a]<-st_reduced[j,i] # Store estimate.
  }
}

# Separate stratum and time betas.
SB_data<-subset(st_data,!grepl('-',st_data$Species)) # Subset strata.
colnames(SB_data)[1]<-'Stratum' # Name parameter group.
TB_data<-subset(st_data,grepl('-',st_data$Source)) # Subset time.
colnames(TB_data)[1]<-'Date' # Name parameter group.

# Format strata df for plotting.
SB_data$Stratum<-as.factor(SB_data$Stratum) # Convert stratum to factor.
SB_data$Site<-c(rep("Gibson Lakes",500*4),rep("Ponds Lake",500*5)) # Set Site.

# Format time df for plotting.
TB_data$Date<-as.Date(TB_data$Date) # Convert to date.
TB_data<-TB_data[order(TB_data$Date),] # Order by date.
TB_data$Date<-as.factor(TB_data$Date) # Convert date to factor.
TB_data$Site<-rep(c(rep("Gibson Lakes",500),rep("Ponds Lake",500)),7) # Set Site.

# Plot MCMC subsample for stratum parameter estimates.
library(ggplot2)
ggplot(data=SB_data,aes(x=Stratum,y=Estimate,fill=Site))+
  geom_hline(yintercept=0,linetype=2)+
  geom_violin(alpha=0.8)+
  ggtitle('Stratum MCMC Subsample (n = 500 each)')+
  ylab('Parameter Estimate\nHigher = Greater Relative Condition')+
  scale_x_discrete(labels=c('5'='1','6'='2','7'='3','8'='4','9'='5'))

# Plot MCMC subsample for time parameter estimates.
ggplot(data=TB_data,aes(x=Date,y=Estimate,fill=Site))+
  geom_hline(yintercept=0,linetype=2)+
  geom_violin(alpha=0.8)+
  ggtitle('Time MCMC Subsample (n = 500 each)')+
  ylab('Parameter Estimate\nHigher = Greater Relative Condition')

# Create function for subsetting 10-fold CV data such that no levels of
# categorical predictors are lost.
protected_CV_subset<-function(CV_data,colnum){
  good<-c() # Create empty storage vector.
  while(length(good)==0 | 'Discard' %in% good){ # Repeat until all subsets are good.
    good<-c() # Clear empty storage vector.
    # Label 10-fold CV subsets.
    CV_data$TenfoldID<-sample(rep(1:10,length=nrow(CV_data)),size=nrow(CV_data),replace=F)
    for(i in 1:10){ # Loop through all 10-fold CV subsets.
      sub<-subset(CV_data,TenfoldID!=i) # Subset by 10-fold CV ID.
      for(j in colnum){ # Loop through specified columns.
        # Ensure that the same number of levels are present in the original and subsetted data.
        good<-c(good,ifelse(length(unique(sub[,j]))==length(unique(CV_data[,j])),'Good','Discard'))
      }
    }
  }
  return(CV_data$TenfoldID) # Return good subset label numbers.
}

# Apply protected 10-fold CV subset function to data.
full$TenfoldID<-protected_CV_subset(CV_data=full,colnum=c(2,3,11,13,14))
# CV_data: data to cross validate on.
# colnum: column numbers for which at least one observation is desired for each level.
# Returns vector of tenfold CV IDs.

# Define CV parameter estimation function.
CV_parameter_estimates<-function(CV_data,number_parameters,variable_names,
                                 data_statement,model_structure){
  gelman_gt_1.05<-c() # Create empty storage vector for Gelman Diagnostic > 1.05.
  gelman_gt_1.10<-c() # Create empty storage vector for Gelman Diagnostic > 1.10.
  point_estimates_df_whole<-as.data.frame(matrix(ncol=number_parameters))
  data_statement<-gsub(CV_data,'sub',data_statement) # Change data statement to work with subsets.
  CV_data<-eval(parse(text=CV_data)) # Get data frame by name.
  for(i in 1:10){
    sub<-subset(CV_data,TenfoldID!=i) # Subset by 10-fold CV ID.
    mod_data<-eval(parse(text=data_statement)) # Get model data.
    mod_structure_eval<-eval(parse(text=model_structure))
    model<-jags.model(mod_structure_eval,data=mod_data,n.chains=3) # Compile model.
    update(model,n.iter=1000) # Initiate model.
    CV_Out<-coda.samples(model=model,variable.names=variable_names,n.iter=15000,thin=3) # Run model.
    # Get count of parameters with Gelman Diagnostic greater than 1.05.
    gelman_gt_1.05<-c(gelman_gt_1.05,sum(gelman.diag(CV_Out,multivariate=F)$psrf[,2]>1.05))
    # Get count of parameters with Gelman Diagnostic greater than 1.10.
    gelman_gt_1.10<-c(gelman_gt_1.10,sum(gelman.diag(CV_Out,multivariate=F)$psrf[,2]>1.10))
    MCMC_Out<-rbind(CV_Out[[1]],CV_Out[[2]],CV_Out[[3]]) # Combine MCMC chains.
    point_estimates<-c() # Create empty storage vector.
    for(j in 1:ncol(MCMC_Out)){ # Loop through parameters.
      point_estimates<-c(point_estimates,median(MCMC_Out[,j])) # Store point estimates.
    }
    point_estimates_df<-t(data.frame(point_estimates)) # Turn point estimates into a dataframe.
    # Store point estimates.
    point_estimates_df_whole<-rbind(point_estimates_df_whole,point_estimates_df)
  }
  
  point_estimates_df_whole<-point_estimates_df_whole[-1,] # Remove top NA row.
  colnames(point_estimates_df_whole)<-colnames(MCMC_Out) # Name columns by parameter names.
  row.names(point_estimates_df_whole)<-1:nrow(point_estimates_df_whole) # Rename rows.
  point_estimates_df_whole$TenfoldIDExcluded<-1:nrow(point_estimates_df_whole) # Give a 10-fold CV ID.
  # Rearrange columns within dataframe. Place the 10-fold CV ID in first column.
  point_estimates_df_whole<-point_estimates_df_whole[,c(ncol(point_estimates_df_whole),
                                                        1:(ncol(point_estimates_df_whole)-1))]
  
  print('##############################')
  print('Gelman Convergence Diagnostic:')
  for(i in 1:10){
    print('--------')
    print(paste0('Fold ',i,':'))
    print('--------')
    print(paste0('>1.05 = ',gelman_gt_1.05[i]))
    print(paste0('>1.10 = ',gelman_gt_1.10[i]))
  }
  
  return(point_estimates_df_whole) # Return point estimates data frame.
  
}

# Run CV parameter estimation function.
Kn_CV_parameter_estimates<-CV_parameter_estimates(
  CV_data='full',
  number_parameters=30,
  variable_names=c("beta0","beta1","beta2","beta3","beta4","beta5","tau","nu"),
  data_statement=dataKn,
  model_structure=modKn)
# CV_data: name as text of data frame to cross validate on. Must have a 'TenfoldID' field.
# number_parameters: number of parameters being estimated.
# variable_names: names of parameters being estimated. For beta4[1], beta4[2], etc., just put beta4.
# data_statement: JAGS data list code stored as text.
# model_structure: JAGS model structure stored as text, not as a text connection.
# Returns data frame of parameter point estimates for each tenfold CV ID.

### The line above performed 10-fold cross validation. The number of parameters whose Gelman
### Convergence Diagnostic is greater than 1.05 and 1.10 is reported for each fold.

# Prepare for prediction.
dataKn<-eval(parse(text=dataKn)) # Execute dataKn. Turns into a list.
dataKn$TenfoldID<-full$TenfoldID # Copy tenfold ID to data list.

# Generate predicted values.
CV_predict<-function(fixed_beta_numbers,hierarchical_beta_numbers,data,CV_par_ests){
  pred_val_vec<-c() # Create empty storage vector for predicted values.
  for(i in 1:data$N){ # Loop through all observations in data list.
    # Subset full parameter data frame to just those parameters relevant to the observation.
    parameter_sub<-subset(CV_par_ests,TenfoldIDExcluded==data$TenfoldID[i])[,-1]
    for(j in fixed_beta_numbers){ # Loop through fixed betas.
      if(j==0){ # Set predicted value to beta0 if first time through loop.
        pred_val<-parameter_sub[1,colnames(parameter_sub)==paste0('beta',j)]
      } else{ # Otherwise, add next beta estimate, if applicable.
        pred_val<-pred_val+parameter_sub[1,colnames(parameter_sub)==paste0('beta',j)]*data[[j+1]][i]
      }
    }
    for(j in hierarchical_beta_numbers){ # Loop through hierarchical betas.
      # Add effect of hierchical beta.
      pred_val<-pred_val+parameter_sub[1,colnames(parameter_sub)==paste0('beta',j,'[',data[[j+1]][i],']')]
    }
    pred_val_vec<-c(pred_val_vec,pred_val) # Store predicted value.
  }
  return(pred_val_vec) # Return vector of predicted values.
}

# Run CV predict function.
Kn_pred_values<-CV_predict(fixed_beta_numbers=c(0:3),
                           hierarchical_beta_numbers=c(4:5),
                           data=dataKn,
                           CV_par_ests=Kn_CV_parameter_estimates)
# fixed_beta_numbers: vector of numbers on end of fixed betas. If beta0 and beta1, then c(0,1).
# hierarchical_beta_numbers: vector of numbers on end of hierarchical betas. 
## For beta4[1], beta4[2], etc., just put 4.
# data: data list to predict on. Not stored as text. Must have TenfoldID included.
# CV_par_ests: data frame of parameter point estimates for each tenfold CV ID.
# Returns vector of predicted values.

full$Predicted<-Kn_pred_values # Add predicted values to full data frame.
full$Actual<-log(full$Kn) # Add actual values to full data frame.

# Visualize CV results.
Kn_pred_v_actual<-lm(Actual~Predicted,data=full)
ggplot(data=full,aes(x=Predicted,y=Actual))+
  geom_abline(intercept=0,slope=1,color='black')+
  geom_point(alpha=0.5,color='blue')+
  geom_smooth(method=lm,fill=NA)+
  ggtitle('Kn Predicted vs. Actual Values')+
  labs(subtitle=paste0('r = ',round(cor(x=full$Predicted,y=full$Actual),3),"; R^2 = ",round(cor(x=full$Predicted,y=full$Actual)^2,3),"; Slope = ",round(Kn_pred_v_actual$coefficients[2],3)))

# Create residual plotting function.
residual_plot<-function(predicted,actual){
  par(mfrow=c(1,2))
  qqnorm(actual-predicted,main="Normal Quantile Plot of Residuals",col="blue",
         panel.first=qqline(actual-predicted,datax=F,distribution=qnorm))
  plot(x=predicted,y=actual-predicted,main='Residuals vs Fitted',
       xlab='Fitted Values',ylab='Residuals',bg=alpha('blue',alpha=0.5),
       pch=21,panel.first=abline(h=0,lty=2))
  lines(lowess(x=predicted,y=actual-predicted),col='red',lwd=1)
  par(mfrow=c(1,1))
}

# Plot residuals.
residual_plot(predicted=full$Predicted,actual=full$Actual)

# Done!
