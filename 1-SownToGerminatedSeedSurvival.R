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
      sown_linear[i,t] <- beta1[species[i],t]+ beta2[source[i],t]+ beta3[wp[i],t] + beta4[temp[i],t] + beta5*seedmass[i] + beta6[cup[i]] + beta7[chamber[i]] + beta8[rep[i]] + beta9*((t*2)-1)

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
  for (chm in 1:(NChambers-1)){ 
    beta7[chm] ~ dnorm(0, pow(nu[2],2)) 
  }#chm
  beta7[NChambers] <- -1 * sum(beta7[1:(NChambers-1)]) # sum-to-zero constraint

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
  nu[2] ~ dt(0, pow(2.5,-2), 1)T(0,)
  nu[3] ~ dt(0, pow(2.5,-2), 1)T(0,)

}")#end model block

##################
### RUN MODELS ###
##################
# read in data- long format
#data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/DummyGerminationData_Trial3_23Apr20.csv',header=T,stringsAsFactors = FALSE)
# read in data- matrix format
germination_matrix <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/GerminationMatrix_FINAL_15DEC2020.csv', header=TRUE, stringsAsFactors=FALSE)
head(germination_matrix)

#### NIMBLE FORMAT ####
{## GOOD WEBSITE FOR NIMBLE TUTORIAL: https://r-nimble.org/html_manual/cha-lightning-intro.html

# set up data
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

# create model
germModel <- nimbleModel(code = SownSurv,
                          name = "Survival to Germination",
                          constants = germData, # just like your usual jags.data
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

}#group Nimble format

#### JAGS FORMAT ####
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
                     NTimes = 37,
                     NChambers = length(unique(germination_matrix$Chamber)))

## Initialize model
germModel <- jags.model(file=SownSurv, data=germData, n.chains=3, n.adapt=1000)
## Run model
update(germModel, n.iter=1000)
out_model <- coda.samples(model=germModel, variable.names = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "kappa1", "tau1", "kappa2", "tau2", "kappa3", "tau3", "kappa4", "tau4", "mu1", "mu2", "sigma1", "sigma2", "nu"), thin=10, n.iter=10000)

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
plot(out_model) # plot chains
#dev.off() # close PDF
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

#### CHECKING COEFFICIENTS THAT ARE NOT CONVERGING ####
## Checking beta0 and beta9 against each other 
## NOTE: Subscripts for sums may need to be changed depending on model parameterization
# pull out coda samples from each chain
chain1 <- out_model[[1]]
chain2 <- out_model[[2]]
chain3 <- out_model[[3]]
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
summary <- summary(out_model)
summary_quantiles <- as.data.frame(summary$quantiles)
rownames(summary_quantiles) <- rownames(summary$quantiles)
chains_merged <- rbind(out_model[[1]], out_model[[2]], out_model[[3]])

#### set indices for betas based on # timesteps ####
## NOTE: These will need to be updated if betas are added or removed!
timesteps = germData$NTimes
beta1_end_index = timesteps * germData$NSpecies
beta2_end_index = beta1_end_index + (timesteps*germData$NSources)
beta3_end_index = beta2_end_index + (timesteps*2)
beta4_end_index = beta3_end_index + (timesteps*3)
beta5_index = beta4_end_index + 1
beta6_end_index = beta5_index + germData$NCups
beta7_end_index = beta6_end_index + germData$NChambers
beta8_end_index = beta7_end_index + 3
beta9_index = beta8_end_index + 1
kappa1_end_index = beta9_index + germData$NSpecies
kappa2_end_index = kappa1_end_index + germData$NSources
kappa3_index = kappa2_end_index + 1
kappa4_index = kappa3_index + 1
mu1_index = kappa4_end_index + 1
mu2_index = mu1_index + 1
nu1_index = mu2_index + 1
nu2_index = nu1_index + 1
nu3_index = nu2_index + 1
sigma1_index = nu3_index + 1
sigma2_index = sigma1_index + 1
tau1_end_index = sigma2_index + germData$NSpecies
tau2_end_index = tau1_end_index + germData$NSources
tau3_index = tau2_end_index + 1
tau4_index = tau3_index + 2

#### mu estimates ####
# boxplot
ggplot(mapping=aes(x="mu1",
                   ymin=summary_quantiles[mu1_index,1],
                   lower=summary_quantiles[mu1_index,2],
                   middle=summary_quantiles[mu1_index,3],
                   upper=summary_quantiles[mu1_index,4],
                   ymax=summary_quantiles[mu1_index,5]))+
  geom_boxplot(stat="identity")+
  ggtitle("Overall species effect")+
  theme_classic()+
  xlab("Mu1")
# violin plots
vioplot(summary_quantiles[mu1_index], col = "gray", names=levels(as.factor(germination_matrix$Species)))

#### Kappa1 estimates ####
# species means (kappa1[spp]) boxplot
ggplot(mapping=aes(x=levels(as.factor(germination_matrix$Species)),
                   ymin=summary_quantiles[(beta9_index+1):kappa1_end_index,1], # 2.5% percentile
                   lower=summary_quantiles[(beta9_index+1):kappa1_end_index,2], # 25% percentile
                   middle=summary_quantiles[(beta9_index+1):kappa1_end_index,3], # median
                   upper=summary_quantiles[(beta9_index+1):kappa1_end_index,4], # 75% percentile
                   ymax=summary_quantiles[(beta9_index+1):kappa1_end_index,5]))+ # 97.5% percentile
  geom_boxplot(stat="identity")+
  ggtitle("Kappa1 estimates: Mean species effects")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Species")+
  geom_abline(aes(slope=0, intercept=0), col="red")
# species means (kappa1[spp]) violin plots
kappa1_1 <- chains_merged[,(beta9_index+1)]
kappa1_2 <- chains_merged[,(beta9_index+2)]
kappa1_3 <- chains_merged[,(beta9_index+3)]
kappa1_4 <- chains_merged[,(beta9_index+4)]
kappa1_5 <- chains_merged[,(beta9_index+5)]
kappa1_6 <- chains_merged[,(beta9_index+6)]
kappa1_7 <- chains_merged[,(beta9_index+7)]
#vioplot(chains_merged[,(beta9_index+1:beta9_index+7)], col="grey", names=levels(as.factor(germination_matrix$Species)))
vioplot(kappa1_1, kappa1_2, kappa1_3, kappa1_4, kappa1_5, kappa1_6, kappa1_7, col = "gray", names=levels(as.factor(germination_matrix$Species)))
mtext("Species", 1, cex = 1.4, line = 2.5)
mtext("Species mean", 2, cex = 1.4, line = 2.)
abline(h = 0, lty = 2)

#### kappa2 estimates ####
# source means (kappa2[src])
## NOTE: missing kappa2[NSources] due to sum-to-zero constraint
ggplot(mapping=aes(x=levels(as.factor(germination_matrix$Source))[1:germData$NSources-1],
                   ymin=summary_quantiles[(kappa1_end_index+1):kappa2_end_index, 1], # 2.5% percentile
                   lower=summary_quantiles[(kappa1_end_index+1):kappa2_end_index, 2], # 25% percentile
                   middle=summary_quantiles[(kappa1_end_index+1):kappa2_end_index, 3], # median
                   upper=summary_quantiles[(kappa1_end_index+1):kappa2_end_index, 4], # 75% percentile
                   ymax=summary_quantiles[(kappa1_end_index+1):kappa2_end_index, 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa2 estimates: Mean source effects")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### kappa3 estimates ####
# wp means
ggplot(mapping=aes(x="WP=1",
                   ymin=summary_quantiles[kappa3_index, 1], # 2.5% percentile
                   lower=summary_quantiles[kappa3_index, 2], # 25% percentile
                   middle=summary_quantiles[kappa3_index, 3], # median
                   upper=summary_quantiles[kappa3_index, 4], # 75% percentile
                   ymax=summary_quantiles[kappa3_index, 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa3 estimate: Mean WP effect for wp=1")+
    theme_classic()+
    #theme(axis.text.x=element_text(angle=90))+
    xlab("Water potential")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### kappa4 estimates ####
# temp means
ggplot(mapping=aes(x=c("28-10", "32-15"),
                   ymin=summary_quantiles[kappa4_index:(kappa4_index+1), 1], # 2.5% percentile
                   lower=summary_quantiles[kappa4_index:(kappa4_index+1), 2], # 25% percentile
                   middle=summary_quantiles[kappa4_index:(kappa4_index+1), 3], # median
                   upper=summary_quantiles[kappa4_index:(kappa4_index+1), 4], # 75% percentile
                   ymax=summary_quantiles[kappa4_index:(kappa4_index+1), 5]))+ # 97.5% percentile
    geom_boxplot(stat="identity")+
    ggtitle("Kappa4 estimates: Mean temp effects")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")+
    geom_abline(aes(slope=0, intercept=0), col="red")

#### beta1 estimates: Species effects ####
# BOMA over time
timesteps=germData$NTimes # number of timesteps
spp_level=1 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=2 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=3 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=4 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=5 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=6 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
timesteps=germData$NTimes # number of timesteps
spp_level=7 # integer for species level
# calculate indices for variables 
total_steps <- 5*7
index = seq(spp_level, total_steps, 7)
ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
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
for (source_level in 1:germData$NSources){
  timesteps=germData$NTimes
  start = beta1_end_index + source_level
  index = seq(start, beta1_end_index, 38)
  source_name = levels(as.factor(germination_matrix$Source))[source_level]
  ggplot(mapping=aes(x=rownames(summary_quantiles)[index],
                     ymin=summary_quantiles[index,1],
                     lower=summary_quantiles[index,2],
                     middle=summary_quantiles[index,3],
                     upper=summary_quantiles[index,4],
                     ymax=summary_quantiles[index,5]))+ 
    geom_boxplot(stat="identity")+
    ggtitle(paste("Beta2 estimates: Source effects for", source_name))+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")
}#for

for (t in 1:germData$NTimes){
  start = (beta1_end_index+1) + 38*(t-1)
  end = start + (germData$NTimes-1)
  ggplot(mapping=aes(x=rownames(summary_quantiles)[start:end],
                     ymin=summary_quantiles[start:end,1],
                     lower=summary_quantiles[start:end,2],
                     middle=summary_quantiles[start:end,3],
                     upper=summary_quantiles[start:end,4],
                     ymax=summary_quantiles[start:end,5]))+
    geom_bosplot(stat="identity")+
    ggtitle(paste("Beta2 estimates: Source effects per timestep"))+
    theme(axis.text.x=element_text(angle=90))+
    xlab("Source")
}#for

#### beta3 estimates: WP ####
ggplot(mapping=aes(x=rownames(summary_quantiles)[(beta2_end_index+1):beta3_end_index],
                   ymin=summary_quantiles[(beta2_end_index+1):beta3_end_index,1],
                   lower=summary_quantiles[(beta2_end_index+1):beta3_end_index,2],
                   middle=summary_quantiles[(beta2_end_index+1):beta3_end_index,3],
                   upper=summary_quantiles[(beta2_end_index+1):beta3_end_index,4],
                   ymax=summary_quantiles[(beta2_end_index+1):beta3_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta3 estimates: Water potential effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Water Potential")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta4 estimate: temp effect ####
ggplot(mapping=aes(x=rownames(summary_quantiles)[(beta3_end_index+1):beta4_end_index],
                   ymin=summary_quantiles[(beta3_end_index+1):beta4_end_index,1],
                   lower=summary_quantiles[(beta3_end_index+1):beta4_end_index,2],
                   middle=summary_quantiles[(beta3_end_index+1):beta4_end_index,3],
                   upper=summary_quantiles[(beta3_end_index+1):beta4_end_index,4],
                   ymax=summary_quantiles[(beta3_end_index+1):beta4_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta3 estimates: Water potential effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Temp effects")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta5 estimate: seed mass ####
ggplot(mapping=aes(x=rownames(summary_quantiles)[beta5_index],
                   ymin=summary_quantiles[beta5_index,1],
                   lower=summary_quantiles[beta5_index,2],
                   middle=summary_quantiles[beta5_index,3],
                   upper=summary_quantiles[beta5_index,4],
                   ymax=summary_quantiles[beta5_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta5 estimates: Seed mass effect")+
  xlab("Beta5")

#### beta6 estimates: cup effect #####
ggplot(mapping=aes(x=rownames(summary$quantiles)[(beta5_index+1):beta6_end_index],
                   ymin=summary$quantiles[(beta5_index+1):beta6_end_index,1],
                   lower=summary$quantiles[(beta5_index+1):beta6_end_index,2],
                   middle=summary$quantiles[(beta5_index+1):beta6_end_index,3],
                   upper=summary$quantiles[(beta5_index+1):beta6_end_index,4],
                   ymax=summary$quantiles[(beta5_index+1):beta6_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  geom_abline(aes(slope=0, intercept=0), col="red")
  ggtitle("Beta6 estimates: Cup effect")+
  xlab("Beta6")

#### beta7 estimates: chamber effect #####
ggplot(mapping=aes(x=rownames(summary_quantiles)[(beta6_end_index+1):beta7_end_index],
                   ymin=summary_quantiles[(beta6_end_index+1):beta7_end_index,1],
                   lower=summary_quantiles[(beta6_end_index+1):beta7_end_index,2],
                   middle=summary_quantiles[(beta6_end_index+1):beta7_end_index,3],
                   upper=summary_quantiles[(beta6_end_index+1):beta7_end_index,4],
                   ymax=summary_quantiles[(beta6_end_index+1):beta7_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta7 estimates: Chamber effects")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Chamber")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### beta8 estimates: rep effect ####
ggplot(mapping=aes(x=rownames(summary$quantiles)[(beta7_end_index-1):beta8_end_index],
                   ymin=summary$quantiles[(beta7_end_index-1):beta8_end_index:1],
                   lower=summary$quantiles[(beta7_end_index-1):beta8_end_index:,2],
                   middle=summary$quantiles[(beta7_end_index-1):beta8_end_index:,3],
                   upper=summary$quantiles[(beta7_end_index-1):beta8_end_index:,4],
                   ymax=summary$quantiles[(beta7_end_index-1):beta8_end_index:,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta8 estimates: Chamber effects")+
  xlab("Chamber * rep number")+
  theme_classic()

#### beta9 estimates ####
ggplot(mapping=aes(x=rownames(summary$quantiles)[beta9_index],
                   ymin=summary$quantiles[beta9_index,1],
                   lower=summary$quantiles[beta9_index,2],
                   middle=summary$quantiles[beta9_index,3],
                   upper=summary$quantiles[beta9_index,4],
                   ymax=summary$quantiles[beta9_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta9 estimate: Day effect")+
  theme_classic()+
  xlab("Parameter")+
  geom_abline(aes(slope=0, intercept=0), col="red")

#### Variance estimates #####
# nu1 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[nu1_index],
                   ymin=summary$quantiles[nu1_index,1],
                   lower=summary$quantiles[nu1_index,2],
                   middle=summary$quantiles[nu1_index,3],
                   upper=summary$quantiles[nu1_index,4],
                   ymax=summary$quantiles[nu1_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu1: variation across cups")+
  xlab("Nu1")

# nu2 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[nu2_index],
                   ymin=summary$quantiles[nu2_index,1],
                   lower=summary$quantiles[nu2_index,2],
                   middle=summary$quantiles[nu2_index,3],
                   upper=summary$quantiles[nu2_index,4],
                   ymax=summary$quantiles[nu2_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu2: variation across chambers")+
  xlab("Nu2")

# nu3 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[nu3_index],
                   ymin=summary$quantiles[nu3_index,1],
                   lower=summary$quantiles[nu3_index,2],
                   middle=summary$quantiles[nu3_index,3],
                   upper=summary$quantiles[nu3_index,4],
                   ymax=summary$quantiles[nu3_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Nu3: variation across reps")+
  xlab("Nu3")

# sigma1 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[sigma1_index],
                   ymin=summary$quantiles[sigma1_index,1],
                   lower=summary$quantiles[sigma1_index,2],
                   middle=summary$quantiles[sigma1_index,3],
                   upper=summary$quantiles[sigma1_index,4],
                   ymax=summary$quantiles[sigma1_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Sigma1: variation across overall species estimates")+
  xlab("Sigma1")

# sigma2 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[sigma2_index],
                   ymin=summary$quantiles[sigma2_index,1],
                   lower=summary$quantiles[sigma2_index,2],
                   middle=summary$quantiles[sigma2_index,3],
                   upper=summary$quantiles[sigma2_index,4],
                   ymax=summary$quantiles[sigma2_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Sigma2: variation across overall source estimates")+
  xlab("Sigma2")

# tau1 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[(sigma2_index+1):tau1_end_index],
                   ymin=summary$quantiles[(sigma2_index+1):tau1_end_index,1],
                   lower=summary$quantiles[(sigma2_index+1):tau1_end_index,2],
                   middle=summary$quantiles[(sigma2_index+1):tau1_end_index,3],
                   upper=summary$quantiles[(sigma2_index+1):tau1_end_index,4],
                   ymax=summary$quantiles[(sigma2_index+1):tau1_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Tau1: variation in species estimates over time")+
  xlab("Tau1")

# tau2 estimates
ggplot(mapping=aes(x=rownames(summary$quantiles)[(tau1_end_index+1):tau2_end_index],
                   ymin=summary$quantiles[(tau1_end_index+1):tau2_end_index,1],
                   lower=summary$quantiles[(tau1_end_index+1):tau2_end_index,2],
                   middle=summary$quantiles[(tau1_end_index+1):tau2_end_index,3],
                   upper=summary$quantiles[(tau1_end_index+1):tau2_end_index,4],
                   ymax=summary$quantiles[(tau1_end_index+1):tau2_end_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Tau2: variation in source estimates over time")+
  xlab("Tau2")

# tau3 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[tau3_index],
                   ymin=summary$quantiles[tau3_index,1],
                   lower=summary$quantiles[tau3_index,2],
                   middle=summary$quantiles[tau3_index,3],
                   upper=summary$quantiles[tau3_index,4],
                   ymax=summary$quantiles[tau3_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Tau3: variation in WP2 estimates over time")+
  xlab("Tau3")

# tau4 estimate
ggplot(mapping=aes(x=rownames(summary$quantiles)[(tau3_index+1):tau4_index],
                   ymin=summary$quantiles[(tau3_index+1):tau4_index,1],
                   lower=summary$quantiles[(tau3_index+1):tau4_index,2],
                   middle=summary$quantiles[(tau3_index+1):tau4_index,3],
                   upper=summary$quantiles[(tau3_index+1):tau4_index,4],
                   ymax=summary$quantiles[(tau3_index+1):tau4_index,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Tau4: variation in Temp estimates over time")+
  xlab("Tau4")

########################
### CROSS-VALIDATION ###
########################
# CODE FROM KENEN- NEEDS TO BE ADJUSTED!!

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
#############################
### SIMULATING POPULATION ###
#############################
#### Individual level ####
## set input values based on factor levels
species <- 6 # options: BOMA, DISP, ELPA, JUBA, PHAU, SCAC, SCAM
source <- 20 # options restricted by species
source_name <- "KIWA"
wp <- 1
temp <- 3
NSown <- 1000
NTimes <- germData$NTimes
seedmass_csv <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedMass_FINAL_Processed.csv', header=TRUE)
names(seedmass_csv)[1] <- "Source"

## pull matrix of coefficients from model
chains_merged <- rbind(out_model[[1]], out_model[[2]], out_model[[3]])
sampled_rows <- sample(1:nrow(chains_merged), NSown, replace=TRUE)
coeff_matrix <- chains_merged[sampled_rows,]

## set indices for betas based on # timesteps
## NOTE: These will need to be updated if betas are added or removed!
timesteps = NTimes
beta1_end_index = timesteps * germData$NSpecies
beta2_end_index = beta1_end_index + (timesteps*germData$NSources)
beta3_end_index = beta2_end_index + (timesteps*2)
beta4_end_index = beta3_end_index + (timesteps*3)
beta5_index = beta4_end_index + 1
beta6_end_index = beta5_index + germData$NCups
beta7_end_index = beta6_end_index + germData$NChambers
beta8_end_index = beta7_end_index + 3
beta9_index = beta8_end_index + 1

## find parameter index for input parameters
## NOTE: Check these using rownames(summary_quantiles)[beta1_cols]
beta1_cols <- seq(species, beta1_end_index, germData$NSpecies) # 238 for 34 timesteps
beta2_cols <- seq(beta1_end_index+source, beta2_end_index, germData$NSources) # 1530 for 34 samples 
beta3_cols <- seq(beta2_end_index+wp, beta3_end_index, 2) #1598 for 34 timesteps
beta4_cols <- seq(beta3_end_index+temp, beta4_end_index, 3) #1700 for 34 timesteps
seedmass <- mean(seedmass_csv$AvgSeedMass[which(seedmass_csv$Site==source_name)])
cups <- sample(1:germData$NCups, NSown, replace=TRUE)
chambers <- sample(1:germData$NChambers, NSown, replace=TRUE)
reps <- sample(1:3, NSown, replace=TRUE)

# set up matrix to hold transition probabilities
prob_sown <- matrix(NA, nrow=NSown, ncol=germData$NTimes)

## simulate by timesteps for each individual
## NOTE: Check indices by using "head(coeff_matrix[,beta1_cols], 1)"
for (i in 1:NSown){
  for (t in 1:germData$NTimes){
    beta1 <- coeff_matrix[i,beta1_cols[t]]
    beta2 <- coeff_matrix[i,beta2_cols[t]]
    beta3 <- coeff_matrix[i,beta3_cols[t]]
    beta4 <- coeff_matrix[i,beta4_cols[t]]
    beta5 <- coeff_matrix[i,beta5_index]*seedmass
    beta6 <- coeff_matrix[i,beta5_index+cups[i]]
    beta7 <- coeff_matrix[i,beta6_end_index+chambers[i]]
    beta8 <- coeff_matrix[i,beta7_end_index+reps[i]]
    beta9 <- coeff_matrix[i,beta9_index]*((t*2)-1)
    linear <- beta1 + beta2 + beta3 + beta4 + beta5 + beta6 + beta7 + beta8 + beta9
    prob_sown[i,t] <- exp(linear) / (1+exp(linear))
  }
}

## multiply transition probabilities through time
germ <- matrix(NA, nrow=NSown, ncol=germData$NTimes)
for (i in 1:NSown){
  germ[i,1] <- rbinom(n=1, size=1, prob=prob_sown[i,1])
  for (t in 2:germData$NTimes){
    ifelse (germ[i,t-1]==0, #statement to test
      germ[i,t] <- 0,  #if true
      germ[i,t] <- rbinom(n=1, size=1, prob=prob_sown[i,t]))  #if false
  }
}

## check mean predictions against mean transitions for data
# calculate predicted mean transition rate per timestep (p)
## NOTE: remember that the real-life transition rate is (1-p)
plot_title <- paste("Species:",species,"Source:",source_name,"WP:",wp,"Temp:",temp)
predicted_means <- c()
for (t in 1:germData$NTimes){ 
  predicted_means[t] <- mean(prob_sown[,t], na.rm=TRUE)
}
# plot predicted means
plot(predicted_means, col="blue", ylim=c(0,1), ylab="1 - Transition rate", xlab="Timestep", main=plot_title)
lines(predicted_means, col="blue")
# calculate means in real data
subset_df <- germination_matrix[which(germination_matrix$Source==source_name & germination_matrix$Temp==temp & germination_matrix$WP==wp),]
real_means <- c()
ones_count <- c()
for (t in 1:germData$NTimes){
  ones_count[t] <- length(which(subset_df[,t]==1))
  ifelse (t==1, # condition
      num_inds_remaining <- nrow(subset_df), # if true
      num_inds_remaining <- ones_count[t-1]) # if false
  real_means[t] <- ones_count[t] / num_inds_remaining
}#t
# add reference lines
abline(h=0, col="darkgrey")
abline(h=1, col="darkgrey")
# plot real data
points(real_means, col="red", lty="dotted")
lines(real_means, col="red", lty="dotted")

#### Population level ####
# set input values based on factor levels
species <- 1 # options: BOMA, DISP, ELPA, JUBA, PHAU, SCAC, SCAM
source <- 10 # options restricted by species
source_name <- "DIST"
wp <- 1
temp <- 2
NSown <- 1000
NTimes <- germData$NTimes # number of timesteps
seedmass_csv <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedMass_FINAL_Processed.csv', header=TRUE)
names(seedmass_csv)[1] <- "Source"
seedmass <- mean(seedmass_csv$AvgSeedMass[which(seedmass_csv$Site==source_name)])

## set indices for betas based on # timesteps
## NOTE: These will need to be updated if betas are added or removed!
timesteps = NTimes
beta1_end_index = timesteps * germData$NSpecies
beta2_end_index = beta1_end_index + (timesteps*germData$NSources)
beta3_end_index = beta2_end_index + (timesteps*2)
beta4_end_index = beta3_end_index + (timesteps*3)
beta5_index = beta4_end_index + 1
beta6_end_index = beta5_index + germData$NCups
beta7_end_index = beta6_end_index + germData$NChambers
beta8_end_index = beta7_end_index + 3
beta9_index = beta8_end_index + 1

# pull point estimates
estimates <- summary(out_model)$statistics[,1]

# find parameter estimates for input parameters
## NOTE: These can be checked by calling the variable
beta1_time <- estimates[seq(1, beta1_end_index, germData$NSpecies)]
beta2_time <- estimates[seq(beta1_end_index+1, beta2_end_index, germData$NSources)]
beta3_time <- estimates[seq(beta2_end_index+1, beta3_end_index, 2)]
beta4_time <- estimates[seq(beta3_end_index+1, beta4_end_index, 3)]
beta5 <- estimates[beta5_index]
# betas 6-8 should have an average value of zero, so left out
beta9 <- estimates[beta9_index]

# set up matrix to hold transition probabilities
prob_sown <- matrix(NA, nrow=1, ncol=germData$NTimes)
  
# simulate by timesteps for each individual
# NOTE: Check indices by using "head(coeff_matrix[,beta1_cols], 1)"
for (t in 1:germData$NTimes){
  beta1 <- beta1_time[t]
  beta2 <- beta2_time[t]
  beta3 <- beta3_time[t]
  beta4 <- beta4_time[t]
  beta5 <- beta5 * seedmass
  beta9 <- beta9 * ((t*2)-1)
  linear <- beta1 + beta2 + beta3 + beta4 + beta5 #+ beta6 + beta7 + beta8 + beta9
  prob_sown[t] <- exp(linear) / (1+exp(linear))
}

# multiply transition probabilities through time
germ <- matrix(NA, nrow=1, ncol=germData$NTimes)
germ[1] <- rbinom(n=1, size=NSown, prob=prob_sown[1])
for (t in 2:germData$NTimes){
  germ[1,t] <- rbinom(n=1, size=germ[t-1], prob=prob_sown[t])
}
