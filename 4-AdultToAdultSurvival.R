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

######################################
### Transition from adult to adult ###
######################################
AdultSurv <- textConnection("model{ # has the same effect as writing all this in a separate text file
                           
  #########################################
  ### RUN PRINCIPAL COMPONENTS ANALYSES ###
  #########################################
  # NOTE: This will need to be run outside of model before incorporating so that
    # we know how many PCs to incorporate
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
                           
  ########################### 
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
  # NAdult = number of adult observations
  # block = greenhouse block
  for (i in 1:NAdult) {
    adult_linear[i] <- beta1[species[i]]+beta2[source[i]]+beta3[wp[i]]+beta4[temp[i]]+beta5*pc1[i]+beta6*pc2[i]+beta7*pc3[i]+beta8*pc4[i]+beta9[block[i]]+eps[i]
  # transformation on logit scale
  prob_adult[i] <- exp(seed_linear) / (1+exp(seed_linear))
  }
                           
  ## Likelihood parameter depends on each observation 
  for (i in 1:NAdult)    { y[i] ~ dbinom(prob_adult[i], length(y)) }
                           
  ## Uninformative prior on overall variability
  for (i in 1:NAdult)   { eps[i] ~ dgamma(0.01, 0.01) } 

  #########?  
  ## Transform to cover variable
  cover = psi*NAdult*prob 

  #############################
  ### PRIOR ON MAIN EFFECTS ###
  #############################
  # prior on species effect - categorical
  for (spp in 1:NSpecies)   { beta1[spp] ~ dnorm(0, nu[1]) }
                           
  # prior on source effect - categorical
  for (src in 1:NSources)   { beta2[src] ~ dnorm(0, nu[2]) }
                           
  # prior on water potential effect - categorical
  for (wp in 1:2)   { beta3[wp] ~ dnorm(0, nu[3]) }
                           
  # prior on temp effect - categorical
  for (tmp in 1:2)   { beta4[tmp] ~ dnorm(0, nu[4]) }
                           
  # prior on PC1 effect - continuous
  beta5 ~ dnorm(0, nu[5])
                           
  # prior on PC2 effect - continuous
  beta6 ~ dnorm(0, nu[6])
                           
  # prior on PC3 effect - continuous
  beta7 ~ dnorm(0, nu[7])
                           
  # prior on PC4 effect - continuous
  beta8 ~ dnorm(0, nu[8])
  
  # prior on greenhouse block effect                            
  for (i in 1:(NBlocks-1))   { beta9 ~ dnorm(0, nu[9]) }
  beta9[NBlocks] <- -1 * sum(beta9[1:(NBlocks-1)]) # sum-to-zero

  #######################################################
  ## model seedling traits or principal components???  ##
  #######################################################
                           
  ############################
  ### Growth rates per cup ###
  ############################
  # relative growth rate
  for (i in 1:Nrgr)  { rgr ~ dnorm(mu1, tau1) }
  for (i in 1:Nrgr) {
    mu1 <- kappa1[species[i]]+kappa2[source[i]]+kappa3[temp[i]]+kappa4[wp[i]]+kappa5[cup[i]]+kappa6[rep[i]]+kappa7[chamber[i]]
    tau1 ~ dgamma(0.01, 0.01)
  }
    for (spp in 1:NSpecies)   { kappa1[spp] ~ dnorm(0, nu[10]) }
    for (src in 1:(NSources-1))   { kappa2[src] ~ dnorm(0, nu[11]) }
    kappa2[NSources] <- -1 * sum(kappa2[1:(NSources-1)])
    for (tmp in 1:2)   { kappa3[tmp] ~ dnorm(0, nu[12]) }
    kappa3[3] <- -1 * sum(kappa3[1:2])
    kappa4[1] ~ dnorm(0, nu[13])
    kappa4[2] <- -1 * sum(kappa4[1])
    for (cup in 1:(NCups-1))   { kappa5[cup]~ dnorm(0, nu[14]) }
    kappa5[NCups] <- -1 * sum(kappa5[1:(NCups-1)])
    for (rep in 1:2)   { kappa6[rep] ~ dnorm(0, nu[15]) }
    kappa6[3] <- -1 * sum(kappa6[1:2])
    for (chm in 1:2)   { kappa7[chm] ~ dnorm(0, nu[16]) }
    kappa7[3] <- -1 * sum(kappa7[1:2])
                           
  # root elongation rate
  for (i in 1:Nrer)  { rer ~ dnorm(mu2, tau2) }
  for (i in 1:Nrer) {
    mu2 <- alpha1[species[i]]+alpha2[source[i]]+alpha3[temp[i]]+alpha4[wp[i]]+alpha5[cup[i]]+alpha6[rep[i]]+alpha7[chamber[i]]
    tau2 ~ dgamma(0.01, 0.01)
  }
    for (spp in 1:NSpecies)   { alpha1[spp] ~ dnorm(0, nu[17]) }
    for (src in 1:(NSources-1))   { alpha2[src] ~ dnorm(0, nu[18]) }
    alpha2[NSources] <- -1 * sum(alpha2[1:(NSources-1)])
    for (tmp in 1:2)   { alpha3[tmp] ~ dnorm(0, nu[19]) }
    alpha3[3] <- -1 * sum(alpha3[1:2])
    alpha4[1] ~ dnorm(0, nu[20])
    alpha4[2] <- -1 * sum(alpha4[2:3])
    for (cup in 1:(NCups-1))   { alpha5[cup]~ dnorm(0, nu[21]) }
    alpha5[NCups] <- -1 * sum(alpha5[1:(NCups-1)])
    for (rep in 1:2)   { alpha6[rep] ~ dnorm(0, nu[22]) }
    alpha6[3] <- -1 * sum(alpha6[1:2])
    for (chm in 1:2)   { alpha7[chm] ~ dnorm(0, nu[23]) }
    alpha7[3] <- -1 * sum(alpha7[1:2])
                           
  ################################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIABILITY ###
  ################################################
  nu[1] ~ dgamma(0.01, 0.01)
  nu[2] ~ dgamma(0.01, 0.01)
  nu[3] ~ dgamma(0.01, 0.01)
  nu[4] ~ dgamma(0.01, 0.01)
  nu[5] ~ dgamma(0.01, 0.01)
  nu[6] ~ dgamma(0.01, 0.01)
  nu[7] ~ dgamma(0.01, 0.01)
  nu[8] ~ dgamma(0.01, 0.01)
  nu[9] ~ dgamma(0.01, 0.01)
  nu[10] ~ dgamma(0.01, 0.01)
  nu[11] ~ dgamma(0.01, 0.01)
  nu[12] ~ dgamma(0.01, 0.01)
  nu[13] ~ dgamma(0.01, 0.01)
  nu[14] ~ dgamma(0.01, 0.01)
  nu[15] ~ dgamma(0.01, 0.01)
  nu[16] ~ dgamma(0.01, 0.01)
  nu[17] ~ dgamma(0.01, 0.01)
  nu[18] ~ dgamma(0.01, 0.01)
  nu[19] ~ dgamma(0.01, 0.01)
  nu[20] ~ dgamma(0.01, 0.01)
  nu[21] ~ dgamma(0.01, 0.01)
  nu[22] ~ dgamma(0.01, 0.01)
  nu[23] ~ dgamma(0.01, 0.01)

  ###########################################
  ### BACK-CALCULATE PRINCIPAL COMPONENTS ###
  ###########################################
  rho1 <- seedling_pca$rotation[1,1] + seedling_pca$rotation[1,2] + seedling_pca$rotation[1,3] + seedling_pca$rotation[1,4]
  rho2 <- seedling_pca$rotation[2,1] + seedling_pca$rotation[2,2] + seedling_pca$ rotation[2,3] + seedling_pca$rotation[2,4]
  rho3 <- seedling_pca$rotation[3,1] + seedling_pca$rotation[3,2] + seedling_pca$ rotation[3,3] + seedling_pca$rotation[3,4]
  rho4 <- seedling_pca$rotation[4,1] + seedling_pca$rotation[4,2] + seedling_pca$rotation[4,3] + seedling_pca$rotation[4,4]
  rho5 <- seedling_pca$rotation[5,1] + seedling_pca$rotation[5,2] + seedling_pca$rotation[5,3] + seedling_pca$rotation[5,4]
  rho6 <- seedling_pca$rotation[6,1] + seedling_pca$rotation[6,2] + seedling_pca$rotation[6,3] + seedling_pca$rotation[6,4]
  rho7 <- seedling_pca$rotation[7,1] + seedling_pca$rotation[7,2] + seedling_pca$rotation[7,3] + seedling_pca$rotation[7,4]
  rho8 <- seedling_pca$rotation[8,1] + seedling_pca$rotation[8,2] + seedling_pca$rotation[8,3] + seedling_pca$rotation[8,4]
  rho9 <- seedling_pca$rotation[9,1] + seedling_pca$rotation[9,2] + seedling_pca$rotation[9,3] + seedling_pca$rotation[9,4]
  rho10 <- seedling_pca$rotation[10,1] + seedling_pca$rotation[10,2] + seedling_pca$rotation[10,3] + seedling_pca$rotation[10,4]
  rho11 <- seedling_pca$rotation[11,1] + seedling_pca$rotation[11,2] + seedling_pca$rotation[11,3] + seedling_pca$rotation[11,4]
  rho12 <- seedling_pca$rotation[12,1] + seedling_pca$rotation[12,2] + seedling_pca$rotation[12,3] + seedling_pca$rotation[12,4]
  rho13 <- seedling_pca$rotation[13,1] + seedling_pca$rotation[13,2] + seedling_pca$rotation[13,3] + seedling_pca$rotation[13,4]
  rho14 <- seedling_pca$rotation[14,1] + seedling_pca$rotation[14,2] + seedling_pca$rotation[14,3] + seedling_pca$rotation[14,4]
  rho15 <- seedling_pca$rotation[15,1] + seedling_pca$rotation[15,2] + seedling_pca$rotation[15,3] + seedling_pca$rotation[15,4]
  rho16 <- seedling_pca$rotation[16,1] + seedling_pca$rotation[16,2] + seedling_pca$rotation[16,3] + seedling_pca$rotation[16,4]
  rho17 <- seedling_pca$rotation[17,1] + seedling_pca$rotation[17,2] + seedling_pca$rotation[17,3] + seedling_pca$rotation[17,4]
  rho18 <- seedling_pca$rotation[18,1] + seedling_pca$rotation[18,2] + seedling_pca$rotation[18,3] + seedling_pca$rotation[18,4]
  rho19 <- seedling_pca$rotation[19,1] + seedling_pca$rotation[19,2] + seedling_pca$rotation[19,3] + seedling_pca$rotation[19,4]
  rho20 <- seedling_pca$rotation[20,1] + seedling_pca$rotation[20,2] + seedling_pca$rotation[20,3] + seedling_pca$rotation[20,4]

}")
