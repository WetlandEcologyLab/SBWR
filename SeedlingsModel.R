#install.packages("rjags")
library(rjags)

# load data


##############################################################
### DEFINE SURVIVAL PROBABILITY HIERARCHICAL PROBABILITIES ###
##############################################################

####################################################
### Transition from sown seed to germinated seed ###
####################################################
SownSurv <- textConnection("model{ # has the same effect as writing all this in a separate text file

  ###########################
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
  # y = vector of all observations
  # beta1 has multiple estimates for each species category
  # beta3 and beta4 have single coefficients and use each observations hydrothermal time and seed mass value as a continuous parameter 
  for (i in 1:NGerm) {
    sown_linear <- beta1[species[i]]+beta2[source[i]]+beta3*htt[i]+beta4*seedmass[i]+beta5*sct[i]+beta6[cup[i]]+beta7[chamber[i]]+beta8[rep[i]]+eps[i]
    # transformation to logit scale
    prob_sown[i] <- exp(sown_linear) / (1+exp(sown_linear))
  }

  ## Likelihood parameter depends on each observation
  for (i in 1:NSown)  { germ[i] <- dbinom(prob_sown[i], length(y)) }         

  ## prior on overall variation
  for (i in 1:NSown)  { eps[i] <- dgamma(0.01, 0.01) }

  ##############################
  ### PRIORS ON MAIN EFFECTS ###
  ########################summary(pca######
  # species effect on germination probability - categorical
  for (spp in 1:NSpecies)   { beta1[spp] ~ dnorm(0, nu[1]) }
  
  # source effect on germination probability - categorical
  for (src in 1:NSources)   { beta2[src] ~ dnorm(0, nu[2]) }
  #beta2[NSources] <- sum(beta2[1:(NSources-1)]) # sum-to-zero constraint

  # hydrothermal time effect on germination - continuous
  beta3 ~ dnorm(0, nu[3])

  # seed mass effect on germination - continuous
  beta4 ~ dnorm(0, nu[4])

  # seed coat thickness effect on germination - continuous
  beta5 ~ dnorm(0, nu[5])

  # cup effect on germination - categorical
  for (cup in 1:NCups)   { beta6[cup] ~ dnorm(0, nu[6]) }

  # chamber effect on germination - categorical
  for (chm in 1:3)   { beta7[chm] ~ dnorm(0, nu[7]) }

  # rep effect on germination - categorical
  for (rep in 1:3)   { beta8[rep] ~ dnorm(0, nu[8]) }
  
  ###############################
  ### HYDROTHERMAL TIME MODEL ###
  ###############################
  # likelihood
  for (i in 1:Nhtt)   { htt[i] ~ dnorm(mu1[i], tau1) }
  # linear model
  for (i in 1:Nhtt) {
    mu1[i] <- eta1[species[i]] + eta2[source[i]] + eta3[cup[i]] + eta4[chamber[i]] + eta5[rep[i]]
  }
    # prior on overall variation
    tau1 ~ dgamma(0.01, 0.01)
    # prior on species effect
    for (spp in 1:(NSpecies))   { eta1[spp] ~ dnorm(0, nu[9]) 
    # prior on source effect
    for (src in 1:NSources)  { eta2[src] ~ dnorm(0, nu[10]) }
    #eta2[NSources] <- -1 * sum(eta2[1:(NSources-1)])
    # prior on cup effect
    for (cup in 1:NCups)  { eta3[cup] ~ dnorm(0, nu[11]) }
    #eta3[NCups] <- -1 * sum(eta3[1:(NCups-1)])
    # prior on chamber effect
    for (chm in 1:3)  { eta4[chm] ~ dnorm(0, nu[12]) }
    #mu3[3] <- -1 * sum(eta4[1:2])
    # prior on rep effect
    for (rep in 1:3)  { eta5[rep] ~ dnorm(0, nu[13]) }
    # eta5[3] <- -1 * sum(eta5[1:2])

  ##############################################
  ### UNINFORMATIVE HYPERPRIORS ON VARIATION ###
  ##############################################
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
}")


###################################################
### Transition from germinated seed to seedling ###
###################################################
GermSurv <- textConnection("model{ # has the same effect as writing all this in a separate text file

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

  ########################### 
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
  for (i in 1:NGerm) {
    germ_linear[i] <- beta1[species[i]]+beta2[source[i]]+beta3[wp[i]]+beta4[temp[i]]+beta5*pc1[i]+beta6*pc2[i]+beta7*pc3[i]+beta8*pc4[i]+beta9[cup[i]]+beta10[chamber[i]]+beta11[rep[i]]+eps[i]
    # transformation on logit scale
    prob_germ[i] <- exp(germ_linear) / (1+exp(germ_linear))
  }
  
  ## Likelihood parameter depends on each observation
  for (i in 1:NGerm)    { y[i] ~ dbinom(prob_germ[i], length(y)) }

  ## Uninformative prior on overall variability
  for (i in 1:NGerm)   { eps[i] ~ dgamma(0.01, 0.01) } 

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
  for (tmp in 1:3)   { beta4[tmp] ~ dnorm(0, nu[4]) }

  # prior on PC1 - continuous
  beta5 ~ dnorm(0, nu[5])

  # prior on PC2 - continuous
  beta6 ~ dnorm(0, nu[6])

  # prior on PC3 - continuous
  beta7 ~ dnorm(0, nu[7])

  # prior on PC4 - continuous
  beta8 ~ dnorm(0, nu[8])

  # prior on cup effect - categorical
  for (cup in 1:NCups)   { beta9[cup] ~ dnorm(0, nu[9]) }

  # prior on chamber effect - categorical
  for (chm in 1:3)   { beta10[chm] ~ dnorm(0, nu[10]) }

  # prior on rep effect - categorical
  for (rep in 1:3)   { beta11[rep] ~ dnorm(0, nu[11]) }

  #######################################################
  ## model seedling traits or principal components???  ##
  #######################################################
  
  ############################
  ### Growth rates per cup ###
  ############################
  # relative growth rate model
  # Nrgr = number of observations of RGR
  for (i in 1:Nrgr)  { rgr ~ dnorm(mu1, tau1) }
  for (i in 1:Nrgr) {
    mu1 <- kappa1[species[i]]+kappa2[source[i]]+kappa3[temp[i]]+kappa4[wp[i]]+kappa5[cup[i]]+kappa6[rep[i]]+kappa7[chamber[i]]
    tau1 ~ dgamma(0.01, 0.01)
  }
    for (spp in 1:NSpecies)   { kappa1[spp] ~ dnorm(0, nu[12]) }
    for (src in 1:(NSources-1))   { kappa2[src] ~ dnorm(0, nu[13]) }
    kappa2[NSources] <- -1 * sum(kappa2[1:(NSources-1)])
    for (tmp in 1:2)   { kappa3[tmp] ~ dnorm(0, nu[14]) }
    kappa3[3] <- -1 * sum(kappa3[1:2])
    kappa4[1] ~ dnorm(0, nu[15])
    kappa4[2] <- -1 * sum(kappa4[1])
    for (cup in 1:(NCups-1))   { kappa5[cup]~ dnorm(0, nu[16]) }
    kappa5[NCups] <- -1 * sum(kappa5[1:(NCups-1)])
    for (rep in 1:2)   { kappa6[rep] ~ dnorm(0, nu[17]) }
    kappa6[3] <- -1 * sum(kappa6[1:2])
    for (chm in 1:2)   { kappa7[chm] ~ dnorm(0, nu[18]) }
    kappa7[3] <- -1 * sum(kappa7[1:2])
  
  # root elongation rate
  # Nrer = number of observations of RER
  for (i in 1:Nrer)  { rer ~ dnorm(mu2, tau2) }
  for (i in 1:Nrer) {
      mu2 <- alpha1[species[i]]+alpha2[source[i]]+alpha3[temp[i]]+alpha4[wp[i]]+alpha5[cup[i]]+alpha6[rep[i]]+alpha7[chamber[i]]
      tau2 ~ dgamma(0.01, 0.01)
    }
    for (spp in 1:NSpecies)   { alpha1[spp] ~ dnorm(0, nu[19]) }
    for (src in 1:(NSources-1))   { alpha2[src] ~ dnorm(0, nu[20]) }
    alpha2[NSources] <- -1 * sum(alpha2[1:(NSources-1)])
    for (tmp in 1:2)   { alpha3[tmp] ~ dnorm(0, nu[21]) }
    alpha3[3] <- -1 * sum(alpha3[1:2])
    alpha4[1] ~ dnorm(0, nu[22])
    alpha4[2] <- -1 * sum(alpha4[2:3])
    for (cup in 1:(NCups-1))   { alpha5[cup]~ dnorm(0, nu[23]) }
    alpha5[NCups] <- -1 * sum(alpha5[1:(NCups-1)])
    for (rep in 1:2)   { alpha6[rep] ~ dnorm(0, nu[24]) }
    alpha6[3] <- -1 * sum(alpha6[1:2])
    for (chm in 1:2)   { alpha7[chm] ~ dnorm(0, nu[25]) }
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
  nu[24] ~ dgamma(0.01, 0.01)
  nu[25] ~ dgamma(0.01, 0.01)

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

#########################################
### Transition from seedling to adult ###
#########################################
SeedSurv <- textConnection("model{ # has the same effect as writing all this in a separate text file
                           
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
                           
  ########################### 
  ### DEFINE LINEAR MODEL ###
  ###########################
  ## Linear model of rate parameter based on each observation
  for (i in 1:NSeed) {
    seed_linear[i] <- beta1[species[i]]+beta2[source[i]]+beta3[wp[i]]+beta4[temp[i]]+beta5*pc1[i]+beta6*pc2[i]+beta7*pc3[i]+beta8*pc4[i]+beta9[cup[i]]+beta10[chamber[i]]+beta11[rep[i]]+eps[i]
    # transformation on logit scale
    prob_seed[i] <- exp(seed_linear) / (1+exp(seed_linear))
  }
                           
  ## Likelihood parameter depends on each observation 
  for (i in 1:NSeed)    { y[i] ~ dbinom(prob_seed[i], length(y)) }
                           
  ## Uninformative prior on overall variability
  for (i in 1:NSeed)   { eps[i] ~ dgamma(0.01, 0.01) } 
                           
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
  for (tmp in 1:3)   { beta4[tmp] ~ dnorm(0, nu[4]) }
                           
  # prior on PC1 effect - continuous
  beta5 ~ dnorm(0, nu[5])
                           
  # prior on PC2 effect - continuous
  beta6 ~ dnorm(0, nu[6])
                           
  # prior on PC3 effect - continuous
  beta7 ~ dnorm(0, nu[7])
                           
  # prior on PC4 effect - continuous
  beta8 ~ dnorm(0, nu[8])
                           
  # prior on cup effect - categorical
  for (cup in 1:NCups)   { beta9[cup] ~ dnorm(0, nu[9]) }
                           
  # prior on chamber effect - categorical
  for (chm in 1:3)   { beta10[chm] ~ dnorm(0, nu[10]) }
                           
  # prior on rep effect - categorical
  for (rep in 1:3)   { beta11[rep] ~ dnorm(0, nu[11]) }
                           
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
    for (spp in 1:NSpecies)   { kappa1[spp] ~ dnorm(0, nu[12]) }
    for (src in 1:(NSources-1))   { kappa2[src] ~ dnorm(0, nu[13]) }
    kappa2[NSources] <- -1 * sum(kappa2[1:(NSources-1)])
    for (tmp in 1:2)   { kappa3[tmp] ~ dnorm(0, nu[14]) }
    kappa3[3] <- -1 * sum(kappa3[1:2])
    kappa4[1] ~ dnorm(0, nu[15])
    kappa4[2] <- -1 * sum(kappa4[1])
    for (cup in 1:(NCups-1))   { kappa5[cup]~ dnorm(0, nu[16]) }
    kappa5[NCups] <- -1 * sum(kappa5[1:(NCups-1)])
    for (rep in 1:2)   { kappa6[rep] ~ dnorm(0, nu[17]) }
    kappa6[3] <- -1 * sum(kappa6[1:2])
    for (chm in 1:2)   { kappa7[chm] ~ dnorm(0, nu[18]) }
    kappa7[3] <- -1 * sum(kappa7[1:2])
                           
  # root elongation rate
  for (i in 1:Nrer)  { rer ~ dnorm(mu2, tau2) }
  for (i in 1:Nrer) {
    mu2 <- alpha1[species[i]]+alpha2[source[i]]+alpha3[temp[i]]+alpha4[wp[i]]+alpha5[cup[i]]+alpha6[rep[i]]+alpha7[chamber[i]]
    tau2 ~ dgamma(0.01, 0.01)
  }
    for (spp in 1:NSpecies)   { alpha1[spp] ~ dnorm(0, nu[19]) }
    for (src in 1:(NSources-1))   { alpha2[src] ~ dnorm(0, nu[20]) }
    alpha2[NSources] <- -1 * sum(alpha2[1:(NSources-1)])
    for (tmp in 1:2)   { alpha3[tmp] ~ dnorm(0, nu[21]) }
    alpha3[3] <- -1 * sum(alpha3[1:2])
    alpha4[1] ~ dnorm(0, nu[22])
    alpha4[2] <- -1 * sum(alpha4[2:3])
    for (cup in 1:(NCups-1))   { alpha5[cup]~ dnorm(0, nu[23]) }
    alpha5[NCups] <- -1 * sum(alpha5[1:(NCups-1)])
    for (rep in 1:2)   { alpha6[rep] ~ dnorm(0, nu[24]) }
    alpha6[3] <- -1 * sum(alpha6[1:2])
    for (chm in 1:2)   { alpha7[chm] ~ dnorm(0, nu[25]) }
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
  nu[24] ~ dgamma(0.01, 0.01)
  nu[25] ~ dgamma(0.01, 0.01)
  
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

#########################################
### Transition from adult to adult ###
#########################################
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


##################
### RUN MODELS ###
##################

# save data as appropriate list
seedlingData <- list(y=data$Length, species=data$Species, source=data$Source, wp=data$WP,
                   temp=data$Temp, cup=data$Cup, 
                   # chamber=data$Chamber, 
                   NSpecies = length(unique(data$Species)), 
                   NSources = length(unique(data$Source)),
                   #NChambers = length(unique(data$Chamber)), 
                   #NReps = length(unique(data$Trial)),
                   NCups = length(unique(data$Cup)))

## Initialize model
seedlingModel <- jags.model(file=DraftModel, data=seedlingData, n.chains = 3)

## Run model
update(seedlingModel, n.iter=200000)
out <- coda.samples(model=seedlingModel, variable.names = c("beta1", "beta2", "beta3", "beta4", "beta6", "nu", "tau"), thin = 100, n.iter = 10000)

## Assess convergence
orig.mar <- par('mar')
par(mar=c(1,1,1,1))
plot(out)
par(mar=orig.mar)
effectiveSize(out)
gelman.diag(out,multivariate=F)

## View outputs
summary(out)
# beta0 = overall mean
# beta1 = species (BOMA, DISP, ELPA, JUBA, PHAU, SCAC, SCAM)
# beta2 = sources
# beta3 = water potential (1, 2)
# beta4 = temperature (1, 2, 3)
# beta5 = chamber (1, 2, ..., 9)
# beta6 = rep (1, 2, 3)
# alpha[1] = average species effect
# alpha[2] = average source effect
# alpha[3] = average wp effect
# alpha[4] = average temp effect
# nu[1] = overall species variance
# nu[2] = overall source variance
# nu[3] = overall wp variance
# nu[4] = overall temp variance
# nu[5] = overall chamber variance
# nu[6] = overall rep variance
# tau = overall variance



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
