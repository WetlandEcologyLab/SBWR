#### View estimates ###
## summarize outputs
library(grid)
summary <- summary(out_model)
write.csv(as.data.frame(summary$statistics), 'C:/Users/Maggie/Desktop/test_seedling_summary.csv')
write.csv(as.data.frame(summary$quantiles), 'C:/Users/Maggie/Desktop/test_quantiles_summary.csv')
pdf(file='C:/Users/Maggie/Desktop/GermModel_11Mar2021.pdf')
plot(out_model) # plot posteriors
dev.off()


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
sum_chain3 <- chain3[,1] + chain3[,473]; traceplot(sum_chain3)
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
# plot all posteriors 
# helpful website: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
#library(bayesplot) can't get the plots to work but they look nice...

##################
## violin plots ##
##################
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
kappa2_end_index = kappa1_end_index + germData$NSources - 1
kappa3_index = kappa2_end_index + 1
kappa4_index = kappa3_index + 1
mu1_index = kappa4_index + 2
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

# open PDF to save plots
library(grid)
pdf(file='C:/Users/Maggie/Desktop/GermModel_24DEC2020.pdf')

#### mu estimates ####
# boxplot
library(ggplot2)
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
library(vioplot)
vioplot(chains_merged[,mu1_index], col = "gray", main="Overall Species Effect")
#mtext("Species", 1, cex = 1.4, line = 2.5)
#mtext("Species effect size", 2, cex = 1.4, line = 2.)
#abline(h = 0, lty = 2)


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
#vioplot(chains_merged[,(beta9_index+1):(beta9_index+7)], col="grey", names=levels(as.factor(germination_matrix$Species)))
vioplot(kappa1_1, kappa1_2, kappa1_3, kappa1_4, kappa1_5, kappa1_6, kappa1_7, col = "gray", names=levels(as.factor(germination_matrix$Species)))
mtext("Species", 1, cex = 1.4, line = 2.5)
mtext("Species effect size", 2, cex = 1.4, line = 2.)
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
kappa2_1 <- chains_merged[,(kappa1_end_index+1)]
kappa2_2 <- chains_merged[,(kappa1_end_index+2)]
kappa2_3 <- chains_merged[,(kappa1_end_index+3)]
kappa2_4 <- chains_merged[,(kappa1_end_index+4)]
kappa2_5 <- chains_merged[,(kappa1_end_index+5)]
kappa2_6 <- chains_merged[,(kappa1_end_index+6)]
kappa2_7 <- chains_merged[,(kappa1_end_index+7)]
kappa2_8 <- chains_merged[,(kappa1_end_index+8)]
kappa2_9 <- chains_merged[,(kappa1_end_index+9)]
kappa2_10 <- chains_merged[,(kappa1_end_index+10)]
kappa2_11 <- chains_merged[,(kappa1_end_index+11)]
kappa2_12 <- chains_merged[,(kappa1_end_index+12)]
kappa2_13 <- chains_merged[,(kappa1_end_index+13)]
kappa2_14 <- chains_merged[,(kappa1_end_index+14)]
kappa2_15 <- chains_merged[,(kappa1_end_index+15)]
kappa2_16 <- chains_merged[,(kappa1_end_index+16)]
kappa2_17 <- chains_merged[,(kappa1_end_index+17)]
kappa2_18 <- chains_merged[,(kappa1_end_index+18)]
kappa2_19 <- chains_merged[,(kappa1_end_index+19)]
kappa2_20 <- chains_merged[,(kappa1_end_index+20)]
kappa2_21 <- chains_merged[,(kappa1_end_index+21)]
kappa2_22 <- chains_merged[,(kappa1_end_index+22)]
kappa2_23 <- chains_merged[,(kappa1_end_index+23)]
kappa2_24 <- chains_merged[,(kappa1_end_index+24)]
kappa2_25 <- chains_merged[,(kappa1_end_index+25)]
kappa2_26 <- chains_merged[,(kappa1_end_index+26)]
kappa2_27 <- chains_merged[,(kappa1_end_index+27)]
kappa2_28 <- chains_merged[,(kappa1_end_index+28)]
kappa2_29 <- chains_merged[,(kappa1_end_index+29)]
kappa2_30 <- chains_merged[,(kappa1_end_index+30)]
kappa2_31 <- chains_merged[,(kappa1_end_index+31)]
kappa2_32 <- chains_merged[,(kappa1_end_index+32)]
kappa2_33 <- chains_merged[,(kappa1_end_index+33)]
kappa2_34 <- chains_merged[,(kappa1_end_index+34)]
kappa2_35 <- chains_merged[,(kappa1_end_index+35)]
kappa2_36 <- chains_merged[,(kappa1_end_index+36)]
kappa2_37 <- chains_merged[,(kappa1_end_index+37)]
# violin plots
vioplot(kappa2_1, kappa2_2, kappa2_3, kappa2_4, kappa2_5, kappa2_6, kappa2_7, kappa2_8, kappa2_9, kappa2_10, kappa2_11, kappa2_12, kappa2_13, kappa2_14, kappa2_15, kappa2_16, kappa2_17, kappa2_18, kappa2_19, kappa2_20, kappa2_21, kappa2_22, kappa2_23, kappa2_24, kappa2_25, kappa2_26, kappa2_27, kappa2_28, kappa2_29, kappa2_30, kappa2_31, kappa2_32, kappa2_33, kappa2_34, kappa2_35, kappa2_36, kappa2_37, col = "gray", names=levels(as.factor(germination_matrix$Source))[1:(germData$NSources-1)], main="Source effect (relative to WASP-UT)")
mtext("Source", 1, cex = 1.4, line = 2.5)
mtext("Source effect size", 2, cex = 1.4, line = 2.)
abline(h = 0, lty = 2)


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
# violin plot
kappa3 <- chains_merged[,kappa3_index]
vioplot(kappa3, col = "gray", names="WP1", main="WP effect (relative to WP2)")
abline(h=0, lty=2)
mtext("Effect size", 2, cex = 1.4, line = 2.)

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
# violin plots
kappa4_1 <- chains_merged[,kappa4_index]
kappa4_2 <- chains_merged[,kappa4_index+1]
vioplot(kappa4_1,kappa4_2, col = "gray", names=c("28-10", "32-15"), main="Temp effect (relative to 36-20)")
abline(h=0, lty=2)
mtext("Effect size", 2, cex = 1.4, line = 2.)


#### beta1 estimates: Species effects ####
# BOMA over time
timesteps=germData$NTimes # number of timesteps
spp_level=1 # integer for species level
# calculate indices for variables 
total_steps <- germData$NTimes*7
index = seq(spp_level, total_steps, 7)
beta1s <- summary_quantiles[index,]
rownames(beta1s) <- as.integer(1:germData$NTimes)
ggplot(mapping=aes(x=factor(rownames(beta1s), levels=as.character(rownames(beta1s))),
                   ymin=beta1s[,1],
                   lower=beta1s[,2],
                   middle=beta1s[,3],
                   upper=beta1s[,4],
                   ymax=beta1s[,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle(paste("Beta1 estimates: BOMA through time", sep=''))+
  xlab("Timestep")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")
# violin plot
beta1_1 <- chains_merged[,index[1]]
beta1_2 <- chains_merged[,index[2]]
beta1_3 <- chains_merged[,index[3]]
beta1_4 <- chains_merged[,index[4]]
beta1_5 <- chains_merged[,index[5]]
beta1_6 <- chains_merged[,index[6]]
beta1_7 <- chains_merged[,index[7]]
beta1_8 <- chains_merged[,index[8]]
beta1_9 <- chains_merged[,index[9]]
beta1_10 <- chains_merged[,index[10]]
beta1_11 <- chains_merged[,index[11]]
beta1_12 <- chains_merged[,index[12]]
beta1_13 <- chains_merged[,index[13]]
beta1_14 <- chains_merged[,index[14]]
beta1_15 <- chains_merged[,index[15]]
beta1_16 <- chains_merged[,index[16]]
beta1_17 <- chains_merged[,index[17]]
beta1_18 <- chains_merged[,index[18]]
beta1_19 <- chains_merged[,index[19]]
beta1_20 <- chains_merged[,index[20]]
beta1_21 <- chains_merged[,index[21]]
beta1_22 <- chains_merged[,index[22]]
beta1_23 <- chains_merged[,index[23]]
beta1_24 <- chains_merged[,index[24]]
beta1_25 <- chains_merged[,index[25]]
beta1_26 <- chains_merged[,index[26]]
beta1_27 <- chains_merged[,index[27]]
beta1_28 <- chains_merged[,index[28]]
beta1_29 <- chains_merged[,index[29]]
beta1_30 <- chains_merged[,index[30]]
beta1_31 <- chains_merged[,index[31]]
beta1_32 <- chains_merged[,index[32]]
beta1_33 <- chains_merged[,index[33]]
beta1_34 <- chains_merged[,index[34]]
beta1_35 <- chains_merged[,index[35]]
beta1_36 <- chains_merged[,index[36]]
# violin plots
vioplot(beta1_1, beta1_2, beta1_3, beta1_4, beta1_5, beta1_6, beta1_7, beta1_8, beta1_9, beta1_10, beta1_11, beta1_12, beta1_13, beta1_14, beta1_15, beta1_16, beta1_17, beta1_18, beta1_19, beta1_20, beta1_21, beta1_22, beta1_23, beta1_24, beta1_25, beta1_26, beta1_27, beta1_28, beta1_29, beta1_30, beta1_31, beta1_32, beta1_33, beta1_34, beta1_35, beta1_36, col = "gray", names=rownames(beta1s), main="Boma effect over time")
abline(h=0, lty=2)
mtext("Effect size", 2, cex = 1.4, line = 2.)


# DISP over time
timesteps=germData$NTimes # number of timesteps
spp_level=2 # integer for species level
# calculate indices for variables 
total_steps <- timesteps*7
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
total_steps <- timesteps*7
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
total_steps <- timesteps*7
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
total_steps <- timesteps*7
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
total_steps <-timesteps*7
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
total_steps <- timesteps*7
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
timesteps=germData$NTimes # number of timesteps
# calculate indices for variables 
total_steps <- timesteps*3
index_temp1 = seq(1+beta3_end_index, total_steps+beta3_end_index, 3)
index_temp2 = seq(2+beta3_end_index, total_steps+beta3_end_index, 3)
# plot temp effect over time for 28-10
ggplot(mapping=aes(x=rownames(summary_quantiles)[index_temp1],
                   ymin=summary_quantiles[index_temp1,1],
                   lower=summary_quantiles[index_temp1,2],
                   middle=summary_quantiles[index_temp1,3],
                   upper=summary_quantiles[index_temp1,4],
                   ymax=summary_quantiles[index_temp1,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta4 estimates: Effect of 29-10 over time")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Timestep -->")+
  theme_classic()+
  geom_abline(aes(slope=0, intercept=0), col="red")
# plot temp effect over time 32-15
ggplot(mapping=aes(x=rownames(summary_quantiles)[index_temp2],
                   ymin=summary_quantiles[index_temp2,1],
                   lower=summary_quantiles[index_temp2,2],
                   middle=summary_quantiles[index_temp2,3],
                   upper=summary_quantiles[index_temp2,4],
                   ymax=summary_quantiles[index_temp2,5]))+ 
  geom_boxplot(stat="identity")+
  ggtitle("Beta4 estimates: Effect of 32-15 over time")+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Timestep -->")+
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
  xlab("Beta5")+
  theme_classic()

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
  geom_abline(aes(slope=0, intercept=0), col="red")+
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
ggplot(mapping=aes(x=rownames(summary$quantiles)[(beta7_end_index+1):beta8_end_index],
                   ymin=summary$quantiles[(beta7_end_index+1):beta8_end_index,1],
                   lower=summary$quantiles[(beta7_end_index+1):beta8_end_index,2],
                   middle=summary$quantiles[(beta7_end_index+1):beta8_end_index,3],
                   upper=summary$quantiles[(beta7_end_index+1):beta8_end_index,4],
                   ymax=summary$quantiles[(beta7_end_index+1):beta8_end_index,5]))+ 
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
  xlab("Timestep")

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

