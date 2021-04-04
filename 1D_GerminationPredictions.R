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
