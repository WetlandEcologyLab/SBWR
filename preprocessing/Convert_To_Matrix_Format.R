######################################
### CONVERT INPUT TO MATRIX FORMAT ###
######################################
data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/DummyGerminationData_Trial3_23Apr20.csv',header=T,stringsAsFactors = FALSE)

n_seeds<- 0
for (i in unique(data$Cup_No)){
  cup_subset <- data[which(data$Cup_No==i & data$RawTime==0),]
  n_seeds <- n_seeds+length(cup_subset$RawTime)
}
# ncol = unique(data$RawTime)
# nrow=n_seeds from above
germination_matrix <- matrix(NA, nrow=n_seeds, ncol=length(unique(data$RawTime)))
species_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
source_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
trial_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
wp_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
temp_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
cupNo_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
SeedlingNo_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
chamber_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
HTT_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
seedmass_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
sct_matrix <- matrix(NA, nrow=n_seeds, ncol=1)
  
row=0
for (i in 1:nrow(data)){
  for (t in 1:34){
    time <- (t-1)*2
    if (data$RawTime[i]==time){
      if (time==0) {
        row=row+1 
        species_matrix[row] <- data$Species[i]
        source_matrix[row] <- data$Source[i]
        trial_matrix[row] <- data$Trial[i]
        wp_matrix[row] <- data$WP[i]
        temp_matrix[row] <- data$Temp[i]
        cupNo_matrix[row] <- data$Cup_No[i]
        SeedlingNo_matrix[row] <- data$Seedling_No[i]
        chamber_matrix[row] <- data$Chamber[i]
        HTT_matrix[row] <- data$HydrothermalTime[i]
        seedmass_matrix[row] <- data$SeedMass[i]
        sct_matrix[row] <- data$SeedCoatThickness[i]
      }
      germination_matrix[row,t] <- data$GERMINATION[i]
    }
  }
}

for (i in 1:length(germination_matrix)){
  if (is.na(germination_matrix[i])) germination_matrix[i] <- 1
}

germination_matrix <- cbind(germination_matrix, species_matrix)
germination_matrix <- cbind(germination_matrix, source_matrix)
germination_matrix <- cbind(germination_matrix, wp_matrix)
germination_matrix <- cbind(germination_matrix, temp_matrix)
germination_matrix <- cbind(germination_matrix, cupNo_matrix)
germination_matrix <- cbind(germination_matrix, SeedlingNo_matrix)
germination_matrix <- cbind(germination_matrix, chamber_matrix)
germination_matrix <- cbind(germination_matrix, HTT_matrix)
germination_matrix <- cbind(germination_matrix, seedmass_matrix)
germination_matrix <- cbind(germination_matrix, sct_matrix)
germination_matrix_subset<- germination_matrix[,5:10]

write.csv(germination_matrix,"C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixGermDataTrial3_30Apr_2020.csv")
