#############################################
### CONVERT TO MATRIX FOR GERMINATION DATA ##
### Sown to Germinated Seed Matrix         ##
### Model 1                                ##
#############################################
## read in data in long format
data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel//DummyData/DummyGerminationData_Trial3_23Apr20.csv', header=T, stringsAsFactors=F)

## calculate total number of seeds
n_seeds<- 0
for (i in unique(data$Cup_No)){
  cup_subset <- data[which(data$Cup_No==i & data$RawTime==0),]
  n_seeds <- n_seeds+length(cup_subset$RawTime)
}

## create blank matrices
# ncol = length(unique(data$RawTime))
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

## populate matrices
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
      }#if(time==0)
      # switch germination observations to match with format in BPA book
      # 1 == not germinated (equivalent to "alive" in book)
      # 0 == germinated (equivalent to "dead" in book)
      if(data$GERMINATION[i]==1) germination_matrix[row,t] <- 0
      if (data$GERMINATION[i]==0) germination_matrix[row,t] <- 1
    }#if(data$RawTime[i]==time)
  }#t
}#i

## take care of NAs - set to 0
for (i in 1:length(germination_matrix)){
  if (is.na(germination_matrix[i])) germination_matrix[i] <- 0
}#i

## combine matrices into one
germination_matrix_full <- cbind(germination_matrix, species_matrix)
germination_matrix_full <- cbind(germination_matrix_full, source_matrix)
germination_matrix_full <- cbind(germination_matrix_full, wp_matrix)
germination_matrix_full <- cbind(germination_matrix_full, temp_matrix)
germination_matrix_full <- cbind(germination_matrix_full, trial_matrix)
germination_matrix_full <- cbind(germination_matrix_full, cupNo_matrix)
germination_matrix_full <- cbind(germination_matrix_full, SeedlingNo_matrix)
germination_matrix_full <- cbind(germination_matrix_full, chamber_matrix)
germination_matrix_full <- cbind(germination_matrix_full, HTT_matrix)
germination_matrix_full <- cbind(germination_matrix_full, seedmass_matrix)
germination_matrix_full <- cbind(germination_matrix_full, sct_matrix)

## save as CSV
write.csv(germination_matrix_full, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixGermDataTrial3_28Apr_2020.csv')

#################################################
### CONVERT TO MATRIX FOR 7-DAY SEEDLING DATA ###
### Germinated Seed to 7-Day Seedling Matrix  ###
### Model 2                                   ###
#################################################
## load dependencies
library(stringr)

## read in data
## populate matrices with pulled seedling data
## NOTE: For this data processing to work, dead seedlings must be included in the dataset with "dead" in the Notes column and estimated age at death in the Age column (if not pulled, otherwise age=age pulled)
seedling_data <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/DummySeedlingsDataset.csv', header=TRUE, stringsAsFactors=FALSE)
germination_matrix <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixGermDataTrial3_28Apr_2020.csv', header=TRUE, stringsAsFactors=FALSE)
seed_mass <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedMass_FINAL_Processed.csv', header=TRUE, stringsAsFactors=FALSE)
names(seed_mass)[1] <- "Source"
## NOTE: Not sure if 21-day seedlings should also be included in here without any predictor data.... for now leaving them in

## Preprocessing data
# add seed mass field
seedling_data$SeedMass <- NA
for (src in unique(seed_mass$Site)){
  seedmass_subset <- seed_mass[which(seed_mass$Site==src),]
  seedmass <- mean(seedmass_subset$AvgSeedMass, na.rm=TRUE)
  seedling_data$SeedMass[which(seedling_data$Source==src)] <- seedmass
}#for
# add seed coat thickness field
seedling_data$SCT <- NA
for (src in unique(seed_dims$Source)){
  seeddims_subset <- seeddims[which(seeddims$Source==src),]
  sct <- mean(seeddims_subset$SCT, na.rm=TRUE)
  seedling_data$SCT[which(seedling_data$Source==src)] <- sct
}#for
# calculate field for "days until early seedling" in seedling dataset
seedling_data$RawTime <- seedling_data$Germ_Date + 8
# calculate pull age for cup dataset
cup_data$PullAge <- cup_data$Pull_Date - cup_data$Germ_Date

## NOTE: If not incorporating 21 day seedlings into this model, this code needs to be changed
## get number of germinated seeds per cup
#ncups = length(unique(germination_matrix$CupNo))
#Ngerm_per_cup <- as.data.frame(matrix(NA, nrow=ncups, ncol=2))
#i = 0 # set counter
#total_germs = 0 # set initial value
#for (cup in unique(seedling_data$Cup_No)){
#  i = i+1
#  germ_cup_subset <- germination_matrix[which(germination_matrix$CupNo==cup),]
#  ngerms <- nrow(germ_cup_subset)
#  Ngerm_per_cup[i,] <- c(cup, ngerms)
#  total_germs <- total_germs + ngerms
#}#for
total_germs <- nrow(seedling_data)

## set up blank matrices
# nrow = total_germs from above
# ncol = max(seedlings_7day$PullDate)/41 - i.e. # of data collection days
seedling7_matrix <- matrix(NA, nrow=total_germs, ncol=41)
species_matrix <- matrix(NA, nrow=total_germs, ncol=1)
source_matrix <- matrix(NA, nrow=total_germs, ncol=1)
trial_matrix <- matrix(NA, nrow=total_germs, ncol=1)
wp_matrix <- matrix(NA, nrow=total_germs, ncol=1)
temp_matrix <- matrix(NA, nrow=total_germs, ncol=1)
cupNo_matrix <- matrix(NA, nrow=total_germs, ncol=1)
SeedlingNo_matrix <- matrix(NA, nrow=total_germs, ncol=1)
chamber_matrix <- matrix(NA, nrow=total_germs, ncol=1)
HTT_matrix <- matrix(NA, nrow=total_germs, ncol=1)
seedmass_matrix <- matrix(NA, nrow=total_germs, ncol=1)
sct_matrix <- matrix(NA, nrow=total_germs, ncol=1)
age_matrix  <- matrix(NA, nrow=total_germs, ncol=1)
rootLength_matrix <- matrix(NA, nrow=total_germs, ncol=1)
shootLength_matrix <- matrix(NA, nrow=total_germs, ncol=1)
rootSA_matrix <- matrix(NA, nrow=total_germs, ncol=1)
shootSA_matrix <- matrix(NA, nrow=total_germs, ncol=1)
totalLength_matrix <- matrix(NA, nrow=total_germs, ncol=1)
totalSA_matrix <- matrix(NA, nrow=total_germs, ncol=1)
SER_matrix <- matrix(NA, nrow=total_germs, ncol=1)
RER_matrix <- matrix(NA, nrow=total_germs, ncol=1)
shootLengthRatio_matrix <- matrix(NA, nrow=total_germs, ncol=1)
AGFM_matrix <- matrix(NA, nrow=total_germs, ncol=1)
BGFM_matrix <- matrix(NA, nrow=total_germs, ncol=1)
AGDM_matrix <- matrix(NA, nrow=total_germs, ncol=1)
BGDM_matrix <- matrix(NA, nrow=total_germs, ncol=1)
totalMass_matrix <- matrix(NA, nrow=total_germs, ncol=1)
shootMassRatio_matrix <- matrix(NA, nrow=total_germs, ncol=1)
germTimePercentile_matrix  <- matrix(NA, nrow=total_germs, ncol=1)
germDate_matrix <- matrix(NA, nrow=total_germs, ncol=1)

# add all seedling data to matrices
NTimes <- 41 # of timesteps (i.e., # of data collection days)
row=0 # set row counter
for (i in 1:nrow(seedling_data)){
  row=row+1
  # add all covariate data to matrices
  species_matrix[row] <- seedling_data$Species[i]
  source_matrix[row] <- seedling_data$Source[i]
  trial_matrix[row] <- seedling_data$Trial[i]
  wp_matrix[row] <- seedling_data$WP[i]
  temp_matrix[row] <- seedling_data$Temp[i]
  cupNo_matrix[row] <- seedling_data$Cup_No[i]
  SeedlingNo_matrix[row] <- seedling_data$Seedling_No[i]
  chamber_matrix[row] <- seedling_data$Chamber[i]
  HTT_matrix[row] <- seedling_data$HydrothermalTime[i]
  seedmass_matrix[row] <- seedling_data$SeedMass[i]
  sct_matrix[row] <- seedling_data$SeedCoatThickness[i]
  age_matrix[row] <- seedling_data$Age[i]
  germDate_matrix[row] <- seedling_data$Germ_Date[i]
  germTimePercentile_matrix[row]  <- seedling_data$GermTimePercentile[i]
  # only add seedling metrics if pulled at early seedling stage
  if (seedling_data$Age[i]==8) {
    rootLength_matrix[row] <- seedling_data$Root_Length[i]
    shootLength_matrix[row] <- seedling_data$Shoot_Length[i]
    rootSA_matrix[row] <- seedling_data$Root_SA[i]
    shootSA_matrix[row] <- seedling_data$Shoot_SA[i]
    totalLength_matrix[row] <- seedling_data$Total_Length[i]
    totalSA_matrix[row] <- seedling_data$Total_SA[i]
    SER_matrix[row] <- seedling_data$ShootElongationRate[i]
    RER_matrix[row] <- seedling_data$RootElongationRate[i]
    shootLengthRatio_matrix[row] <- seedling_data$ShootLengthRatio[i]
    AGFM_matrix[row] <- seedling_data$Above_Fresh[i]
    BGFM_matrix[row] <- seedling_data$Below_Fresh[i]
    AGDM_matrix[row] <- seedling_data$Above_Dry[i]
    BGDM_matrix[row] <- seedling_data$Below_Dry[i]
    totalMass_matrix[row] <- seedling_data$Total_Dry[i]
    shootMassRatio_matrix[row] <- seedling_data$ShootMassRatio[i]
  # set seedling metrics to NA for any seedlings not harvested at this stage
  } else {
    rootLength_matrix[row] <- NA
    shootLength_matrix[row] <- NA
    rootSA_matrix[row] <- NA
    shootSA_matrix[row] <- NA
    totalLength_matrix[row] <- NA
    totalSA_matrix[row] <- NA
    SER_matrix[row] <- NA
    RER_matrix[row] <- NA
    shootLengthRatio_matrix[row] <- NA
    AGFM_matrix[row] <- NA
    BGFM_matrix[row] <- NA
    AGDM_matrix[row] <- NA
    BGDM_matrix[row] <- NA
    totalMass_matrix[row] <- NA
    shootMassRatio_matrix[row] <- NA
  }#ifelse
  # if "dead" in notes, set all to 1 (never transitioned)
  ## NOTE: This assumes that any dead seedling are flagged with "dead" in the notes section!!!
  if (str_detect(seedling_data$Notes[i], "dead")){ 
    # if seedling was dead at or before age 8, then set to never transitioned
    if (seedling_data$Age[i]<=8){
      seedling7_matrix[row,] <- 1
    # if died after age 8, set transition at germDate + 8
    } else {
      t <- seedlingdata$GermDate[i] + 8
      seedling7_matrix[row, 1:(t-1)] <- 1
      seedling7_matrix[row, t:NTimes]
    }
  # else set 0/1 transition responses
  } else {
    for (t in 1:NTimes){
      time <- (t-1)*2
      if (seedling_data$RawTime[i]==time){
        # set all timesteps up to transition (i.e. RawTime) to 1
        seedling7_matrix[row, 1:(t-1)] <- 1
        # set all timesteps at and after transition to 0
        seedling7_matrix[row, t:NTimes] <- 0
      }#if
    }#for
  }#else
}#for

# combine matrices and save as CSV
seedling7_matrix_full <- cbind(seedling7_matrix, species_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, source_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, wp_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, temp_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, trial_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, cupNo_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, SeedlingNo_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, chamber_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, HTT_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, seedmass_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, sct_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, age_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, germDate_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, germTimePercentile_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, rootLength_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, shootLength_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, rootSA_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, shootSA_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, totalLength_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, totalSA_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, SER_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, RER_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, shootLengthRatio_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, AGFM_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, BGFM_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, AGDM_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, BGDM_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, totalMass_matrix)
seedling7_matrix_full <- cbind(seedling7_matrix_full, shootMassRatio_matrix)

## save to file
seedling7_data_frame <- as.data.frame(seedling7_matrix_full)
names(seedling7_data_frame) <- c(seq(2, NTimes*2, 2), "Species", "Source", "WP", "Temp", "Trial", "CupNo", "SeedlingNo", "Chamber", "HTT", "SeedMass", "SCT", "Age", "GermDate", "GermTimePercentile", "RootLength", "ShootLength", "RootSA", "ShootSA", "TotalLength", "TotalSA", "SER", "RER", "ShootLengthRatio", "AboveGroundFresh", "BelowGroundFresh", "AboveGroundDry", "BelowGroundDry", "TotalDryMass", "ShootMassRatio")
write.csv(seedling7_data_frame, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/DummyData/MatrixEarlySeedlingData_21May2020.csv')

