# load in germination data
germination_matrix <- read.csv('C:/Users/Maggie/Desktop/FullGerminationData_FINAL.csv', header=TRUE, stringsAsFactors=FALSE)
names(germination_matrix)[1] <- "Trial"

# add chamber field to dataset
germination_matrix$Chamber <- rep(NA, nrow(germination_matrix))
for (i in 1:nrow(germination_matrix)){
  if(germination_matrix$Trial[i]==3 & germination_matrix$Temp[i]==1) germination_matrix$Chamber[i] <- 1
  if(germination_matrix$Trial[i]==3 & germination_matrix$Temp[i]==2) germination_matrix$Chamber[i] <- 2
  if(germination_matrix$Trial[i]==3 & germination_matrix$Temp[i]==3) germination_matrix$Chamber[i] <- 3
  if(germination_matrix$Trial[i]==4 & germination_matrix$Temp[i]==1) germination_matrix$Chamber[i] <- 4
  if(germination_matrix$Trial[i]==4 & germination_matrix$Temp[i]==2) germination_matrix$Chamber[i] <- 5
  if(germination_matrix$Trial[i]==4 & germination_matrix$Temp[i]==3) germination_matrix$Chamber[i] <- 6
  if(germination_matrix$Trial[i]==5 & germination_matrix$Temp[i]==1) germination_matrix$Chamber[i] <- 7
  if(germination_matrix$Trial[i]==5 & germination_matrix$Temp[i]==2) germination_matrix$Chamber[i] <- 8
  if(germination_matrix$Trial[i]==5 & germination_matrix$Temp[i]==3) germination_matrix$Chamber[i] <- 9
}#for

# make cup numbers unique
for (i in 1:nrow(germination_matrix)){
  if(germination_matrix$Trial[i]==4) germination_matrix$CupNo[i] <- germination_matrix$CupNo[i] + 234
  if(germination_matrix$Trial[i]==5) germination_matrix$CupNo[i] <- germination_matrix$CupNo[i] + 468
}#for

# clean up seed mass data
seedmass_csv <- read.csv('C:/Users/Maggie/Downloads/SeedMass_FINAL.csv', header=TRUE, stringsAsFactors=FALSE)
names(seedmass_csv)[1] <- "Source" 
# find out which source names do not match with the master dataset
sort(unique(seedmass_csv$Site[which(!seedmass_csv$Site%in%germination_matrix$Source)]))
sort(unique(germination_matrix$Source))
# fix source names
for (i in 1:220){
  ifelse(seedmass_csv$Source[i]=="Keith" & seedmass_csv$Species[i]=="JUBA", seedmass_csv$Site[i] <- "JUNB",
    ifelse(!is.na(seedmass_csv$Site[i]),
    ifelse(seedmass_csv$Site[i]=="JUNBAL-15", seedmass_csv$Site[i] <- "JUNB",
    ifelse(seedmass_csv$Site[i]=="RRVA-BW" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "RRVABW1",
    ifelse(seedmass_csv$Site[i]=="RRVA-BW" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "RRVABW2",
    ifelse(seedmass_csv$Site[i]=="RRVA-BW" & seedmass_csv$Species[i]=="SCAM", seedmass_csv$Site[i] <- "RRVABW3",
    ifelse(seedmass_csv$Site[i]=="PAHR" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "PAHR1",
    ifelse(seedmass_csv$Site[i]=="PAHR" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "PAHR2",
    ifelse(seedmass_csv$Site[i]=="FROU" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "FROU1",
    ifelse(seedmass_csv$Site[i]=="FROU" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "FROU2",
    ifelse(seedmass_csv$Site[i]=="WASP-MT", seedmass_csv$Site[i] <- "WASPMT",
    ifelse(seedmass_csv$Site[i]=="SACR" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "SACR1",
    ifelse(seedmass_csv$Site[i]=="SACR" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "SACR2",
    ifelse(seedmass_csv$Site[i]=="BERI" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "BERI1",
    ifelse(seedmass_csv$Site[i]=="BERI" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "BERI2",
    ifelse(seedmass_csv$Site[i]=="BERI" & seedmass_csv$Species[i]=="PHAU", seedmass_csv$Site[i] <- "BERI3",
    ifelse(seedmass_csv$Site[i]=="FABA" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "FABA1",
    ifelse(seedmass_csv$Site[i]=="FABA" & seedmass_csv$Species[i]=="PHAU", seedmass_csv$Site[i] <- "FABA3",
    ifelse(seedmass_csv$Site[i]=="FISP" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "FISP1",
    ifelse(seedmass_csv$Site[i]=="FISP" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "FISP2",
    ifelse(seedmass_csv$Site[i]=="FISP" & seedmass_csv$Species[i]=="SCAM", seedmass_csv$Site[i] <- "FISP3",
    ifelse(seedmass_csv$Site[i]=="CLLA" & seedmass_csv$Species[i]=="BOMA", seedmass_csv$Site[i] <- "CLLA1",
    ifelse(seedmass_csv$Site[i]=="CLLA" & seedmass_csv$Species[i]=="SCAC", seedmass_csv$Site[i] <- "CLLA2",
    ifelse(seedmass_csv$Site[i]=="WASP-UT", seedmass_csv$Site[i] <- "WASPUT",
    ifelse(seedmass_csv$Site[i]=="SCHAME-113", seedmass_csv$Site[i] <- "SHAM",
    ifelse(seedmass_csv$Site[i]=="SCHACU-115", seedmass_csv$Site[i] <- "SHAC",
    ifelse(seedmass_csv$Site[i]=="DIST-BE", seedmass_csv$Site[i] <- "DIST",
    ifelse(seedmass_csv$Species[i]=="ELPA", seedmass_csv$Site[i] <- "SHANE",
           print(seedmass_csv$Site[i])))))))))))))))))))))))))))))
}#for
# check that all sources that now match with the master dataset (all should except for unused Keith sources and Typha)
sort(unique(seedmass_csv$Site[which(!seedmass_csv$Site %in% germination_matrix$Source)]))
# Save edited dataset
write.csv(seedmass_csv, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedMass_ModelEdits_15DEC2020.csv', row.names=FALSE)

# Add seed mass values to master dataset based on source
germination_matrix$SeedMass <- rep(NA, nrow(germination_matrix))
for (src in unique(seedmass_csv$Site)){
  seedmass_subset <- seedmass_csv[which(seedmass_csv$Site==src),]
  seedmass <- mean(seedmass_subset$AvgSeedMass, na.rm=TRUE)
  germination_matrix$SeedMass[which(germination_matrix$Source==src)] <- seedmass
}#for
# check for missing values- there should not be any
germination_matrix$SeedMass[which(is.na(germination_matrix$SeedMass))]

# Clean up seed coat thickness dataset
sct_csv <- read.csv('C:/Users/Maggie/Downloads/RAW_DATA_2020.09_SCT_Dimensions.csv', header=TRUE, stringsAsFactors=FALSE)
names(sct_csv)[1] <- "Species"
# convert calculated average seed coat thickness to a numeric and check that missing values match with "divide by zero" error
sct_csv$AVG_SCT <- as.numeric(sct_csv$AVG_SCT)
which(is.na(sct_csv$AVG_SCT[1:1065])) #all good
# check for source names that do not match with master dataset
sort(unique(sct_csv$Source[which(!sct_csv$Source%in%germination_matrix$Source)]))
sort(unique(germination_matrix$Source))
# fix source names
for (i in 1:1065){
  ifelse(sct_csv$Source[i]=="Keith" & sct_csv$Species[i]=="JUBA", sct_csv$Source[i] <- "JUNB",
    ifelse(!is.na(sct_csv$Source[i]),
    ifelse(sct_csv$Source[i]=="JUNBAL-115", sct_csv$Source[i] <- "JUNB",
    ifelse(sct_csv$Source[i]=="RRVA-BW" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "RRVABW1",
    ifelse(sct_csv$Source[i]=="RRVA-BW" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "RRVABW2",
    ifelse(sct_csv$Source[i]=="RRVA-BW" & sct_csv$Species[i]=="SCAM", sct_csv$Source[i] <- "RRVABW3",
    ifelse(sct_csv$Source[i]=="PAHR" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "PAHR1",
    ifelse(sct_csv$Source[i]=="PAHR" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "PAHR2",
    ifelse(sct_csv$Source[i]=="FROU" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "FROU1",
    ifelse(sct_csv$Source[i]=="FROU" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "FROU2",
    ifelse(sct_csv$Source[i]=="WASP-MT", sct_csv$Source[i] <- "WASPMT",
    ifelse(sct_csv$Source[i]=="SACR" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "SACR1",
    ifelse(sct_csv$Source[i]=="SACR" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "SACR2",
    ifelse(sct_csv$Source[i]=="BERI" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "BERI1",
    ifelse(sct_csv$Source[i]=="BERI" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "BERI2",
    ifelse(sct_csv$Source[i]=="BERI" & sct_csv$Species[i]=="PHAU", sct_csv$Source[i] <- "BERI3",
    ifelse(sct_csv$Source[i]=="FABA" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "FABA1",
    ifelse(sct_csv$Source[i]=="FABA" & sct_csv$Species[i]=="PHAU", sct_csv$Source[i] <- "FABA3",
    ifelse(sct_csv$Source[i]=="FISP" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "FISP1",
    ifelse(sct_csv$Source[i]=="FISP" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "FISP2",
    ifelse(sct_csv$Source[i]=="FISP" & sct_csv$Species[i]=="SCAM", sct_csv$Source[i] <- "FISP3",
    ifelse(sct_csv$Source[i]=="CLLA" & sct_csv$Species[i]=="BOMA", sct_csv$Source[i] <- "CLLA1",
    ifelse(sct_csv$Source[i]=="CLLA" & sct_csv$Species[i]=="SCAC", sct_csv$Source[i] <- "CLLA2",
    ifelse(sct_csv$Source[i]=="WASP-UT", sct_csv$Source[i] <- "WASPUT",
    ifelse(sct_csv$Source[i]=="SCHAME-113", sct_csv$Source[i] <- "SHAM",
    ifelse(sct_csv$Source[i]=="SCHACU-115", sct_csv$Source[i] <- "SHAC",
    ifelse(sct_csv$Source[i]=="DIST-BE", sct_csv$Source[i] <- "DIST",
    ifelse(sct_csv$Species[i]=="ELPA", sct_csv$Source[i] <- "SHANE",
           print(sct_csv$Source[i])))))))))))))))))))))))))))))
}#for
# check that all sources that now match with the master dataset (all should except for unused Keith sources and Typha)
sort(unique(sct_csv$Source[which(!sct_csv$Source %in% germination_matrix$Source)]))
# Save edited dataset
write.csv(sct_csv, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/SeedDimensionsSCT_ModelEdits_15DEC2020.csv', row.names=FALSE)

# Add seed coat thickness to master dataset based on source
germination_matrix$SCT <- rep(NA, nrow(germination_matrix))
for (src in unique(sct_csv$Source)){
  sct_subset <- sct_csv[which(sct_csv$Source==src),]
  sct <- mean(sct_subset$AVG_SCT, na.rm=TRUE)
  germination_matrix$SCT[which(germination_matrix$Source==src)] <- sct
}#for
# check for missing values
germination_matrix[which(is.na(germination_matrix$SCT)),]#missing values for JUBA and SHAC (Keith's SCAC source)
# convert all missing values to same type of NA
germination_matrix$SCT[which(is.na(germination_matrix$SCT))] <- NA

# Save edited master germination dataset
write.csv(germination_matrix, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/GerminationData_FINAL_Processed_15DEC2020.csv', row.names=FALSE)

###################################################
# Convert Original Data to Individual Seed Format #
###################################################
# load cleaned up germination data
in_file <- 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/GerminationData_FINAL_Processed_15DEC2020.csv'
raw_germ_data <- read.csv(in_file, header=TRUE, stringsAsFactors=FALSE)

# set up empty data frame to save data to
expanded_df <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(raw_germ_data)+1))
names(expanded_df) <- c(names(raw_germ_data[1:6]), "Time0", names(raw_germ_data[7:53]))

# loop through each row of time to event data and expand counts to represent each seed
for (i in 1:nrow(raw_germ_data)){
  # set up data frame for cup
  row_df <- as.data.frame(matrix(1, nrow=raw_germ_data$NoSown[i], ncol=ncol(raw_germ_data)+1))
  names(row_df) <- names(expanded_df)
  # fill common values for cup
  row_df[,1:6] <- raw_germ_data[i,1:6]
  row_df[,44:54] <- raw_germ_data[i, 43:53]
  # add rows for each germination
  j=0 # start row counter
  for (t in 7:42){
    ngerm <- raw_germ_data[i,t]
    if (!is.na(ngerm) & ngerm > 0){
      k = j + ngerm
      row_df[(j+1):k,(t+1):43] <- 0 # transitioned = 0(modeled after survival/mortality analysis where 0=transitioned to death)
      j = k
    }#if
  ## The rest of the remaining rows will be all 1's: i.e., seeds that never transitioned to germination
  }#t
  # append data for cup to full data frame
  expanded_df <- rbind(expanded_df, row_df)
}#i

# save expanded germination matrix to new file
nrow(expanded_df); sum(raw_germ_data$NoSown) # CHECK: gained 5 "seeds" somewhere...
# find extra values...
expanded_df[which(is.na(expanded_df$Trial)),]
# remove extra values
expanded_df <- expanded_df[-which(is.na(expanded_df$Trial)),]
# save final matrix
write.csv(expanded_df, 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/GerminationMatrix_FINAL_15DEC2020.csv', row.names=FALSE)

# check for any missing values
names(expanded_df)
which(is.na(expanded_df$Trial))
which(is.na(expanded_df$CupNo))
which(is.na(expanded_df$Species))
which(is.na(expanded_df$Source))
which(is.na(expanded_df$WP))
which(is.na(expanded_df$Temp))
which(is.na(expanded_df[,7:43]))
which(is.na(expanded_df$SeedMass))
# these would need to be filled if SCT were added back into the model...
which(is.na(expanded_df$SCT)); expanded_df$Species[which(is.na(expanded_df$SCT))]
