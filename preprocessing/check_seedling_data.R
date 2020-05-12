#####################
## RUN DATA CHECKS ##
#####################

# load dependencies
library(pracma)

# specify file paths
in_file <- 'C:/Users/Maggie/Documents/WetlandEcology/ModelData/Processed_CSVs/Trial3_Germination_Processed.csv'
lookup_file <- 'C:/Users/Maggie/Documents/WetlandEcology/ModelData/Processed_CSVs/Seedling_Lookup_Table_CupNos.csv'
checked_file <-'C:/Users/Maggie/Documents/WetlandEcology/ModelData/Processed_CSVs/Trial3_Germination_Checked.csv'

# Run species and source categories check
## NOTE: This function assumes field names of 'Species','Source','WP', and 'Temp'
## It also assumes that 'WP' and 'Temp' are in numeric format
sources_checked <- check_field_values(in_file)
# make corrections as needed and rerun until no warnings show...
sources_checked$Source[141] <- 'THNA'
sources_checked[177,]
sources_checked$Source[177] <- 'FISP2'
sources_checked[179,]
sources_checked$Source[179] <- 'NIPI'
names(sources_checked) <- c("Trial", "CupNo", "Species", "Source", "WP", "Temp", "10/6/2019", "10/8/2019", "10/10/2019", "10/12/2019", "10/14/2019", "10/16/2019", "10/18/2019", "10/20/2019", "10/22/2019", "10/24/2019", "10/26/2019", "10/28/2019", "10/30/2019", "11/1/2019", "11/3/2019", "11/5/2019", "11/7/2019", "11/9/2019", "11/11/2019", "11/13/2019", "11/15/2019", "11/17/2019", "11/19/2019", "11/21/2019", "11/23/2019", "11/25/2019", "11/27/2019", "11/29/2019", "12/1/2019", "12/3/2019", "12/5/2019", "12/7/2019", "12/9/2019", "12/11/2019", "12/13/2019", "12/14/2019", "12/16/2019", "NoGermDate", "NoSown")
write.csv(sources_checked, checked_file)

# Run cup numbers check
## Note: This function assummes field names of 'Species', 'Source', 'CupNo', 'WP', and 'Temp'. It also assumes that 'Temp' is numeric 1-3 and 'WP' is numeric 1-2.
## Any blank or NA values in the in_file will cause this function to fail.
## Sources should already be corrected via check_field values before running this.
check_cup_numbers(checked_file, lookup_file)
## Once any discrepancies are fixed within the in_file, proceed...


#############################
### FIELDS CHECK FUNCTION ###
#############################

check_field_values <- function(in_file){
  # read in data and define species and source
  data <- read.csv(in_file, header=TRUE,stringsAsFactors = F)
  data$Species <- as.character(data$Species)
  data$Source <- as.character(data$Source)
  
  # define species and sources
  known_species <- c("BOMA", "DISP","ELPA", "JUBA", "PHAU", "SCAC", "SCAM", "TYSP")
  boma_sources <- c("ALK2", "BENLA", "BERI1", "BLHO", "CLLA1", "FABA1", "FISP1", "FROU1", "PAHR1", "RRVABW1", "SACR1", "WASPUT")
  scac_sources <- c("BERI2", "CLLA2", "CUMA", "FISP2", "FROU2", "KIWA", "MULA", "NIPI", "PAHR2", "PRBA", "RRVABW2", "SACR2", "SHAC", "THNA", "WASPMT")
  scam_sources <- c("HACR", "FISP3", "RRVABW3", "SASP", "SHAM")
  phau_sources <- c("BERI3", "FABA3", "OGBA")
  elpa_sources <- c("SHANE")
  disp_sources <- c("DIST")
  juba_sources <- c("JUNB")
  tysp_sources <- c("TATIANA")
  
  # check species
  for (i in 1:length(data$Species)){
    if(! data$Species[i] %in% known_species){
      print(paste("Species not in defined species: row ", as.character(i), " = '", data$Species[i], "'"))
    }# end of if
  }# end of for loop
  
  # remove "-" in any source names
  for(i in 1:nrow(data)) data$Source <- gsub("-", "", data$Source)
  
  # rename sources for each species
  for (i in 1:nrow(data)){
    # check & rename BOMA sources
    if (strcmp(data$Species[i],"BOMA")){
      if (strcmp(data$Source[i],"KEITH")) data$Source[i] <- "ALK2"
      if (strcmp(data$Source[i],"BERI")) data$Source[i] <- "BERI1"
      if (strcmp(data$Source[i],"CLLA")) data$Source[i] <- "CLLA1"
      if (strcmp(data$Source[i],"FABA")) data$Source[i] <- "FABA1"
      if (strcmp(data$Source[i],"FISP")) data$Source[i] <- "FISP1"
      if (strcmp(data$Source[i],"FROU")) data$Source[i] <- "FROU1"
      if (strcmp(data$Source[i],"PAHR")) data$Source[i] <- "PAHR1"
      if (strcmp(data$Source[i],"RRVABW")) data$Source[i] <- "RRVABW1"
      if (strcmp(data$Source[i],"SACR")) data$Source[i] <- "SACR1"
      if (strcmp(data$Source[i],"WASPMT") && data$CupNo[i] %in% c(60,71,100,116,194,196)) data$Source[i] <- "WASPUT"
      if (!data$Source[i] %in% boma_sources) print(paste("Source not in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename SCAC sources
    if (strcmp(data$Species[i],"SCAC")){
      if (strcmp(data$Source[i],"BERI")) data$Source[i] <- "BERI2"
      if (strcmp(data$Source[i],"CLLA")) data$Source[i] <- "CLLA2"
      if (strcmp(data$Source[i],"FISP")) data$Source[i] <- "FISP2"
      if (strcmp(data$Source[i],"FROU")) data$Source[i] <- "FROU2"
      if (strcmp(data$Source[i],"PAHR")) data$Source[i] <- "PAHR2"
      if (strcmp(data$Source[i],"RRVABW")) data$Source[i] <- "RRVABW2"
      if (strcmp(data$Source[i],"SACR")) data$Source[i] <- "SACR2"
      if (strcmp(data$Source[i],"KEITH")) data$Source[i] <- "SHAC"
      if (strcmp(data$Source[i],"WASPUT") && data$CupNo[i] %in% c(10,29,97,110,184,211)) data$Source[i] <- "WASPMT"
      if (!data$Source[i] %in% scac_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    } #end of if
    # check & rename SCAM sources
    if (strcmp(data$Species[i],"SCAM")){
      if (strcmp(data$Source[i],"FISP")) data$Source[i] <- "FISP3"
      if (strcmp(data$Source[i],"RRVABW")) data$Source[i] <- "RRVABW3"
      if (strcmp(data$Source[i],"KEITH")) data$Source[i] <- "SHAM"
      if (!data$Source[i] %in% scam_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename PHAU sources
    if (strcmp(data$Species[i],"PHAU")){
      if (strcmp(data$Source[i],"BERI")) data$Source[i] <- "BERI3"
      if (strcmp(data$Source[i],"FABA")) data$Source[i] <- "FABA3"
      if (!data$Source[i] %in% phau_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename ELPA sources
    if (strcmp(data$Species[i],"ELPA")){
      if (strcmp(data$Source[i],"THNA")) data$Source[i] <- "SHANE"
      if (!data$Source[i] %in% elpa_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename DISP sources
    if (strcmp(data$Species[i],"DISP")){
      if (strcmp(data$Source[i],"KEITH")) data$Source[i] <- "DIST"
      if (!data$Source[i] %in% disp_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename JUBA sources
    if (data$Species[i]=="JUBA"){
      if (data$Source[i]=="KEITH") data$Source[i] <- "JUNB"
      if (!data$Source[i] %in% juba_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    # check & rename TYSP sources
    if (data$Species[i]=="TYSP"){
      if (data$Source[i]=="-") data$Source[i] <- "TATIANA"
      if (!data$Source[i] %in% tysp_sources) print(paste("Source not found in defined list: row", as.character(i), "= '", data$Source[i], "'"))
    }# end of if
    
  # check WP data
    if(!data$WP[i] %in% c(1,2)) print(paste("Check water potential value: row",as.character(i)))
    
    # check temp data
    if(!data$Temp[i] %in% c(1,2,3)) print(paste("Check water potential value: row",as.character(i)))
  }#end of for loop
  return(data)
}#end of function


#############################
## CHECK CUP NUMBERS MATCH ##
## WITH LOOKUP TABLE       ##
#############################

check_cup_numbers <- function(data_file, lookup_file){
  # read in files
  lookup_table <- read.csv(lookup_file, header=T,stringsAsFactors = F)
  data <- read.csv(data_file, header=T,stringsAsFactors = F)
  names(lookup_table) <- c("CupNo", "Species", "Source", "Temp", "WP")
  
  # set field types
  lookup_table$CupNo <- as.integer(lookup_table$CupNo)
  lookup_table$Species <- as.character(lookup_table$Species)
  lookup_table$Source <- as.character(lookup_table$Source)
  lookup_table$Temp <- as.integer(lookup_table$Temp)
  lookup_table$WP <- as.integer(lookup_table$WP)
  data$CupNo <- as.integer(data$CupNo)
  data$Species <- as.character(data$Species)
  data$Source <- as.character(data$Source)
  data$Temp <- as.integer(data$Temp)
  data$WP <- as.integer(data$WP)

  # check that all values match for cup number
  for (i in 1:nrow(data)){
    cupNo <- data$CupNo[i]
    lookup <- lookup_table[which(lookup_table$CupNo==cupNo),]
    if (!strcmp(data$Species[i],lookup$Species)) print(paste("Species doesn't match for row", as.character(i)))
    if (!strcmp(data$Source[i],lookup$Source)) print(paste("Source doesn't match for row", as.character(i)))
    if (data$Temp[i]!=lookup$Temp) print(paste("Temp doesn't match for row", as.character(i)))
    if (data$WP[i]!=lookup$WP) print(paste("WP doesn't match for row", as.character(i)))
  }#end for loop
} #end function

