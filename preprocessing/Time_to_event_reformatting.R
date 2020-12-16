##############################################
## Reformats germination/mortality data in  ##
## long format to short tabular format with ##
## one date-count record per line           ##
##############################################

# load dependencies
library(timeDate)

####################################
## Run time to event reformatting ##
####################################
in_file <- 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs/GerminationData_FINAL_Processed_15DEC2020.csv'
out_file <- 'C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ModelData/Processed_CSVs//TimeToEvent_Germination_FINAL_15DEC2020.csv'
reformat_germ_TTE_data(in_file, out_file)

##############
## FUNCTION ##
##############
reformat_germ_TTE_data <- function(in_file, out_file){
  #NOTE: ASSUMES THAT THE FIRST ENTRY FOR EACH CUP NUMBER HAS THE APPROPRIATE 
  #SPECIES, SOURCE, TEMP, WP, AND GERMINATION DATA ASSOCIATED. 
  # Checks for unique entries per cup
  #param in_file: filepath for formatted input germination data
  #param out_fle: filepath for formatted output reformatted time-to-event germination data
  #returns germination data formatted in time-to-event style for HTT model entry
  
  # read in original dataset and set standard names
  germ_original <- read.csv(in_file, header=T, stringsAsFactors=FALSE)
  
  # rename columns
  ncols = ncol(germ_original)
  names(germ_original)[1:6] <- c("Trial", "CupNo", "Species", "Source", "WP", "Temp")
  names(germ_original)[(ncols-1):ncols] <- c("NoGermDate", "NoSown")

  # remove any empty rows
  for (i in 1:length(germ_original[,1])){
    if(is.na(germ_original$CupNo[i])) germ_original <- germ_original[-i,]
  }
  
  # Check that there are no repeat cup numbers
  if(length(unique(germ_original$CupNo)) != length(germ_original$CupNo))
  {stop("ERROR! Cup numbers are not unique. Double check cup numbers in input file.")}

  # Check that all cup numbers are represented
  tysp_numbers <- c(34, 76, 120, 144, 207, 221)
  cup_numbers <- germ_original$CupNo
  for (i in 1:702){
    if (! i %in% cup_numbers){
      if (i %in% c(34, 76, 120, 144, 207, 221))
      {print(paste("WARNING: Cup number is missing from data set: #", as.character(i), "- Typha cup"))}
      else {print(paste("WARNING: Cup number is missing from data set: #", as.character(i)))}
    }#close if
  }#close for loop
    
  # set data types for each field
  germ_original$CupNo <- as.character(germ_original$CupNo)
  germ_original$Species <- as.character(germ_original$Species)
  germ_original$Source <- as.character(germ_original$Source)
  germ_original$WP <- as.character(germ_original$WP)
  germ_original$Temp <- as.character(germ_original$Temp)
  germ_original$NoGermDate <- as.integer(germ_original$NoGermDate)
  germ_original$NoSown <- as.numeric(germ_original$NoSown)
  germ_original$Trial <- as.character(germ_original$Trial)
  
  # rename date fields based on timestep
  first_time_col <- 7
  last_time_col <- ncols-11
  for (i in 1:(last_time_col-first_time_col+1)){
    field_name <- paste("Time", as.character(i*2), sep="")
    names(germ_original)[i+6] <- field_name
    suppressWarnings(germ_original[,i+6] <- as.numeric(germ_original[,i+6]))
  }

  # find the timeframe
  timesteps <- (last_time_col - first_time_col)
  # create new dataset and set column names
  new_germ_df <- as.data.frame(matrix(NA, ncol=14, nrow=((timesteps+1)*702)))
  names(new_germ_df) <- c("Trial", "CupNo", "Species", "Source", "Temp", "WP", "Chamber", "BeginTime", "EndTime", "SeedMass", "SCT", "NoGerm", "PropGerm", "CumPropGerm")
  # set beginning and end times for each interval
  new_germ_df$BeginTime <- rep(seq(0, timesteps*2, by=2), 234)
  new_germ_df$EndTime <- rep(seq(2, (timesteps*2)+2, by=2), 234)
  new_germ_df$EndTime[which(new_germ_df$EndTime==(timesteps*2)+2)] <- Inf

  # calculate total germinations
  for (i in 1:nrow(germ_original)){
    germ_numeric <- gsub("END", "", germ_original[i,first_time_col:last_time_col])
    suppressWarnings(germ_numeric <- as.numeric(germ_numeric))
    germ_original$TotalGerm[i] <- sum(germ_numeric,na.rm=T) + germ_original$NoGermDate[i]
  }#end for loop
  
  # copy over fields attached to cup numbers
  for (i in 1:702){
    # get entry for cup number from original data frame
    cup_data <- germ_original[which(germ_original$CupNo==i),]
    # if cup_data not found, skip
    ifelse(nrow(cup_data)==0,
           print(paste("COULD NOT FIND DATA FOR CUP NUMBER",as.character(i))),
           {
    # get start and end row for each cup
    start = (timesteps+1)*(i-1)+1
    end = (timesteps+1)*i
    
    # set cup number
    new_germ_df$CupNo[start:end] <- i
    # set species
    cup_spec <- cup_data$Species
    new_germ_df$Species[start:end] <- cup_spec
    # set source
    cup_source <- cup_data$Source
    new_germ_df$Source[start:end] <- cup_source
    # set water potential
    cup_wp <- cup_data$WP
    new_germ_df$WP[start:end] <- cup_wp
    # set temp
    cup_temp <- cup_data$Temp
    new_germ_df$Temp[start:end]  <- cup_temp
    # set chamber
    cup_chm <- cup_data$Chamber
    new_germ_df$Chamber[start:end] <- cup_chm
    # set seed mass
    cup_seedmass <- cup_data$SeedMass
    new_germ_df$SeedMass[start:end] <- cup_seedmass
    # set sct
    cup_sct <- cup_data$SCT
    new_germ_df$SCT[start:end] <- cup_sct
    
    ## find and copy over germination number data for each cup number and time step
    # extract cup germ data from cup entry
    cup_germ_data <- cup_data[first_time_col:last_time_col]
    
    # find end timstep of cup data collection
    end_timestep <- which(cup_germ_data=="END")
    # if no end date, set to last date
    if (length(end_timestep)==0) end_timestep <- length(cup_germ_data)
    # calculate which row will have the end date in new data frame
    new_row_end <- end_timestep+start-1
    
    # reset end time based on data collection end time
    new_germ_df$EndTime[new_row_end] <- Inf

    # set germ numbers for all timesteps before data collection end for cup
    new_germ_df$NoGerm[start:(new_row_end-1)] <- cup_germ_data[1:(end_timestep-1)]
    # set any NAs to 0
    new_germ_df$NoGerm[start:(new_row_end-1)] <- as.numeric(new_germ_df$NoGerm[start:(new_row_end-1)])
    for (j in start:(new_row_end-1)){
      if (is.na(new_germ_df$NoGerm[j])) new_germ_df$NoGerm[j] <- 0
    }#for
    
    # set germ number to remaining seeds for end date
    new_germ_df$NoGerm[new_row_end] <- (cup_data$NoSown-cup_data$NoGermDate) - sum(as.data.frame(new_germ_df$NoGerm[start:(new_row_end-1)]),na.rm=T)
    # any rows above end date of data collection will have "NA" value by default
    
    # set germination proportions per day
    new_germ_df$PropGerm[start:(new_row_end-1)] <- as.data.frame(new_germ_df$NoGerm[start:(new_row_end-1)]) / (cup_data$NoSown-cup_data$NoGermDate)
  
    }) # close of ifelse statement
  } # close of for loop
  
  ## CALCULATE CUMULATIVE PROPORTIONS
  for (i in 1:702){
    # get start and end row for each cup
    start = (timesteps+1)*(i-1)+1
    end = (timesteps+1)*i

    # find end timstep of cup data collection
    end_timestep <- which(cup_germ_data=="END")
    # if no end date, set to last date
    if (length(end_timestep)==0) end_timestep <- length(cup_germ_data)
    # calculate which row will have the end date in new data frame
    new_row_end <- end_timestep+start-1
    
    # first cumulative proportion is equal to first proportion
    new_germ_df$CumPropGerm[start] <- as.data.frame(new_germ_df$PropGerm[start])
    # rest of cumulative proportions are sum of previous plus current row
    for (j in (start+1):(new_row_end-1)){
      previous_cumulative <- as.data.frame(new_germ_df$CumPropGerm[j-1])
      current_raw_prop <- as.data.frame(new_germ_df$PropGerm[j])
      new_germ_df$CumPropGerm[j] <- previous_cumulative+current_raw_prop  
    }
    # end of data collection and onward is NA by default
  } # end of for loop
  
  # delete any data beyond end date
  new_germ_df <- new_germ_df[!is.na(new_germ_df$NoGerm),]

  # double check all fields are vectors and not lists
  new_germ_df$NoGerm <- as.numeric(new_germ_df$NoGerm)
  new_germ_df$PropGerm <- as.numeric(new_germ_df$PropGerm)
  new_germ_df$CumPropGerm <- as.numeric(new_germ_df$CumPropGerm)
  
  # save new data frame to specified output file
  write.csv(new_germ_df, out_file, row.names = FALSE)
  
# if you want the function to also return the data frame,
# for example for troubleshooting, uncomment the following line:
  #return(new_germ_df) 
} # end of function


