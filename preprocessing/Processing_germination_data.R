##########################################################################
##   PROCESSES GERMINATION DATA                                         ##
## Adds germ synchrony, mean germ time, and germination proportion      ##
##########################################################################


##########
## Run ###
##########
# load dependencies
library(dplyr)
library(GerminaR)

# set file paths - must be CSV
in_file <- 'C:/Users/Maggie/Documents/WetlandEcology/Seedling_Trial_2_Germinations.csv'
out_file <- 'C:/Users/Maggie/Documents/WetlandEcology/Trial2_Processed_GermData.csv'
process_germ_data(in_file, out_file)



##############
## Function ##
##############
process_germ_data <- function(in_file, out_fle){
  # read in file
  germ <- read.csv(in_file,header=T,stringsAsFactors = F)

  # set column breaks
  # rename columns
  ncols = ncol(germ)
  names(germ)[1:6] <- c("Trial", "CupNo", "Species", "Source", "WP", "Temp")
  names(germ)[(ncols-1):ncols] <- c("NoGermDate", "NoSown")
  
  # set data types for each field
  germ$CupNo <- as.factor(germ$CupNo)
  germ$Species <- as.factor(germ$Species)
  germ$Source <- as.factor(germ$Source)
  germ$WP <- as.factor(germ$WP)
  germ$Temp <- as.factor(germ$Temp)
  germ$NoGermDate <- as.numeric(germ$NoGermDate)
  germ$NoSown <- as.numeric(germ$NoSown)
  germ$Trial <- as.factor(germ$Trial)
  
  # remove any empty rows
  for (i in 1:length(germ[,1])){
    if(is.na(germ$CupNo[i])) germ <- germ[-i,]
  }#end for loop
  
  # Check that there are no repeat cup numbers
  if(length(unique(germ$CupNo)) != length(germ$CupNo))
  {stop("ERROR! Cup numbers are not unique. Double check cup numbers in input file.")}
  
  # Check that all cup numbers are represented
  tysp_numbers <- c(34, 76, 120, 144, 207, 221)
  cup_numbers <- germ$CupNo
  for (i in 1:234){
    if (! i %in% cup_numbers){
      if (i %in% c(34, 76, 120, 144, 207, 221))
      {warning(paste("WARNING: Cup number is missing from data set: #", as.character(i), "- Typha cup"))}
      else {warning(paste("WARNING: Cup number is missing from data set: #", as.character(i)))}
    }#close if
  }#close for loop
  
  # rename date fields based on timestep
  first_time_col <- 7
  last_time_col <- ncols-2
  for (i in 1:(last_time_col-first_time_col+1)){
    field_name <- paste("Time", as.character(i*2), sep="")
    names(germ)[i+6] <- field_name
    germ[,i+6] <- gsub("END", 0, germ[,i+6])
    germ[,i+6] <- as.numeric(germ[,i+6])
  }#end for loop
  
  
  # calculate total germinations
  for (i in 1:nrow(germ)){
    germ_numeric <- germ[i,first_time_col:last_time_col]
    ifelse(is.na(germ$NoGermDate[i]), missing<-0, missing<-germ$NoGermDate[i])
    germ$TotalGerm[i] <- sum(germ_numeric, na.rm=T) + missing
  }#end for loop
  
  # calculate mean germination rate
  germ$GRP <- germ$TotalGerm / germ$NoSown
  
  # calculate germination synchrony
  for (i in 1:nrow(germ)){
    # make list of germination days for cup
    germ_row <- germ[i,first_time_col:last_time_col]
    list <- c()
    for(j in 1:length(germ_row)){
      ifelse(is.na(germ_row[j]), count<-0, count<-germ_row[j])
      reps <- rep((j*2-1), count)
      if (length(reps)>0) list <- c(list, reps)
    }#end for loop
    # calculate percentiles
    q <- quantile(list, c(0.1,0.9))
    # synchrony: # days between when 10% and 90% of seeds have germinated
    germ$GermSynch[i] <- q[2] - q[1]
  }#end for loop
  
  # calcualate mean germination time
  timesteps <- last_time_col-first_time_col
  # set number of days for each data collection data: sequence from 1, 3, 5, 7...
  multipliers <- as.matrix(seq(from=1, to=(timesteps*2)+1, by=2),dim(37,1))
  germ_days <- as.matrix(germ[,first_time_col:last_time_col])
  germ_days[is.na(germ_days)] <- 0 #set NAs to zero
  # matrix math: # of germinations multiplied by germination day, all summed together per cup
  sum_days <- germ_days %*% multipliers
  germ$MeanGermTime <- sum_days / germ$TotalGerm

}#end function