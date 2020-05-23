##### Dry Weight Regression analysis
#### Before uploading your data sets, prepare the data by standardizing the column names and verifying
#### NA, <0.001, and 0 cells.Remove any extra spaces from the source column.
#### Important date and time formatting. I ran into the problem where all dates from 2019 were recorded
#### as 2020, so this must be changed. The dates must be formatted mm-dd-yyyy. Time must be hh:mm.
#### Imporant header names to change before importing the data are "cup", "species", "wp", "temp",  "weigh_time",
#### "pull_date", "pull_time", "ag_fresh_mass_mg", "ag_dry_mass_mg", "bg_fresh_mass_mg", and "bg_dry_mass_mg".

### Load needed packages. 
library (tidyr)
library(lubridate)
library(tidyverse)
library(ggplot2)
library(ggiraphExtra)

### Set working directory. (Change to the location of your files.)
setwd("C:/Users/Maggie/Downloads")
getwd()

### Read dry weight, fresh weight, and cup number datasets.(Change to the names of your files.)
dry_mass <- read.csv("Trial_4_Dry_Mass_Data_Sheet_Entry_edited4_3_2020.csv", header = T)
fresh_mass <- read.csv("T4_Fresh_Mass_Data_Entry_dlwd5_1.csv", header = T)
cup_number <- read.csv("CupNo.csv", header = T)

#### DATA MANIPULATION ####
### Create a unique identifier for each row for both fresh and dry mass. The selected columns to combine
### should be cup, source, and seedling number which is represented by 3:5 in the dry mass dataset and 
### 2:3,5 in the fresh mass dataset.(Change to the correct columns in your dataset if necessary.)
dry_mass <- unite(dry_mass, "seedling_id", 3:5, sep = "", remove = F)
fresh_mass <- unite(fresh_mass, "seedling_id", 2:3,5, sep = "", remove = F)

### Add species, wp, and temp variables to fresh mass dataset from the cup number dataset. 
### (Change column names if different.)
fresh_mass <- merge(fresh_mass, cup_number[, c("cup", "species", "wp", "temp")], by = "cup", all.x = T)

### Calculate the difference in minutes between the the pull and weigh time. (You will need to check 
### that both times are on the same date... i.e. if weighing went until 1 am the next day, you will 
### need to create a weigh date column. This will most likely need to be done when all fresh weight 
### data is used. Date also must be formatte as month/day/year.)

## Create a variable that combines date and time for both weigh and pull time.
pull_datetime <- with(fresh_mass, mdy(pull_date) + hm(pull_time))
weigh_datetime <- with(fresh_mass, mdy(pull_date) + hm(weigh_time))

## Change the time zone to Denver.
pull_datetime <- force_tz(pull_datetime, tzone = "America/Denver")
weigh_datetime <- force_tz(weigh_datetime, tzone = "America/Denver")

## Calculate difference between weigh and pull time.
minutes <- weigh_datetime - pull_datetime

## Check there are no negative values.Fix or delete seedlings with incorrect pull times in original
## csv document. You Will need to run script again if you made any changes to the csv.
minutes
na_minutes <- which(is.na(minutes))
minutes <- minutes[-na_minutes]

## Since this variable type is now difftime, ggplot will not be able to use it,
## so convert it to numeric.
minutes <- as.numeric(minutes,units = "mins")
## Check it is numeric
str(minutes)

### Add all variables to the merged data frame.
fresh_mass <- cbind(fresh_mass, pull_datetime, weigh_datetime, minutes)

### Only keep the observations with a positive time variable
fresh_mass <- subset(fresh_mass, minutes > 0)

### Combine fresh mass data with dry mass data, keeping all dry mass and fresh mass rows that match
### and all columns from both data sets.(You should be able to change all.x to TRUE once all fresh 
### weight data has been entered and corrected.I still have it false because there is still some missing
### fresh mass data.) 
dry_fresh_merge <- merge(dry_mass, fresh_mass, by = "seedling_id", all.x = F, no.dups = T)


#### Plot the data to see visual patterns for aboveground data ####
### fresh and minutes compared to dry
ag_plot_min <- ggplot(dry_fresh_merge,aes(y=ag_dry_mass_mg,x=ag_fresh_mass_mg,
                color=minutes))+geom_point()+stat_smooth(method="lm",
                se=F)+ggtitle("Above Ground"); ag_plot_min
### fresh and species compared to dry
ag_plot_species <- ggplot(dry_fresh_merge,aes(y=ag_dry_mass_mg,x=ag_fresh_mass_mg,
                color=species))+geom_point()+stat_smooth(method="lm",
                se=F)+ggtitle("Above Ground"); ag_plot_species
### fresh and temperature compared to dry
ag_plot_temp <- ggplot(dry_fresh_merge,aes(y=ag_dry_mass_mg,x=ag_fresh_mass_mg,
                color=temp))+geom_point()+stat_smooth(method="lm",
                se=F)+ggtitle("Above Ground"); ag_plot_temp
### fresh and wp compared to dry
ag_plot_wp <- ggplot(dry_fresh_merge,aes(y=ag_dry_mass_mg,x=ag_fresh_mass_mg,
               color=wp))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Above Ground"); ag_plot_wp
### fresh and source compared to dry
ag_plot_source.x <- ggplot(dry_fresh_merge,aes(y=ag_dry_mass_mg,x=ag_fresh_mass_mg,
               color=source.x))+geom_point()+stat_smooth(method="lm",
                se=F)+ggtitle("Above Ground"); ag_plot_source.x

#### Plot the data to see visual patterns for belowground data
### fresh and minutes compared to dry
bg_plot_min <- ggplot(dry_fresh_merge,aes(y=bg_dry_mass_mg,x=bg_fresh_mass_mg,
               color=minutes))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Below Ground"); bg_plot_min
### fresh and species compared to dry
bg_plot_species <- ggplot(dry_fresh_merge,aes(y=bg_dry_mass_mg,x=bg_fresh_mass_mg,
               color=species))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Below Ground"); bg_plot_species
### fresh and temp compared to dry
bg_plot_temp <- ggplot(dry_fresh_merge,aes(y=bg_dry_mass_mg,x=bg_fresh_mass_mg,
               color=temp))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Below Ground"); bg_plot_temp
### fresh and wp compared to dry
bg_plot_wp <- ggplot(dry_fresh_merge,aes(y=bg_dry_mass_mg,x=bg_fresh_mass_mg,
               color=wp))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Below Ground"); bg_plot_wp
### fresh and source compared to dry
bg_plot_source.x <- ggplot(dry_fresh_merge,aes(y=bg_dry_mass_mg,x=bg_fresh_mass_mg,
               color=source.x))+geom_point()+stat_smooth(method="lm",
               se=F)+ggtitle("Below Ground"); bg_plot_source.x

#### MODEL SELECTION ####
#### Above Ground multiple linear regression (all combinations of speceies,temp, wp, and source.x)
#### source.x has this name because during previous merge, there were two source columns.
### According to se, R-squared, and AIC values, I think model ag2 is the best
"""
ag1 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + temp + wp + source.x, 
         data = dry_fresh_merge); ag1
summary(ag1)

ag2 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + temp + wp, 
                  data = dry_fresh_merge); ag2
summary(ag2) ### residual standard error = 0.38888; multiple R-squared: 0.8873

ag3 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + temp + source.x, 
          data = dry_fresh_merge); ag3
summary(ag3)

ag4 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + wp + source.x, 
          data = dry_fresh_merge); ag4
summary(ag4)

ag5 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + temp + wp + source.x, 
          data = dry_fresh_merge); ag5
summary(ag5)

ag6 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + temp, 
          data = dry_fresh_merge); ag6
summary(ag6)

ag7 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + wp, 
          data = dry_fresh_merge); ag7
summary(ag7)

ag8 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + source.x, 
          data = dry_fresh_merge); ag8
summary(ag8)

ag9 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + temp + wp, 
          data = dry_fresh_merge); ag9
summary(ag9)

ag10 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + temp + source.x, 
          data = dry_fresh_merge); ag10
summary(ag10)

ag11 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + wp + source.x, 
          data = dry_fresh_merge); ag11
summary(ag11)

ag12 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species, 
          data = dry_fresh_merge); ag12
summary(ag12)

ag13 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + temp, 
          data = dry_fresh_merge); ag13
summary(ag13)

ag14 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + wp, 
          data = dry_fresh_merge); ag14
summary(ag14)

ag15 <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + source.x, 
          data = dry_fresh_merge); ag15
summary(ag15)

### Determinie which model has the better fit. ag_nosource is the better fit
AIC(ag1) 
AIC(ag2) ## Best fit
AIC(ag3)
AIC(ag4)
AIC(ag5)
AIC(ag6)
AIC(ag7)
AIC(ag8)
AIC(ag9)
AIC(ag10)
AIC(ag11)
AIC(ag12)
AIC(ag13)
AIC(ag14)
AIC(ag15)

#### Below ground multiple linear regression (all combinations of speceies,temp, wp, and source)
### According to se, R-squared, and AIC values, I think model ag2 is the best

bg1 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + temp + wp + source.x, 
          data = dry_fresh_merge); bg1
summary(bg1)

bg2 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + temp + wp, 
          data = dry_fresh_merge); bg2
summary(bg2) ### Residual standard error = 0.1618; Multiple R-squared = 0.8069

bg3 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + temp + source.x, 
          data = dry_fresh_merge); bg3
summary(bg3)

bg4 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + wp + source.x, 
          data = dry_fresh_merge); bg4
summary(bg4)

bg5 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp + wp + source.x, 
          data = dry_fresh_merge); bg5
summary(bg5)

bg6 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + temp, 
          data = dry_fresh_merge); bg6
summary(bg6)

bg7 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + wp, 
          data = dry_fresh_merge); bg7
summary(bg7)

bg8 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + source.x, 
          data = dry_fresh_merge); bg8
summary(bg8)

bg9 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp + wp, 
          data = dry_fresh_merge); bg9
summary(bg9)

bg10 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp + source.x, 
           data = dry_fresh_merge); bg10
summary(bg10)

bg11 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + wp + source.x, 
           data = dry_fresh_merge); bg11
summary(bg11)

bg12 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species, 
           data = dry_fresh_merge); bg12
summary(bg12)

bg13 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp, 
           data = dry_fresh_merge); bg13
summary(bg13)

bg14 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + wp, 
           data = dry_fresh_merge); bg14
summary(bg14)

bg15 <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + source.x, 
           data = dry_fresh_merge); bg15
summary(bg15)

### Determinie which model has the better fit. bg_nosource is the better fit
AIC(bg1) 
AIC(bg2) ## best fit 
AIC(bg3)
AIC(bg4)
AIC(bg5)
AIC(bg6)
AIC(bg7)
AIC(bg8)
AIC(bg9)
AIC(bg10)
AIC(bg11)
AIC(bg12)
AIC(bg13)
AIC(bg14)
AIC(bg15)
"""

#### BACKWARD SELECTION
# backward selection - aboveground
ag_global_model <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + temp + wp + source.x, data = dry_fresh_merge)
step(ag_global_model) #best: ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + wp
# backward selection - aboveground
bg_global_model <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + species + temp + wp + source.x, data = dry_fresh_merge[-c(212, 140),])
step(bg_global_model) #best: bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp + wp

#### REGRESSION DIAGNOSTICS ####
# run models & diagnostics
ag_model <- lm(ag_dry_mass_mg ~ ag_fresh_mass_mg + minutes + species + wp, data=dry_fresh_merge)
bg_model <- lm(bg_dry_mass_mg ~ bg_fresh_mass_mg + minutes + temp + wp, data=dry_fresh_merge[-c(212,140),])
# basic residual diagnostic plots
plot(ag_model)
plot(bg_model) # based on qq-plot, looks overdispersed but I'm not sure how to handle- maybe predictions will shake out better with more data?
# check for outliers
leveragePlots(ag_model)
leveragePlots(bg_model)
influencePlot(ag_model)
influencePlot(bg_model) # rerun without observations 212 and 140
# test homoscedasticity
ncvTest(ag_model)
ncvTest(bg_model)
spreadLevelPlot(ag_model)
spreadLevelPlot(bg_model)
# check nonlinearity
crPlots(ag_model)
crPlots(bg_model)
ceresPlots(ag_model)
ceresPlots(bg_model)
# check for nonindependence
durbinWatsonTest(ag_model)
durbinWatsonTest(bg_model)

#### PREDICTED VS. OBSERVED ####
plot(x=predict(ag_model), y=dry_fresh_merge$ag_dry_mass_mg, xlab="predicted", ylab="actual"); abline(a=0, b=1, col="red")
plot(x=predict(bg_model), y=dry_fresh_merge$bg_dry_mass_mg[-c(212,140, which(is.na(dry_fresh_merge$bg_fresh_mass_mg)))], xlab="predicted", ylab="actual"); abline(a=0, b=1, col="red")

#### Predict the dry weights based on the best-fittting model. The data you use must have the correct
#### date (m/d/y) and time (h:m) format to be prepared correctly. For whatever data you are analyzing,
#### the data frame must include these fields with these exact column names to be compatible with the 
#### ag2 and bg2 models: "ag_fresh_mass_mg", "bg_fresh_mass_mg", "minutes", "species", "temp", "wp"
#### "pull_date", "pull_time", and "weigh_time". 

### Import data with fresh weights to be predicted.You will need the cup_number data from the beginning 
### of the analysis as well.
new_fresh_mass <- read.csv("T5_Fresh_Mass_Data_Entry_edited5_18_20.csv", header = T)

#### Add the "minutes", "species", "temp", and "wp" columns to the data.
### Add species, wp, and temp variables to fresh mass dataset from the cup number dataset. 
### (Change column names if different.)
new_fresh_mass <- merge(new_fresh_mass, cup_number[, c("cup", "species", "wp", "temp")], by = "cup", all.x = T)

### Calculate the difference in minutes between the the pull and weigh time. (You will need to check 
### that both times are on the same date... i.e. if weighing went until 1 am the next day, you will 
### need to create a weigh date column. This will most likely need to be done when all fresh weight 
### data is used. Date also must be formatte as month/day/year.)

## Create a variable that combines date and time for both weigh and pull time.
pull_datetime <- with(new_fresh_mass, mdy(pull_date) + hm(pull_time))
weigh_datetime <- with(new_fresh_mass, mdy(pull_date) + hm(weigh_time))

## Change the time zone to Denver.
pull_datetime <- force_tz(pull_datetime, tzone = "America/Denver")
weigh_datetime <- force_tz(weigh_datetime, tzone = "America/Denver")

## Calculate difference between weigh and pull time.
minutes <- weigh_datetime - pull_datetime

## Check there are no negative values.Fix or delete seedlings with incorrect pull times in original
## csv document. You Will need to run script again if you made any changes to the csv.
minutes

## Since this variable type is now difftime, convert it to numeric.
minutes <- as.numeric(minutes,units = "mins")
## Check it is numeric
str(minutes)

### Add all variables to the merged data frame.
new_fresh_mass <- cbind(new_fresh_mass, pull_datetime, weigh_datetime, minutes)

### Only keep the observations with a positive time variable
fresh_mass <- subset(fresh_mass, minutes > 0)

#### Predict above and below ground dry weights based on the variables this new data frame.

### Above ground predicted dry weights with standard error and confidence interval level 95%.
ag_predict <- predict(ag2, new_fresh_mass, se.fit=T); ag_predict
### Below ground predicted dry weights with standard error and confidence interval level 95%.
bg_predict <- predict(bg2, new_fresh_mass, se.fit=T); bg_predict

### Add these predicted values and the se to the new data frame.
new_fresh_mass$predict_ag_dry <- ag_predict$fit
new_fresh_mass$se_ag_dry <- ag_predict$se.fit
new_fresh_mass$predict_bg_dry <- bg_predict$fit
new_fresh_mass$se_bg_dry <- bg_predict$se.fit

### Plot the predicted dry values against the fresh mass with error bars (se).
## AG plot
predicted_ag_dry_plot <- ggplot(new_fresh_mass ,aes(y=predict_ag_dry,x=ag_fresh_mass_mg))+
  geom_point()+
  geom_errorbar(aes(ymin=predict_ag_dry-se_ag_dry, ymax=predict_ag_dry+se_ag_dry), width=1,
                position=position_dodge(0.05))+
  ggtitle("Above Ground Predicted dry"); predicted_ag_dry_plot

## BG plot
predicted_bg_dry_plot <- ggplot(new_fresh_mass ,aes(y=predict_bg_dry,x=bg_fresh_mass_mg))+
  geom_point()+
  geom_errorbar(aes(ymin=predict_bg_dry-se_bg_dry, ymax=predict_bg_dry+se_bg_dry), width=1,
                position=position_dodge(0.05))+
  ggtitle("Below Ground Predicted dry"); predicted_bg_dry_plot

### Export data (Change the file location for your computer.)
write.csv(new_fresh_mass,"C:/Users/audre/Documents/WetlandEcologyLab/EMILY/Dry_vs_wet_weight\\T5FreshAndPredictedDryWeights.csv", row.names = FALSE)
