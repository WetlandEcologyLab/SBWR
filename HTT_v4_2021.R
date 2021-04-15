#### Hydrothermal time model for seed germination
#### Water potentials = 0 and -0.6
#### Temperatures = 28/10, 32/15, 36/20
#### See this website as a good reference: https://rpubs.com/mbmesgaran/495385
#rm(list=ls())

### Load packages
library(devtools)
library(tidyverse)
library(drc)
library(ggplot2)
library(gridExtra)
library(drc) #library(drcSeedGerm)
library(reshape2)
library(cowplot)
library(timeDate)
library(knitr)

### Set working directory
#setwd("C:/Users/Maggie/Downloads/2021_03_TimeTevent_Germ_reformat.csv")

### Import data & set data types
germ <- read.csv("C:/Users/Maggie/Downloads/2021_03_TimeTevent_Germ_reformat.csv", as.is=T)
germ$Trial <- as.factor(germ$Trial)
germ$Chamber <- as.factor(germ$Chamber)
germ$CupNo <- as.factor(germ$CupNo)
germ$Species <- as.factor(germ$Species)
germ$Source <- as.factor(germ$Source)
str(germ)

### Change temp and WP to original numeric values
germ <-  germ %>%
  mutate(WP = ifelse(WP=="1",0,-0.6))
with(germ,table(WP))

germ <-  germ %>%
  mutate(Temp = ifelse(Temp=="1",28,
                        ifelse(Temp=="2",32,36)))

### Inspect data
with(germ,table(Chamber))
with(germ,table(CupNo))
with(germ,table(Species))
with(germ,table(Source))
with(germ,table(Temp))
with(germ,table(WP))

ggplot(filter(germ, Species=="DISP"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)
ggplot(filter(germ, Species=="PHAU"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)
ggplot(filter(germ, Species=="DISP"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)
ggplot(filter(germ, Species=="SCAC"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)
ggplot(filter(germ, Species=="SCAM"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)
ggplot(filter(germ, Species=="JUBA"), aes(x=Temp, y=PropGerm)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter() +
  facet_wrap(~WP)

#### Subset out by species ####
htt.phau <- subset(germ, Species=="PHAU")
htt.phau$Species <- factor(htt.phau$Species)
htt.phau$Source <- factor(htt.phau$Source)
htt.phau$CupNo <- factor(htt.phau$CupNo)
htt.phau <- htt.phau[!is.na(htt.phau$CumPropGerm),]
str(htt.phau)

htt.disp <- subset(germ, Species=="DISP")
htt.disp$Species <- factor(htt.disp$Species)
htt.disp$Source <- factor(htt.disp$Source)
htt.disp$CupNo <- factor(htt.disp$CupNo)
htt.disp <- htt.disp[!is.na(htt.disp$CumPropGerm),]
str(htt.disp)

htt.boma <- subset(germ, Species=="BOMA")
htt.boma$Species <- factor(htt.boma$Species)
htt.boma$Source <- factor(htt.boma$Source)
htt.boma$CupNo <- factor(htt.boma$CupNo)
htt.boma <- htt.boma[!is.na(htt.boma$CumPropGerm),]
str(htt.boma)

htt.scac <- subset(germ, Species=="SCAC")
htt.scac$Species <- factor(htt.scac$Species)
htt.scac$Source <- factor(htt.scac$Source)
htt.scac$CupNo <- factor(htt.scac$CupNo)
htt.scac <- htt.scac[!is.na(htt.scac$CumPropGerm),]
str(htt.scac)

htt.scam <- subset(germ, Species=="SCAM")
htt.scam$Species <- factor(htt.scam$Species)
htt.scam$Source <- factor(htt.scam$Source)
htt.scam$CupNo <- factor(htt.scam$CupNo)
htt.scam <- htt.scam[!is.na(htt.scam$CumPropGerm),]
str(htt.scam)

htt.juba <- subset(germ, Species=="JUBA")
htt.juba$Species <- factor(htt.juba$Species)
htt.juba$Source <- factor(htt.juba$Source)
htt.juba$CupNo <- factor(htt.juba$CupNo)
htt.juba <- htt.juba[!is.na(htt.juba$CumPropGerm),]
str(htt.juba)

htt.elpa <- subset(germ, Species=="ELPA")
htt.elpa$Species <- factor(htt.elpa$Species)
htt.elpa$Source <- factor(htt.elpa$Source)
htt.elpa$CupNo <- factor(htt.elpa$CupNo)
htt.elpa <- htt.elpa[!is.na(htt.elpa$CumPropGerm),]
str(htt.elpa)
#####

#### THERMAL TIME ANALYSIS ####
## Thermal time model; logistic function
# G(t) = Gmax / (1+exp(b(log(t) - log(t50))))
# G(t) = cumulative germination over time (t)
# Gmax = maximum germination as t approaches infinity
# b = slope around the inflection point
# t50 = time at which germination is half the Gmax

## Extract the sub-optimal germination temperatures
## Set sub-optimal at 28/10 (temp 28)
phau.sub <- subset(htt.phau, Temp==28)
phau.sub$Chamber <- factor(phau.sub$Chamber)
phau.sub$CupNo <- factor(phau.sub$CupNo)
str(phau.sub)

# Visualize raw data (not averaged over replicate cups)
ggplot(phau.sub,aes(x= BeginTime, y=CumPropGerm*100, color=WP)) +
  geom_point() +
  facet_wrap(~Source)

# Average over water potential (focus only on thermal time)
phau.sub.nowp <- aggregate(CumPropGerm ~ BeginTime + Temp + Source + Chamber + CupNo, mean, data=phau.sub)

## Fit nonlinear regression model
# LL.3 = 3 parameter log-logistic function where lower limit is equal to 0
# see getMeanFunctions() for possible drm distributions
reg_sub_model_phau <- drc::drm(CumPropGerm ~ BeginTime, curveid=Source, data=phau.sub.nowp, fct=drc::LL.3())
summary(reg_sub_model_phau) # b = slope around inflection point, d = Gmax, and e = t50

# Visualize model fit
fits <- fitted(reg_sub_model_phau)
phau.sub.nowp$fits <- fits

ggplot(phau.sub.nowp) +
  geom_point(aes(x= BeginTime, y=CumPropGerm*100, color=Source)) +
  geom_line(aes(x= BeginTime, y=fits*100, color=Source))
#####


#### CALCULATE GR, WB_50, SIG_WB, AND TB FROM DATA ####
## Subset to only one source& less than optimal temp
ogba.sub <- subset(htt.phau, Temp<36 & Source=="OGBA")
ogba.sub$Chamber <- factor(ogba.sub$Chamber)
ogba.sub$CupNo <- factor(ogba.sub$CupNo)
str(ogba.sub)

# Average over water potential (focus only on thermal time)
ogba.sub.nowp <- aggregate(CumPropGerm ~ BeginTime + Temp + Source + Chamber + CupNo, mean, data=ogba.sub)

# set temp as factor
ogba.sub.nowp$Temp <- as.factor(ogba.sub.nowp$Temp)

## Run nonlinear regression model for germination based on temp
reg_sub_model_ogba <- drc::drm(CumPropGerm ~ BeginTime, curveid=Temp, data=ogba.sub.nowp, fct=drc::LL.3(), na.action=na.omit)
summary(reg_sub_model_ogba)

# Visualize model fit
fits <- fitted(reg_sub_model_ogba)
ogba.sub.nowp$fits <- fits
ggplot(ogba.sub.nowp) +
  geom_point(aes(x= BeginTime, y=CumPropGerm*100, color=Temp)) +
  geom_line(aes(x= BeginTime, y=fits*100, color=Temp))

## Calculate time to germination using effective dose function
drc::ED(reg_sub_model_ogba, 0.5, type = "absolute")
## TODO : Decide what to do when model never predicts 50% of seeds germinating
## Possible alternative: when 50% of total seeds that are going to germinate do?
# to do this use type="relative"
Tg_ogba <- drc::ED(reg_sub_model_ogba, seq(0.05,0.9,0.05), type="absolute")
plot(seq(0.05,0.9,0.05), Tg_ogba[1:18,1], type='l', xlab="Percentile", ylab="Time to germ", main="Temp = 32C")
plot(seq(0.05,0.9,0.05), Tg_ogba[19:36,1], type='l', xlab="Percentile", ylab="Time to germ", main="Temp = 28C")

## Calculate germination rate GR based on temperatures
GR_ogba <- 1 / Tg_ogba
head(GR_ogba)

## Explore GR-Temp relationships
# Add parameter names as fields
GR_ogba <- as.data.frame(GR_ogba)
rn <- stringr::str_split_fixed(rownames(GR_ogba), ":", 3) # split rownames
GR_ogba$Temp <- as.numeric(rn[,2]) # add temps from rownames
GR_ogba$g <- as.factor(rn[,3]) # add time from rownames
GR_ogba <- na.omit(GR_ogba) # remove NAs
timesteps <- length(unique(GR_ogba$g))

# Plot germination rate by temperature
ggplot(GR_ogba, aes(x=Temp, y=Estimate, color=g)) +
  geom_point()+
  geom_line()+
  labs(x= "Temperature (C)" , 
       y="Germination rate (1/Tg)", color="Percentile" )+
  xlim(c(10,35))

## Estimate the base temperature Tb
GR_model_ogba <- nls(Estimate~b[g]*(Temp-Tb), start=list(b=rep(0.10,timesteps), Tb=0), data=GR_ogba)
sum <- summary(GR_model_ogba)
tb_coeff <- sum$coefficients["Tb",1]
b_coeff <- mean(sum$coefficients[c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15","b16"),1], na.rm=T)

# visualize model fit
GR_fits <- fitted(GR_model_ogba)
GR_ogba$fits <- GR_fits
ggplot(GR_ogba)+
  geom_point(aes(x=Temp, y=Estimate, color=g))+
  geom_line(aes(x=Temp, y=fits, color=g))+
  labs(x= "Temperature (C)" , y="Germination rate (1/tg)", color="Percentile" ) 

# test base temp with predictions from model
X <- expand.grid(Temp = seq(tb_coeff, 36, 2), g = unique(GR_ogba$g))
X$fits <- predict(GR_model_ogba, newdata=X)
ggplot(GR_ogba, aes(x=Temp, y=Estimate, color=g))+
  geom_point()+
  geom_line(data=X,aes(x= Temp, y=fits, color=g))+
  labs(x="Temperature (C)", y="Germination rate (1/tg)", color="Percentile" ) 

# test other model option...
GR_model2 <- nls(Estimate~b*(Temp-Tb[g]), start=list(Tb=rep(tb_coeff,timesteps),b=b_coeff), data=GR_ogba)
sum2 <- summary(GR_model2)
tb_index <- paste0("Tb", as.character(timesteps/2))
b_coeff <- sum2$coefficients["b",1]
tb_coeff_2 <- sum2$coefficients[tb_index,1]
tb_coeffs_2 <- sum2$coefficients[c("Tb1","Tb2","Tb3","Tb4","Tb5","Tb6","Tb7","Tb8","Tb9","Tb10","Tb11","Tb12","Tb13","Tb14","Tb15","Tb16"),1]

# plot fit of model 2
X2 <- expand.grid(Temp = seq(tb_coeff_2, 36, 1), g=unique(GR_ogba$g))
X2$fits <- predict(GR_model2, newdata=X2)
ggplot(GR_ogba, aes(x=Temp, y=Estimate, color=g))+
  geom_point()+
  geom_line(data=X2,aes(x= Temp, y=fits, color=g))+
  ylim(0,3)+
  labs(x= "Temperature (C)" , y="Germination rate (1/tg)", color="Percentile" )

# test model 3
## can't get this one running....
GR_model3 <- nls(Estimate~b[g]*(Temp-Tb[g]), start=list(Tb=rep(20,timesteps), b=rep(0.03,timesteps)), data=GR_ogba, trace=TRUE)
sum3 <- summary(GR_model3)
tb_coeff_2 <- sum3$coefficients[tb_index,1]

# plot fit of model 3
X3 <- expand.grid(Temp = seq(tb_index, 36, 0.1), g=unique(GR_ogba$g))
X3$fits <- predict(GR_model3, newdata=X3)
ggplot(GR_ogba, aes(x=Temp, y=Estimate, color=g)) +
  geom_point()+
  geom_line(data=X3, aes(x=Temp, y=fits, color=g)) +
  ylim(0,0.7)+
  labs(x="Temperature (C)" , y="Germination rate (1/tg)", color="Percentile" ) 

# compare fit of models
AIC(GR_model_ogba, GR_model2)
anova(GR_model_ogba, GR_model2)
# choose model 1 since lower AIC and lower residuals
# that means our base temp is:
sum$coefficients["Tb",]

## Relating b to g
# plot g vs. b
parms <- coefficients(GR_model_ogba)
b <- parms[1:timesteps]
knitr::kable(b, col.names="b")
gdd_g <- data.frame(b=b, g=unique(GR_ogba$g))
gdd_g$g <- as.numeric(as.character(gdd_g$g))
ggplot(gdd_g, aes(x=g, y=b)) +
  geom_point()+
  labs(x= "Germination fraction (g)" , y="Slope (b)") 

# fit models
# normal model
b_normal <- nls(b~qnorm(1-g,mu,sigma), start=list(mu=0.04,sigma=0.05), data=gdd_g)
summary(b_normal)
#logistic model
b_logistic <- nls(b~b50*(-g/(g-1))^(1/c), start=list(b50=0.04, c=-2), data=gdd_g)
summary(b_logistic)

# compare models
AIC(b_normal,b_logistic)

X_gdd <- expand.grid(g=seq(0.09,0.91,0.01))
X_gdd$fits <- predict(b_normal, newdata=X_gdd)

X_gdd2 <- expand.grid(g=seq(0.09,0.91,0.01))
X_gdd2$fits <- predict(b_logistic, newdata=X_gdd2)

ggplot() +
  geom_line(data=X_gdd, aes(x=g,y=fits), color="#D7690A",linetype="longdash" )+
  geom_line(data=X_gdd2, aes(x=g,y=fits), color="#046C9A")+
  geom_point(data=gdd_g, aes(x=g,y=b), color="#398882")+
  labs(x= "Germination percentile (g)" , y="Slope (b)") 
# logistic regression fits the best by far

## Calculate Wb_50 & sig_50 using hydrotime model
ogba.sub.28 <- subset(ogba.sub, Temp==28)
ogba28.sub.mean <- aggregate(CumPropGerm ~ BeginTime + Temp + WP + Source + Chamber + CupNo, mean, data=ogba.sub.28)

# model GR
ht_reg <- drc::drm(CumPropGerm~BeginTime, as.factor(WP), fct=drc::LL.3(), data=ogba28.sub.mean)
summary(ht_reg)

# plot fit
ogba28.sub.mean$reg_fit <- fitted(ht_reg)
ggplot(ogba28.sub.mean)+
  geom_point(aes(x=BeginTime, y=CumPropGerm*100, color=as.factor(WP)))+
  geom_line(aes(x=BeginTime, y=reg_fit*100, color=as.factor(WP)))+
  labs(x= "Time (day)" , y="Cumulative germination (%)", color="Water potential (MPa)")

# estimate GR values by water potential
g <- seq(0.05, 0.9, 0.05)
gr_ht <- 1/drc::ED(ht_reg, g, type="absolute", display=FALSE)

gr_ht <- as.data.frame(gr_ht)
rn_ht <- stringr::str_split_fixed(rownames(gr_ht), ":", 3) 
gr_ht$water <- as.numeric(rn_ht[,2])
gr_ht$g <- as.factor(rn_ht[,3])
gr_ht <- na.omit(gr_ht)

# plot GR vs. water potential
ggplot(gr_ht)+
  geom_line(aes(x=water, y=Estimate, color=as.factor(g)))+
  geom_point(aes(x=water, y=Estimate, color=as.factor(g)))+
  labs(x= "Water potential (MPa)" , y="Germination rate (1/time)", color="Germination percentile (g)")

# fit hydrotime model
ht_ln <- length(unique(gr_ht$g))
ht_gr_mod <- nls(Estimate~b*(water-wb[g]), start=list(b=0.1, wb=rep(-0.5,ht_ln)), data=gr_ht)
sum_ht <- summary(ht_gr_mod)
ht_coeffs1 <- sum_ht$coefficients[c("wb1","wb2","wb3","wb4","wb5","wb6","wb7","wb8","wb9","wb10","wb11","wb12","wb13","wb14"),1]

# visually check fit
X_wt <- expand.grid(water=seq(-2.3,0,0.01), g = unique(gr_ht$g))
X_wt$fits <- predict(ht_gr_mod, newdata=X_wt)
ggplot(gr_ht) +
  geom_point(aes(x= water, y=Estimate, color=as.factor(g)))+
  geom_line(data=X_wt, aes(x= water, y=fits, color=as.factor(g)))+
  ylim(0, 2.5)+
  labs(x="Water potential (MPa)" , y="Germination rate (1/time)", color="Germination percentile (g)")

# check other options...
ht_gr_mod2 <- nls(Estimate~b[g]*(water-wb), start=list(wb=-2.5,b=rep(0.18,ht_ln)), data=gr_ht, trace=TRUE)
sum2 <- summary(ht_gr_mod2)
tb_index <- paste0("Tb", as.character(timesteps/2))
#tb_coeffs_2 <- sum2$coefficients[c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14"),1]
#b_coeff <- sum2$coefficients["b",1]

# check fit
X_wt <- expand.grid(water=seq(-2.3,0,0.01), g = unique(gr_ht$g))
X_wt$fits <- predict(ht_gr_mod2, newdata=X_wt)
ggplot(gr_ht)+
  geom_point(aes(x= water, y=Estimate, color=as.factor(g)))+
  geom_line(data=X_wt, aes(x= water, y=fits, color=as.factor(g)))+
  ylim(0, 2.5)+
  xlim(-1.5,0)+
  labs(x="Water potential (MPa)" , y="Germination rate (1/time)", color="Germination percentile (g)")

# best model is obviously model 2- based on this, base WP is -1.5

# HT model
# testing starting values
test <- data.frame(theta=2, wb_50=-0.2, sig_wb=0.8)
germ_test_wp1 <- pnorm(0-(test$theta/ogba28.sub.mean$BeginTime), test$wb_50, test$sig_wb)
germ_test_wp2 <- pnorm(-0.6-(test$theta/ogba28.sub.mean$BeginTime), test$wb_50, test$sig_wb)
ggplot()+
  geom_line(aes(ogba28.sub.mean$BeginTime, germ_test_wp1), col="red")+
  geom_point(data=ogba28.sub.mean[which(ogba28.sub.mean$WP==0),], aes(BeginTime, CumPropGerm))

# run HT model
## lower and upper set parameter boundaries- only works with algorithm "port" 
ht_cumGerm <- nls(CumPropGerm ~ pnorm(WP-(theta/BeginTime), wb_50, sig_wb), start = list(theta=2, wb_50=-0.2, sig_wb=0.8), data=ogba28.sub.mean, trace=T)#, algorithm="port", lower=c(1,-0.6,0.01), upper=c(100,1,2))
summary(ht_cumGerm)
phal16_mean$ht_fits<-fitted(ht_cumGerm)

# plot HT model
g1<- ggplot(phal16_mean)+
  geom_point(aes(x=timeBef, y=propCum*100, color=as.factor(water)))+
  geom_line(aes(x=timeBef, y=ht_fits*100, color=as.factor(water)))+
  labs(x= "Time (h)" , y="Cumulative germination (%)", color="Water potential (MPa)")

# check distribution of wb(g)
water <- seq(-4,2,0.05)
fit_gmax<-dnorm(water,coef(ht_cumGerm)["wb_50"],coef(ht_cumGerm)["sig_wb"])
ggplot()+
    geom_line(aes(x=water, y=fit_gmax), color="#00998a")+
    geom_area(aes(x=water[water>=0], y=fit_gmax[water>=0]), fill="#00998a", alpha =0.5)+
  labs(x= "Water potential (MPa)" , y="Probability density")

# calculate proportion of dormant seeds
Q <- 1-pnorm(0,coef(ht_cumGerm)["wb_50"],coef(ht_cumGerm)["sig_wb"])
Q


#### Hydrothermal Time ####
## Create HTT function for the suboptimal temperature range
sub_HTT_fun <- function(WP, Temp, tg, theta_ht, Tb, wb_50, sigma_wb){
  wbg <- WP-(theta_ht/((Temp-Tb)*tg))
  CumGerm <- pnorm(wbg, wb_50, sigma_wb)
  return(CumGerm)
}#function

## PHAU analysis
# Fit HTT model 
sub_HTT<- nls(CumPropGerm ~ sub_HTT_fun(WP, Temp, BeginTime, theta_ht, Tb, wb_50, sigma_wb), 
              start = c(theta_ht = 50, Tb = 0, wb_50 = 0.3, sigma_wb = 0.15),
              data=phau.sub, trace=T) 
summary(sub_HTT)

sub_HTT2<- nls(CumPropGerm ~ sub_HTT_fun(WP, Temp, BeginTime, theta_ht, Tb, wb_50, sigma_wb), 
              start = c(theta_ht = 50, Tb = 0, wb_50 = 0.3, sigma_wb = 0.15),
              data=htt.phau, trace=T) 
summary(sub_HTT2)

summary(phau.sub)
summary(htt.phau)
str(htt.phau)
dim(phau.sub)
str(phalaris)
str(htt.phau)

## HTT model for sub & supra-optimal temperature ranges
HTT_AB_fun <- function(WP, Temp, tg, theta_ht, Tb, To, wb_o, k, sigma_wb) {
  wb_50 <- ifelse(Temp < To, wb_o, wb_o+k*(Temp-To))
  wb_g<- WP-(theta_ht/((Temp-Tb)*tg))
  cumgerm <- pnorm(wb_g, wb_50,sigma_wb)
  return(cumgerm)
}

# Fit model to mean germination data
phau.mean <- aggregate(CumPropGerm ~ Source + Temp + WP + BeginTime, mean, data=htt.phau)

HTT_phau.all <- nls(CumPropGerm ~ HTT_AB_fun(WP, Temp, BeginTime, theta_ht, Tb, To, wb_o, k, sigma_wb), 
                    start=c(theta_ht = 200, Tb=20, To=32, wb_o=-1, k=0.01, sigma_wb=0.4), 
                    data=htt.phau, trace=T) 
summary(HTT_AB)

summary(phau.mean)
