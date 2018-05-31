## Exploration Data Analysis

rm(list = ls())
library(ggplot2)
library(tidyverse)
# library(lme4)
library(arm)
library(GGally)
library(ggpubr)
library(nlme)

# Load data ---------------------------------------------------------------
list.files(pattern = 'Rdata')
load("lnfsForagingTrips_ShelfVars.Rdata")
ft$year <- lubridate::year(ft$date)
ft$floc <- factor(ft$floc, levels = c('shelf', 'shelf break', 'oceanic'))
ft$season <- factor(ft$season, levels = c('Summer', 'Autumn', 'Winter', 'Spring'))

# remove non-lactating seals
ft2 <- ft %>% filter(id != '307', id != '317', id != '324',id!= '340') %>% 
  dplyr::select(-StartDate, -LastDate, -Dur, -colonytime, -colony, -z, -grp, -grp2, -pmean, -psd)
ft2$id <- factor(ft2$id)
ft2$year <- factor(ft2$year)
ft2 <- ft2 %>% mutate(mSSTa = mSSTa * -1, maxSSTa = maxSSTa * -1, minSSTa = minSSTa * -1)


# Check collinearity ------------------------------------------------------
# results show that ASSTa (count data) should be squareroot transformed,
# proportional data xSSTa1 should be arcsine transform, many of the other
# variables are bimodal
ggpairs(ft2[,-c(1:4)]) + 
  theme( text = element_text(size = 8))

## try conditional plot on bimodal variables e.g. all ASL CTD data, xSSTa data
# conditioning on year reduces bimodality so should add year and interactions into model later
ft2 %>% select(floc, tsd, tmean, salsd, salmean, mldsd, mldmean) %>% 
  ggpairs(aes(colour = year, alpha = 0.4)) + 
  theme( text = element_text(size = 8))

# ft2 %>% select(-c(1:3), -ASSTa, -sdSSTa1, mSSTa1, minSSTa1, maxSSTa1) %>% 
#   ggpairs(aes(colour = season, alpha = 0.4), size = 0.5) + 
#   theme( text = element_text(size = 8))
# 
# ft2 %>% select(-c(1:3), -ASSTa, -sdSSTa1, mSSTa1, minSSTa1, maxSSTa1) %>% 
#   ggpairs(aes(colour = floc, alpha = 0.4), size = 0.5) + 
#   theme( text = element_text(size = 8))

## Check for possibility of binning mean, max, min SSTa
# ft2 %>% select(sdSSTa, mSSTa, maxSSTa, minSSTa) %>% gather() %>% filter(value >0) %>% 
#   ggplot(aes(value)) + geom_histogram(bins = 6) + facet_wrap(~key)
# 
# ft2$maxSSTa1 <- cut(ft2$maxSSTa, c(-1, 0.1, 0.5, 1, 2), labels = c('none', 'small', 'med', 'high'))
# ft2$minSSTa1 <- cut(ft2$minSSTa, c(-1, 0.1, 0.5, 1, 2), labels = c('none', 'small', 'med', 'high'))
# ft2$mSSTa1 <- cut(ft2$mSSTa, c(-1, 0.1, 0.5, 1, 2), labels = c('none', 'small', 'med', 'high'))
# ft2$sdSSTa1 <- cut(ft2$sdSSTa, c(-1, 0.01, 0.05, 1), labels = c('none', 'small', 'med'))
#                                                                                 

## Do transformations
# arcsine for proportional data
# asinTransform <- function(p) { asin(sqrt(p)) }
# a <- ft2 %>% select(mSSTa1, sdSSTa1, maxSSTa1, minSSTa1) %>% map_df(~asinTransform(.x))

# create transformed dataset
# sqrt for count data, log for right skewed 
# ft3 <- ft2 %>% mutate(ASSTa = sqrt(ASSTa)) %>% select(-sdSSTa, -mSSTa, -maxSSTa, -minSSTa)
# logt <- ft2 %>% select(sdSSTa, mSSTa, maxSSTa, minSSTa) %>% map_df(~log(.x + 0.001))
# ft3 <- bind_cols(ft3, logt)
# ft3 <- ft3[,colnames(ft2)]
ft3 <- ft2 %>% mutate(ASSTa = sqrt(ASSTa)) %>% dplyr::select(-minLat, -Lon_minLat)
# ggpairs(ft3[,-c(1:3)]) + 
#   theme( text = element_text(size = 8))

# no strong correlation = good 
ft2[,c('tmean', 'salmean', 'mldmean')] %>% ggpairs() + 
  theme( text = element_text(size = 8))

# Check for outliers ------------------------------------------------------
ft4 <- ft3 %>% dplyr::select(-id, -floc, -season, -year, -X, -date) %>% gather(key = 'x', value = 'value')
ft4$x <- factor(ft4$x)

# x vars with outliers = ASSTa, mldsd, maxSSTa, salmean, sam, sdSSTa1, tmean
ft4 %>% 
ggplot(aes(x = x, y = value)) + 
  geom_boxplot()

ft4_out <- ft3 %>% select(mldsd, maxSSTa, salmean, sam, tmean) %>% 
  gather(key = 'x', value = 'value') %>% 
  mutate(x = factor(x))

# The x var with the most potential outliers is mldsd
ft4_out %>% ggplot(aes(x = x, y = value)) + 
  geom_boxplot()

# squareroot and log transform seems to work well on mldsd to get rid of extreme values
ft4_out %>% mutate(value = log(value)) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_boxplot()

# check density plot for sqrt and log transform on mldsd - LOG is better for normality as well.
ft3 %>% select(mldsd,year) %>% mutate(mldsd = log(mldsd)) %>% 
  ggplot(aes(x = mldsd, fill = year )) + 
  geom_density(alpha = .3)

## Create new transformed dataset
# ft5 <- ft3 %>% mutate(mldsd = log(mldsd + 0.01))
# ggpairs(ft5[,-c(1:3)]) + 
  # theme( text = element_text(size = 8))
ft5 <- ft3

# model preprocessing -----------------------------------------------------
## Remove collinear x vars - vif > 3
source("HighstatLib.r")
# remove all non-continuous vars
# keep removing variables 1 by 1 with vif > 3 until no more vars have vif > 3
# ft6 <- ft5 %>% select(-season, -year, -id, -floc, -minLat, -Lon_minLat, -mSSTa1, -minSSTa1,
#                       -minSSTa, -sdSSTa, -tsd, -mSSTa)
# corvif(ft6)
# ft6 <- bind_cols(ft6, ft5[,c('season', 'year', 'id', 'floc')])

ft6 <- ft5 %>% dplyr::select(-season, -year, -id, -floc, -X, -date,
                             -minSSTa, -mSSTa, -tsd)
corvif(ft6)
ft6 <- bind_cols(ft6, ft5[,c('season', 'year', 'id', 'floc', 'X','date')])

## scale and centre, remove 0 correlation variables
ft7 <- ft6 %>% mutate(tmean = (tmean^2), salmean = (salmean^2), mldmean = log(mldmean))
# ft7 <- ft6


ft7_pp <- caret::preProcess(ft7, method = c('center', 'scale'))
ftpp <- predict(ft7_pp, ft7)

# Modeling 2 --------------------------------------------------------------
# library(qpcR)
## compare only shelf vs oceanic (ie remove shelfbreak)
ftpp <- ftpp %>% filter(floc != 'shelf break')
ftpp <- ftpp %>% mutate(floc = ifelse(floc == 'shelf', 1, 0))
ftpp$floc <- factor(ftpp$floc)
data <- ftpp[!is.na(ftpp$tmean),]
data$X <- as.integer(data$X)

mod1 <- glmer(floc ~ tmean + salsd + salmean + mldsd + mldmean + sam + ASSTa + sdSSTa + maxSSTa + year + X+ (1 | id), data = data, family = binomial)
update(mod1, corr=corAR1())
summary(mod1, corr = F)
drop1(mod1,test="Chisq")
mod2 <- update(mod1, .~.-mldsd)
drop1(mod2,test="Chisq")
mod3 <- update(mod1, .~.-mldsd -maxSSTa)
drop1(mod3,test="Chisq")
mod4 <- update(mod1, .~.-mldsd -maxSSTa -sam)
drop1(mod4,test="Chisq")
mod5 <- update(mod1, .~.-mldsd -maxSSTa -sam -year)
drop1(mod5,test="Chisq")
mod6 <- update(mod1, .~.-mldsd -maxSSTa -sam -year -X)
drop1(mod6,test="Chisq")
mod7 <- update(mod1, .~.-mldsd -maxSSTa -sam -year -X -salsd)
drop1(mod7,test="Chisq")
mod8 <- update(mod1, .~.-mldsd -maxSSTa -sam -year -X -salsd -sdSSTa)
drop1(mod8,test="Chisq")
mod9 <- update(mod1, .~.-mldsd -maxSSTa -sam -year -X -salsd -sdSSTa -ASSTa)
drop1(mod9,test="Chisq")
summary(mod9)
mod10 <- update(mod1, .~.-mldsd -maxSSTa -sam -year -X -salsd -sdSSTa + ASSTa*salmean)
drop1(mod10,test="Chisq")

mod11 <- glmer(floc ~ tmean + salmean*ASSTa + mldmean + (1 | id/X), data = data, family = binomial)
mod12 <- glmer(floc ~ tmean + salmean*ASSTa + mldmean + (1 | id) + (1|X), data = data, family = binomial)
mod13 <- glm(floc ~ tmean + salmean*ASSTa + mldmean, data = data, family = binomial)



# Diagnostics of best model -----------------------------------------------
bbmle::AICctab(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)
# bbmle::AICtab(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)
best <- mod10
plot(best)
p1 <- plot(best,id=0.05,idLabels=~.obs)
p2 <- plot(best,ylim=c(-1.5,1),type=c("p","smooth"))
gridExtra::grid.arrange(p1,p2,nrow=1)
plot(best,id~resid(.))
# i guess here we're looking at residuals within the confidence intervals is
# good
binnedplot(fitted(best),resid(best, type = 'response')) 
# check overdispersion:
# theoretically residual deviance = residual df. Overdispersion may be an
# indicator that important predictor variables have not been considered, and
# variability in these factors may result in our data being biased
source('overdispersion.R')
overdisp_fun(best) # p > 0.05 therefore model not overdispersed
# fixed effects estimates
fixef(best)
# random effect intercept estimates
lattice::dotplot(ranef(best,condVar=TRUE)) # seems like not many unique ids.
## partial effects
library(remef)
# tmean
y_partial <- remef(best, fix = c("salmean", "mldmean"), ran = "all")
p3 <- ggplot(data = NULL, aes(x = data$tmean, y = y_partial)) + 
  geom_smooth() +
  geom_point()
# salmean
y_partial <- remef(best, fix = c("tmean", "mldmean"), ran = "all")
p4 <- ggplot(data = NULL, aes(x = data$salmean, y = y_partial)) + 
  geom_smooth() +
  geom_point()
# mldmean
y_partial <- remef(best, fix = c("tmean", "salmean"), ran = "all")
p5 <- ggplot(data = NULL, aes(x = data$mldmean, y = y_partial)) + 
  geom_smooth() +
  geom_point()
ggarrange(p3,p4,p5, nrow = 2, ncol=2)

qqnorm(resid(best))
qqline(resid(best))
summary(best)
# residuals to be approximately symmetrically distributed, ie if the min is -1
# but the max is 8 something is wrong. The residuals should also be centered
# around 0, ie a median more or less close to 0.

library(DHARMa)
set.seed(123)
simulationOutput <- simulateResiduals(fittedModel = best, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput) # you want p value > 0.05
# plotResiduals(data$tmean, simulationOutput$scaledResiduals) 
# test uniform residuals you want p value > 0.05, it also gives some indication of overdispersion (if p < 0.05)
testUniformity(simulationOutput = simulationOutput) 
# a more powerful test for overdispersion, you want p value > 0.05
testOverdispersion(simulationOutput = simulationOutput)
testTemporalAutocorrelation(simulationOutput = simulationOutput, time = data$date) # no temporal autocorrelation since p  > 0.05

# interpretation ----------------------------------------------------------

library(effects)
plot(allEffects(best))
with(ft2,table(floc, year, season))
## shelf (1) => higher tmean, lower mldmean
ft7 %>% ggplot(aes(x = floc, y = salmean)) +
  geom_boxplot(aes(fill = season)) + 
  facet_grid(~year)

data %>% ggplot(aes(x = season, y = salmean)) +
  geom_boxplot() +
  facet_wrap(floc~year)

data %>% ggplot(aes(x = season, y = salmean)) +
  geom_boxplot() +
  facet_wrap(floc~year)



# Modeling: consider shelf break as oceanic --------------------------
ftpp <- predict(ft7_pp, ft7)
ftpp <- ftpp %>% mutate(floc = ifelse(floc == 'shelf', 1, 0))
ftpp$floc <- factor(ftpp$floc)
data <- ftpp[!is.na(ftpp$tmean),]
data$X <- as.integer(data$X)

mod1 <- glmer(floc ~ tmean + salsd + salmean + mldsd + mldmean + sam + ASSTa + sdSSTa + maxSSTa + year + X+ (1 | id), data = data, family = binomial)
summary(mod1, corr = F)
drop1(mod1,test="Chisq")
mod2 <- update(mod1, .~.-mldsd)
drop1(mod2,test="Chisq")
mod3 <- update(mod1, .~.-mldsd -X)
drop1(mod3,test="Chisq")
mod4 <- update(mod1, .~.-mldsd -X -maxSSTa)
drop1(mod4,test="Chisq")
mod5 <- update(mod1, .~.-mldsd -X -maxSSTa -year)
drop1(mod5,test="Chisq")
mod6 <- update(mod1, .~.-mldsd -X -maxSSTa -year-sam)
drop1(mod6,test="Chisq")
summary(mod6)
mod7 <- update(mod1, .~.-mldsd -X -maxSSTa -year-sam + ASSTa*salmean)
drop1(mod7,test="Chisq")
mod8 <- update(mod1, .~.-mldsd -X -maxSSTa -year -sam + ASSTa*salmean -salsd)
drop1(mod8,test="Chisq")
mod9 <- update(mod1, .~.-mldsd -X -maxSSTa -year -sam + ASSTa*salmean -salsd -sdSSTa)
drop1(mod9,test="Chisq")

mod11 <- glmer(floc ~ tmean + salmean*ASSTa + mldmean + (1 | id/X), data = data, family = binomial)
mod12 <- glmer(floc ~ tmean + salmean*ASSTa + mldmean + (1 | id) + (1|X), data = data, family = binomial)
mod13 <- glm(floc ~ tmean + salmean*ASSTa + mldmean, data = data, family = binomial)
mod14 <- glm(floc ~ tmean + salmean*ASSTa + mldmean +id, data = data, family = binomial)

## Diagnostics
bbmle::AICctab(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod11,mod12,mod13,mod14)
best <- mod9
set.seed(123)
simulationOutput <- simulateResiduals(fittedModel = best, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput) # you want p value > 0.05
# plotResiduals(data$tmean, simulationOutput$scaledResiduals) 
# test uniform residuals you want p value > 0.05, it also gives some indication of overdispersion (if p < 0.05)
testUniformity(simulationOutput = simulationOutput) 
# a more powerful test for overdispersion, you want p value > 0.05
testOverdispersion(simulationOutput = simulationOutput)
testTemporalAutocorrelation(simulationOutput = simulationOutput, time = data$date) # no temporal autocorrelation since p  > 0.05
plot(allEffects(best))

# not sure ---------------------------------
# how to compare glm vs glmm AIC? to test random effect. 


# because not many winter and spring trips (and winter 2017 doesn't have IMOS data)


## remove non-lactating seals
ftl <- ft %>% filter(id != '307', id != '317', id != '324',id!= '340')
ftl$id <- factor(ftl$id)
ftl$year <- factor(ftl$year)
ftl <- ftl %>% mutate(mSSTa = mSSTa * -1, maxSSTa = maxSSTa * -1, minSSTa = minSSTa * -1)

## make histo/density plot of duration times
# see if bimodal to determine threshold between short and long trips. oceanic
# may not always mean long trips because they may only go as far as past shelf
# break. and there's a degree of accuracy of gls tracks

ftl %>% ggplot(aes(x = Dur)) +
  geom_histogram(bins = 2) + 
  facet_wrap(~year)

ftl %>% ggplot(aes(x = Dur)) +
  geom_density() + 
  facet_wrap(~id)


# based individual trip duration sd, we could set a threshold to determine seals
# that have an alternate long/short trips vs seal that only had one type of trip
# - although I would choose it based on the smallest trip sd for 2016 seals - if
# that's the case only 3 seals do not show alternating foraging strategy
ftl %>% group_by(id) %>% summarise(tripsd = sd(Dur)) %>% 
  ggplot(aes(tripsd)) + geom_histogram(bins = 2)

# out of the seals that show foraging strategy, how do we determine or definte
# the switch in foraging strategy. e.g. 2016 seals pretty obvious when seals
# make their first significant long trip after lots of shorter ones. but how to
# define what's significantly long enough? maybe use the histogram again and
# find the break for 2 bins of each individual i.e. each individual has their
# own unique threshold for difference foraging trip duration

ftl %>% ggplot(aes(Dur)) + geom_histogram(bins = 2) + facet_wrap(~id)
# actually seems like the histogram with 2 bins on trip duration can also show
# which seals have 2 distinct foraging trip lengths

ftl %>% filter(id != '305', id != '315', id!= '353') %>% 
  ggplot(aes(x = date, y = Dur)) + 
  # geom_line() + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  geom_line(aes(color = id), alpha = .8) + 
  geom_point(aes(color = id), size = .3)


# find the biggest jump in trip duration 
ftl %>% filter(id != '305', id != '315', id!= '353') %>% 
  group_by(id) %>% 
  mutate(DurDiff = Dur - lag(Dur)) %>% 
  filter(DurDiff == max(DurDiff, na.rm = T)) %>% 
  ggplot(aes(x = date, y = DurDiff)) + 
  # geom_line() + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  geom_point(aes(color = id), alpha = .8)


# find transition date by shelf vs oceanic locations - seems like the best 
# method the transition in a switch of foraging strategy can be seen in the big 
# jumps, and generally from a shelf/shelf break to consecutive oceanic foraging 
# oceanic foraging locations near the shelf break are not considered the time of
# transition as geolocation locations have up to +-63km (~0.5degree) accuracy
# (based on validation with gps) so the seals could have been on the shelf at
# these locations too.
ftl %>% 
  # filter(id != '305', id != '315', id!= '353') %>% 
  group_by(id) %>% 
  ggplot(aes(x = as.Date(format(date, '2016-%m-%d')), y = minLat)) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  geom_line(aes(group = id, color = floc), alpha = .8) +
  geom_point( aes(color = floc), size = 0.8) + 
  facet_wrap(year~id) +
  geom_text(aes(label = 1:nrow(ftl)), size = 3) +
  labs(x = 'Date')

ftl %>% filter(id != '305', id != '315', id!= '353') %>% 
  group_by(id) %>% 
  mutate(minLatDiff = minLat - lag(minLat)) %>% 
  ggplot(aes(x = as.Date(format(date, '2016-%m-%d')), y = minLatDiff)) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_line(aes(group = id, color = floc), alpha = .8) +
  geom_point( aes(color = floc), size = .8) + 
  facet_wrap(year~id) +
  labs(x = 'Date')

fts <- ftl[c(13, 34, 56, 81, 101, 113, 121, 151, 159, 170, 176),]

## Find variability
sdates <- as.Date(format(fts$date, '2016-%m-%d'))
sd(sdates)
range(sdates)
fts %>% group_by(year) %>% summarise(sdatesd = sd(as.numeric(date)), mindate = min(date), maxdate = max(date))
  
