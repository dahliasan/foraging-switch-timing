## Create dataframe of foraging trips classifying each as either shelf or
## oceanic foraging and associated environmental predictors to put into model

rm(list = ls())
library(tidyverse)
library(RcppRoll)

## Import foraging trips
# combine into single dataframe
ff <- list.files('./foraging trip analysis', pattern = 'foragingtrips', full.names = T)
ft <- data.frame()
for(i in seq_along(ff)){
  f <- read.csv(ff[i])
  f$StartDate <- as.POSIXct(f$StartDate, 'GMT')
  f$LastDate <- as.POSIXct(f$LastDate, 'GMT')
  f$id <- strsplit(str_remove(ff[i], './foraging trip analysis/'), split = '_')[[1]][1]
  ft <- rbind(ft, f)
}
ft <- tbl_df(ft)
# ft <- ft %>% select(X, StartDate, LastDate, Dur, colonytime, colony, id)


# Classify trips into shelf, shelf break or oceanic types ----------------------------------------------

## Import foraging tracks
load("glsTracksAll.RData")
tracks <- tbl_df(tracks)
colnames(tracks)[1] <- 'Date'

## For each foraging trip get the max south Latitude and corresponding Lon 
x1 <- split(ft, ft$id)[[1]]
y1 <- split(tracks, tracks$id)[[1]]

# iterate over each seal
tmp <- Map(function(trac, trip){
  # interate over each trip
  Map(function(t1, t2){ 
    a <- min(trac$Lat[trac$Date >= t1 & trac$Date <= t2])
    b <- trac$Lon[trac$Lat == a]
    # d <- trac$Date[trac$Lat == a]
    return(data.frame(a,b))
    # return(a)
    }, as.list(trip$StartDate), as.list(trip$LastDate))
}, split(tracks, tracks$id), split(ft, ft$id))
tmp <- lapply(tmp, function(x) do.call(rbind,x))
tmp <- do.call(rbind,tmp)

ft$minLat <- tmp$a
ft$Lon_minLat <- tmp$b

## Get shelf break data
library(marmap)
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000 & z <= -1000)

# get bathy for max south location of each foraging trip
# Do this by finding the minimum euclidean distance
ft$z <- sapply(split(ft, 1:nrow(ft)), function(x){
  bathy_map$z[which.min(sqrt((x$minLat - bathy_map$y)^2) + sqrt((x$Lon_minLat - bathy_map$x)^2))] })

## classify foraging locations as shelf, shelfbreak, oceanic
ft$floc <- 'oceanic'
ft$floc[ft$z <= -1000 & ft$z >= -2000] <- 'shelf break'
ft$floc[ft$z > -1000] <- 'shelf'

# check if you did it right - yeah i did!
ggplot() +
  geom_path(data = b2000, aes(x,y)) +
  geom_point(data = ft, aes(Lon_minLat, minLat, color = floc), size = 0.5)


# Extract shelf enviro parameters -----------------------------------------
load('./enviro data/bonneyupwelling_ssta_2016-2017.Rdata')
load('./enviro data/shelf(iso2000)_ssta_2016-2017.Rdata')
load('./enviro data/SSHA_shelf(isobath2000)_2016-17.Rdata')
load("./enviro data/bonneycoast_10m_upwelling_wind_index_1997-2017.RData")
she <- tbl_df(she)
ssha_shelf <- tbl_df(ssha_shelf)
w <- uw %>% rename(wind_index = upwell_wind)

## Bonney Upwelling SSTa area -----
# Calculate daily BU SSTa area (each cell is 0.25 degrees resolution)
bu <- bud %>% group_by(date) %>% summarise(SSTa_ncells = sum(!is.na(v1)))

## Shelf SSTA, SSHA -----
# Calculate daily mean, sd, max and min of shelf SST anomaly (SSTA) and shelf SSH anomaly (SSHA)
a <- she %>% filter(!is.na(v1)) %>% group_by(date) %>% summarise(SSTA_SD = sd(v1, na.rm = T), SSTA_mean = mean(v1, na.rm = T), SSTA_max = max(v1, na.rm = T), SSTA_min = min(v1, na.rm = T))
b <- ssha_shelf %>% filter(!is.na(v1)) %>% group_by(date) %>% summarise(SSHA_SD = sd(v1, na.rm = T), SSHA_mean = mean(v1, na.rm = T), SSHA_max = max(v1, na.rm = T), SSHA_min = min(v1, na.rm = T))
shelf <- left_join(bu, a, by = 'date')
shelf <- left_join(shelf, b, by = 'date')

# Calculating rolling 5-day means (past 5 days) 
shelf5d <- shelf %>% ungroup() %>% map_df(~ roll_mean(.x, n = 5, align = 'right', na.rm = T))
shelf5d$date <- shelf$date[-c(1:4)]

## Alongshore wind stress -----
# calc 5-day rolling mean and sd of windstress
w5d<- w %>% 
  ungroup() %>% 
  dplyr::select(date, wind_index) %>% 
  slice(-1:-4) %>% 
  mutate(wind_mean = roll_mean(w$wind_index, n = 5, align = 'right', na.rm = T), 
         wind_SD = roll_sd(w$wind_index, n= 5, align = 'right', na.rm = T),
         wind_min = roll_min(w$wind_index, n= 5, align = 'right', na.rm = T),
         wind_max = roll_max(w$wind_index, n= 5, align = 'right', na.rm = T))

shelf5d <- left_join(shelf5d, w5d, by = 'date')

## IMOS ASL CTD data on shelf: MLD, TEMP, SALINITY -----
imos <- tbl_df(read.csv("./enviro data/IMOS_CTD_ASL_FEB2016.csv")) %>% dplyr::select(timestamp, lon, lat, profile_id, pressure, temp_vals, sal_corrected_vals)
imosb <- tbl_df(read.csv("./enviro data/IMOS_CTD_ASL_MAR2017.csv", skip = 251) %>% dplyr::select(timestamp, lon, lat, profile_id,  pressure, temperature, salinity))
colnames(imos) <- colnames(imosb)
imos <- rbind(imos, imosb)
imos$timestamp <- as.POSIXct(strptime(imos$timestamp, format = "%Y-%m-%dT%H:%M:%S"), tz = "GMT")
anytime::anydate(sapply(split(imos$timestamp, lubridate::year(imos$timestamp)), range))

# write calculate MLD function
library(oce)
calculate.mld <- function(sigma, z, deltaSigma = 0.125) { 
  # check that the number of sigma and z values are equal 
  if (length(sigma) != length(z)) 
  { 
    stop('sigma and z vectors must have equal length') 
  } 
  # remove the na's 
  keepers <- !(is.na(z) | is.na(sigma)) 
  z <- z[keepers] 
  sigma <- sigma[keepers] 
  # return an NA and a warning if there aren't at least two 
  # numeric values 
  if (length(sigma) < 2) 
  { 
    warning('fewer than two valid sigma and depth-value 
            combinations entered, 
            NA returned') 
    NA 
  } else { 
    # I use negative depths to be consistent with the Scripps database 
    if (all(z >= 0)) 
    { 
      z <- z * -1 
      pos <- TRUE 
    } else { 
      pos <- FALSE 
    } 
    # be sure the data are sorted by depths 
    ord <- order(z, decreasing = TRUE) 
    z <- z[ord] 
    sigma <- sigma[ord] 
    #z <- sort(z, decreasing = TRUE) 
    #sigma <- sigma[order(z, decreasing = TRUE)] 
    # Manuscript uses a z of 10 m as the initial sigmaRef, but we will 
    # take the closest value in case the 10-m measurement is missing. 
    minDepth <- which(abs(z + 10) == min(abs(z + 10))) 
    minDepth <- ifelse(length(minDepth) > 1, minDepth[2], minDepth) 
    sigmaRef <- sigma[minDepth] 
    sigma <- sigma[minDepth:length(sigma)] 
    z <- z[minDepth:length(z)] 
    diffs <- abs(sigma - sigmaRef) 
    # if sigma never changes by at least deltaSigma, return the lowest depth 
    # Otherwise proceed
    if (max(diffs, na.rm = TRUE) >= deltaSigma) 
    { 
      # the uniform region, if present, occurs where the change between any
      # two points is <= deltaSigma * 1/10, and the change in sigma over the
      # profile has not yet exceeded deltaSigma
      uniformRegion <- (abs(sigma[2:length(sigma)] - 
                              sigma[1:(length(sigma) - 1)]) <= 
                          (deltaSigma / 10)) & 
        (diffs[2:length(diffs)] < deltaSigma) 
      if (any(uniformRegion)) 
      { 
        sigmaRefPos <- max(which(uniformRegion)) 
        # change sigmaRef to the base of the uniform region 
        sigmaRef <- sigma[sigmaRefPos] 
        # calculate change from the new reference 
        reachedDeltaSigma <- 
          abs(sigma[(sigmaRefPos + 1):length(sigma)] - 
                sigmaRef) >= 
          deltaSigma 
        # if any deeper measurements of sigma reach or exceed deltaSigma, 
        # linearly interpolate between the nearest points to find mixed-layer
        # depth
        if (any(reachedDeltaSigma)) 
        { 
          pair <- min(which(reachedDeltaSigma)) + sigmaRefPos - 1 
          linmod <- lm(z[c(pair, pair + 1)] ~ sigma[c(pair, 
                                                      pair + 1)]) 
          mld <- as.vector(linmod$coefficients[1] + 
                             linmod$coefficients[2] * 
                             (sigmaRef + deltaSigma)) 
        } else { 
          # otherwise, mixed-layer depth is the deepest point 
          mld <- min(z) 
        } 
      } else { 
        # if there is no uniform region, just linearly interpolate mld
        pair <- min(which(diffs >= deltaSigma)) - 1 
        linmod <- lm(z[c(pair, pair + 1)] ~ sigma[c(pair, pair + 1)]) 
        mld <- as.vector(linmod$coefficients[1] + 
                           linmod$coefficients[2] * 
                           (sigmaRef + deltaSigma)) 
      } 
    } else { 
      mld <- min(z) 
    } 
    if (pos) mld <- abs(mld) 
    mld 
  } 
} 

# calculate MLD 
# imos_1profile <- imos %>% filter(profile_id == 6154029)
tmp_sigma <- lapply(split(imos, imos$profile_id), function(x) return(with(x, swSigma0(salinity, temperature, pressure, lon, lat))))
tmp_mld <- map2(tmp_sigma, split(imos$pressure, imos$profile_id), calculate.mld)
tmp_mld <- data.frame(profile_id = as.integer(as.character(names(tmp_mld))), mld = unlist(tmp_mld))
imos <- left_join(imos, tmp_mld, by = 'profile_id')

## Summary statistics for MLD, TEMP, SALNITY per dive profile
# check which profile id only has 1 observation
unique(table(imos$profile_id))
r <- names(table(imos$profile_id))[ table(imos$profile_id) < 20]
# remove profile ids with only 1 observation
imos2 <- imos %>% filter(!profile_id %in% r) %>% mutate(date = as.Date(timestamp)) %>% 
  group_by(date) %>%
  summarise(t_SD = sd(temperature, na.rm = T), t_mean = mean(temperature, na.rm = T),
            t_max = max(temperature, na.rm = T), t_min = min(temperature, na.rm = T),
            sal_SD = sd(salinity, na.rm = T), sal_mean = mean(salinity, na.rm = T),
            p_SD = sd(pressure, na.rm = T), p_mean = mean(pressure, na.rm = T),
            mld_SD = sd(mld, na.rm = T), mld_mean = mean(mld, na.rm = T))

# save(imos2, file = 'ASLCTD_IMOS_2016-2017.RData')

# Get SAM data ------------------------------------------------------------
# download.file('ftp://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.aao.index.b790101.current.ascii',
#               destfile = 'sam.ascii')

sam <- readr::read_table('./enviro data/sam.ascii', col_names = c('year', 'month', 'day', 'sam'))
sam <- sam %>% mutate(date = as.Date(sprintf('%s-%s-%s', year, month, day)))
# save(sam, file = 'sam.Rdata')

# Add all shelf parameters to foraging trips ------------------------------
ft$date <- as.Date(ft$StartDate)

## interpolate ASL IMOS data to foraging trip dates
tmp <- map_df(imos2[,-1], function(y) approx(x = imos2$date, y = y, xout = ft$date)$y)
# ft <- left_join(ft, imos2, by = 'date')
ft <- bind_cols(ft, tmp)

## combine the rest of the shelf variables
ft <- left_join(ft, sam[,c('sam', 'date')], 'date' )
ft <- left_join(ft, shelf5d, 'date')

## Add season factor
getSeason <- function(DATES) {
  wi <- as.Date("2012-06-01", format = "%Y-%m-%d") # Winter 
  sp <- as.Date("2012-09-01",  format = "%Y-%m-%d") # Spring
  su <- as.Date("2012-12-01",  format = "%Y-%m-%d") # Summer
  au <- as.Date("2012-03-01",  format = "%Y-%m-%d") # Fall 
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= wi & d < sp, "Winter",
          ifelse (d >= sp & d < su, "Spring",
                  ifelse (d >= au & d < wi, "Autumn", "Summer")))
}

ft$season <- getSeason(ft$date)
ft$season <- factor(ft$season)

ft$floc <- factor(ft$floc)
save(ft, file=  './extract shelf enviro/lnfsForagingTrips_ShelfVars(iso2000).Rdata')

