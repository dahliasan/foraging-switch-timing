### PLOT GLS TEMPERATURE AND WET/DRY DATA ###
rm(list=ls()) 

library(ggplot2)
library(reshape)
library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)

setwd("~/Dropbox/phd/lnfs/lnfs processed data/gls deg")

# create get deg file names
degf <- dir(pattern = ".deg")
locf <- dir(pattern = "estelle")

# read all deg file and merge to one dataframe
# t = internal temperature, wets = wet count (0-480), cond = conductivity (0-127)
df <- data.frame()
gl <- list()

for (i in seq_along(degf)){
  print (i)
  d <- read.table(degf[i], sep = "", skip = 20)
  id <- as.factor(str_split(degf[i], "_")[[1]][1])
  colnames(d) <- c("date", "time", "tmin", "tmax", "tmean", "wetcount", "cond")
  d$id <- id
  
  # combine date and time column
  gmt <-  as.POSIXct(strptime((paste(d$date, d$time)), format = "%d/%m/%Y %H:%M:%S"), tz = "GMT")
  d$gmt <- gmt
  d <- d[,-c(1:2)]
  
  ## get trip durations 
  d2 <- tbl_df(d)
  d2$wetok <- 0
  d2$wetok[d2$wetcount != 0] <- 1
  d2 <- d2 %>% mutate(wetokdiff = abs(wetok - lag(wetok)))
  d2$wetokdiff[1] <- 1
  d2 <- d2 %>% mutate(tripnum = cumsum(wetokdiff))
  
  ## summarise into per foraging trip
  tr <- d2 %>% group_by(tripnum) %>% summarise(id = first(id), start = first(gmt), end = last(gmt),
                                         dur = difftime(last(gmt), first(gmt), units = "days"),
                                         sea = first(wetok),
                                         meantmean = mean(tmean),
                                         meantmin = mean(tmin),
                                         meantmax = mean(tmax))
  
  ## subset out at sea foraging trips greater than 0.5 days
  foragetr <- tr %>% filter(sea == 1 & dur > 0.5)
  foragetr_gather <- foragetr %>% gather(key = variable, value = value, dur, meantmean, meantmax, meantmin)
  
  ## read in estelle csv (from processing light data using SGAT and BASTAG) lat lons 
  loc <- tbl_df(read.csv(locf[i]))
  loc$time <- as.POSIXct(strptime(loc$time, format = "%Y-%m-%d %H:%M:%S"), tz = "GMT")
  loc$id <- id
  colnames(loc)[1] <- "start"
  
  ## read IMOS ASL CTD data
  imos <- read.csv("IMOS_CTD_ASL_FEB2016.csv")
  imos$timestamp <- as.POSIXct(strptime(imos$timestamp, format = "%Y-%m-%dT%H:%M:%S"), tz = "GMT")
  imos <- tbl_df(imos)
  imos_temp <- imos %>% group_by(timestamp) %>% summarise(start = first(timestamp),
                                                          tmax = max(temp_vals),
                                                          tmin = min(temp_vals),
                                                          tmean = mean(temp_vals))
  imos_temp_gather <- imos_temp %>% gather(key = variable, value = value, tmean, tmax, tmin)
  
  
  ## make plots!
  xmin = min(c(foragetr$start, loc$start))
  xmax = max(c(foragetr$start, loc$start))
  
  g1 <- ggplot(data = NULL, aes(x = start, y = value, color = variable)) + 
    geom_line(data = foragetr_gather) +
    geom_point(data = foragetr_gather, size = 0.15) + 
    geom_point(data = imos_temp_gather, size = 0.1) + 
    facet_grid(id~ .) + 
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b", limits = c(xmin, xmax)) + 
    labs(x = "date", y = " value") + 
    theme(legend.position = "top") +
    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold")) + 
    theme(legend.text = element_text(size = 8))
  
  g2 <- ggplot(data = loc, aes(x = start, y = Lat)) + 
    geom_line() +
    facet_grid(id~ .) + 
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b", limits = c(xmin, xmax)) + 
    labs(x = "date", y = " lat") + 
    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"))
  
  legend <- get_legend(g1)
  
  g1 <- g1 + theme(legend.position = "none")
  
  g1 <- ggplot_gtable(ggplot_build(g1))
  g2 <- ggplot_gtable(ggplot_build(g2))
  
  maxWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3])
  
  g1$widths[2:3] <- maxWidth
  g2$widths[2:3] <- maxWidth
  
  # grid.arrange(g1, g2, heights = c(3, 2))
  
  ## merge into main dataframes
  df <- rbind(df, d2)
  gl <- c(gl, list(g1,g2))
  
  
}

## plot all g plots
grid.arrange(legend, gl[[1]], gl[[2]], gl[[3]], gl[[4]], gl[[5]],
             gl[[6]], gl[[7]], gl[[8]], gl[[9]], gl[[10]], ncol = 2, nrow = 7,
             layout_matrix = rbind(c(1,1), c(2,8), c(3,9), c(4,10), c(5,11), c(6,NA), c(7, NA)),
             heights = c(0.1, rep(c(0.3,0.2), 3)), widths = c(1,1))

## summarise into per foraging trip
tr <- df %>% group_by(id, tripnum) %>% summarise(seal = first(id), start = first(gmt), end = last(gmt),
                                             dur = difftime(last(gmt), first(gmt), units = "days"),
                                             sea = first(wetok),
                                             meantmean = mean(tmean),
                                             meantmin = mean(tmin),
                                             meantmax = mean(tmax))

foragetr <- tr %>% filter(sea == 1 & dur > 0.5)
foragetr_gather <- foragetr %>% gather(key = variable, value = value, dur, meantmean, meantmin, meantmax)

                                                        
## read IMOS ASL CTD data
imos <- read.csv("IMOS_CTD_ASL_FEB2016.csv")
imos$timestamp <- as.POSIXct(strptime(imos$timestamp, format = "%Y-%m-%dT%H:%M:%S"), tz = "GMT")
imos <- tbl_df(imos)
imos_temp <- imos %>% group_by(timestamp) %>% summarise(start = first(timestamp),
                                                                  tmax = max(temp_vals),
                                                                  tmin = min(temp_vals),
                                                                  tmean = mean(temp_vals))
imos_temp_gather <- imos_temp %>% gather(key = variable, value = value, tmean, tmax, tmin)

g <- ggplot(data = NULL, aes(x = start, y = value, color = variable)) + 
  geom_line(data = foragetr_gather) +
  geom_point(data = imos_temp_gather, size = 0.1) + 
  facet_grid(id~ .) + 
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") + 
  ggtitle("LNFS foraging trip duration and sea temperature") + 
  xlab("Date")


# write.csv(foragetr, file = "lnfs_gls_foraging_trips_2016.csv")
