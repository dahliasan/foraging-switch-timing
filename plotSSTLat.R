## plot gls sst against ASL IMOS sst. with latitude and foraging trip duration
## from estelle tracks
rm(list = ls())
library(tidyverse)
library(lubridate)
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(RColorBrewer)

load('sstAll2017(outliersRemoved).Rdata')
load('degAll2016(outliersRemoved).Rdata')

# Import IMOS SST ---------------------------------------------------------
dir(pattern = 'IMOS')
#2016
imos1 <- tbl_df(read.csv("IMOS_CTD_ASL_FEB2016.csv"))
imos1$timestamp <- as.POSIXct(strptime(imos1$timestamp, format = "%Y-%m-%dT%H:%M:%S"), tz = "GMT")
imos1 <- imos1 %>% group_by(timestamp) %>% 
  summarise(tmax = max(temp_vals), tmin = min(temp_vals), tmean = mean(temp_vals))
colnames(imos1)[1] <- "Date"
#2017
imos2 <- tbl_df(read.csv("IMOS_CTD_ASL_MAR2017.csv", skip = 251))
imos2$timestamp <- as.POSIXct(strptime(imos2$timestamp, format = "%Y-%m-%dT%H:%M:%S"), tz = "GMT")
imos2 <- imos2 %>% group_by(timestamp) %>% 
  summarise(tmax = max(temperature), tmin = min(temperature), tmean = mean(temperature))
colnames(imos2)[1] <- "Date"

imos1$date <- as.POSIXct(strftime(imos1$Date, format = "2016-%m-%d %H:%M:%S", tz = "GMT"))
imos2$date <- as.POSIXct(strftime(imos2$Date, format = "2016-%m-%d %H:%M:%S", tz = "GMT"))
  

# Import Estelle tracks ---------------------------------------------------
load("glsTracksAll.RData")
tracks <- tbl_df(tracks)
colnames(tracks)[1] <- 'Date'

## Import foraging trips
ff <- list.files(pattern = 'foragingtrips')
# ft <- data.frame()
# for(i in seq_along(ff)){
#   f <- read.csv(ff[i])
#   f$StartDate <- as.POSIXct(f$StartDate, 'GMT')
#   f$LastDate <- as.POSIXct(f$LastDate, 'GMT')
#   f$id <- strsplit(ff[i], split = '_')[[1]][1]
#   ft <- rbind(ft, f)
# }
# ft <- tbl_df(ft)

ft <- map(ff, function(x){
  f <- read_csv(x)
  f$StartDate <- as.POSIXct(f$StartDate, 'GMT')
  f$LastDate <- as.POSIXct(f$LastDate, 'GMT')
  f$id <- strsplit(x, split = '_')[[1]][1]
  return(f)
}) %>% reduce(bind_rows)


# At-colony dates ---------------------------------------------------------
## Set seal to be at colony in between foraging trips

tmp <- Map(function(trac, trip){
  Map(function(t1, t2){which(trac$Date >= t1 & trac$Date <= t2)},
      as.list(trip$LastDate), as.list(lead(trip$StartDate)))
}, split(tracks, tracks$id), split(ft, ft$id))
tmp <- lapply(tmp, unlist) # get index of rows for each seal (as a list) that are to become at the colony
tracks <- Map(function(x, y){
  x$colony <- 0
  x$colony[y] <- 1
  return(x)
}, split(tracks, tracks$id), tmp)
tracks <- do.call(rbind, tracks)
tracks$Lon[tracks$id != '307' & tracks$id != '317' & tracks$id != '324' & tracks$id != '340' & tracks$colony == 1] <- 137.4680
tracks$Lat[tracks$id != '307' & tracks$id != '317' & tracks$id != '324' & tracks$id != '340' & tracks$colony == 1] <- -36.06900


## add trip id to locations 
tracks$trip <- NA
tracks2 <- map(split(ft, 1:nrow(ft)), function(x){
  start <- x$StartDate
  end <- x$LastDate
  trip <- x$X1
  ID <- x$id
  y <- with(tracks, tracks[Date >= start & Date <= end & id == ID,])
  y$trip <- trip
  return(y)
}) %>%  reduce(bind_rows)

tracks3 <- tracks[!tracks$Date %in% tracks2$Date,]
tracks2 <- bind_rows(tracks2, tracks3) %>% arrange(id, Date)
tracks <- tracks2
# save(tracks, file = "glsTracksAll_cleaned.RData" )

# Import Bonny Upwelling SSTA 16-17 ---------------------------------------
dir(pattern = 'bonney')
load('bonneyupwelling_ssta_2016-2017.Rdata')
load("bonneyupwelling_chlanom_2016-2017.Rdata")
# bud$date <- as.POSIXct(strptime(bud$date, "%Y-%m-%d"), tz = "GMT")
# bud <- buchld
bud$date <- as.POSIXct(strptime(bud$date, "%Y-%m-%d"), tz = "GMT")


# Plot glsSST x IMOS x Latitude ------------------------------------------
### Do for all animals - combine multiple plots in one. t1 = sst data, t2 =
### latitude data
## function for plotting sst vs imos vs lat
sstLatPlot <- function(t1, t2, imos, xmin, xmax, bu){
  maxt <- max(t1$temperature, na.rm = T)
  id <- t1$id[1]
  #animal temp vs imos temp vs sstanom 
  # E.g. To set Temperature 0 – 30 and Precipitation 0 – 400, scale Precipitation by multiplying 30 / 400 to fit range of Temperature
  g1 <- ggplot() +
    geom_line(data = t1, aes(x = Date, y = temperature, color = variable), alpha = 0.5, size = 0.5) +
    geom_point(data = imos, aes(x = Date, y = temperature, color = variable), alpha = 0.7, size = 0.5)+
    scale_color_brewer(palette = "Paired") +
    scale_x_datetime(date_labels = "%b",date_breaks = '1 month', limits = c(xmin, xmax)) +
    labs(y = "Temp") +
    geom_label(aes(x = xmax - 60*60*24*5, y = 26, label = id), nudge_y = -3, size = 3 ) +
    geom_point(data = bu, aes(x = date, y = v1 * -10/1), alpha = 0.7, size = 0.5) +
    scale_y_continuous(sec.axis = sec_axis(~./-10, name = "SSTa"), limits = c(10, 26)) +
    theme(
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(1, 1, 0, 1), "mm"),
      text = element_text(size = 7),
      panel.background = element_rect(fill = "white", colour = 'black')
    )
  
  ## to include chl anomaly 
  # g1 <- ggplot() + 
  #   geom_line(data = t1, aes(x = Date, y = temperature, color = variable), alpha = 0.5, size = 0.5) +
  #   geom_point(data = imos, aes(x = Date, y = temperature, color = variable), alpha = 0.7, size = 0.5)+
  #   scale_color_brewer(palette = "Paired") +
  #   scale_x_datetime(date_labels = "%b",date_breaks = '1 month', limits = c(xmin, xmax)) +
  #   labs(y = "Temp") +
  #   geom_label(aes(x = xmax - 60*60*24*5, y = 26, label = id), nudge_y = -3, size = 3 ) + 
  #   geom_point(data = bu, aes(x = date, y = v1 /0.6), color = 'darkgrey', alpha = 0.6, size = 0.3) + 
  #   scale_y_continuous(sec.axis = sec_axis(~.*0.6, name = "aCHL"), limits = c(0, 26)) + 
  #   theme(
  #     legend.position = 'none',
  #     axis.title.x = element_blank(),
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     plot.margin = unit(c(1, 1, 0, 1), "mm"),
  #     text = element_text(size = 7),
  #     panel.background = element_rect(fill = "white", colour = 'black')
  #   )
  
  #latitude plot
  g2 <- ggplot() +
    geom_path(data = t2, aes(x = Date, y = Lat)) + 
    scale_x_datetime(date_labels = "%b",date_breaks = '1 month', date_minor_breaks = '5 days', limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(-46, -36), breaks = seq(-46, -36, by = 3)) +
    theme(
      plot.title = element_blank(), 
      plot.margin = unit(c(0, 1, 1, 1), "mm"),
      text = element_text(size = 7),
      panel.background = element_rect(fill = "white", colour = 'black'))
  
  return(list(g1,g2))
}

# ymin <- min(sst$MinSST, na.rm = T)
# ymax <- max(sst$MaxSST, na.rm = T)
## first for 2017 seals
sst$id <- factor(sst$id, levels = c(305, 311, 315, 319, 322, 326, 351, 353, 307, 317, 324, 340))
sstg <- gather(sst, variable, temperature, MinSST, MaxSST, SST)
sstg$variable <- as.factor(sstg$variable)
imos2g <- gather(imos2, variable, temperature, tmin, tmax, tmean)
imos2g$variable <- as.factor(imos2g$variable)
xmin <- as.POSIXct('2017-02-01 00:00:00', tz = "GMT")
xmax <- as.POSIXct('2017-08-10 00:00:00', tz = "GMT")
d1 <- split(sstg, sstg$id)
d2 <- droplevels.data.frame(tracks %>% filter(year(Date) == 2017))
d2$id <- factor(d2$id, levels = c(305, 311, 315, 319, 322, 326, 351, 353, 307, 317, 324, 340))
d2 <- split(d2, d2$id)
bu <- bud %>% filter(lubridate::year(date) == 2017)
plots <- mapply(sstLatPlot, t1 = d1, t2 = d2, list(imos2g), xmax = xmax, xmin = xmin, bu = list(bu))
plots[1]
a <- ggarrange(plotlist = plots[1:12], heights = c(1, 1),
          ncol = 1, nrow = 12, align = 'v')
b <- ggarrange(plotlist = plots[13:24], heights = c(1, 1),
               ncol = 1, nrow = 12, align = 'v')
# ggarrange(a,b,ncol = 2, nrow = 1)

## next for 2016 seals
degg <- gather(deg, variable, temperature, MinSST, MaxSST, SST)
degg$variable <- as.factor(degg$variable)
imos1g <- gather(imos1, variable, temperature, tmin, tmax, tmean)
imos1g$variable <- as.factor(imos1g$variable)
# imos1g$variable <- factor(imos1g$variable, levels = c('tmax', 'tmean', 'tmin'),labels = c('Shelf Max', 'Shelf Mean','Shelf Min'), ordered = FALSE)
xmin <- as.POSIXct('2016-01-31 00:00:00', tz = "GMT")
xmax <- as.POSIXct('2016-09-30 00:00:00', tz = "GMT")
d1 <- split(degg, degg$id)
d2 <- droplevels.data.frame(tracks %>% filter(year(Date) == 2016))
d2 <- split(d2, d2$id)
bu <- bud %>% filter(lubridate::year(date) == 2016)
plots2 <- mapply(sstLatPlot, t1 = d1, t2 = d2, list(imos1g), xmax = xmax, xmin = xmin, bu = list(bu))
plots2[1]

c <- ggarrange(plotlist = plots2, heights = c(1, 1), ncol = 1, nrow = 10, align = 'v')
## Get legend
mypal <- brewer.pal(6, 'Paired')
# display.brewer.pal(6, 'Paired')
p <- ggplot() + 
  geom_line(data = d1[[1]], aes(x = Date, y = temperature, color = variable), alpha = 0.5, size = 0.5) +
  geom_point(data = imos1g, aes(x = Date, y = temperature, color = variable), alpha = 0.7, size = 0.5)+
  scale_color_manual(values = mypal, name = 'SST',
                     breaks = c('MaxSST', 'tmax',  'SST',  'tmean', 'MinSST', 'tmin'),
                       labels = c('MaxSST', 'Shelf Max', 'SST', 'Shelf Mean', 'MinSST', 'Shelf Min' )) +
  scale_x_datetime(date_labels = "%b",date_breaks = '1 month', limits = c(xmin, xmax)) +
  labs(y = "Temp") +
  theme(
    legend.position = 'top',
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 1, 0, 1), "mm"),
    text = element_text(size = 8)
  )

legend <- get_legend(p)
plots2[[11]] <- legend

## Combine 2016 and 2017 seal tracks
d <- ggarrange(c, legend, heights = c(10, 2), ncol = 1, nrow = 2)

tiff(filename = 'SSTxShelfSSTxLATxSSTa.tiff',  width=10, height=7, units= "in", res = 300)
ggarrange(a,b,d, ncol = 3, nrow = 1)
dev.off()


# LEF Grant Plot ----------------------------------------------------------
# only for 2017 seals - ASSTa + CTD + Latitude
# CTD plot
# imos temp vs sstanom 
# E.g. To set Temperature 0 – 30 and Precipitation 0 – 400, scale Precipitation by multiplying 30 / 400 to fit range of Temperature
imos2g <- gather(imos2, variable, temperature, tmin, tmax, tmean)
imos2g$variable <- as.factor(imos2g$variable)
bu <- bud %>% filter(lubridate::year(date) == 2017)
xmin <- as.POSIXct('2017-02-01 00:00:00', tz = "GMT")
xmax <- as.POSIXct('2017-08-10 00:00:00', tz = "GMT")
p1 <- ggplot() +
  geom_point(data = imos2g, aes(x = Date, y = temperature, color = variable), alpha = 0.7, size = 1)+
  # scale_color_brewer(palette = "Paired") +
  scale_colour_manual(values = wesanderson::wes_palette('Darjeeling', 3)) + 
  scale_x_datetime(date_labels = "%b",date_breaks = '1 month', limits = c(xmin, xmax)) +
  labs(y = "Temp") +
  geom_point(data = bu, aes(x = date, y = v1 * -10/1), alpha = 0.7, size = 1) +
  geom_vline(aes(xintercept = as.POSIXct('2017-04-03')), colour = 'black') + 
  scale_y_continuous(sec.axis = sec_axis(~./-10, name = "SSTa")) +
  theme(
    legend.position = c(0.75, 0.75),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 9),
    panel.background = element_rect(fill = "white", colour = 'black'),
    panel.border = element_rect(linetype = "solid", fill = NA),
    legend.key.size = unit(.25,'cm'),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.background = element_rect(fill = alpha('white', 0), colour = NA)
  ) 

# latitude plots
LatPlot <- function(t2, xmin, xmax){
  id <- t2$id[1]
  #latitude plot
  g2 <- ggplot() +
    geom_path(data = t2, aes(x = Date, y = Lat)) + 
    scale_x_datetime(date_labels = "%b",date_breaks = '1 month', date_minor_breaks = '5 days', limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(-44, -36), breaks = seq(-44, -36, by = 3)) +
    geom_label(aes(x = xmax - 60*60*24*5, y = -42, label = id), size = 3 ) +
    theme(
      plot.title = element_blank(),
      text = element_text(size = 9),
      panel.background = element_rect(fill = "white", colour = 'black'),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  return(list(g2))
}

d2 <- droplevels.data.frame(tracks %>% filter(year(Date) == 2017, 
                                              id != '307', id != '317', id != '324',id!= '340', id != '315', id!='353'))
d2$id <- factor(d2$id, levels = c(305, 311, 319, 322, 326, 351))
d2 <- split(d2, d2$id)
latplots <- mapply(LatPlot, t2 = d2, xmax = xmax, xmin = xmin)
plotlist <- list(p1)
plotlist <- c(plotlist, latplots)
a <- ggarrange(plotlist = plotlist, ncol = 1, nrow = 7, align = 'v', 
               heights = c(1.5,1,1,1,1,1,1) )
# ylab <-  as_ggplot(text_grob("Temp (°C)", rot=90, size = 10, vjust = 5)) 
# ylab2 <-  as_ggplot(text_grob("Latitude", rot=90, size = 10, vjust = 5)) # makes common ylab
# llabs <- ggarrange(ylab, ylab2, ncol = 1, nrow = 2, align = 'v', heights = c(1,7))
# ylab3 <- as_ggplot(text_grob(expression("Shelf SST"[a]*' (°C)'), rot=-90, size = 10, vjust = -10))
# rlab <- ggarrange(ylab3, ncol = 1, nrow = 2, align = 'v', heights = c(1,4))
# b <- ggarrange(llabs, a, rlab, ncol = 3,  align = 'h', widths = c(1,8,1))
# c <- ggarrange(b, rlab, ncol = 2,  align = 'h', widths = c(8,1))
# xlab <- as_ggplot(text_grob("Date", size = 10)) # makes common xlab

# combine ctd and lat plots
jpeg(filename = 'CTDxSSTaxLAT_2017 (LEF Grant).jpeg',  width = 7, height=10, units= "in", res = 300)
ggarrange(plotlist = plotlist, ncol = 1, nrow = 7, align = 'v', 
          heights = c(1.5,1,1,1,1,1,1) )
dev.off()
