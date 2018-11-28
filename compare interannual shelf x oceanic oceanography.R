## Compare inter-annual oceanographic properties on shelf / oceanic regions 

library(tidyverse)
library(lubridate)
library(marmap)
library(mapdata)
library(purrrlyr)

rm(list = ls())
# load data ---------------------------------------------------------------
setwd('~/Dropbox/OneDrive - University of Tasmania/phd/LNFS research/LNFS processed data/foraging switch timing/enviro data/')
## SSTA data
load("SSTAnom_131.5-43E_36-44S_2000-01_2005.Rdata")
a <- ssta
load("SSTAnom_131.5-43E_36-44S_2016-17.Rdata" )
b <- ssta
c <- tbl_df(rbind(a,b))
ssta <- c %>% 
  mutate(year = year(date))

## SST data
load("SST_131.5-43E_36-44S_2000-01_2005.Rdata")
a <- sst
load("SST_131.5-43E_36-44S_2016-17.Rdata" )
b <- sst
c <- tbl_df(rbind(a,b))
sst <- c %>% 
  mutate(year = year(date))

## SSTA shelf data
load("SSTAnom_shelf(isobath2000)_2000-01_2005.Rdata")
a <- ssta_shelf
load("SSTAnom_shelf(isobath2000)_2016-17.Rdata" )
b <- ssta_shelf
c <- tbl_df(rbind(a,b))
ssta_shelf <- c %>% 
  mutate(year = year(date))

## SSTA oceanic data
load("SSTAnom_oceanic(isobath2000)_2000-01_2005.Rdata")
a <- ssta_oce
load("SSTAnom_oceanic(isobath2000)_2016-17.Rdata" )
b <- ssta_oce
c <- tbl_df(rbind(a,b))
ssta_oce <- c %>% 
  mutate(year = year(date))

## SSHA shelf data
load("SSHA_shelf(isobath2000)_2000-01_2005.Rdata")
a <- ssha_shelf
load("SSHA_shelf(isobath2000)_2016-17.Rdata" )
b <- ssha_shelf
c <- tbl_df(rbind(a,b))
ssha_shelf <- c %>% 
  mutate(year = year(date))

## SSHA oceanic data
load("SSHA_oceanic(isobath2000)_2000-01_2005.Rdata")
a <- ssha_oce
load("SSHA_oceanic(isobath2000)_2016-17.Rdata" )
b <- ssha_oce
c <- tbl_df(rbind(a,b))
ssha_oce <- c %>% 
  mutate(year = year(date))

## SST shelf data
load("SST_shelf(isobath2000)_2000-01_2005.Rdata")
a <- sst_shelf
load("SST_shelf(isobath2000)_2016-17.Rdata" )
b <- sst_shelf
c <- tbl_df(rbind(a,b))
sst_shelf <- c %>% 
  mutate(year = year(date))

## SST oceanic data
load("SST_oceanic(isobath2000)_2000-01_2005.Rdata")
a <- sst_oce
load("SST_oceanic(isobath2000)_2016-17.Rdata" )
b <- sst_oce
c <- tbl_df(rbind(a,b))
sst_oce <- c %>% 
  mutate(year = year(date))


# plots -------------------------------------------------------------------
## get map 
map <- map_data("worldHires")

## stf
stf <- sst %>% group_by(x, y, year) %>% 
  summarise(v1 = mean(v1)) %>% 
  filter(v1 >= 12 & v1 < 12+1) 

# get bathy data
bathy <- getNOAA.bathy(130,143,-44, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
sb <- b2000 %>% group_by(x) %>% summarise(y = min(y))


ssta %>% group_by(x, y, year) %>% 
  summarise(mean = mean(v1), sd = sd(v1)) %>% 
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ year) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.2, fill = 'black') +
  scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'SSTA')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-44, -36) +
  xlim(131.5, 143) +
  labs(y = 'Lat', x = 'Lon') +
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5)
  


# shelf x oceanic summary -----------------------------------------------------------
## combine shelf and oceanic datasets
ssta_oce <- ssta_oce %>% mutate(loc = 'oceanic')
ssta_shelf <- ssta_shelf %>% mutate(loc = 'shelf')
ssta_so <- rbind(ssta_shelf, ssta_oce)
ssha_oce <- ssha_oce %>% mutate(loc = 'oceanic')
ssha_shelf <- ssha_shelf %>% mutate(loc = 'shelf')
ssha_so <- rbind(ssha_shelf, ssha_oce)
sst_oce <- sst_oce %>% mutate(loc = 'oceanic')
sst_shelf <- sst_shelf %>% mutate(loc = 'shelf')
sst_so <- rbind(sst_shelf, sst_oce)

## plots by year
setwd('~/Dropbox/OneDrive - University of Tasmania/phd/LNFS research/LNFS processed data/foraging switch timing/')
mytheme <- theme(
  text = element_text(size = 8.5),
  panel.background = element_rect(fill = "white", colour = 'black'),
  panel.border = element_rect(linetype = "solid", fill = NA),
  axis.text = element_text(size = 8.5)
)

ssta_so %>% 
  mutate(loc = factor(loc, levels = c('shelf', 'oceanic'))) %>% 
  group_by(year, date, loc) %>% 
  filter(v1 <= -1) %>%
  summarise(ncell = n(), mean = mean(v1)) %>% 
  mutate(dm = as.Date(format(date, '2000-%m-%d'))) %>% 
  ggplot(aes(x = dm, y = mean)) +
  geom_point(aes(size = ncell, colour = loc), alpha = 0.4) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') + 
  facet_wrap(~year) +
  labs(y = 'Mean SSTc \n (°C)', x = 'month') +
  mytheme + 
  scale_colour_brewer(palette = 'Set1', name = 'Location') +
  scale_size(name = 'No. cells (scaled)') +
  theme(legend.position = c(0.85, 0.3),
        legend.background = element_rect(fill = alpha('white', 0), colour = NA),
        legend.key.size = unit(.25,'cm'), 
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        legend.direction = 'horizontal',
        legend.spacing = unit(-.5, 'lines')) -> p1

png(filename = './plots/interannual_SSTc_2001-17.png', units = 'in', width = 10, height = 8, res = 300)
p1
dev.off()


ssha_so %>% 
  filter(loc == 'shelf') %>% 
  group_by(year, date) %>% 
  summarise( min = min(v1)) %>% 
  mutate(dm = as.Date(format(date, '2000-%m-%d'))) %>% 
  ggplot(aes(x = dm, y = min)) +
  geom_point(size = 0.5) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') + 
  facet_wrap(~year) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'red') +
  labs(y = 'min SSHA (m)') +
  mytheme -> p2

png(filename = './plots/interannual_minSSHA_2001-17.png', units = 'in', width = 10, height = 8, res = 300)
p2
dev.off()


## t-tests (using wilcox.. t test because non-parametric)
### is ssta shelf and oceanic different between 2016 and 2017?
ssta_so %>% filter(year %in% c(2016, 2017)) %>% 
  slice_rows('loc') %>% 
  by_slice(partial(wilcox.test, v1~year)) -> test

test$.out #shows results, arranged by oceanic, shelf. both are significantly different between years. 

ssta_so %>% group_by(year, loc) %>%
  ggplot(aes(x = factor(year), y = v1)) +
  geom_boxplot(aes(fill = loc), position = 'dodge')

ssha_so %>% group_by(year, loc) %>%
  summarise(mean = mean(v1), sd = sd(v1)) %>% 
  ggplot(aes(x = factor(year), y = mean)) +
  geom_col(aes(fill = loc), position = 'dodge') +
  geom_point(aes(y = sd, colour = loc))

sst_so %>% group_by(year, loc) %>%
  summarise(mean = mean(v1), sd = sd(v1)) %>% 
  ggplot(aes(x = factor(year), y = mean)) +
  geom_col(aes(fill = loc), position = 'dodge') +
  geom_point(aes(y = sd, colour = loc))


