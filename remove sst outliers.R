## Removing outlier SST from .deg and .sst files. Script used before plotting
## SST + Lat data.
library(dplyr)

## Load merged data
# load('degAll2016.Rdata')
# load('sstAll2017.Rdata')


## Remove specific seals
# sst <- sst %>% filter(id != "356", id != "357", id != "358", id != "359", id != "316")
# sst <- droplevels.data.frame(sst)

## Because SST in .deg files include temp when tag was dry, there are some
## unrealistic temp values, so we remove them by finding the outliers
# source("http://goo.gl/UUyEzD")
# par(mar=c(1,1,1,1))
# outlierKD(deg, MaxSST)
# outlierKD(deg, MinSST)
# outlierKD(deg, SST)
# save(deg, file = 'degAll2016(outliersRemoved).Rdata')

# outlierKD(sst, MaxSST)
# outlierKD(sst, MinSST)
# outlierKD(sst, SST)
# save(sst, file = "sstAll2017(outliersRemoved).Rdata")


# Remove temperatures when tag completely dry for 4h (tag records every 4h)
deg <- deg %>% filter(Wet > 0)