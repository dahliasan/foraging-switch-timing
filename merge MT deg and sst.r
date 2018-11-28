## Merge all .deg .sst files ##

library(BAStag)
library(dplyr)
library(stringr)

## Import deg/sst files
degfiles <- dir(pattern = ".deg")
deg <- lapply(degfiles, function(x){
 d <- tbl_df(readMTdeg(x))
 d$id <- as.factor(str_split(x, pattern = '_')[[1]][1])
 return(d)
})

deg <- do.call(rbind,deg)

sstfiles <- dir(pattern = '.sst')
sst <- lapply(sstfiles, function(x){
 d <- tbl_df(readMTsst(x))
 d$id <- as.factor(str_split(x, pattern = '_')[[1]][1])
 return(d)
})

sst <- do.call(rbind, sst)


save(deg, file = 'degAll2016.Rdata')
save(sst, file = 'sstAll2017.Rdata')
