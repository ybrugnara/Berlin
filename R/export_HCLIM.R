# Export HCLIM data to one SEF file

library(dataresqc)

setwd("../.")

x <- read.csv("HCLIM-data.tsv", sep="\t")
meta <- read.csv("HCLIM-meta.tsv", sep="\t")

id <- 3747 # can be a vector of ids (data will be merged into a single SEF file)

out <- x[x$STATION_ID%in%id & x$PARAMETER=="PRES",4:16]
out <- data.frame(YEAR=rep(out$YEAR,each=12), MONTH=rep(1:12,nrow(out)), P=c(t(out[,2:13]),recursive=T))
out$DAY <- NA
out$HOUR <- NA
out$MINUTE <- NA
out <- out[,c("YEAR","MONTH","DAY","HOUR","MINUTE","P")]

lat <- meta$LATITUDE[meta$STATION_ID==id[1]]
lon <- meta$LONGITUDE[meta$STATION_ID==id[1]]
ele <- meta$ELEVATION[meta$STATION_ID==id[1]]
nam <- meta$NAME[meta$STATION_ID==id[1]]

write_sef(out,"data","p",nam,nam,lat,lon,ele,"HCLIM",units="hpa",stat="mean",period="month",keep_na=T)
