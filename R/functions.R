library(dataresqc)


# Read one data file (inpath/f) and calculate daily means
read_daily <- function(inpath, f) {
  
  col_names <- c("Year","Month","Day","Value")
  extension <- strsplit(f, "[.]")[[1]][2]
  if (extension == "tsv") {
    x <- read_sef(paste(inpath,f,sep="/"))
  } else {
    x <- read.csv(paste(inpath,f,sep="/"), sep=";")
    x$Year <- as.integer(substr(x$MESS_DATUM,1,4))
    x$Month <- as.integer(substr(x$MESS_DATUM,5,6))
    x$Day <- as.integer(substr(x$MESS_DATUM,7,8))
    x$Value <- x$PP_TER
  }
  x <- x[, col_names]
  x$Value[x$Value < 0] <- NA
  x$Day[is.na(x$Day)] <- 1
  xday <- aggregate(x$Value, list(x$Year,x$Month,x$Day), mean)
  names(xday) <- c("Year", "Month", "Day", "Value")
  xday <- xday[order(ISOdate(xday$Year,xday$Month,xday$Day)), ]
  xday$Value <- round(xday$Value, 1)
  
  return(xday)
  
}


# Calculate monthly means for months with at least nmin non-missing days
monthly_means <- function(X, nmin=20) {
  
  xmon <- aggregate(X$Value, list(X$Year,X$Month), mean, na.rm=TRUE)
  n <- aggregate(X$Value, list(X$Year,X$Month), function(x) sum(!is.na(x)))
  xmon$x[which(n$x<20)] <- NA
  names(xmon) <- c("Year", "Month", "Value")
  xmon$Time = ISOdate(xmon[,1], xmon[,2], 1)
  xmon <- xmon[, c("Time", "Year", "Month", "Value")]
  xmon <- xmon[order(xmon$Time), ]
  
  return(xmon)
  
}


# Read all files in the inpath folder with pattern patt and return list of monthly series
read_monthly <- function(inpath, patt=NULL) {

  files <- list.files(inpath, pattern=patt)
  out <- list()
  
  # Read files and calculate monthly means
  for (f in files) {
    xday <- read_daily(inpath, f)
    if (grepl("HCLIM|Stefan", f)) {
      xday$Time = ISOdate(xday[,1], xday[,2], 1)
      xmon <- xday[, c("Time", "Year", "Month", "Value")]
      out[[f]] <- xmon
    } else {
      out[[f]] <- monthly_means(xday)
    }
  }
  
  return(out)

}