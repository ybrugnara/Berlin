## Homogenize each target_series and then merge them together by adjusting them with respect to the most recent one
## Breakpoints and reference series for each breakpoint are read from the breakpoints.csv file
## NB1: When two breakpoints are too close to each other (< 6 months) the data inbetween is deleted instead of adjusted
## NB2: When there is only one reference series, a constant adjustment is applied (exception: auxiliary reference series)

library(lubridate)
library(RColorBrewer)
source("functions.R")
load("data.RData")

inpath <- "../data"
outpath <- "../homogenized"
plotpath <- "../homogenized/adjustments"
reanalysis <- "../20CR.txt"
bp_file <- "breakpoints.csv"
startyear <- 1720
endyear <- 2023
nmin <- 24*30 # minimum number of overlapping days between two segments under which the auxiliary reference series is used (if available, otherwise use constant adjustment)
aux_span <- 10*365*24*3600 # maximum length before/after breakpoint of auxiliary series to use (in seconds)
colors <- brewer.pal(7, "Set1")


target_series <- c(
  "DWD_BerlinDahlem_18760101_20231231_00403.txt",
  "Jena_18240101-19401231_p_daily.tsv",
  reanalysis,
  "Brazdil_Zitenice_18001101-18181231_p.tsv",
  "PALAEO-RA_Berlin_Gronau_17750101-18011015_p.tsv",
  "PALAEO-RA_Palatine-Society_Berlin_17810701-17881231_p_daily.tsv",
  "WEAR_Berlin_Brand_17601201_17941227_p_daily.tsv",
  "PALAEO-RA_Berlin_Lambert_17691216-17770925_p.tsv",
  "PALAEO-RA_Europe_Berlin_1_17210123-17740430_p_daily.tsv",
  "PALAEO-RA_Kanold_Loebau_17200101-17300630_p_daily.tsv",
  "PALAEO-RA_Europe_Wittenberg_1_17280320-17290321_p_daily.tsv",
  "PALAEO-RA_Berlin_Jablonski_17250112-17270111_p.tsv"
  )

aux_series <- "HCLIM_Prague_178801-200412_p.tsv" # series to use as reference in the merging when the overlap between two segments is too short
super_series <- c(
  "PALAEO-RA_Palatine-Society_Berlin_17810701-17881231_p_daily",
  "Brazdil_Zitenice_18001101-18181231_p"
  ) # These series overwrite the preceding ones during the merging process


## Read breakpoints
breakpoints_target <- read.csv(bp_file, stringsAsFactors=FALSE,
                               colClasses=c("character","Date","character","character","Date","Date"))


## Remove all files in output folder
file.remove(list.files(outpath, full.names=TRUE, pattern="[.]hom"))


## Read in daily series to be homogenized (fill in list x of homogeneous target segments)
x <- breakpoints <- list()
for (i in 1:length(target_series)) {
  if (target_series[i] == reanalysis) {
    xtmp <- read.table(reanalysis, sep="\t", header=TRUE)
    xtmp <- xtmp[, c("YR","M","D","Berlin")]
    xtmp <- xtmp[xtmp$YR<1830, ]
    xtmp$Berlin <- round(xtmp$Berlin / 100, 1)
    names(xtmp) <- c("Year", "Month", "Day", "Value")
  } else {
    xtmp <- read_daily(inpath, list.files(inpath,pattern=target_series[i]))
    xtmp <- xtmp[which(!is.na(xtmp$Value)), ]
  }
  xtmp$Time <- ISOdate(xtmp$Year, xtmp$Month, xtmp$Day)
  xtmp$julian <- as.integer(format(xtmp$Time,"%j"))
  xtmp$station <- sub("../", "", target_series[i], fixed=TRUE)
  xtmp$Adjustment <- 0
  
  ## Read breakpoints
  breakpoints[[target_series[i]]] <- breakpoints_target[breakpoints_target$Filename==target_series[i], ]
  bps <- sort(breakpoints[[target_series[i]]]$Date)
  bps <- c(min(xtmp$Time), bps, max(xtmp$Time))
  bps <- bps[bps>=bps[1] & bps<=bps[length(bps)] & !duplicated(bps)]
  #if (i == 1) bps <- bps[bps < max(xtmp$Time)]
  
  ## Loop on breakpoints
  j <- 0 # segment index
  if (length(bps) > 2) {
    for (i_bp in length(bps):2) {
      j <- j + 1
      x[[target_series[i]]][[j]] <- xtmp[xtmp$Time<bps[i_bp]&xtmp$Time>=bps[i_bp-1], ]
    }
  } else {
    j <- j + 1
    x[[sub("../","",target_series[i],fixed=TRUE)]][[j]] <- xtmp
  }
}

## Read in monthly series to be homogenized (fill in list xmon)
xmon <- list()
for (nam in target_series) {
  if (nam != sub("../","",reanalysis,fixed=TRUE)) {
    xmon[[nam]] <- Data[[nam]]
  }
}

## Read in reference series (fill in list y)
y <- list() # list of homogeneous reference segments
Refs <- Data
#for (nam in target_series) Refs[[nam]] <- NULL
  
for (ref in names(Refs)) {
  ytmp <- Refs[[ref]]
  y[[ref]] <- ytmp[which(!is.na(ytmp$Value)), ]
}


# Homogenization
months <- 1:12
days <- 1:365
trigon <- cbind(cos(2*pi*(months/12)), sin(2*pi*(months/12)), 
                cos(4*pi*(months/12)), sin(4*pi*(months/12)))
trigon_daily <- cbind(cos(2*pi*(days/365)), sin(2*pi*(days/365)), 
                      cos(4*pi*(days/365)), sin(4*pi*(days/365)))
xhom <- list()

for (stat in sub("../","",target_series,fixed=TRUE)) {
  
  if (length(x[[stat]]) > 1) {
    for (i in 2:length(x[[stat]])) {
      
      message("WORKING ON ", stat, ": ", 
              min(x[[stat]][[i]]$Year), "-", max(x[[stat]][[i]]$Year))
      
      ## Read reference series to use
      ref_string <- gsub(" ", "", breakpoints[[stat]]$Ref_series[i-1]) # pipe-separated station names
      ref_series <- strsplit(ref_string, "[|]")[[1]]
      n_refs <- length(ref_series)
      
      ## Read min and max dates to use
      dmin <- breakpoints[[stat]]$Ref_start[i-1]
      dmax <- breakpoints[[stat]]$Ref_end[i-1]
      
      ## Merge all segments homogenized so far into one and calculate monthly means
      xhom[[stat]] <- Reduce(rbind, x[[stat]][1:(i-1)])
      xhom[[stat]] <- xhom[[stat]][!duplicated(xhom[[stat]]$Time), ]
      xhom[[stat]] <- monthly_means(xhom[[stat]])
      
      k <- grep(ref_string, names(y), ignore.case=TRUE)
      if (length(k) > 0) {
        ref_adjustments <- array(dim=c(length(k),12))
        for (ik in 1:length(k)) {
          xref1 <- y[[k[ik]]][y[[k[ik]]]$Time%in%xhom[[stat]]$Time & y[[k[ik]]]$Time>max(x[[stat]][[i]]$Time) & y[[k[ik]]]$Time<=dmax, ]
          xref2 <- y[[k[ik]]][y[[k[ik]]]$Time%in%x[[stat]][[i]]$Time & y[[k[ik]]]$Time>=dmin, ]
          x1 <- xhom[[stat]][xhom[[stat]]$Time%in%xref1$Time, ] # series to be homogenized against
          x2 <- xmon[[stat]][xmon[[stat]]$Time%in%xref2$Time, ] # series to homogenize
          for (im in months) {
            ref_adjustments[ik,im] <- mean(x1$Value[x1$Month==im]-xref1$Value[xref1$Month==im],na.rm=TRUE) -
              mean(x2$Value[x2$Month==im]-xref2$Value[xref2$Month==im],na.rm=TRUE)
          }
        }
        rownames(ref_adjustments) <- names(y)[k]
        adjustments <- round(ref_adjustments, 1)
      } else {
        stop("Reference series not found")
      }
      # }
      print(adjustments)
      
      ## Take median and smooth adjustments
      if (length(k) > 1) {
        #if (!is.null(dim(adjustments))) {
        adjustments_med <- apply(adjustments, 2, median, na.rm=TRUE)
        #}
        params <- lm(adjustments_med ~ trigon)$coefficients
        adjustments_med <- params[1] + params[2]*trigon[,1] + params[3]*trigon[,2] + 
          params[4]*trigon[,3] + params[5]*trigon[,4]
        adjustments_daily <- params[1] + params[2]*trigon_daily[,1] + params[3]*trigon_daily[,2] + 
          params[4]*trigon_daily[,3] + params[5]*trigon_daily[,4]
      } else if (sum(!is.na(colMeans(adjustments))) >= 6) { # adjustments based at least on 6 months
        adjustments_med <- rep(mean(adjustments, na.rm=TRUE), 12)
        adjustments_daily <- rep(mean(adjustments, na.rm=TRUE), 365)
      } else {
        adjustments_med <- rep(NA, 12)
        adjustments_daily <- rep(NA, 365)
      }
      print(adjustments_med)
      
      ## Apply daily adjustments
      adjustments_daily[366] <- adjustments_daily[365]
      for (j in 1:366) {
        k <- which(x[[stat]][[i]]$julian==j)
        x[[stat]][[i]]$Value[k] <- x[[stat]][[i]]$Value[k] + round(adjustments_daily[j], 1)
        x[[stat]][[i]]$Adjustment[k] <- round(adjustments_daily[j], 1)
      }
      
      ## Plot adjustments
      if (sum(!is.na(adjustments_med)) > 0) {
        plotname <- paste0(plotpath, "/", stat, "_",
                           min(x[[stat]][[i]]$Year), "-", max(x[[stat]][[i]]$Year), ".pdf")
        pdf(plotname, width=5, height=4)
        par(mar=c(2.5,2.5,2,1), mgp=c(1.5,0.7,0))
        plot(1:12, adjustments_med, type="l", lwd=3, 
             ylim=c(min(adjustments,na.rm=TRUE),max(adjustments,na.rm=TRUE)),
             xlab="Month", ylab="hPa", 
             main=paste0(min(x[[stat]][[i]]$Year), " - ", max(x[[stat]][[i]]$Year)))
        for (i_ref in 1:nrow(adjustments)) {
          lines(1:12, adjustments[i_ref,], col=colors[i_ref])
        }
        lines(seq(0.5,12.5,length.out=366), adjustments_daily, lwd=2, col="royalblue")
        grid(lwd=0.3)
        legend("topleft", c(rownames(adjustments),"Daily","Median"), cex=0.5,
               lwd=c(rep(1,i_ref),2,3), col=c(colors[1:i_ref],"royalblue","black"))
        dev.off()
      }
      
    }
  }
  
}


## Write homogenized data
i <- 0
for (stat in sub("../","",target_series,fixed=TRUE)) {
  i <- i + 1
  out <- Reduce(rbind, x[[stat]])
  out <- out[!is.na(out$Value), ]
  out <- out[order(out$Time), ]
  filename <- paste(formatC(i, width=2, flag=0), sub("tsv|txt", "hom", stat), sep="_")
  out <- out[out$station==stat, c("Year","Month","Day","Value","Adjustment")]
  write.csv(out, file=paste(outpath,filename,sep="/"), row.names=FALSE)
}


## Merge homogeneous segments
message("Merging data...")
for (stat in list.files(outpath, pattern="[.]hom") ) {
  
  if (substr(stat, 1, 2) == "01") {
    out_merged <- read.csv(paste(outpath, stat, sep="/"))
    out_merged$Source <- stat
    out_merged$Time <- ISOdate(out_merged$Year, out_merged$Month, out_merged$Day)
  
  } else {
      
    ## Calculate adjustments
    tmp <- read.csv(paste(outpath, stat, sep="/"))
    tmp$Source <- stat
    tmp$Time <- ISOdate(tmp$Year, tmp$Month, tmp$Day)
    tmp$julian <- as.integer(format(tmp$Time,"%j"))
    cand <- tmp[tmp$Time%in%out_merged$Time, ]
    if (nrow(cand) >= nmin | max(tmp$Time) < min(y[[aux_series]]$Time)) {
      ref <- out_merged[out_merged$Time%in%tmp$Time, ]
      adjustments <- aggregate(ref$Value - cand$Value, list(cand$Month), mean, na.rm=TRUE)
      if (nrow(adjustments) < 12) {
        for (i in 1:12) {
          if (!i %in% adjustments[,1]) {
            adjustments <- rbind(adjustments, c(i, NA))
          }
        }
      }
      adjustments <- adjustments[order(adjustments[,1]), "x"]
      label <- "Using overlap"
    } else {
      adjustments <- rep(NA, 12)
      tmpmon <- monthly_means(tmp)
      outmon <- monthly_means(out_merged)
      ref <- y[[aux_series]][y[[aux_series]]$Time>=min(outmon$Time)-aux_span & 
                               y[[aux_series]]$Time<=min(outmon$Time)+aux_span, ]
      xref1 <- ref[ref$Time%in%outmon$Time, ]
      xref2 <- ref[ref$Time%in%tmpmon$Time, ]
      x1 <- outmon[outmon$Time%in%xref1$Time, ] # series to be homogenized against
      x2 <- tmpmon[tmpmon$Time%in%xref2$Time, ] # series to homogenize
      for (im in months) {
        adjustments[im] <- mean(x1$Value[x1$Month==im]-xref1$Value[xref1$Month==im],na.rm=TRUE) -
          mean(x2$Value[x2$Month==im]-xref2$Value[xref2$Month==im],na.rm=TRUE)
      }
      label <- paste("Using auxiliary series", aux_series)
    }
    
    ## Smooth adjustments
    if (nrow(cand) >= nmin) {
      params <- lm(adjustments ~ trigon)$coefficients
      adjustments_daily <- params[1] + params[2]*trigon_daily[,1] + params[3]*trigon_daily[,2] + 
        params[4]*trigon_daily[,3] + params[5]*trigon_daily[,4]
    } else {
      adjustments_daily <- rep(mean(adjustments,na.rm=TRUE), 365)
    }
    
    ## Apply daily adjustments
    adjustments_daily[366] <- adjustments_daily[365]
    for (j in 1:366) {
      k <- which(tmp$julian==j)
      tmp$Value[k] <- tmp$Value[k] + round(adjustments_daily[j], 1)
      tmp$Adjustment[k] <- tmp$Adjustment[k] + round(adjustments_daily[j], 1)
    }
    out_merged <- rbind(out_merged, tmp[,1:7])
    if (grepl(paste(super_series,collapse="|"), stat)) {
      out_merged <- out_merged[!duplicated(out_merged$Time, fromLast=TRUE), ]
    } else {
      out_merged <- out_merged[!duplicated(out_merged$Time), ]
    }
    
    ## Plot adjustments
    plotname <- paste0(plotpath, "/", stat, ".pdf")
    pdf(plotname, width=5, height=4)
    par(mar=c(2.5,2.5,2,1), mgp=c(1.5,0.7,0))
    plot(1:12, adjustments, type="l", lwd=3, 
         ylim=c(min(adjustments,na.rm=TRUE),max(adjustments,na.rm=TRUE)),
         xlab="Month", ylab="hPa")
    lines(seq(0.5,12.5,length.out=366), adjustments_daily, lwd=2, col="royalblue")
    mtext(label)
    grid(lwd=0.3)
    legend("topleft", c("Daily","Monthly"), cex=0.5, lwd=c(2,3), col=c("royalblue","black"))
    dev.off()
  }
  
  last_stat <- stat
  
}

## Add NAs
all_dates <- seq(ISOdate(startyear,1,1), ISOdate(endyear,12,31), by="day")
df_nas <- data.frame(
  Year = year(all_dates),
  Month = month(all_dates),
  Day = day(all_dates),
  Value = NA,
  Adjustment = NA,
  Source = NA,
  Time = all_dates
)
out_merged <- rbind(out_merged, df_nas)
out_merged <- out_merged[out_merged$Year>=startyear & out_merged$Year<=endyear, ]
out_merged <- out_merged[!duplicated(out_merged$Time), ]

## Write merged data
out_merged <- out_merged[order(out_merged$Time), 1:6]
filename <- "Berlin_p_1720-2023.csv"
write.csv(out_merged, file=paste(outpath,filename,sep="/"), row.names=FALSE)

## Print statistics
for (stat in list.files(outpath, pattern="[.]hom")) {
  message(paste(stat, sum(!is.na(out_merged$Value[which(out_merged$Source==stat)]))))
}

## Plot comparison with Mod-ERA and standard deviation
#source("compare_ModERA.R")
#source("plot_stdev.R")
