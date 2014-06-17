


###############################################################################
###############################################################################

# You should be able to use this by editing this section only.

setwd("C:/Users/Work/AAA/Programming/ProgramOutput/Nino")


lat.range <- 13:14
lon.range <- 142:143

firstyear <- 1957
lastyear <- 1958

outputfilename <- paste0("Scotland-", firstyear, "-", lastyear, ".txt")


###############################################################################
###############################################################################

#                  Explanation


# 1. Use setwd() to set the working directory to the one 
# containing the .nc files such as air.sig995.1951.nc.
# Example:
# setwd("C:/Users/Work/AAA/Programming/ProgramOutput/Nino")

# 2. Supply the latitude and longitude range.  The NOAA data is 
# every 2.5 degrees. The ranges are supplied as the number of steps of this 
# size. For latitude, 1 means North Pole, 73 means South Pole. For longitude,
# 1 means 0 degrees East, 37 is 90E, 73 is 180, 109 is 90W or 270E, 144 is 2.5W.

# These roughly cover Scotland.
# lat.range <- 13:14
# lon.range <- 142:143

# These are the area used by Ludescher et al, 2013. It is 27x69 points
# which then subsampled to 9 by 23.
#lat.range <- 24:50
#lon.range <- 48:116 

# 3. Supply the years
# firstyear <- 1950
# lastyear <- 1952

# 4. Supply the output name as a text string. paste0() concatenates strings
# which you may find handy:
# outputfilename <- paste0("Pacific-", firstyear, "-", lastyear, ".txt")


###############################################################################
###############################################################################


#                      Example of output


#  S013E142 S013E143 S014E142 S014E143
#  Y1950P001 281.60000272654 281.570002727211 281.60000272654 280.970002740622
#  Y1950P002 280.740002745762 280.270002756268 281.070002738386 280.49000275135
#  Y1950P003 280.100002760068 278.820002788678 281.120002737269 280.070002760738
#  Y1950P004 281.070002738386 279.420002775267 281.620002726093 280.640002747998
#  ...
#  Y1950P193 285.450002640486 285.290002644062 285.720002634451 285.75000263378
#  Y1950P194 285.570002637804 285.640002636239 286.070002626628 286.570002615452
#  Y1950P195 285.92000262998 286.220002623275 286.200002623722 286.620002614334
#  ...
#  Y1950P364 276.100002849475 275.350002866238 276.37000284344 275.200002869591
#  Y1950P365 276.990002829581 275.820002855733 276.020002851263 274.72000288032
#  Y1951P001 278.220002802089 277.470002818853 276.700002836064 275.870002854615
#  Y1951P002 277.750002812594 276.890002831817 276.650002837181 275.520002862439
#  ...
#  Y1952P365 280.35000275448 280.120002759621 280.370002754033 279.390002775937

# There is one row for each day, and 365 days in each year (leap days are 
# omitted). In each row, you have temperatures in Kelvin for each grid 
# point in a rectangle.

# S13E142 means 13 steps South from the North Pole and 142 steps East from 
# Greenwich. The points are in reading order, starting at the top-left 
# (Northmost, Westmost) and going along the top row first.

# Y1950P001 means year 1950, day 1. (P because longer periods might be used 
# later.)

###############################################################################
###############################################################################

library(RNetCDF)


n.lat <- length(lat.range)
n.lon <- length(lon.range)
n.points <- n.lat * n.lon
smallest.lat <- lat.range[1]
smallest.lon <- lon.range[1]
biggest.lat <- lat.range[n.lat]
biggest.lon <- lon.range[n.lon]


# extract the data as a vector (not a matrix) at each time
# These functions do the translation
index.from.latlon <- function(i, j) { (i-lat.range[1]) * n.lon + (j-lon.range[1]) + 1 }
lon.from.index <- function(idx) { (idx-1) %% n.lon  +  lon.range[1] }
lat.from.index <- function(idx) { floor((idx-1) / n.lon)  +  lat.range[1] }

index.table <- matrix(0, nrow=biggest.lat, ncol=biggest.lon)
for (lat in lat.range) {
  for (lon in lon.range) {
    index.table[lat, lon] <- index.from.latlon(lat, lon)
  }
}


point.names <- function() {
  pointnames <- rep("", n.points)
  for (idx in 1:n.points) {
    lon <- lon.from.index(idx)
    lat <- lat.from.index(idx)
    pointnames[idx] <- paste0("S", sprintf("%03d", lat), "E", sprintf("%03d", lon))
    stopifnot(idx == index.from.latlon(lat, lon))
  }
  pointnames
}


make.Kvals.for.year <- function(year) { 
  onc <- open.nc(paste0("air.sig995.", year, ".nc"))
  rnc <- read.nc(onc)
  close.nc(onc)
  yearof.Kvals <- matrix(0, nrow=365, ncol=n.points)
  for (lat in lat.range) {
    for (lon in lon.range) { 
      yearof.Kvals[, index.table[lat, lon] ] <- rnc$air[ lon, lat, 1:365 ] # ignore leap days
    }
  }
  
  colnames(yearof.Kvals) <- point.names()
  rownames(yearof.Kvals) <- paste0("Y", year, "P", sprintf("%03d", 1:365))
  yearof.Kvals
}


setof.Kvals <- NULL
for (i in firstyear:lastyear) {
  setof.Kvals <- rbind(setof.Kvals, make.Kvals.for.year(i))
}

write.table(x=round(setof.Kvals, digits=2), file=outputfilename, quote=FALSE)
