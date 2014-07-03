setwd("/Users/macbookpro/home/elnino/data")

# functions by Graham Jones:

lat.from.label <- function(label) {
  as.integer(substr(label,2,4))
}

lon.from.label <- function(label) {
  as.integer(substr(label,6,8))
}

lat.range.from.colnames <- function(cnames) {
  lat.from.label(cnames[1]) : lat.from.label(cnames[length(cnames)])
}

lon.range.from.colnames <- function(cnames) {
  lon.from.label(cnames[1]) : lon.from.label(cnames[length(cnames)])
}  

textual.longitude.from.3D.lon <- function(lon, cnames) {
  lon.range <- lon.range.from.colnames(cnames)
  x <- lon.range[1] + lon - 2
  x <- 2.5 * x
  if (x < 179.99) {
    txt <- paste0(x, "E")
  } else if (x > 180.01) {
    txt <- paste0(360-x, "W")
  } else {
    txt <- "180"
  }
  txt
} 

textual.latitude.from.3D.lat <- function(lat, cnames) {
  lat.range <- lat.range.from.colnames(cnames)
  y <- lat.range[1] + lat - 2
  y <- 2.5 * y
  if (y < 89.99) {
    txt <- paste0(90-y, "N")
  } else if (y > 90.01) {
    txt <- paste0(y-90, "S")
  } else {
    txt <- "0"
  }
  txt
} 

data.to.3D <- function(vals) {
  lat.range <- lat.range.from.colnames(colnames(vals))
  lon.range <- lon.range.from.colnames(colnames(vals))
  n.times <- dim(vals)[1]
  n.lat <- length(lat.range)
  n.lon <- length(lon.range)
  vals3D <- array(0, dim=c(n.times, n.lat, n.lon))
  for (tim in 1:n.times) {
    for (lat in 1:n.lat) {
      vals3D[tim, lat, ] <- vals[tim, ((lat-1) * n.lon) +  (1:n.lon)]
    }
  }
  vals3D
}

seasonally.adjust <- function(vals) {
  stopifnot(dim(vals)[1] %% 365 == 0)
  n.years <- as.integer(dim(vals)[1]/365)
  offsets <- (0:(n.years-1))*365
  y.means <- array(0, dim=c(365, dim(vals)[2], ncols=dim(vals)[3]))
  for (d in 1:365) {
    for (lat in 1:dim(vals)[2]) {
      for (lon in 1:dim(vals)[3]) {
        y.means[d, lat, lon] <- mean(vals[d + offsets, lat, lon])
      }                   
    }
  }  
  for (y in 1:n.years ) {
    for (d in 1:365) {
      tim <- (y-1)*365 + d
      vals[tim, , ] <- vals[tim, , ] - y.means[d, , ]
    }
  }
  vals
}

# untested. aimed at values calculated eg every 14 days. 14*26=364, miss one day
seasonally.adjust.periods <- function(vals, nppy) {
  # nppy = number of periods per year
  stopifnot(dim(vals)[1] %% nppy == 0)
  n.years <- as.integer(dim(vals)[1]/nppy)
  offsets <- (0:(n.years-1))*nppy
  y.means <- array(0, dim=c(nppy, dim(vals)[2], ncols=dim(vals)[3]))
  for (p in 1:nppy) {
    for (lat in 1:dim(vals)[2]) {
      for (lon in 1:dim(vals)[3]) {
        y.means[p, lat, lon] <- mean(vals[p + offsets, lat, lon])
      }                   
    }
  }  
  for (y in 1:n.years ) {
    for (p in 1:nppy) {
      tim <- (y-1)*nppy + p
      vals[tim, , ] <- vals[tim, , ] - y.means[p, , ]
    }
  }
  vals
}

print("reading file...")

#Kvals <- as.matrix(read.table(file="Scotland-1950-1952.txt", header=TRUE))
Kvals <- as.matrix(read.table(file="Pacific-1950-1979.txt", header=TRUE))
#Kvals <- as.matrix(read.table(file="Pacific-1950-1952.txt", header=TRUE))
Kvals.cnames <- colnames(Kvals)

print("converting to 3D...")

Kvals.3D <- data.to.3D(Kvals)

print("seasonally adjusting...")

SAvals.3D <- seasonally.adjust(Kvals.3D)

# find covariances at time (day) p between each point x in region and a point
# offset dS steps South, dE steps East, and with time lag dp, and
# using a period length covp for the covariances. Positive
# dp means the points in the region are at time p, and
# the offset points are earlier, so is how the points
# in the region are influenced by others. (Negative dp means
# the offset point are at time p, not yet implemented.) 
covariances <- function(vals, region, p, dS, dE, dp, covp) {
  S0 <- region$S0
  S1 <- region$S1
  E0 <- region$E0
  E1 <- region$E1
  stopifnot(S1 + dS <=  dim(vals)[2])
  stopifnot(E1 + dE <=  dim(vals)[3])
  stopifnot(p - dp - covp >= 1)  
  stopifnot(dp >= 0)  
  covs <- matrix(0, nrow=S1-S0+1, ncol=E1-E0+1)
  hrange <- p + ((-covp+1):0)
  orange <- p - dp + ((-covp+1):0)
  for (lat in S0:S1) {
    for (lon in E0:E1) {
      covs[lat-S0, lon-E0] <- cov(vals[orange, lat+dS, lon+dE], vals[hrange, lat, lon])
    }
  }                                                                             
  covs
}

# functions by Dave Tanzer:

compute.covariance <- function(vals, t, cov.days, lon1, lat1, lon2, lat2) {

  range <- t + (0:(cov.days-1)) 

  cov(vals[range, lat1, lon1], vals[range, lat2, lon2]) 
}

nlat <- dim(SAvals.3D)[2]
nlon <- dim(SAvals.3D)[3]

band.size <- 1

max.radius = min(nlat,nlon) - 1

numBands <- (max.radius + 1) / band.size 

print("numBands")
print(numBands)

observations.needed <- 10000 

check.done <- function(all.bands) {
  for (i in 1:numBands) {
    if (length(all.bands[[i]]) < observations.needed) {
      return(FALSE)
    }
  }
  TRUE 
}

compute.bands <- function(start.day, cov.days) {

  all.bands = list()
  for (i in 1:numBands) {
    all.bands[[i]] <- vector() 
  } 
  
  repeat {
     x1 <- sample(1:nlon, 1) 
     y1 <- sample(1:nlat, 1) 
  
     x2 <- sample(1:nlon, 1) 
     y2 <- sample(1:nlat, 1)
  
     line.length <- sqrt((x2 - x1)^2 + (y2 - y1)^2) 
  
     if (line.length > max.radius) {
        next
     }
  
     band.index <- 1 + floor(line.length / band.size) 
  
     point.count <- length(all.bands[[band.index]])
  
     if (point.count >= observations.needed) {
       next
     }
  
     all.bands[[band.index]][point.count + 1] <- compute.covariance(SAvals.3D, start.day, cov.days, x1, y1, x2, y2) 
  
     if (point.count + 1 == observations.needed) {
        if (check.done(all.bands)) {
           break
        }
     }
  }
  
  all.bands
}

print.year.data <- function(year) {

  all.bands <- compute.bands((year - 1950) * 365, 365)
  
  cat(sprintf("%d, ", year))
  
  for (i in 1:numBands) {
    median.val <- median(all.bands[[i]])
  
    if (i < numBands) {
      cat(sprintf("%f, ", median.val))
    } else {
      cat(sprintf("%f\n", median.val))
    }
  }
  
}

for (year in 1950:1979) {
  print.year.data(year)
}
