### R code from vignette source 'examplecc.Rnw'

###################################################
### code chunk number 1: examplecc.Rnw:14-15
###################################################
setwd("C:/Users/Work/AAA/Programming/misc/Azi-Nino")


###################################################
### code chunk number 2: examplecc.Rnw:18-33
###################################################
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


###################################################
### code chunk number 3: examplecc.Rnw:36-49
###################################################
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


###################################################
### code chunk number 4: examplecc.Rnw:52-65
###################################################
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


###################################################
### code chunk number 5: examplecc.Rnw:68-82
###################################################
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


###################################################
### code chunk number 6: examplecc.Rnw:85-105
###################################################
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


###################################################
### code chunk number 7: examplecc.Rnw:108-130
###################################################
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


###################################################
### code chunk number 8: examplecc.Rnw:133-140
###################################################
#Kvals <- as.matrix(read.table(file="Scotland-1950-1952.txt", header=TRUE))
Kvals <- as.matrix(read.table(file="Pacific-1950-1979.txt", header=TRUE))
Kvals.cnames <- colnames(Kvals)

Kvals.3D <- data.to.3D(Kvals)

SAvals.3D <- seasonally.adjust(Kvals.3D)


###################################################
### code chunk number 9: examplecc.Rnw:143-169
###################################################
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


###################################################
### code chunk number 10: examplecc.Rnw:172-187
###################################################
region <- list(S0=12,S1=16,E0=20,E1=60)
pvals <- seq(from=366, to=30*365, by=5)

mediancovs1 <- matrix(0, nrow=8, ncol=length(pvals))
mediancovs5 <- matrix(0, nrow=8, ncol=length(pvals))
for (i in 1:length(pvals)) {
  p <- pvals[i]
  for (dE in (0:7)) {
    covs <- covariances(SAvals.3D, region, p, dS=0, dE=dE, dp=1, covp=183)
    mediancovs1[dE+1, i] <- median(c(covs))
    covs <- covariances(SAvals.3D, region, p, dS=0, dE=dE, dp=5, covp=183)
    mediancovs5[dE+1, i] <- median(c(covs))
  }
}
timeinyears <- 1951 + (0:(length(pvals)-1))/73


###################################################
### code chunk number 11: examplecc.Rnw:190-197
###################################################
plot(timeinyears, mediancovs1[1,], ylim=c(min(mediancovs1), max(mediancovs1)), type='l')
for (dE in 1:7) {
  lines(timeinyears, mediancovs1[dE,], type='l', col=rgb(dE/8,dE/8,dE/8,1))
}
for (y in 1951:1980)  {
  lines(c(y,y), c(0,max(mediancovs1)))
}


###################################################
### code chunk number 12: examplecc.Rnw:202-209
###################################################
plot(timeinyears, mediancovs5[1,], ylim=c(min(mediancovs5), max(mediancovs5)), type='l')
for (dE in 1:7) {
  lines(timeinyears, mediancovs5[dE,], type='l', col=rgb(dE/8,dE/8,dE/8,1))
}
for (y in 1951:1980)  {
  lines(c(y,y), c(0,max(mediancovs5)))
}

