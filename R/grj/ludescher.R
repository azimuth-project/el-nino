
# replicating Ludescher et al


setwd("C:/Users/Work/AAA/Programming/ProgramOutput/Nino")

options(warn=2)

###################################################################
############ For reading in data

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

year.from.label <- function(label) {
  as.integer(substr(label,2,5))
}


year.range.from.rownames <- function(rnames) {
  year.from.label(rnames[1]) : year.from.label(rnames[length(rnames)])
}



#####################################################################
########## Conversion, Seasonal adjustment, subsampling

# Converts from a vector per day to a 2D array per day
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

# subtracts the climatological seasonal cycle (mean over years for each grid point, each day-in-year)
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

# the data per day is reduced from eg 27x69 to 9x23. 
subsample.3x3 <- function(vals) {
  stopifnot(dim(vals)[2] %% 3 == 0)
  stopifnot(dim(vals)[3] %% 3 == 0)
  n.sslats <- dim(vals)[2]/3
  n.sslons <- dim(vals)[3]/3
  ssvals <- array(0, dim=c(dim(vals)[1], n.sslats, n.sslons))  
  for (d in 1:dim(vals)[1]) {
    for (slat in 1:n.sslats) {
      for (slon in 1:n.sslons) {
        ssvals[d, slat, slon] <- mean(vals[d, (3*slat-2):(3*slat), (3*slon-2):(3*slon)])
      }
    }
  }
  ssvals
}



#####################################################################
############ Making the time-delayed cross correlations - arbitrary points and days


# input: 3D array vals of temperatures vals for s days, a running mean length n (=365)

# output: 3D array of the same dimensions as vals containing running means
make.runningmeans <- function(vals) {
  rmeans <- array(0, dim=dim(vals))
  for (lat in 1:dim(vals)[2]) {
    for (lon in 1:dim(vals)[3]) { 
      rmeans[, lat, lon] <- as.vector(filter(vals[, lat, lon], rep(1/n,n), sides=1))
    }
  }
  rmeans
}




convolve.filter.nextpow2 <- function(delayed, current) {
  nm <- length(delayed)
  n <- length(current)
  m <- nm - n
  v <- nm+n-1
  N <- nextn(v, 2)
  d <- N-v
  convolve(c(delayed,rep(0,d)), current, type="filter")[1:(m+1)]
}




# input: two vectors x and y of length s, corresponding running means
# of period n, two numbers n and m with
# n + m < s, and a day d in (n+m):s. 

# output: an array convs of length (2*m+1) columns containing 
# covariances between certain subsets of x and y of length n.
# convs[m+1-j] = cov( x[(d-j-n+1):(d-j)],  y[(d-n+1):d] )
# convs[m+1+j] = cov( x[(d-n+1):d],  y[(d-j-n+1):(d-j)] )
# where j is in 0:m.
convolve.covs <- function(x, y, xrm, yrm, n, m, d) {
  s <- length(x)
  stopifnot(s == length(y))
  stopifnot(n + m <= s)
  
  convs <- rep(0, 2*m+1)
  
  yn <- y[(d-n+1):d] - yrm[d]
  xnm <-  x[(d-n-m+1):d]    
  convs[1:(m+1)] <- convolve.filter.nextpow2(xnm, yn)   
  
  xn <- x[(d-n+1):d] - xrm[d]
  ynm <-  y[(d-n-m+1):d]
  convs[(2*m+1):(m+1)] <- convolve.filter.nextpow2(ynm, xn)   

  convs <- convs/n
}





# input: 3D array of temperatures vals for s days, a convolution-length n (eg 365)
# a maximum time-delay m (eg 200), and a day d in (n+m):s. 
# output: 3D array of time-delayed standard deviations for day d. 
# The first dimension is the time delay -m:m, second is lat, third lon
make.standarddevs <- function(vals, n, m, d) {
  sds <- array(0, dim=c(2*m+1, dim(vals)[2], dim(vals)[3]))
  for (lat in 1:dim(vals)[2]) {
    for (lon in 1:dim(vals)[3]) {
      recentx <- vals[(d-n-m+1):d , lat, lon]
      xrm <- as.vector(filter(recentx, rep(1/n,n), sides=1))
      x2rm <- as.vector(filter(recentx^2, rep(1/n,n), sides=1))      
      fsds <- sqrt(x2rm - xrm^2)*sqrt(n/(n-1))
      fsds <- fsds[n:(n+m)]
      sds[1:(m+1) ,lat, lon] <- fsds
      sds[(m+1):(2*m+1) ,lat, lon] <- rev(fsds)
    }
  }
  sds
}



# Not used for Ludescher replication
# input: 3D array of temperatures vals for s days, corresponding array
# rmeans of running means of length n, a grid point xlat, xlon,
# a convolution-length n (eg 365)
# a maximum time-delay m (eg 200), and a day d in (n+m):s. 
# output: 3D array of time-delayed covariances for point (xlat, xlon) for day d. 
# The first dimension is the time delay -m:m, second is lat, third lon
make.covs <- function(vals, rmeans, xlat, xlon, n, m, d) {
  covs <- array(0, dim=c(2*m+1, dim(vals)[2], dim(vals)[3]))
  x <- vals[ , xlat, xlon]
  xrm <- rmeans[, xlat, xlon]           
  for (ylat in 1:dim(vals)[2]) {#
   for (ylon in 1:dim(vals)[3]) {
     y <- vals[ , ylat, ylon]
     yrm <- rmeans[ , ylat, ylon]
     covs[ , ylat, ylon] <- convolve.covs(x, y, xrm, yrm, n, m, d)
   }
  }
  covs
}




# input: 3D array of temperatures vals for s days, corresponding array
# rmeans of running means of length n, an array sds from make.standarddevs()
# for the same day d, a grid point (xlat, xlon),
# a convolution-length n (eg 365)
# a maximum time-delay m (eg 200), and a day d in (n+m):s. 
# output: 3D array of time-delayed correlations for point (xlat, xlon) for day d. 
# The first dimension is the time delay -m:m, second is lat, third lon
make.cors <- function(vals, rmeans, sds, xlat, xlon, n, m, d) {
  cors <- array(0, dim=c(2*m+1, dim(vals)[2], dim(vals)[3]))
  x <- vals[ , xlat, xlon]
  xrm <- rmeans[, xlat, xlon]           
  for (ylat in 1:dim(vals)[2]) {#
    for (ylon in 1:dim(vals)[3]) {
      y <- vals[ , ylat, ylon]
      yrm <- rmeans[ , ylat, ylon]
      covs <- convolve.covs(x, y, xrm, yrm, n, m, d)
      cors[ , ylat, ylon] <- covs / (sds[ , ylat, ylon] * sds[ , xlat, xlon])      
    }
  }
  cors
}




######################################################################
############ Making the link strengths - basin vs rest



basin <- function() {
  lats <- c( 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6)
  lons <- c(11,12,13,14,15,16,17,18,19,20,21,22,16,22)
  stopifnot(length(lats) == length(lons))
  list(lats=lats,lons=lons)
}


# returns TRUE if (lat,lon) is in 'focus'
in.focus <- function(fa, lat, lon) {
  w <- which(fa$lats==lat)  
  (length(which(fa$lons[w]==lon)) > 0)
}



# for each point (lat,lon) in grid, calculates an average link strength between
#  (lat,lon) and the points in Ludescher et al's "El Nino basin" 
signalstrength <- function(linktype, vals, rmeans, sds, n, m, d) {
  S <- 0
  ba <- basin()
  nba <- length(ba$lats)
  nlats <- dim(vals)[2]
  nlons <- dim(vals)[3]
  for (b in 1:nba) {
    latb <- ba$lats[b]
    lonb <- ba$lons[b]
    if (linktype == "correlations") {
      links <- make.cors(vals, rmeans, sds, latb, lonb, n, m, d)
    } else if (linktype == "covariances") {
      links <- make.covs(vals, rmeans, latb, lonb, n, m, d)
    } else {
      stop()
    }
    
    for (lat in 1:nlats) {
      for (lon in 1:nlons) {
        if (!in.focus(ba, lat, lon)) {
          tdccs <- links[ , lat, lon]
          sig <- (max(tdccs) - mean(tdccs)) / sd(tdccs)
          S <- S + sig
        }
      }
    } 
  }
  S / nba / (nlats*nlons - nba)
}



#########################################################################

plot.nino.3.4.background.rectangle <- function(mint, maxt, col) {
  polygon(x=c(1,13,13,1,1), y=c(maxt,maxt,mint,mint,maxt), border = NA, col=col)
}


plot.nino.3.4 <- function(firstyear, lastyear, miny, maxy, clr) {
  nini <- read.table("nino34-anoms.txt", skip=1, header=TRUE)
  nini <- as.matrix(nini)
  w <- which((nini[,"YR"] >= firstyear) & (nini[,"YR"] <= lastyear))
  stopifnot((length(w) %% 12) == 0)
  yrnini <- nini[w,"ANOM"]
  offset <- min(yrnini) 
  scaling <- (maxy - miny) / (max(yrnini) - min(yrnini))
  yrnini <- miny + scaling * (yrnini - offset)
  zp5 <- miny + scaling * (0.5 - offset)
  time.axis <- firstyear + (0:(length(w)-1))/12
  # plot the index
  lines(time.axis, yrnini, col=clr)
  # plot the 0.5 line
  lines(c(time.axis[1], time.axis[length(time.axis)]), c(zp5,zp5), col=clr)
  
}



############################################################################

Kvals <- as.matrix(read.table(file="Pacific-1950-1979.txt", header=TRUE))
Kvals.cnames <- colnames(Kvals)

Kvals.3D <- data.to.3D(Kvals)

SAvals.3D <- seasonally.adjust(Kvals.3D)

SAvals.3D.3x3 <- subsample.3x3(SAvals.3D)

n <- 365
m <- 200
w <- seq(from = 2*365, to=dim(SAvals.3D.3x3)[1], by=73)

rmeans <- make.runningmeans(SAvals.3D.3x3)

S <- rep(0, length(w))
for (i in 1:length(w)) {
  d <- w[i]
  sds <- make.standarddevs(SAvals.3D.3x3, n, m, d)  
  S[i] <- signalstrength("correlations", SAvals.3D.3x3, rmeans, sds, n, m, d)
  cat("done day", d, "S(d)=", S[i], "\n")
  
}


firstyear <- 1952
lastyear <- 1979
plot(firstyear+(1:length(S))/5, S, type='n', xlab="years")
for (yr in firstyear:(lastyear+1)) {
  lines(c(yr,yr), c(min(S),max(S)), col="grey80")
}
lines(firstyear+(0:(length(S)-1))/5, S)
plot.nino.3.4(firstyear, lastyear, min(S), max(S), "red")
lines(c(firstyear,(lastyear+1)), rep(2.3,2))

############################################################################






























