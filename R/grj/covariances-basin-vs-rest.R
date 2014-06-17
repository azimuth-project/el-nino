


setwd("C:/Users/Work/AAA/Programming/ProgramOutput/Nino")





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




# find covariances at time (day) p between point (latx,lonx) and
# (laty,lony), with time lag dp, and
# using a period length covp for the covariances. Positive
# dp means (latx,lonx) is at time p, and (laty,lony) is
# earlier by |dp|. Negative dp means
# the (laty,lony) point is at time p, and (latx, lonx)
# earlier by |dp|.
covariance.of.2pts <- function(vals, latx, lonx, laty, lony, p, dp, covp) {
  stopifnot(latx >= 1  &&  latx <=  dim(vals)[2])
  stopifnot(lonx >= 1  &&  lonx <=  dim(vals)[3])
  stopifnot(laty >= 1  &&  laty <=  dim(vals)[2])
  stopifnot(lony >= 1  &&  lony <=  dim(vals)[3])  
  stopifnot(p - dp - covp >= 1)  
  if (dp >= 0) { 
    px <- p + ((-covp+1):0)
    py <- p - dp + ((-covp+1):0)
  } else {
    px <- p + dp + ((-covp+1):0)
    py <- p + ((-covp+1):0)    
  }
  cov(vals[px, latx, lonx], vals[py, laty, lony])
}


# returns TRUE if (lat,lon) is in Ludescher et al's "El Nino basin"
in.basin <- function(lat, lon) {
  if (lat==5  && lon >= 11  && lon <= 22) {
    return (TRUE)
  }
  if (lat == 6 && lon == 16) {
    return (TRUE)
  }
  if (lat == 6 && lon == 22) {
    return (TRUE)
  }
  return (FALSE)
}

# returns lon value for i'th point in Ludescher et al's "El Nino basin"
basin.point.lon <- function(i) {
  stopifnot(i >= 1  && i <= 14)
  if (i <= 12) {
    return (i+10)
  } else {
    return (16 + 6*(i-13))
  }
}

# returns lon value for i'th point in Ludescher et al's "El Nino basin"
basin.point.lat <- function(i) {
  stopifnot(i >= 1  && i <= 14)
  if (i <= 12) {
    return (5)
  } else {
    return (6)
  }
}


# for each point (lat,lon) in grid, calculates an average covariance between
#  (lat,lon) and the points in Ludescher et al's "El Nino basin" 
cov.basin.vs.rest <- function(vals, d, covp) {
  bcovs <- matrix(0, nrow=dim(vals)[2], ncol=dim(vals)[3]) 
  for (lat in 1:dim(vals)[2]) {
    for (lon in 1:dim(vals)[3]) {
      if (!in.basin(lat, lon)) {
        bcov <- 0
        for (b in 1:14) {
          latb <- basin.point.lat(b)
          lonb <- basin.point.lon(b)
          cov0 <- covariance.of.2pts(vals, lat, lon, latb, lonb, d, 0, covp)
          #cov.p1 <- covariance.of.2pts(vals, lat, lon, latb, lonb, d,  1, covp)
          #cov.n1 <- covariance.of.2pts(vals, lat, lon, latb, lonb, d, -1, covp)          
          #bcov <- bcov + mean(cov0, cov.n1, cov.p1)
          bcov <- bcov + cov0
        }
        bcovs[lat,lon] <- bcov
      }     
    }    
  }
  bcovs
}

#########################################################################

plot.nino.3.4.background.rectangle <- function(mint, maxt, col) {
  polygon(x=c(1,13,13,1,1), y=c(maxt,maxt,mint,mint,maxt), border = NA, col=col)
}


plot.nino.3.4 <- function(year) {
  nini <- read.table("nino34-anoms.txt", skip=1, header=TRUE)
  nini <- as.matrix(nini)
  w <- which(nini[,"YR"] == year)
  stopifnot(length(w)==12)
  if (year+1 <= nini[nrow(nini), "YR"]) {
    w <- c(w,w[12]+1)
  } else {
    w <- c(w,w[12])
  }
  yrnini <- nini[w,"ANOM"]

  time.axis <- 1:13
  # empty plot, so text and 'guide' lines can go underneath
  plot(time.axis, yrnini, type='n', ylim=c(-2.5,2.5), yaxt='n', xaxt='n')

  plot.nino.3.4.background.rectangle(1.5, 2.5, col="#ff9999ff")
  plot.nino.3.4.background.rectangle(1, 1.5, col="#ffbbbbff")
  plot.nino.3.4.background.rectangle(.5, 1, col="#ffddddff")
  
  plot.nino.3.4.background.rectangle(-1, -.5, col="#ddddffff")
  plot.nino.3.4.background.rectangle(-1.5, -1, col="#bbbbffff")
  plot.nino.3.4.background.rectangle(-2.5, -1.5, col="#9999ffff")
  
  text(time.axis[7], 0, paste0(year), col="grey80", cex=2)
  
  
  for (y in c(1,4,7,10,13)) {
    lines(c(y, y), c(-2,2), col="grey70")
  }
  
  # plot the index
  lines(time.axis, yrnini, type='l', ylim=c(-1.5,2.5), col="black", lwd=1, yaxt='n', xaxt='n')
  
}



plot.bcovs.image <- function(bcovs) {
  # image() has inconvenient conventions for x- and y- directions 
  # for current purpose. Transpose, then reverse columns.
  x <-  t(bcovs)
  for (r in 1:nrow(x)) {
    x[r,] <- rev(x[r,])
  }
  # squash big positive and negative values
  x <- sign(x) * sqrt(abs(x))
  # make colours. Black is zero, red negative, green positive
  clrs <- rep(rgb(1,1,1,1),21)
  clrs[1] <- rgb(1,0,0,1)
  for (i in 0:9) {
    clrs[i+2] <- rgb(1-.08*i,.9-.1*i,.9-.1*i,1)
  }
  clrs[12] <- rgb(0,0,0,1)
  for (i in 0:9) {
    clrs[22-i] <- rgb(.9-.1*i,1-.08*i,.9-.1*i,1)
  }  
  clrs[23] <- rgb(0,1,0,1)
  brks <- c(-99, seq(from=-3,to=3,length.out=22), 99)
  image(1:nrow(x), 1:ncol(x), x, breaks=brks, col=clrs, yaxt='n', xaxt='n')
  
}

############################################################################

Kvals <- as.matrix(read.table(file="Pacific-1950-1979.txt", header=TRUE))
Kvals.cnames <- colnames(Kvals)

Kvals.3D <- data.to.3D(Kvals)

SAvals.3D <- seasonally.adjust(Kvals.3D)

SAvals.3D.3x3 <- subsample.3x3(SAvals.3D)

png(filename="covs-b-vs-r-maps.png", width=800, height=2000)
oldpar <- par(mfrow=c(29,5), mar=rep(.5,4))
firstyear <- year.range.from.rownames(rownames(Kvals))[1]
for (y in 1951:1979) {
  for (p in c(45,136,227,318)) {
    d <- (y - firstyear)*365 + p
    bcovs <- cov.basin.vs.rest(SAvals.3D.3x3, d, 365)
    plot.bcovs.image(bcovs)    
  }
  plot.nino.3.4(y)
}
par(oldpar)
dev.off()

############################################################################





