
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
make.runningmeans <- function(vals, n) {
  rmeans <- array(0, dim=dim(vals))
  for (lat in 1:dim(vals)[2]) {
    for (lon in 1:dim(vals)[3]) { 
      rmeans[, lat, lon] <- as.vector(filter(vals[, lat, lon], rep(1/n,n), sides=1))
    }
  }
  rmeans
}




convolve.filter.nextpow2and3 <- function(delayed, current) {
  nm <- length(delayed)
  n <- length(current)
  m <- nm - n
  v <- nm+n-1
  N <- nextn(v, c(2,3))
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
  convs[1:(m+1)] <- convolve.filter.nextpow2and3(xnm, yn)   
  
  xn <- x[(d-n+1):d] - xrm[d]
  ynm <-  y[(d-n-m+1):d]
  convs[(2*m+1):(m+1)] <- convolve.filter.nextpow2and3(ynm, xn)   

  convs <- convs/n
}


# FOR TESTING
make.standarddevs.for.one.point <- function(vals, lat, lon, n, m, d) {
  sds <- rep(0, 2*m+1)
  recentx <- vals[(d-n-m+1):d , lat, lon]
  xrm <- as.vector(filter(recentx, rep(1/n,n), sides=1))
  x2rm <- as.vector(filter(recentx^2, rep(1/n,n), sides=1))      
  fsds <- sqrt(x2rm - xrm^2)*sqrt(n/(n-1))
  fsds <- fsds[n:(n+m)]
  sds[1:(m+1)] <- fsds
  sds[(m+1):(2*m+1)] <- rev(fsds)
  sds
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


# FOR TESTING
make.covs.for.one.pair <- function(vals, rmeans, xlat, xlon, ylat, ylon, n, m, d) {
  x <- vals[ , xlat, xlon]
  xrm <- rmeans[, xlat, xlon]           
  y <- vals[ , ylat, ylon]
  yrm <- rmeans[ , ylat, ylon]
  convolve.covs(x, y, xrm, yrm, n, m, d)
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



ludescher.basin <- function() {
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



make.links.to.one.point <- function(linktype, lat, lon, vals, rmeans, sds, n, m, d) {
  if (linktype == "abs_corrs") {
    links <- abs(make.cors(vals, rmeans, sds, lat, lon, n, m, d))
  } else if (linktype == "sgn_corrs") {
    links <- make.cors(vals, rmeans, sds, lat, lon, n, m, d)
  } else if (linktype == "abs_covs") {
    links <- abs(make.covs(vals, rmeans, lat, lon, n, m, d))
  } else if (linktype == "sgn_covs") {
    links <- make.covs(vals, rmeans, lat, lon, n, m, d)
  } else if (linktype == "stddevs") {
    links <- sds    # not really a link; just a convenient implementation
  } else {
    stop("Unknown link type")
  }
  links 
}



# for each point (lat,lon) in grid, calculates an average link strength between
#  (lat,lon) and the points in Ludescher et al's "El Nino basin" 
signalstrength <- function(linktype, sigtype, focus, vals, rmeans, sds, n, m, d) {
  S <- 0
  nba <- length(focus$lats)
  nlats <- dim(vals)[2]
  nlons <- dim(vals)[3]
  for (b in 1:nba) {
    links <- make.links.to.one.point(linktype, focus$lats[b], focus$lons[b], 
                                     vals, rmeans, sds, n, m, d)
    for (lat in 1:nlats) {
      for (lon in 1:nlons) {
        if (!in.focus(focus, lat, lon)) {
          tdccs <- links[ , lat, lon]
          if (sigtype == "ludescher") {
            sig <- (max(tdccs) - mean(tdccs)) / sd(tdccs)
          } else if (sigtype == "smalldelay") {
            mid <- (length(tdccs) + 1) / 2
            rnge <- min(10, mid-1)
            mx <- max(tdccs[(mid-rnge):(mid+rnge)])
            #noiz <- mean(abs(tdccs[-1] - tdccs[-length(tdccs)]))
            #sig <- mx / noiz
            sig <- mx
          }
          
          S <- S + sig
        }
      }
    } 
  }
  S / nba / (nlats*nlons - nba)
}


make.signal.strength <- function(linktype, sigtype, vals, rmeans, params) {
  firstyear <- params$firstyear
  lastyear <- params$lastyear
  n <- params$n
  m <- params$m
  step <- params$step
  w <- seq(from = (firstyear-1950)*365, to=dim(vals)[1], by=step)
  S <- rep(0, length(w))
  for (i in 1:length(w)) {
    d <- w[i]
    sds <- make.standarddevs(vals, n, m, d)  
    S[i] <- signalstrength(linktype, sigtype, ludescher.basin(), vals, rmeans, sds, n, m, d)
    cat("done day", d, "S(d)=", S[i], "\n")  
  }
  S
}



##############################################################################
############## Plotting results


############# For plotting Nino index

plot.nino.3.4.background.rectangle <- function(mint, maxt, col) {
  polygon(x=c(1,13,13,1,1), y=c(maxt,maxt,mint,mint,maxt), border = NA, col=col)
}


find.nino.plotting.info <- function(firstyear, lastyear, miny, maxy) {
  nini <- read.table("nino34-anoms.txt", skip=1, header=TRUE)
  nini <- as.matrix(nini)
  w <- which((nini[,"YR"] >= firstyear) & (nini[,"YR"] <= lastyear))
  stopifnot((length(w) %% 12) == 0)
  yrnini <- nini[w,"ANOM"]
  offset <- min(yrnini) 
  scaling <- (maxy - miny) / (max(yrnini) - min(yrnini))
  yrnini <- miny + scaling * (yrnini - offset)
  zp5 <- miny + scaling * (0.5 - offset)
  labels <- c("-2", "-1", "0", "+1", "+2")
  ticks <- miny + scaling * (c(-2,-1,0,1,2) - offset)
  time.axis <- firstyear + (0:(length(w)-1))/12
  list(time.axis=time.axis, yrnini=yrnini, zp5=zp5, ticks=ticks, labels=labels,
       firstyear=firstyear, lastyear=lastyear, miny=miny, maxy=maxy)
}

plot.nino.zp5.rect <- function(plotinfo, col) {
  time.axis <- plotinfo$time.axis
  minx <- time.axis[1]
  maxx <- time.axis[length(time.axis)]
  miny <- plotinfo$miny
  zp5 <- plotinfo$zp5
  polygon(x=c(minx,maxx,maxx,minx,minx), y=c(zp5,zp5,miny,miny,zp5), border = NA, col=col)
  
}



plot.nino.3.4 <- function(plotinfo, col) {
  lines(plotinfo$time.axis, plotinfo$yrnini, col=col)
}



############################################################################
# Pages of 5 by 12 cross-correlation graphs.

make.graphs.of.links <- function(linktype, description, vals, rmeans, lat, lon, d, params) {
  firstyear <- params$firstyear
  lastyear <- params$lastyear
  n <- params$n
  m <- params$m
  step <- params$step
  sds <- make.standarddevs(vals, n, m, d)  
  links <- make.links.to.one.point(linktype, lat, lon, vals, rmeans, sds, n, m, d)
  year <- 1950+floor((d-1)/365)
  dayinyear <- 1+(d-1)%%365
  png(paste0(linktype, "_", year, "_", dayinyear, "_", lat, "_", lon, "_", d, ".png"), width=1200, height=600, pointsize=18)  
  nfmat <- matrix(c((1:60),rep(61,12)), nrow=6, ncol=12, byrow=TRUE)
  nf <- layout(nfmat)
  layout.show(nf)
  oldpar <- par(mar=rep(.1,4))
  for (laty in seq(from=1, to=9, by=2)) {
    for (lony in seq(from=1, to=23, by=2)) {
      plot((-m:m), abs(links[,laty,lony]), type='l', ylim=c(min(links),max(links)), 
           col="grey", xaxt='n', yaxt='n', ann=FALSE)  
      lines((-m:m), links[,laty,lony], ylim=c(min(links),max(links)))
      lines(c(-m,m), c(0,0)) 
      cex=1
      col = "grey"
      if (laty==lat && lony==lon) {
        cex <- 1.5
      }
      if (laty==5 && lony >= 11  && lony <= 21) {
        col="red"
      }
      text(x=-m, y=max(links), paste0(laty, ",", lony), adj=c(0,1), col=col, cex=cex)
    }  
  }
  plot(c(-10,10), c(-2,2), type='n', xaxt='n', yaxt='n', ann=FALSE)
  if (linktype == "abs_corrs") {
    text(0,1, paste0("Delayed cross-abs(correlations) between grid point (", lat, ",", lon, 
                     ") and others in the Pacific",
                     " for year ", year, ", day ", dayinyear, ". "))    
  } else if (linktype == "sgn_corrs") {
    text(0,1, paste0("Delayed cross-correlations between grid point (", lat, ",", lon, 
                     ") and others in the Pacific",
                     " for year ", year, ", day ", dayinyear, ". "))    
    text(0,0, "The black wiggles show the value for tau=-200 to 200, the grey wiggles the absolute value (where it differs)")    
  } else if (linktype == "abs_covs") {
    text(0,1, paste0("Delayed cross-abs(covariances) between grid point (", lat, ",", lon, 
                     ") and others in the Pacific",
                     " for year ", year, ", day ", dayinyear, ". "))        
  } else if (linktype == "sgn_covs") {
    text(0,1, paste0("Delayed cross-covariances between grid point (", lat, ",", lon, 
                     ") and others in the Pacific",
                     " for year ", year, ", day ", dayinyear, ". "))    
    text(0,0, "The black wiggles show the value for tau=-200 to 200, the grey wiggles the absolute value (where it differs)")      
  } else if (linktype == "stddevs") {
    text(0,1, paste0("Delayed stddevs at each grid point in the Pacific",
                     " for year ", year, ", day ", dayinyear, ". "))      
    text(0,0, "The black wiggles show the value for tau=-200 to 200")      
  } 
  text(0,-1,paste0("The y-axis range is ", round(min(links), digits=2), " to ", round(max(links), digits=2), "."))
  par(oldpar)
  dev.off()
}



make.setof.graphs.of.links <- function(linktype, description, vals, rmeans, params) {
  make.graphs.of.links(linktype, description, vals, rmeans, 3, 5, 7*365, params)
  make.graphs.of.links(linktype, description, vals, rmeans, 3, 5, 8*365, params)
  make.graphs.of.links(linktype, description, vals, rmeans, 3, 5, 9*365, params)
  make.graphs.of.links(linktype, description, vals, rmeans, 5, 13, 7*365, params)
  make.graphs.of.links(linktype, description, vals, rmeans, 5, 13, 8*365, params)
  make.graphs.of.links(linktype, description, vals, rmeans, 5, 13, 9*365, params)
}


# main plot 
plot.signalstrength.vs.nino <- function(S, params) {
  firstyear <- params$firstyear
  lastyear <- params$lastyear
  n <- params$n
  m <- params$m
  step <- params$step
  time.axis <- firstyear+(0:(length(S)-1)) * step / 365
  par(mar=c(5, 4, 4, 5))
  plot(time.axis, S, type='n', xlab="Years", ylab="Signal strength S", 
       main="S and theta in red. NINO index in blue, below 0.5C shaded")
  ninoplotinfo <- find.nino.plotting.info(firstyear, lastyear, min(S), max(S))
  plot.nino.zp5.rect(ninoplotinfo, "#eeeeffff")
  for (yr in firstyear:(lastyear+1)) {
    lines(c(yr,yr), c(min(S),max(S)), col="grey80")
  }
  lines(time.axis, S, col="red")
  plot.nino.3.4(ninoplotinfo, "blue")
  lines(c(firstyear,(lastyear+1)), rep(20,2), col="red")
  axis(side=4, at=ninoplotinfo$ticks, labels=ninoplotinfo$labels)
  mtext(text="NINO 3.4 index", side = 4, line = 3)
}



############################################################################
#################### Main analysis

Kvals <- as.matrix(read.table(file="Pacific-1950-1979.txt", header=TRUE))
analysis.params <- list(firstyear=1952, lastyear=1979, n=365, m=200, step=100)

Kvals.3D <- data.to.3D(Kvals)
rm(Kvals)
SAvals.3D <- seasonally.adjust(Kvals.3D)
rm(Kvals.3D)
SAvals.3D.3x3 <- subsample.3x3(SAvals.3D)
SAvals.3D.3x3.rm <- make.runningmeans(SAvals.3D.3x3, analysis.params$n)


make.setof.graphs.of.links("stddevs", "stddevs", SAvals.3D.3x3, SAvals.3D.3x3.rm, analysis.params)
make.setof.graphs.of.links("sgn_covs", "covariances", SAvals.3D.3x3, SAvals.3D.3x3.rm, analysis.params)
make.setof.graphs.of.links("sgn_corrs", "correlations", SAvals.3D.3x3, SAvals.3D.3x3.rm, analysis.params)


#S <- make.signal.strength("abs_covs", "smalldelay", SAvals.3D.3x3, SAvals.3D.3x3.rm, analysis.params)
#plot.signalstrength.vs.nino(S, analysis.params)
































