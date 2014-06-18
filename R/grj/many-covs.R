# For efficient Ludescher calculation

# input: two vectors x and y of length s, two numbers n and m with
# n + m < s, and a subset w of (n+m):s. 

# output: a matrix M with length(w) rows and (2*m+1) columns containing 
# covariances between  certain subsets of x and y of length N.
# M[i, m+1-j] = cov(x[(w[i]-j-n+1):(w[i]-j)], y[(w[i]-n+1):(w[i])])
# M[i, m+1+j] = cov(x[(w[i]-n+1):(w[i])], y[(w[i]-j-n+1):(w[i]-j)])
# where j is in 0:m.

# Eg s=10, n=3, m=2. w = c(6,10)

# M[1, 1] = cov(x[(6-2-3+1):(6-2)], y[(6-3+1):6]) = cov(x[2,3,4], y[4,5,6])
# M[1, 2] = cov(x[(6-1-3+1):(6-1)], y[(6-3+1):6]) = cov(x[3,4,5], y[4,5,6])
# M[1, 3] = cov(x[(6-0-3+1):(6-0)], y[(6-3+1):6]) = cov(x[4,5,6], y[4,5,6])
# M[1, 3] = cov(x[(6-3+1):6], y[(6-0-3+1):(6-0)]) = cov(x[4,5,6], y[4,5,6])
# M[1, 4] = cov(x[(6-3+1):6], y[(6-1-3+1):(6-1)]) = cov(x[4,5,6], y[3,4,5])
# M[1, 5] = cov(x[(6-3+1):6], y[(6-2-3+1):(6-2)]) = cov(x[4,5,6], y[2,3,4])

# M[2, 1] = cov(x[(10-2-3+1):(10-2)], y[(10-3+1):10]) = cov(x[6,7,8], y[8,9,10])
# M[2, 2] = cov(x[(10-1-3+1):(10-1)], y[(10-3+1):10]) = cov(x[7,8,9], y[8,9,10])
# M[2, 3] = cov(x[(10-0-3+1):(10-0)], y[(10-3+1):10]) = cov(x[8,9,10], y[8,9,10])
# M[2, 3] = cov(x[(10-3+1):10], y[(10-0-3+1):(10-0)]) = cov(x[8,9,10], y[8,9,10])
# M[2, 4] = cov(x[(10-3+1):10], y[(10-1-3+1):(10-1)]) = cov(x[8,9,10], y[7,8,9])
# M[1, 5] = cov(x[(10-3+1):10], y[(10-2-3+1):(10-2)]) = cov(x[8,9,10], y[6,7,8])

# The repetition of M[1, 3] and M[2, 3] is deliberate: I think it makes things clearer.

# Algorithm

# Construct xrm, yrm (rm=running mean) of length s.
# xrm[j] = mean(x[(j-n+1):j])
# yrm[j] = mean(y[(j-n+1):j])
# (valid for n <= j <= s)


# Construct matrix xy with s rows and m+1 columns
# xy[i, j] = x[i-j] * y[i]
# (valid for m+1 <= i <=  s and  0 <= j <= m )

# Construct matrix xyrm with s rows and m+1 columns
# xyrm[i, j+1] = mean( x[(i-j-n+1):(i-j)] * y[(i-n+1):(i)] )
#              = mean( xy[(i-n+1):i, j] )
# (valid for n+m  <= i <=  s, 0 <= j <= m)

# xyrm[5, 0+1] = mean(x[(5-0-3+1):(5-0)] * y[(5-3+1):(5)]) = mean(x[3,4,5] * y[3,4,5])
# xyrm[5, 1+1] = mean(x[(5-1-3+1):(5-1)] * y[(5-3+1):(5)]) = mean(x[2,3,4] * y[3,4,5])
# xyrm[5, 2+1] = mean(x[(5-2-3+1):(5-2)] * y[(5-3+1):(5)]) = mean(x[1,2,3] * y[3,4,5])
# ...
# xyrm[10, 0+1] = mean(x[(10-0-3+1):(10-0)] * y[(10-3+1):(10)]) = mean(x[8,9,10] * y[8,9,10])
# xyrm[10, 1+1] = mean(x[(10-1-3+1):(10-1)] * y[(10-3+1):(10)]) = mean(x[7,8,9]  * y[8,9,10])
# xyrm[10, 2+1] = mean(x[(10-2-3+1):(10-2)] * y[(10-3+1):(10)]) = mean(x[6,7,8]  * y[8,9,10])

# Also a similar one yxrm with x and y swapped.

# From xrm, yrm, xyrm, and yxrm, the matrix M can be found. Eg

# M[1, 1] =  cov(x[2,3,4], y[4,5,6]) 
#         = mean(x[2,3,4]*y[4,5,6]) - mean(x[2,3,4]) * mean(y[4,5,6])
#         = xyrm[6, 2+1] - xrm[4]*yrm[6]


# M[i, m+1-j] = cov(x[(w[i]-j-n+1):(w[i]-j)], y[(w[i]-n+1:(w[i])])
#             = mean( x[(w[i]-j-n+1):(w[i]-j)] * y[(w[i]-n+1:(w[i])] )  -  mean( x[(w[i]-j-n+1):(w[i]-j)] ) * mean( y[(w[i]-n+1:(w[i])] )
#             = xyrm[w[i], j+1]                                         -  xrm[w[i]-j] * yrm[w[i]]
            




simple.covs <- function(x, y, n, m, w) {
  #M[i, m+1-j] = cov(x[(w[i]-j-n+1):(w[i]-j)], y[(w[i]-n+1:(w[i])])
  # M[i, m+1+j] = cov(x[(w[i]-n+1):(w[i])], y[(w[i]-j-n+1):(w[i]-j)]) 
  M <- matrix(0, nrow=length(w), ncol=2*m+1)
  for (i in 1:length(w)) {
    for (j in 0:m) {
      M[i, m+1-j] = cov( x[(w[i]-j-n+1):(w[i]-j)],   y[(w[i]-n+1):(w[i])] )
      M[i, m+1+j] = cov( x[(w[i]-n+1):(w[i])],   y[(w[i]-j-n+1):(w[i]-j)] ) 
    }
  }
  M
}


fast.covs <- function(x, y, n, m, w) {
  s <- length(x)
  stopifnot(s == length(y))
  stopifnot(n + m <= s)
  stopifnot(min(w) >= n + m)  
  stopifnot(max(w) <= s) 
  
  # this is to improve numerical accuracy
  x <- x - mean(x)
  y <- y - mean(y)
  
  xcs <- c(rep(0,n), cumsum(x))
  xrm <- (xcs[-(1:n)] - xcs[-((s+1):(s+n))]) / n
  ycs <- c(rep(0,n), cumsum(y))
  yrm <- (ycs[-(1:n)] - ycs[-((s+1):(s+n))]) / n
  
  xy <- matrix(0, nrow=s, ncol=m+1)
  for (i in (m+1):s) {
    xy[i, ] = x[i-(0:m)] * y[i]    
  }
  
  xyrm <- matrix(0, nrow=s, ncol=m+1)
  for (i in 1:length(w)) {
    xyrm[w[i], ] <- colMeans( xy[(w[i]-n+1):w[i], ] ) 
  }
  
  yx <- matrix(0, nrow=s, ncol=m+1)
  for (i in (m+1):s) {
    yx[i, ] = y[i-(0:m)] * x[i]    
  }
  
  yxrm <- matrix(0, nrow=s, ncol=m+1)
  for (i in 1:length(w)) {
    yxrm[w[i], ] <- colMeans( yx[(w[i]-n+1):w[i], ] ) 
  }
  
  M <- matrix(0, nrow=length(w), ncol=2*m+1)
  for (i in 1:length(w)) {
    M[i, 1:(m+1)] <- rev(xyrm[w[i], 1:(m+1)]   -  xrm[w[i]-(0:m)] * yrm[w[i]])   
  }
  for (i in 1:length(w)) {
    M[i, (m+1):(2*m+1)] <- yxrm[w[i], 1:(m+1)]   -  yrm[w[i]-(0:m)] * xrm[w[i]]    
  }
  
  # to agree with R version of cov
  M*n/(n-1)
}




s <- 365*30

x <- 1:s  + rnorm(s)
y <- 1:s + rnorm(s)
n <- 365
m <- 200
w <- c(seq(from=1000,to=s,by=10))



print( system.time( Mf <- simple.covs(x, y, n, m, w) ) )
print( system.time( Ms <- fast.covs(x, y, n, m, w) ) )
max(abs(Ms-Mf))
  
  



