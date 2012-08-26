##create a data matrix
t <- seq(0,4*pi,len=50)
X.org <- cbind(cos(t),sin(2*t)) %*% matrix(rnorm(20),2,10)

##add some normal errors
X <- X.org + .25*rnorm(length(X.org))
##and mark some data as missing
X[runif(length(X))<.25] <- NA

##Ensure that we have complet columns/rows
while( any(rowSums(is.na(X))==dim(X)[2]) || any(colSums(is.na(X))==dim(X)[1]) ){
  X <- X.org + .25*rnorm(length(X.org))
  X[runif(length(X))<.25] <- NA
}

##compute two smooth basis functions
res <- SVDsmooth(X, n.basis=2, niter=100)

##plot the two smooth basis functions
par(mfcol=c(3,2), mar=c(4,4,.5,.5))
plot(t, res[,1], ylim=range(res), type="l")
lines(t, res[,2], col=2)
##and some of the data fitted to the smooths
for(i in 1:5){
  plot(t, X[,i])
  lines(t, predict.lm(lm(X[,i]~res), data.frame(res)) )
  lines(t, X.org[,i], col=2)
}

##compute cross-validation for 1 to 4 basis functions
res.cv <- SVDsmoothCV(X, n.basis=1:4, niter=100)

##study cross-validation results
print(res.cv)

##plot cross-validation statistics
plot(res.cv)
##plot the BIC for each column
plot(res.cv, pairs=TRUE)
