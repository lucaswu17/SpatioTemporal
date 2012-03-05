############################################
# FUNCTIONS THAT COMPUTE MISSING DATA SVDs #
############################################
#missing data SVD
SVD.miss <- function(X, niter = 25, ncomp = 4, conv.reldiff = 0.001)
{
  if( all(!is.na(X)) ){ #if X has no missing data this is simple
    svd0 <- svd(X)
    XF <- X
    i <- diff <- reldiff <- 0
  }else{ #X has missing data, use iterative method
    # Iterative svd calculation with missing data.
    # Initial first element of U matrix is average curve.
    u1 <- apply(X, 1, mean, na.rm = TRUE)
    XM <- matrix(1, nrow(X), ncol(X))
    XM[is.na(X)] <- 0
    XZ <- X
    # v1 is proportional to X'u1/(u1'u1), but with calculations
    # filling in zeros for missing values in the sums.
    XZ[is.na(X)] <- 0.
    # Fill in missing values for initial complete SVD calculation.
    # Then iterate  using only complete data.
    v1 <- diag(t(XZ) %*% (XM * u1))/diag(t(XM * u1) %*% (XM * u1))
    XF <- X
    XF[is.na(X)] <- (matrix(u1, ncol = 1) %*% matrix(v1, nrow = 1))[is.na(X)]
    if( any(is.na(XF)) )
      stop("In 'SVD.miss': Unable to complete matrix, too much missing data")
    reldiff <- conv.reldiff+1
    i <- 0
    while(i<niter && reldiff>conv.reldiff){
      svd0 <- svd(XF)
      Xnew <- X
      Xnew[is.na(X)] <- (svd0$u[, 1:ncomp] %*%
                         diag(svd0$d[1:ncomp],nrow=length(svd0$d[1:ncomp])) %*%
                         t(svd0$v[, 1:ncomp]))[is.na(X)]
      #Xnew[is.na(X)] <- (svd0$u %*% diag(svd0$d,nrow=length(svd0$d)) %*%
      #t(svd0$v))[is.na(X)]
      diff <- max(abs(Xnew - XF))
      reldiff <- diff/max(abs(XF[is.na(X)]))
      XF <- Xnew
      i <- i+1
    }
  }#if( all(!is.na(X)) ) ... else ...
  final.diff <- c(diff,reldiff,i)
  names(final.diff) <- c("diff","rel.diff","n.iter")
  return(list(svd=svd0, Xfill=XF, status=final.diff))
}

#smoother for the missing data SVD, uses cross-validated smooth splines.
SVD.smooth <- function(data, n.basis, date.ind=NA, scale.data=TRUE,
                       niter=100, conv.reldiff=0.001, df=NULL, spar=NULL){
  if( any(is.na(date.ind)) ){
    date.ind <- try(as.Date(rownames(data)),silent=TRUE)
    if( class(date.ind)=="try-error" )
      date.ind <- 1:dim(data)[1]
  }
  if( length(date.ind)!=dim(data)[1] )
    stop("In 'SVD.smooth': length(date.ind)!=dim(data)[1]")
  if( max(n.basis)>dim(data)[2] )
    stop("In 'SVD.smooth': Number of basis functions cannot exceed the number of columns in the data matrix")
  if( scale.data ) #scale data to zero mean and unit variance
    data <- scale(data)
  #missing data SVD
  data.svd <- SVD.miss(data, niter=niter, ncomp=n.basis,
                       conv.reldiff=conv.reldiff)
  #calculate basis functions
  data.comps <- matrix(NA, length(date.ind), n.basis)
  for(j in 1:n.basis){
    if( is.null(df) && is.null(spar) ){
      data.comps[,j] <- smooth.spline(date.ind,
                                      data.svd$svd$u[,j])$y
    }else{
      data.comps[,j] <- smooth.spline(date.ind, data.svd$svd$u[,j],
                                      df=df, spar=spar)$y
    }
    #scale components to unit variance and zero mean
    data.comps[,j] <- scale(data.comps[,j])
    #Ensure that components have alternating sign
    data.comps[,j] <- (-1)^j*data.comps[,j]*sign(data.comps[1,j])
  }
  rownames(data.comps) <- as.character(date.ind)
  colnames(data.comps) <- paste("V",1:n.basis,sep="")
  return(data.comps)
}

#cross-validation for the SVD smoother, either does cv for all the n.basis
#specified
SVD.smooth.cv <- function(data, n.basis, date.ind=NA, scale.data=TRUE,
                          niter=100, conv.reldiff=0.001, df=NULL, spar=NULL){
  if( max(n.basis)>(dim(data)[2]-1) )
    stop("In 'SVD.smooth.cv': Number of basis functions cannot exceed the number of columns in the data matrix minus one when doing CV.")
  #matrix holding CV statistics
  CV.stat <- matrix(NA,length(n.basis),3)
  colnames(CV.stat) <- c("RMSE","R2","BIC")
  rownames(CV.stat) <- paste("n.basis",as.character(n.basis),sep=".")
  #matrix that holds the individual BIC:S
  BIC <- matrix(NA,length(n.basis),dim(data)[2])
  colnames(BIC) <- colnames(data)
  rownames(BIC) <- rownames(CV.stat)
  RSS <- matrix(NA,length(n.basis),dim(data)[2])
  #loop over the number of basis functions requested
  trend <- list()
  for(i in 1:length(n.basis)){
    #for each basis functions do leave one out CV
    trend[[i]] <- array(NA,c(dim(data)[1],n.basis[i],dim(data)[2]))
    err <- matrix(NA,dim(data)[1],dim(data)[2])
    for(j in 1:dim(data)[2]){
      trend[[i]][,,j] <- SVD.smooth(data[,-j,drop=FALSE], n.basis[i],
                                    date.ind=date.ind, scale.data=scale.data,
                                    niter=niter, conv.reldiff=conv.reldiff,
                                    df=df, spar=spar)
      I <- !is.na(data[,j])
      #compute CV-error
      err[I,j] <- lm(data[,j]~trend[[i]][,,j])$residuals
      #compute BIC
      RSS[i,j] <- sum(err[I,j]^2)
      BIC[i,j] <- sum(I)*log(RSS[i,j]/sum(I)) + n.basis[i]*log(sum(I))
    }
    #compute CV statistics
    CV.stat[i,"RMSE"] <- sqrt(mean(err^2,na.rm=TRUE))
    CV.stat[i,"R2"] <- 1 - CV.stat[i,"RMSE"]^2/var(c(data),na.rm=TRUE)
    dimnames(trend[[i]]) <- list(rownames(data),
                                 paste("V",1:n.basis[i],sep=""),
                                 colnames(data))
  }
  n.all <- sum(!is.na(data))
  RSS.all <- rowSums(RSS)
  p.all <- (1+n.basis)*dim(data)[2]
  CV.stat[,"BIC"] <- n.all*log(RSS.all/n.all) + p.all*log(n.all)
  names(trend) <- rownames(CV.stat)
  return( list(CV.stat=as.data.frame(CV.stat),
               BIC.all=as.data.frame(t(BIC)),smooth.SVD=trend) )
}

#wraper functions that creates smoothed SVD:s from data structures.
#takes an option for also calculating leave one site out SVDs.
calc.smooth.trends <- function(mesa.data=NA, obs=mesa.data$obs$obs,
                               date=mesa.data$obs$date, ID=mesa.data$obs$ID,
                               subset=NA, n.basis=2, cv=FALSE, niter=100,
                               conv.reldiff=0.001, df=NULL, spar=NULL){
  data <- create.data.matrix(mesa.data=NA, obs=obs, date=date,
                             ID=ID, subset=subset)
  ##drop rows with no observations
  if( any(apply(data,1,function(x){all(is.na(x))})) )
    warning("In 'calc.smooth.trends': dates with no observations are being dropped, The choice of 'subset' may be to agressive.")
  data <- data[!apply(data,1,function(x){all(is.na(x))}),,drop=FALSE]
  ##now let's do SVD
  data.comps <- SVD.smooth(data, n.basis, scale.data=TRUE, niter=niter,
                           conv.reldiff=conv.reldiff, df=df, spar=spar)
  data.comps <- as.data.frame(data.comps)
  data.comps$date <- as.Date(rownames(data.comps))
  if(cv){
    svd.cv <- SVD.smooth.cv(data, n.basis, scale.data=TRUE, niter=niter,
                            conv.reldiff=conv.reldiff, df=df, spar=spar)
    svd.cv <- svd.cv$smooth.SVD[[1]]
    svd.tmp <- list()
    for(i in 1:dim(svd.cv)[3]){
      svd.tmp[[i]] <- as.data.frame(svd.cv[,,i])
      svd.tmp[[i]]$date <- as.Date(rownames(svd.cv))
    }
    names(svd.tmp) <- dimnames(svd.cv)[[3]]
  }else{
    svd.tmp <- NA
  }
  return( list(svd=data.comps,svd.cv=svd.tmp) )
}
