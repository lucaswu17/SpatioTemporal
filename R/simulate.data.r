################################################
## FUNCTION THAT SIMULATE DATA FROM THE MODEL ##
################################################
##simulate data from a MESA data set using the locations and
##covariates given in mesa.data
##parameters are given by x (same parameters as for the likelihood calculation)
simulateMesaData <- function(x, mesa.data.model, mesa.data=NA, rep=1,
                             combine.data=FALSE){
  if( combine.data==TRUE && all(is.na(mesa.data)) ){
    combine.data <- FALSE
    warning("In 'cond.expectation': No data to combine, using mesa.data.model",
            immediate. = TRUE)
  }

  ##first figure out a bunch of dimensions
  dimensions <- loglike.dim(mesa.data.model)
  ##extract and check parameters
  if(length(x)==dimensions$nparam.cov){
    ##compute alpha and gamma as cond. exp. given data
    tmp <- cond.expectation(x, mesa.data.model, only.pars=TRUE,
                            type="p")$pars
    x <- c(tmp$gamma.E,tmp$alpha.E, x)
  }
  if(length(x)!=dimensions$nparam){
    stop(sprintf("In 'simulateMesaData': 'x' should have %d or %d elements",
                 dimensions$nparam.cov,dimensions$nparam))
  }
  ##add names to the parameters for future reference
  names(x) <- loglike.var.names(mesa.data.model)

  ##now see if we need to merge or change mesa.data.model
  if( all(!is.na(mesa.data)) ){
    if(combine.data){
      ##combine the two datasets
      mesa.data.model <- combineMesaData(mesa.data.model, mesa.data)
    }else{
      if( is.null(mesa.data$trend) ){
        warning("In 'simulateMesaData': mesa.data lacking a trend element, using trend from mesa.data.model.", immediate.=TRUE)
        mesa.data$trend <- mesa.data.model$trend
      }
      ##create a mesa.data.model from mesa.data
      mesa.data.model <- create.data.model(mesa.data, mesa.data.model$LUR.list,
                                           mesa.data.model$ST.Ind, strip=FALSE)
    }
  }
  ##mesa.data.model might have changed due to merging, recompute dimensions
  dimensions <- loglike.dim(mesa.data.model)

  ##timepoints (create.data.model has already
  ##combined dates from obs and trend)
  T <- mesa.data.model$trend$date
  
  ##create new time vectors F, one for each observation time.
  ##construct a matrix holding the temporal basis for each observation
  ##(including the constant)
  F <- matrix(0,length(T),dimensions$m)
  F[,1] <- 1
  if(dimensions$m>1){
    tmp.trend <- mesa.data.model$trend
    tmp.trend$date <- NULL
    F[,c(2:dimensions$m)] <- as.matrix(tmp.trend)
  }
  F <- as.data.frame(F)
  if(dimensions$m>1)
    names(F)[2:dimensions$m] <- names(tmp.trend)
  names(F)[1] = "const"
  F <- data.matrix(F)
  
  ##now lets extract parameters from x
  tmp <- get.params(x,dimensions)
  gamma <- tmp$gamma
  alpha <- tmp$alpha
  trend.range <- tmp$range
  trend.sill <- tmp$sill
  phi.residuals <- tmp$phi.nu

  ##calculate the land use regression for the temporal trends
  mu.B <- calc.mu.B(dimensions, mesa.data.model$X, alpha)
  ##stack into one long vector
  dim(mu.B) <- c(dimensions$n * dimensions$m, 1)
  
  ##create covariance matrices
  sigma.B <- make.sigma.B(trend.sill, trend.range, mesa.data.model$dist)
  ##parameter order is sill, nugget, range
  ##just create residual matrix for all locations but one time point,
  ##Due to independence we just sample repeatedly from this.
  sigma.nu <- make.sigma.nu(phi.residuals[2], phi.residuals[3],
                            phi.residuals[1],
                            block.sizes=c(dimensions$n),
                            loc.index=c(1:dimensions$n),
                            dists=mesa.data.model$dist)
  
  ##calculate block cholesky factor of the matrices
  ##if any distance==0 we know this will fail for the beta fields
  if( any(mesa.data.model$dist[upper.tri(mesa.data.model$dist)]==0) ){
    Rsigma.B <- NA
    class(Rsigma.B) <- "try-error"
  }else{
    Rsigma.B <- try(makeCholBlock(sigma.B, n.blocks=dimensions$m),
                    silent=TRUE)
  }
  Rsigma.nu <- try(makeCholBlock(sigma.nu, n.blocks=1,
                                 block.sizes=dimensions$n), silent=TRUE)
  B <- array(0,c(dimensions$n,dimensions$m,rep),list(NULL,NULL,NULL))
  dimnames(B)[[1]] <- as.character(mesa.data.model$location$ID)
  dimnames(B)[[2]] <- names(mesa.data.model$X)
  dimnames(B)[[3]] <- paste("rep",as.character(1:rep),sep=".")
  ##simulate data from B
  if( class(Rsigma.B)=="try-error" ){
    for(j in 1:dimensions$m){
      Ind <- (1:dimensions$n) + (j-1)*dimensions$n
      B[,j,] <- matrix(mu.B[Ind],dimensions$n,rep) + 
        t(matrix(mvrnorm(n=rep,mu=rep(0,length(Ind)),Sigma=sigma.B[Ind,Ind]),
                 rep,length(Ind)))
    }
  }else{
    for(j in 1:dimensions$m){
      Ind <- (1:dimensions$n) + (j-1)*dimensions$n
      R <- t(Rsigma.B[Ind,Ind])
      B[,j,] <- matrix(mu.B[Ind],dimensions$n,rep) + 
        R %*% matrix(rnorm(dim(R)[1]*rep),dim(R)[1],rep)
    }
  }
  ##create a time points by observations points matrix
  ##to hold the simulated data.
  X <- array(0,c(dim(F)[1], dimensions$n,rep),list(NULL,NULL,NULL))
  dimnames(X)[[1]] <- as.character(as.Date(T))
  dimnames(X)[[2]] <- as.character(mesa.data.model$location$ID)
  dimnames(X)[[3]] <- paste("rep",as.character(1:rep),sep=".")
  ##if needed add spatio-temporal covariate(s)
  if(dimensions$L!=0){
    ##extract relevant spatio-temporal covariates
    ST <- mesa.data.model$SpatioTemp.all
    ##extract relevant locations and timepoints
    IND <- match(dimnames(X)[[1]],rownames(ST))
    ST <- ST[IND,,,drop=FALSE]
    IND <- match(dimnames(X)[[2]],colnames(ST))
    ST <- ST[,IND,,drop=FALSE]
   
    for(i in 1:dim(ST)[3])
      X[,,1] <- X[,,1] + gamma[i]*ST[,,i]
    ##and replicate for the other X:s
    if(rep>1)
      for(j in 2:rep)
        X[,,j] <- X[,,1]
  }
    
  ##now combine simulated B, with the temporal trends and simulated
  ##data from the residuals field
  for(k in 1:rep)
    X[,,k] <- X[,,k] + F %*% t(B[,,k])
  ##simulate data from B
  if( class(Rsigma.nu)=="try-error" ){
    e <- t(mvrnorm(n=prod(dim(X)[-2]),mu=rep(0,dim(X)[2]),
                   Sigma=sigma.nu))
  }else{
    ##simulate from the residuals
    e <- rnorm( prod(dim(X)) )
    dim(e) <- c(dim(X)[2],prod(dim(X)[-2]))
    e <- t(Rsigma.nu) %*% e
  }
  e <- array(e,dim(X)[c(2,1,3)])
  e <- aperm(e, c(2,1,3))
  X <- X + e

  ##lets order the columns to match locations in mesa.data
  if( combine.data==FALSE && all(!is.na(mesa.data)) )
    X <- X[,as.character(mesa.data$location$ID),,drop=FALSE]

  ##extract a vector matching the observed locations (i.e. matching the
  ##obs vector in data.obs)
  if( length(mesa.data.model$obs$obs)==0 ){
    obs = NULL
  }else{
    I <- (match(as.character(mesa.data.model$location$ID[mesa.data.model$obs$idx]),dimnames(X)[[2]])-1)*dim(X)[1] + 
      match(mesa.data.model$obs$date,as.Date(dimnames(X)[[1]]))
    obs <- list()
    for(i in 1:rep){
      obs[[i]] <- mesa.data.model$obs
      tmp <- X[,,i]
      obs[[i]]$obs <- tmp[I]
      ##sort the observations (should be done since mesa.data.model$obs
      ##are sorted, but better safe than sorry)
      IND <- order(obs[[i]]$date,obs[[i]]$idx)
      obs[[i]] <- obs[[i]][IND,,drop=FALSE]
    }
  }
  return( list(param=x,B=B,X=X,obs=obs) )
}##simulateMesaData
