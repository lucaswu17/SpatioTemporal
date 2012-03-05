#################################################
## FILE CONTAINING THE LOGLIKELIHOOD FUNCTIONS ##
#################################################
loglike <- function(x=NA, mesa.data.model, type="p"){
  
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'loglike': Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")
  ##if x is given as NA (done to find out the ordering/names of input
  ##arguments) return expected parameter names
  if( any(is.na(x)) ){
    ##return the expected variable names
    return( loglike.var.names(mesa.data.model, all=(type=="f")) )
  }
  ##else, calculate loglikelihood
  ##first figure out a bunch of dimensions
  dimensions <- loglike.dim(mesa.data.model)
  
  ##check parameter sizes
  if((type=="f" && length(x)!=dimensions$nparam) ||
     (type!="f" && length(x)!=dimensions$nparam.cov))
    stop("In 'loglike': Size missmatch for x, number of parameters is incorrect")
  
  ##extract parameters from x
  tmp <- get.params(x,dimensions)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  trend.range <- tmp$range
  trend.sill <- tmp$sill
  phi.residuals <- tmp$phi.nu

  ##extract the observations
  Y <- mesa.data.model$obs$obs
  if(type=="f"){
    ##calculate the land use regression for the temporal trends
    mu.B <- calc.mu.B(dimensions, mesa.data.model$X, alpha)
    ##stack into one long vector
    dim(mu.B) <- c(dimensions$n * dimensions$m, 1)
    if(dimensions$L!=0){ ##we have spatio temporal covariates, subtract these
      Y <- Y - mesa.data.model$SpatioTemp %*% gamma
    }
  }

  ##create covariance matrices
  sigma.B <- make.sigma.B(trend.sill, trend.range, mesa.data.model$dist)
  ##parameter order is sill, nugget, range
  sigma.nu <- make.sigma.nu(phi.residuals[2], phi.residuals[3],
                              phi.residuals[1],
                              block.sizes=mesa.data.model$nt,
                              loc.index=mesa.data.model$obs$idx,
                              dists=mesa.data.model$dist)
  ##calculate block cholesky factor of the matrices
  ##storing in-place to conserve memory
  sigma.B <- try(makeCholBlock(sigma.B, n.blocks=dimensions$m),
                 silent=TRUE)
  if(class(sigma.B)=="try-error")
    return(-.Machine$double.xmax)
  sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks=dimensions$T,
                                block.sizes=mesa.data.model$nt), silent=TRUE)
  if(class(sigma.nu)=="try-error")
    return(-.Machine$double.xmax)

  ##loglikelihood calculations follow:
  ##-log(det(sigma.B)^.5)
  l <- -sumLogDiag( sigma.B )
  ##-log(det(sigma.nu)^.5)
  l <- l -sumLogDiag( sigma.nu )

  ##invert the matrices
  ##calculate inverse of sigma.B
  sigma.B <- invCholBlock(sigma.B, n.blocks=dimensions$m)
  ##calculate inverse of sigma.nu
  sigma.nu <- invCholBlock(sigma.nu, n.blocks=dimensions$T,
                              block.sizes=mesa.data.model$nt)
  ##inv(sigma.nu) * Y
  i.sR.Y <- block.mult(sigma.nu, Y, n.x=1,
                       n.blocks=dimensions$T,
                       block.sizes=mesa.data.model$nt)
  if(type=="f"){
    ##calculate inv(matrix)*mu
    ##inv(sigma.B) * mu.B
    i.sB.mu.B <- block.mult(sigma.B, mu.B, n.blocks=dimensions$m)
    ##-1/2 mu.B' * inv(sigma.B) * mu.B
    l <- l - dot.prod(i.sB.mu.B,mu.B)/2
    ##-1/2 Y' * inv(sigma.nu) * Y
    l <- l - dot.prod(i.sR.Y,Y)/2
  }else{
    ##inv(sigma.B) * X
    iS.X <- calc.iS.X(mesa.data.model$X, sigma.B, dimensions)
    if(dimensions$L!=0){
      ##inv(sigma.nu) * M
      i.sR.M <- block.mult(sigma.nu, mesa.data.model$SpatioTemp,
                           n.blocks=dimensions$T,
                           block.sizes=mesa.data.model$nt)
    }
    ##F'*inv(sigma.nu)*Y
    F.i.sR.Y <- calc.tF.times.mat(i.sR.Y, mesa.data.model$F,
                                  mesa.data.model$obs$idx, n.x=1,
                                  n.loc=dimensions$n)
    ##F'*inv(sigma.nu)*M
    if(dimensions$L!=0)
      F.i.sR.M <- calc.tF.times.mat(i.sR.M, mesa.data.model$F,
                                    mesa.data.model$obs$idx,
                                    n.loc=dimensions$n)
  }
  ##F'*inv(sigma.nu)*F
  sigma.B.Y <- calc.tF.mat.F(sigma.nu, mesa.data.model$F,
                             mesa.data.model$obs$idx,
                             n.blocks=dimensions$T,
                             block.sizes=mesa.data.model$nt,
                             n.loc=dimensions$n)
  ##inv(sigma.B|Y) = F'*inv(sigma.nu)*F + inv(sigma.B)
  sigma.B.Y <- sigma.B.Y + sigma.B
  ##calculate cholesky factor of inv(sigma.B|Y)
  sigma.B.Y <- try(makeCholBlock(sigma.B.Y, n.blocks=1), silent=TRUE)
  if(class(sigma.B.Y)=="try-error")
    return(-.Machine$double.xmax)

  ##-log(det( inv(sigma.B|Y) )^.5)
  l <- l - sumLogDiag( sigma.B.Y )

  if(type=="f"){
    ##mu.B.Y = inv(sigma.B)*mu.B + F'*inv(sigma.nu)*Y
    mu.B.Y <- i.sB.mu.B + calc.tF.times.mat(i.sR.Y, mesa.data.model$F,
                                            mesa.data.model$obs$idx,
                                            n.loc=dimensions$n)
    ##calculate inv(R)'*mu.B.Y
    mu.B.Y <- solveTriBlock(sigma.B.Y, mu.B.Y, transpose=TRUE)
    ##+1/2 mu.B.Y * inv(i.sigma.B.Y) * mu.B.Y
    l <- l + norm2(mu.B.Y)/2
  }else{
    ##chol(inv(sigma.B.Y))^-T * F' * inv(sigma.nu) * Y
    Y.hat <- solveTriBlock(sigma.B.Y, F.i.sR.Y, transpose=TRUE)
    ##chol(inv(sigma.B.Y))^-T * F' * inv(sigma.nu) * M
    if(dimensions$L!=0)
      M.hat <- solveTriBlock(sigma.B.Y, F.i.sR.M, transpose=TRUE)
    ##calculate inverse of inv(sigma.B|Y) (keep in place)
    sigma.B.Y <- invCholBlock(sigma.B.Y, n.blocks=1)
    ##sigma.B.Y * inv(sigma.B) * X
    sigma.B.Y.iS.X <- sigma.B.Y %*% iS.X
    ##inv(sigma.alpha|Y) = X'*inv(sigma.B)*X -
    ##       X'*inv(sigma.B)*sigma.B|Y*inv(sigma.B)*X
    i.sigma.alpha.Y <- -t(iS.X) %*% sigma.B.Y.iS.X +
      calc.X.iS.X(mesa.data.model$X, iS.X, dimensions)
    ##calculate cholesky factor of inv(sigma.alpha|Y)
    i.sigma.alpha.Y <- try(makeCholBlock(i.sigma.alpha.Y, n.blocks=1),
                           silent=TRUE)
    if(class(i.sigma.alpha.Y)=="try-error")
      return(-.Machine$double.xmax)
    if(type=="r")
      l <- l - sumLogDiag( i.sigma.alpha.Y )
    ## X' * inv(sigma.B) * sigma.B.Y * F' * inv(sigma.nu) * Y
    Y.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.Y
    ## chol(inv(sigma.alpha.Y))^-T * Y.hat.2
    Y.hat.2 <- solveTriBlock(i.sigma.alpha.Y, Y.hat.2, transpose=TRUE)
    ##Y'*sigma_hat*Y
    Y.sigma.hat.Y <- dot.prod(Y,i.sR.Y) -
      norm2(Y.hat) - norm2(Y.hat.2)
    l <- l - Y.sigma.hat.Y/2
  
    ##parts relevant only with spatio-temporal model
    if(dimensions$L!=0){
      ## X' * inv(sigma.B) * sigma.B.Y * F' * inv(sigma.nu) * M
      M.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.M
      ## chol(inv(sigma.alpha.Y))^-T * M.hat.2
      M.hat.2 <- solveTriBlock(i.sigma.alpha.Y, M.hat.2, transpose=TRUE)
    
      Y.sigma.hat.M <- t(Y) %*% i.sR.M -
        t(Y.hat) %*% M.hat - t(Y.hat.2) %*% M.hat.2
      Y.sigma.hat.M <- t(Y.sigma.hat.M)
      M.sigma.hat.M <- t(mesa.data.model$SpatioTemp) %*% i.sR.M -
        t(M.hat) %*% M.hat - t(M.hat.2) %*% M.hat.2
      ##calculate cholesky factor of M.sigma.hat.M
      M.sigma.hat.M <- try(makeCholBlock(M.sigma.hat.M, n.blocks=1),
                           silent=TRUE)
      if(class(M.sigma.hat.M)=="try-error")
        return(-.Machine$double.xmax)
      if(type=="r") ##-log(det( M.sigma.hat.M )^.5)
        l <- l - sumLogDiag( M.sigma.hat.M )

      ## chol(inv(sigma.alpha.Y))^-T * Y.hat.2
      Y.sigma.hat.M <- solveTriBlock(M.sigma.hat.M, Y.sigma.hat.M, transpose=TRUE)
      l <- l + norm2(Y.sigma.hat.M)/2
    }##if(dimensions$L!=0)
  }##if(type=="f") ... else ...

  ##ensure that l is treated as a double value (and not a sparse matrix)
  l <- as.double(l)
  ##add safe guard (infinite or nan results almost always due to
  ##matrix inprecision, lets return a very small value)
  if( !is.finite(l) )
    l <- -.Machine$double.xmax
  return(l)
}

###############################################
## Log-likelihood using the full formulation ##
###############################################
loglike.naive <- function(x=NA, mesa.data.model, type="p"){
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'loglike.naive': Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")
  ##if x is given as NA (done to find out the ordering/names of input
  ##arguments) return expected parameter names
  if( any(is.na(x)) ){
    ##return the expected variable names
    return( loglike.var.names(mesa.data.model, all=(type=="f")) )
  }
  ##else, calculate loglikelihood
  
  ##first figure out a bunch of dimensions
  dimensions <- loglike.dim(mesa.data.model)
  ##check parameter sizes
  if((type=="f" && length(x)!=dimensions$nparam) ||
     (type!="f" && length(x)!=dimensions$nparam.cov))
    stop("In 'loglike.naive': Size missmatch for x, number of parameters is incorrect")
    
  ##extract parameters from x
  tmp <- get.params(x,dimensions)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  trend.range <- tmp$range
  trend.sill <- tmp$sill
  phi.residuals <- tmp$phi.nu

  ##extract the observations
  Y <- mesa.data.model$obs$obs
  if(type=="f"){
    ##calculate the land use regression for the temporal trends
    mu.B <- calc.mu.B(dimensions, mesa.data.model$X, alpha)
    ##multiply the beta fields with appropriate trends and sum to get the
    ##expectation
    mean.val <- 0
    for(i in c(1:dimensions$m) ) ##loop over diff. trends
      mean.val <- mean.val + mu.B[mesa.data.model$obs$idx,i]*mesa.data.model$F[,i]
    ##subtract mean value from observations
    Y <- Y - mean.val
    if(dimensions$L!=0) ##also subtract spatio-termporal covariate
      Y <- Y - mesa.data.model$SpatioTemp %*% gamma
  }
  
  ##create covariance matrices
  sigma.B.full <- make.sigma.B.full(trend.sill, trend.range,
                                    loc.ind1=mesa.data.model$obs$idx, 
                                    F1=mesa.data.model$F,
                                    dists=mesa.data.model$dist)
  ##parameter order is sill, nugget, range
  sigma.nu <- make.sigma.nu(phi.residuals[2], phi.residuals[3],
                            phi.residuals[1],
                            block.sizes=mesa.data.model$nt,
                            loc.index=mesa.data.model$obs$idx,
                            dists=mesa.data.model$dist)
  ##Total covariance matrix
  sigma.nu <- sigma.nu + sigma.B.full
  ##calculate (block) cholesky factor of the matrices
  ##storing in-place to conserve memory
  sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks=1),silent=TRUE)
  if(class(sigma.nu)=="try-error")
    return(-.Machine$double.xmax)

  ##loglikelihood calculations follow:
  ##-log(det(sigma.nu)^.5)
  l <- -sumLogDiag( sigma.nu )
  ##calculate if(type=="f") inv(R)'*(Y-mean.val) else inv(R)'*Y
  Y <- solveTriBlock(sigma.nu, Y, n.x=1, transpose=TRUE)
  ##if(type=="f")  -1/2 (Y-mean)' * inv(sigma.nu) * (Y-mean)
  ##     else      -1/2 Y' * inv(sigma.nu) * Y
  l <- l - norm2(Y)/2
  if(type!="f"){
    ##Create the Ftmp = [FX M] matrix
    Ftmp <- calc.F.times.X(mesa.data.model$X, mesa.data.model$F,
                           mesa.data.model$obs$idx)
    ##Add the spatio-temporal covariate (if it exists)
    if( dimensions$L!=0 )
      Ftmp <- cbind(Ftmp,mesa.data.model$SpatioTemp)
  
    ##calculate inv(R)'*Ftmp
    Ftmp <- solveTriBlock(sigma.nu, Ftmp, transpose=TRUE)
    ##calculate Ftmp'*inv(Sigma)*Y
    FY <- t(Ftmp) %*% Y
    ##calculate [FX M]*invSigma*[FX M]'
    sigma.alt <- t(Ftmp) %*% Ftmp
    ##calculate cholesky factor (storing in-place to conserve memory)
    sigma.alt <- try(makeCholBlock(sigma.alt, n.blocks=1),silent=TRUE)
    if(class(sigma.alt)=="try-error")
      return(-.Machine$double.xmax)
    
    if(type=="r"){
      ##-log(det(sigma.alt)^.5)
      l <- l - sumLogDiag( sigma.alt )
    }
    ##calculate inv(R)'*Y
    FY <- solveTriBlock(sigma.alt, FY, n.x=1,transpose=TRUE)
    ##+1/2 FY' * inv(sigma.alt) * FY
    l <- l + norm2(FY)/2
  }##if(type!="f")
  
  ##ensure that l is treated as a double value (and not a sparse matrix)
  l <- as.double(l)
  ##add safe guard (infinite or nan results almost always due to
  ##matrix inprecision, lets return a very small value)
  if( !is.finite(l) )
    l <- -.Machine$double.xmax
  return(l)
}
