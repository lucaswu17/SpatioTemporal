############################################
## FUNCTIONS THAT DO ESTIMATION/INFERENCE ##
############################################

##################################
## Wraper for parameter fitting ##
##################################
##function that attempts to fit the parameters of the model
##if x is a matrix it uses each of the columns in x as starting values
##returning the best result with a positive definite matrix.
fit.mesa.model <- function(x, mesa.data.model, type="p",
                           h=1e-3, diff.type=1, lower=-15, upper=15,
                           hessian.all=FALSE,
                           control=list(trace=3, maxit=1000)){
  ##ensure lower case
  type <- tolower(type) 
  ##first check that type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'fit.mesa.model': Unknown option for type. Should be: (f)ull, (r)eml, (p)rofile")

  ##get size of the models
  dimensions <- loglike.dim(mesa.data.model)
  
  ##Second check the input x
  if(is.vector(x))
    x <- as.matrix(x)
  if(length(x)==1){
    ##since x has length one we need to construct a matrix
    ##of initialvalues
    x <- matrix(NA, dimensions$nparam.cov, x)
    init <- seq(-5, 5, length.out=dim(x)[2])
    for(i in 1:dim(x)[2])
      x[,i] <- init[i]
  }

  ##see if size of x matches the expected number of parameters
  ##if not we might want to either truncate or expand
  if(dim(x)[1] != dimensions$nparam && dim(x)[1] != dimensions$nparam.cov)
    stop("In 'fit.mesa.model': Size missmatch for x, number of parameters is incorrect")
  if(dim(x)[1] == dimensions$nparam && type!="f"){
    ##requested REML or profile but provided full parameter
    ##starting points -> truncate
    x <- x[(dimensions$nparam-dimensions$nparam.cov+1):dim(x)[1],]
  }else if(dim(x)[1] == dimensions$nparam.cov && type=="f"){
    ##requested full but provided parameters for REML or profile -> expand
    x.old <- x
    x <- matrix(NA,dimensions$nparam,dim(x)[2])
    for(i in 1:dim(x)[2]){
      tmp <- cond.expectation(x, mesa.data.model, only.pars=TRUE,
                              type="p")$pars
      x[,i] <- c(tmp$gamma.E, tmp$alpha.E, x.old[,i])
    }
  }

  ##define lokal version of the loglikelihood and gradient functions
  loglike.loc <- function(x){ loglike(x, mesa.data.model, type) }
  loglike.grad.loc <- function(x){
    loglike.grad(x, mesa.data.model, type, h=h, diff.type=diff.type)
  }

  ##attempt to fit the model for each of the provided starting values.
  res <- as.list(rep(NA,dim(x)[2]))
  ##vector with convergence information and optimal value for
  ##all the starting points
  conv <- rep(FALSE,dim(x)[2])
  value <- rep(NA,dim(x)[2])
  ##ensure that we are doing maximisation
  control$fnscale=-1
  for(i in 1:dim(x)[2]){
    if(!is.null(control$trace) && control$trace!=0){
      cat( sprintf("In 'fit.mesa.model': Optimisation using starting value %d of %d\n",
                     i,dim(x)[2]) )
    }
    try(res[[i]] <- optim(x[,i], loglike.loc, gr=loglike.grad.loc,
                          method = "L-BFGS-B", control=control,
                          hessian=TRUE, lower=lower, upper=upper),
        silent=TRUE)
    if(all(!is.na(res[[i]]))){
      ##if this is a valid optim point, then compute convergence criteria
      conv[i] <- (res[[i]]$convergence==0 &&
                  all(eigen(res[[i]]$hessian)$value < -1e-10))
      ##extract ML-value
      value[i] <- res[[i]]$value
      ##add convergence and initial parameters
      res[[i]]$conv <- conv[i]
      res[[i]]$par.init <- x[,i]
      ##Add expanded parameters and parameter names
      if(type!="f"){
        ##for profile or REML, we need to expand the parameters
        tmp <- cond.expectation(res[[i]]$par, mesa.data.model, only.pars=TRUE,
                                type=type)$pars
        res[[i]]$par.all <- c(tmp$gamma.E, tmp$alpha.E, res[[i]]$par)
      }else{
        ##replicate all parameters so output is consistent
        res[[i]]$par.all <- res[[i]]$par
      }
      ##add names to the variables.
      names(res[[i]]$par) <- loglike.var.names(mesa.data.model, all=(type=="f"))
      names(res[[i]]$par.init) <- names(res[[i]]$par)
      names(res[[i]]$par.all) <- loglike.var.names(mesa.data.model, all=TRUE)
    }##if(all(!is.na(res[[i]])))
  }##for(i in 1:dim(x)[2])
  if(all(is.na(res)))
    stop("In 'fit.mesa.model': All optimisations failed")
  ##pick out the converged option with the best value
  Ind.overall <- which.max(value)
  if(any(conv==TRUE)){
    ##mark no-converged values as NA to avoid picking these
    value[!conv] <- NA
  }
  ##extract the best value
  Ind <- which.max(value)
  res.best <- res[[Ind]]

  ##figure out status of optimisation
  optim.status <- sprintf("Optimisation: %d converged, %d not converged, %d failed. Best result found for starting point %d (optimisation has %sconverged).%s",
                          sum(conv), sum(!conv), sum(is.na(res)), Ind.overall,
                          switch(conv[Ind.overall]+1,"not ",""),
                          switch(conv[Ind]+1," No optimisations converged.",
                                 sprintf(" Best converged result for starting point %d.",Ind)))
  
  if(hessian.all==TRUE){
    if(type!="f"){
      ##find hessian for all the parameters
      res.best$hessian.all <- loglike.hessian(res.best$par.all, mesa.data.model,
                                              type="f",h=h)
    }else{
      ##replicate all parameters so output is consistent
      res.best$hessian.all <- res.best$hessian
    }
  }##if(hessian.all==TRUE)
  
  ##return result
  return(list(res.best=res.best,res.all=res,message=optim.status))
}##fit.mesa.model

################################################################
## conditional expectation for all locations given parameters ##
################################################################
##input is:
## estimated parameters
## data structure used for estimation
## covariates for the site where we want predictions, should
##   include the equivalent of the
##   $covar, $LUR, and $trend variabels in mesa.data
## the maximum number of sites for which to do predictions per batch.
##   Reasonalbe values depend on available memory, defaults to 5000
## if we should do conditional expectation with or with out nugget
##   (no.nugget=TRUE implies smoothing of observed values, no.nugget=FALSE
##   the values at observed locations will be equal to the observations.
cond.expectation <- function(x, mesa.data.model, mesa.data=NA, Nmax=1000,
                             only.pars=FALSE, compute.beta=FALSE, no.nugget=FALSE,
                             only.obs=FALSE, pred.var=TRUE, pred.covar=FALSE,
                             combine.data=FALSE, type="p"){
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'cond.expectation': Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")
  ##check if only.obs is valid
  if( only.obs && (all(is.na(mesa.data)) || length(mesa.data$obs$obs)==0) )
    stop("In 'cond.expectation': only.obs==TRUE makes no-sense when length(mesa.data$obs$obs)==0")
  ##check for inocompatible options
  if(only.pars && type=="f")
    warning("In 'cond.expectation': Incompatible options: only.pars=TRUE and type=(f)ull.")
  if( combine.data==TRUE && all(is.na(mesa.data)) ){
    combine.data <- FALSE
    warning("In 'cond.expectation': No data to combine, using mesa.data.model",
            immediate. = TRUE)
  }
  if(combine.data==TRUE && only.obs==TRUE)
    stop("In 'cond.expectation': Incompatible options: combine.data and only.obs cannot both be TRUE.")
  if(pred.covar){
    if(only.obs)
      stop("In 'cond.expectation': Incompatible options: pred.covar and only.obs cannot both be TRUE.")
    if(!pred.var)
      warning("In 'cond.expectation': pred.covar=TRUE implies pred.var=TRUE.",
              immediate. = TRUE)
    pred.var <- TRUE
  }
    
  ##create the covar.aux used for predictions
  if( all(is.na(mesa.data)) ){
    covar.aux <- mesa.data.model
  }else{
    if(compute.beta || combine.data)
      ##combine the two datasets
      mesa.data.model.all <- combineMesaData(mesa.data.model, mesa.data)
    if(compute.beta){
      mesa.data.model <- mesa.data.model.all
    }
    if(combine.data){
      covar.aux <- mesa.data.model.all
    }else{
      if( is.null(mesa.data$trend) ){
    warning("In 'cond.expectation': mesa.data lacking a trend element, using trend from mesa.data.model.", immediate.=TRUE)
        mesa.data$trend <- mesa.data.model$trend
      }
      ##extract the land-use-regresion for the unobserved locations
      covar.aux <- create.data.model(mesa.data, LUR=mesa.data.model$LUR.list,
                                     ST.Ind=mesa.data.model$ST.Ind, strip=FALSE)
    }
  }
  ##figure out a bunch of dimensions
  dimensions <- loglike.dim(mesa.data.model)
  ##check if inputs are of reasonable size.
  if(length(x) != dimensions$nparam && length(x) != dimensions$nparam.cov)
    stop("In 'cond.expectation': Size missmatch for x, number of parameters is incorrect")
  
  ##first lets extract parameters from x
  tmp <- get.params(x,dimensions)
  trend.range <- tmp$range
  trend.sill <- tmp$sill
  phi.residuals <- tmp$phi.nu
  
  ##Create the Xtilde = [M FX] matrix
  Xtilde <- calc.F.times.X(mesa.data.model$X, mesa.data.model$F,
                           mesa.data.model$obs$idx)
  ##Add the spatio-temporal covariate (if it exists)
  if( dimensions$L!=0 )
    Xtilde <- cbind(mesa.data.model$SpatioTemp,Xtilde)
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
  sigma.B <- makeCholBlock(sigma.B, n.blocks=dimensions$m)
  sigma.nu <- makeCholBlock(sigma.nu, n.blocks=dimensions$T,
                               block.sizes=mesa.data.model$nt)
  ##invert the matrices
  ##storing in-place to conserve memory
  sigma.B <- invCholBlock(sigma.B, n.blocks=dimensions$m)
  sigma.nu <- invCholBlock(sigma.nu, n.blocks=dimensions$T,
                              block.sizes=mesa.data.model$nt)
  ##F'*inv(sigma.nu)*F
  sigma.B.Y <- calc.tF.mat.F(sigma.nu, mesa.data.model$F,
                             mesa.data.model$obs$idx,
                             n.blocks=dimensions$T,
                             block.sizes=mesa.data.model$nt,
                             n.loc=dimensions$n)
  ##inv(sigma.B|Y) = F'*inv(sigma.nu)*F + inv(sigma.B)
  sigma.B.Y <- sigma.B.Y + sigma.B
  ##calculate cholesky factor of inv(sigma.B|Y)
  sigma.B.Y <- makeCholBlock(sigma.B.Y, n.blocks=1)

  ##First lets compute the regression parameters (alpha and gamma)
  if(length(x)==dimensions$nparam && type=="f"){
    ##regression parameters given
    gamma.E <- c(tmp$gamma)
    alpha.E <- unlist(tmp$alpha)
    ##known parameters, set variances to zero
    gamma.V <- matrix(0,length(gamma.E),length(gamma.E))
    alpha.V <- matrix(0,length(alpha.E),length(alpha.E))
    gamma.alpha.C <- matrix(0,length(gamma.E),length(alpha.E))
    ##create a combined vector
    gamma.alpha <- as.matrix(c(gamma.E,alpha.E))
  }else{
    ##compute iSigma.nu*Xtilde
    iS.nu.X <- block.mult(sigma.nu, Xtilde, n.blocks=dimensions$T,
                          block.sizes=mesa.data.model$nt)
    ##compute F'*iSigma.nu*Xtilde
    tF.iS.nu.X <- calc.tF.times.mat(iS.nu.X, mesa.data.model$F, 
                                    mesa.data.model$obs$idx,
                                    n.loc=dimensions$n)
    ##compute inv(sigma.B|Y)*F'*iSigma.nu*Xtilde
    iSBY.tF.iS.X <- solveTriBlock(sigma.B.Y, tF.iS.nu.X,
                                  n.blocks=1, tr=TRUE)
    iSBY.tF.iS.X <- solveTriBlock(sigma.B.Y, iSBY.tF.iS.X,
                                  n.blocks=1, tr=FALSE)
    ##compute inv(Xtilde'*Sigma^-1*Xtilde)
    i.XSX <- t(Xtilde) %*% iS.nu.X
    i.XSX <- i.XSX - t(tF.iS.nu.X) %*% iSBY.tF.iS.X
  
    ##compute iSigma.nu*Y
    iS.nu.Y <- block.mult(sigma.nu, mesa.data.model$obs$obs,
                          n.x=1, n.blocks=dimensions$T,
                          block.sizes=mesa.data.model$nt)
    ##compute F'*iSigma.nu*Y
    tF.iS.nu.Y <- calc.tF.times.mat(iS.nu.Y, mesa.data.model$F, 
                                    mesa.data.model$obs$idx,
                                    n.loc=dimensions$n)
    ##compute t(Xtilde'*Sigma^-1*Y)
    tmp2 <- mesa.data.model$obs$obs %*% iS.nu.X
    tmp2 <- tmp2 - t(tF.iS.nu.Y) %*% iSBY.tF.iS.X
    ##compute alpha and gamma
    gamma.alpha <- solve(i.XSX, t(tmp2))
    ##remove parts not needed
    rm(iS.nu.X, tF.iS.nu.X, iSBY.tF.iS.X, iS.nu.Y, tF.iS.nu.Y)
    ##extract computed parameters
    if(dimensions$L!=0){
      gamma.E <- c(gamma.alpha[1:dimensions$L])
      alpha.E <- c(gamma.alpha[-(1:dimensions$L)])
    }else{
      gamma.E <- double(0)
      alpha.E <- c(gamma.alpha)
    }
    ##compute the variances
    i.XSX <- solve(i.XSX)
    if(dimensions$L!=0){
      gamma.V <- i.XSX[1:dimensions$L,1:dimensions$L,drop=FALSE]
      alpha.V <- i.XSX[-c(1:dimensions$L),-c(1:dimensions$L),drop=FALSE]
      gamma.alpha.C <- i.XSX[1:dimensions$L,-c(1:dimensions$L),drop=FALSE]
    }else{
      gamma.V <- matrix(0,0,0)
      alpha.V <- i.XSX
      gamma.alpha.C <- matrix(0,length(gamma.E),length(alpha.E))
    }
  }
  ##force to matrices
  gamma.E <- as.matrix(gamma.E)
  alpha.E <- as.matrix(alpha.E)
  ##add names
  names.tmp <- loglike.var.names(mesa.data.model)
  names.tmp <- names.tmp[1:(dimensions$nparam-dimensions$nparam.cov)]
  if(dimensions$L!=0){
    ##first for gamma
    rownames(gamma.E) <- names.tmp[1:dimensions$L]
    colnames(gamma.V) <- rownames(gamma.V) <- rownames(gamma.E)
    rownames(gamma.alpha.C) <- rownames(gamma.E)
    rownames(alpha.E) <- names.tmp[-(1:dimensions$L)]
  }else{
    rownames(alpha.E) <- names.tmp
  }
  ##and then for alpha
  colnames(alpha.V) <- rownames(alpha.V) <- rownames(alpha.E)
  colnames(gamma.alpha.C) <- rownames(alpha.E)
  ##collect the results
  pars <- list(gamma.E=gamma.E, alpha.E=alpha.E, 
               gamma.V=gamma.V, alpha.V=alpha.V,
               gamma.alpha.C=gamma.alpha.C)
  if(only.pars){
    ##return the parameters
    return(list(pars=pars, EX.beta.mu=NULL, EX.beta=NULL, VX.beta=NULL,
                EX.mu=NULL, EX.mu.beta=NULL, EX=NULL, VX=NULL,
                VX.full=NULL, I=NULL))
  }
  ##Compute the beta fields
  if(compute.beta){
    Y <- mesa.data.model$obs$obs
    ##subtract the spatio-temporal component
    if( dimensions$L!=0 )
      Y <- Y - mesa.data.model$SpatioTemp %*% gamma.E
    ##compute iSigma.nu * (Y-M*gamma)
    Y <- block.mult(sigma.nu, Y, n.x=1, n.blocks=dimensions$T,
                    block.sizes=mesa.data.model$nt)
    ##compute F'*iSigma.nu * (Y-M*gamma)
    Y <- calc.tF.times.mat(Y, mesa.data.model$F, 
                           mesa.data.model$obs$idx,
                           n.loc=dimensions$n)
    ##compute iSigma.B * mu.B + F'*iSigma.nu*(Y-M*gamma)
    alpha <- vector("list",dimensions$m)
    offset <- 1
    for(i in 1:length(alpha)){
      alpha[[i]] <- alpha.E[offset:sum(dimensions$p[1:i])]
      offset <- sum(dimensions$p[1:i])+1
    }
    E.beta.mu <- calc.mu.B(dimensions, mesa.data.model$X, alpha)
    Y <- sigma.B %*% c(E.beta.mu) + Y
    ##compute Sigma.B.Y %*% (the above)
    E.beta <- solveTriBlock(sigma.B.Y, Y, n.blocks=1, tr=TRUE)
    E.beta <- solveTriBlock(sigma.B.Y, E.beta, n.blocks=1, tr=FALSE)
    ##reshape result and extract the relevant parts
    E.beta <- matrix(E.beta,ncol=dimensions$m)
    colnames(E.beta.mu) <- colnames(E.beta) <- colnames(mesa.data.model$F)
    rownames(E.beta.mu) <- rownames(E.beta) <- mesa.data.model$location$ID
    ##subset to location in covar.aux
    IND <- match(covar.aux$location$ID,rownames(E.beta))
    E.beta <- E.beta[IND,,drop=FALSE]
    ##compute variance of beta
    if(pred.var){
      V.beta <- invCholBlock(sigma.B.Y, n.blocks=1)
      if(pred.covar){
        V.beta.tmp <- array(NA,c(dimensions$n,dimensions$n,dimensions$m))
        for(i in 1:dimensions$m){
          Ind <- (1:dimensions$n) + (i-1)*dimensions$n
          V.beta.tmp[,,i] <- V.beta[Ind,Ind]
        }
        V.beta <- V.beta.tmp
        dimnames(V.beta) <- list(mesa.data.model$location$ID,
                                 mesa.data.model$location$ID, colnames(E.beta))
        V.beta <- V.beta[IND,IND,,drop=FALSE]
      }else{
        V.beta <- matrix(diag(V.beta),ncol=dimensions$m)
        colnames(V.beta) <- colnames(E.beta)
        rownames(V.beta) <- mesa.data.model$location$ID
        V.beta <- V.beta[IND,,drop=FALSE]
      }
    }else{
      V.beta <- NULL
    }
  }else{
    E.beta.mu <- NULL
    E.beta <- NULL
    V.beta <- NULL
  }##if(compute.beta) else 
    
  ##observations minus the mean value of the field
  C.minus.mu <- mesa.data.model$obs$obs - (Xtilde %*% gamma.alpha)

  ##compute iSigma.nu*(C-mu)
  iS.nu.C <- block.mult(sigma.nu, C.minus.mu, n.blocks=dimensions$T,
                        n.x=1, block.sizes=mesa.data.model$nt)
  ##compute F'*iSigma.nu*(C-mu)
  tF.iS.nu.C <- calc.tF.times.mat(iS.nu.C, mesa.data.model$F, 
                                  mesa.data.model$obs$idx, n.x=1,
                                  n.loc=dimensions$n)
  ##compute inv(sigma.B|Y)*F'*iSigma.nu*(C-mu)
  iSBY.tF.iS.C <- solveTriBlock(sigma.B.Y, tF.iS.nu.C,
                                n.blocks=1, tr=TRUE)
  iSBY.tF.iS.C <- solveTriBlock(sigma.B.Y, iSBY.tF.iS.C,
                                n.blocks=1, tr=FALSE)
  ##compute F'*iSigma.nu
  tF.iS <- calc.tF.times.mat(sigma.nu, mesa.data.model$F, 
                             mesa.data.model$obs$idx,
                             n.loc=dimensions$n)
  ##compute iSigma.tilde*(C-mu)
  ##this is the same for all the conditionals...
  obs.mu <- iS.nu.C - t(tF.iS) %*% iSBY.tF.iS.C

  ##size of the unobserved locations
  N.unobs <- dim(covar.aux$location)[1]
  T.unobs <- dim(covar.aux$trend)[1]

  ##create matrices with site index and temporal trends for all
  if(only.obs){
    ##unobserved points
    idx.unobs <- covar.aux$obs$idx
    ##time of observation for each site
    T1 <- covar.aux$obs$date
    ##trend functions for all points matrix
    F <- covar.aux$F
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- covar.aux$SpatioTemp
    }
  }else{
    ##unobserved points
    idx.unobs <- rep(1:N.unobs,each=T.unobs)
    ##time of observation for each site
    T1 <- rep(covar.aux$trend$date,N.unobs)
    ##trend functions for all points matrix
    F <- matrix(1, (T.unobs*N.unobs), dimensions$m)
    for(i in (1:dimensions$m)){
      ##is not constant: find matching LUR coefficients
      ##(o.w. do nothing -> keep the 1)
      if(colnames(covar.aux$F)[i] != "const")
        F[,i] <- rep(covar.aux$trend[,colnames(covar.aux$F)[i]],N.unobs)
    }##for(i in (1:dimensions$m))
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- matrix(NA, (T.unobs*N.unobs), dimensions$L)
      ST.tmp <- covar.aux$SpatioTemp.all
      ##extract relevant locations and timepoints
      IND <- match(covar.aux$location$ID,colnames(ST.tmp))
      ST.tmp <- ST.tmp[,IND,,drop=FALSE]
      IND <- match(as.character(covar.aux$trend$date),rownames(ST.tmp))
      ST.tmp <- ST.tmp[IND,,,drop=FALSE]
      for(i in (1:dimensions$L))
        ST.unobs[,i] <- c( ST.tmp[,,i] )
    }##if( dimensions$L!=0 )
    if(pred.covar)
      Nmax <- T.unobs
  }##if(only.obs) ... else ...

  if(compute.beta){
    E.beta.list <- list()
    for(i in 1:dim(E.beta)[2]){
      E.beta.list[[i]] <- E.beta[,i,drop=FALSE]
    }
    ##compute beta-fields times the temporal trends
    EX.mu.beta <- calc.F.times.X(E.beta.list, F, idx.unobs)
    ##add the contributions from the different temporal trends
    EX.mu.beta <- matrix(rowSums(EX.mu.beta))
    ##Add the spatio-temporal covariates, if any
    if( dimensions$L!=0 )
      EX.mu.beta <- EX.mu.beta + (ST.unobs %*% gamma.E)

    ##reshape EX.mu.beta to match T-by-N
    if(!only.obs)
      dim(EX.mu.beta) <- c(T.unobs, N.unobs)
  }else{
    EX.mu.beta <- NULL
  }##if(compute.beta) ... else ...
  ##compute LUR times the temporal trends
  Xtilde.unobs <- calc.F.times.X(covar.aux$X, F, idx.unobs)
  
  if( dimensions$L!=0 )
    Xtilde.unobs <- cbind(ST.unobs, Xtilde.unobs)
  ##mean value of the field
  EX.unobs <- Xtilde.unobs %*% gamma.alpha
  ##reshape EX.unobs to match T-by-N
  if(!only.obs)
    dim(EX.unobs) <- c(T.unobs, N.unobs)
  ##lets save the contribution from the mean.
  EX.mu <- EX.unobs
  ##Create a matrix containing pointwise variance
  if(pred.var){
    VX.unobs <- EX.unobs ##right size
    VX.unobs[] <- NA ##set elements to NA
  }else{
    VX.unobs <- NULL
  }
  if(pred.covar){
    VX.full <- list()
  }else{
    VX.full <- NULL
  }
  
  
  ##create distance matrix between unobserved and observed elements
  ##can't find an R function so let's do it by hand
  ##the locations
  loc.unobs <- cbind(covar.aux$location$x,
                     covar.aux$location$y)
  loc.obs <- cbind(mesa.data.model$location$x,
                   mesa.data.model$location$y)
  ##squared distance
  cross.dist <- (rep(loc.unobs[,1],dim(loc.obs)[1]) -
                 rep(loc.obs[,1],each=dim(loc.unobs)[1]) )^2 +
                   (rep(loc.unobs[,2],dim(loc.obs)[1]) -
                    rep(loc.obs[,2],each=dim(loc.unobs)[1]))^2
  ##take the square root and reshape into a matrix of correct size.
  cross.dist <- matrix(sqrt(cross.dist), dim(loc.unobs)[1], dim(loc.obs)[1])
  
  ##We also have that in the stripped data the indecies may not match so we need a
  ##second vector giving which idx in striped matches original idx.
  Ind.2.1 <- match(mesa.data.model$location$ID, covar.aux$location$ID)
  ##determine the nugget for the computations, no.nugget implies smoothing at
  ##observed locations.
  if(no.nugget){
    nugget <- 0
  }else{
    nugget <- phi.residuals[3]
  }
  ##now we need to split the conditional expectation into parts to
  ##reduce memory footprint
  for(i in c(1:ceiling(length(EX.unobs)/Nmax))){
    ##index of the points we want to get conditional expectations for
    Ind <- c((1+(i-1)*Nmax):min(i*Nmax,length(EX.unobs)))
    
    ##create full cross-covariance matrix for all the observations
    sigma.B.full.C <- make.sigma.B.full(trend.sill, trend.range,
                                        loc.ind1=idx.unobs[Ind],
                                        loc.ind2=mesa.data.model$obs$idx,
                                        F1=F[Ind,,drop=FALSE],
                                        F2=mesa.data.model$F,
                                        dists=cross.dist)
    ##parameter order is sill, nugget, range
    sigma.nu.C <- make.sigma.nu.cross.cov(sill=phi.residuals[2],
                                            nugget=nugget,
                                            range=phi.residuals[1],
                                            loc.ind1=idx.unobs[Ind],
                                            loc.ind2=mesa.data.model$obs$idx,
                                            loc.ind2.to.1=Ind.2.1,
                                            T1=T1[Ind],
                                            T2=mesa.data.model$obs$date,
                                            dists=cross.dist)
    sigma.nu.C <- sigma.nu.C + sigma.B.full.C
    ##calculate conditional expectation
    EX.unobs[Ind] <- EX.unobs[Ind] + sigma.nu.C %*% obs.mu
    if(pred.var){
      ##calculate pointwise variance
      ##first the unobserved covariance matrix
      V.uu <- make.sigma.B.full(trend.sill, trend.range,
                                loc.ind1=idx.unobs[Ind],
                                F1=F[Ind,,drop=FALSE],
                                dists=covar.aux$dist)
      V.uu <- V.uu + make.sigma.nu.cross.cov(phi.residuals[2], nugget,
                                             phi.residuals[1],
                                             loc.ind1=idx.unobs[Ind],
                                             T1=T1[Ind],
                                             dists=covar.aux$dist)
      ##compute iSigma.nu*Sigma.ou
      iS.Sou <- block.mult(sigma.nu, t(sigma.nu.C), n.blocks=dimensions$T,
                              block.sizes=mesa.data.model$nt)
      ##compute F'*iSigma.nu*Sigma.ou
      tF.iS.Sou <- tF.iS %*% t(sigma.nu.C)
      ##compute chol(inv(sigma.B.Y))' * (F'*iSigma.nu*Sigma.ou)
      Sby.tF.iS.Sou <- solveTriBlock(sigma.B.Y, tF.iS.Sou, tr=TRUE)
      if(pred.covar){
        ##full matrix
        V.cond <- V.uu - sigma.nu.C %*% iS.Sou +
          t(Sby.tF.iS.Sou) %*% Sby.tF.iS.Sou
      }else{
        ##only diagonal elements
        V.cond <- diag(V.uu) - colSums(t(sigma.nu.C) * iS.Sou) +
          colSums(Sby.tF.iS.Sou * Sby.tF.iS.Sou)
      }
      if(type=="r"){
        stop("NOT YET IMPLEMENTED. Use type=\"p\" instead.")
#        tmp <- Xtilde.unobs[Ind,,drop=FALSE] - t(V.ou) %*% Xtilde
#        if(pred.covar){
#          ##full matrix
#          V.cond <- V.cond + (tmp %*% i.XSX) %*% t(tmp)
#        }else{
#          ##only diagonal elements
#          V.cond <- V.cond + rowSums( (tmp %*% i.XSX) * tmp)
#        }
      }
      ##extract diagonal elements
      if(pred.covar){
        VX.unobs[Ind] <- diag(V.cond)
        VX.full[[i]] <- V.cond
      }else{
        VX.unobs[Ind] <- V.cond
      }
    }##if(pred.var)
  }##for(i in c(1:ceiling(length(EX.unobs)/Nmax)))
  ##ensure positive variances
  if(pred.var)
    VX.unobs <- pmax(VX.unobs,0)
  if(only.obs){
    I <- NULL
  }else{
    ##add names to the EX and VX matrices
    colnames(EX.mu) <- colnames(EX.unobs) <- as.character(covar.aux$location$ID)
    rownames(EX.mu) <- rownames(EX.unobs) <- as.character(covar.aux$trend$date)
    if(compute.beta){
      colnames(EX.mu.beta) <- colnames(EX.unobs)
      rownames(EX.mu.beta) <- rownames(EX.unobs)
    }
    if(pred.var){
      colnames(VX.unobs) <- colnames(EX.unobs)
      rownames(VX.unobs) <- rownames(EX.unobs)
    }
    ##add names to full prediction variances
    if(pred.covar){
      names(VX.full) <- colnames(EX.unobs)
      for(i in 1:length(VX.full))
        colnames(VX.full[[i]]) <- rownames(VX.full[[i]]) <- rownames(EX.unobs)
    }
    ##index for the unobserved values into the predicted values
    if( length(covar.aux$obs$obs)==0 ){
      I <- NULL
    }else{
      I <- (covar.aux$obs$idx-1)*dim(EX.unobs)[1] +
        match(covar.aux$obs$date, covar.aux$trend$date)
    }
  }
  return(list(pars=pars, EX.beta.mu=E.beta.mu, EX.beta=E.beta,
              VX.beta=V.beta, EX.mu=EX.mu, EX.mu.beta=EX.mu.beta,
              EX=EX.unobs, VX=VX.unobs, VX.full=VX.full, I=I))
}##cond.expectation
