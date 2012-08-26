############################################
## S3-METHOD THAT PREDICTS FROM A STmodel ##
############################################
##Functions in this file:
## predict.STmodel         EX:ok
## print.predictSTmode     EX:ok
## plot.predictSTmodel     EX:ok

##' Compute the conditional expectations (i.e. predictions) at the unobserved
##' space-time locations. Predictions are computed for the space-time locations in
##' \code{object} and/or \code{STdata}, conditional on the observations in
##' \code{object} and parameters given in \code{x}.
##' 
##' In addition to computing the conditional expectation at a number of
##' space-time locations the function also computes predictions based on only
##' the regression part of the model as well as the latent beta-fields.
##'
##' Prediction are computed as the conditional expectation of a latent field
##' given observations. This implies that \code{E(X_i| Y_i) != Y_i}, with the
##' difference being due to smoothing over the nugget. Further two possible
##' variance are given, \code{V(X_i|Y_i)} and \code{V(X_i|Y_i)+nugget_i}. Here
##' the nugget for unobserved locations needs to be specified as an additional
##' argument \code{nugget.nobs}. The two variances correspond, losely, to
##' confidence and prediction intervals.
##' 
##' Predictions variances can also be computed. If \code{pred.var=TRUE}
##' point-wise variances for the predictions (and the latent beta-fields) are
##' computed. If instead \code{pred.covar=TRUE} the full covariance matrices for
##' each predicted time series is computed; this implies that the covariances between
##' temporal predictions at the same location are calculated but \emph{not}, due
##' to memory restrictions, any covariances between locations.
##' \code{beta.covar=TRUE} gives the full covariance matrices for the latent
##' beta-fields.
##' 
##' @title Computes Conditional Expectation at Unobserved Locations
##' 
##' @param object \code{STmodel} object for which to compute predictions.
##' @param x Model parameters for which to compute the conditional
##'   expectation. Either as a vector/matrix or an \code{estimateSTmodel} from
##'   \code{\link{estimate.STmodel}}.
##' @param STdata \code{STdata}/\code{STmodel} object with locations/times for
##'   which to predict. If not given predictions are computed for locations/times
##'   in \code{object}
##' @param Nmax Limits the size of matrices constructed when computing
##'   expectations. Use a smaller value if memory becomes a problem.
##' @param only.pars Compute only the regression parameters (using GLS) along
##'   with the related variance.
##' @param nugget.unobs Value of nugget at unonserved locations, either a scalar
##'   or a vector with one element per unobserved site. \strong{NOTE:} All sites in
##'   \code{STdata} are considered unobserved!
##' @param only.obs Compute predictions at only locations specified by
##'   observations in \code{STdata}. Used to limit computations when doing
##'   cross-validation.
##'   \code{only.obs=TRUE} \emph{implies} \code{pred.covar=FALSE} and
##'   \code{combine.data=FALSE}.
##'   Further \code{\link{createSTmodel}} will be called on any \code{STdata}
##'   input, possibly \emph{reordering the observations.}
##' @param pred.var,pred.covar Compute point-wise prediction variances; or
##'   compute covariance matrices for the predicted time series at each location.
##'   \code{pred.covar=TRUE} \emph{implies} \code{pred.var=TRUE} and sets
##'   \code{Nmax} equal to the number of timepoints.
##' @param beta.covar Compute the full covariance matrix for the latent
##'  beta-fields, otherwise only the diagonal elements of V(beta|obs) are
##'  computed. 
##' @param combine.data Combine \code{object} and \code{STdata} and predict for
##'   the joint set of points, see \code{\link{c.STmodel}}.
##' @param type A single character indicating the type of log-likelihood to
##'   compute. Valid options are "f", "p", and "r", for \emph{full},
##'   \emph{profile} or \emph{restricted maximum likelihood} (REML).
##' @param ... Ignored additional arguments.
##' 
##' @return The function returns a list containing (objects not computed
##'   will be missing):
##'   \item{opts}{Copy of options used in the function call.}
##'   \item{pars}{A list with regression parameters and related variances.
##'               \code{pars} contain \code{gamma.E} and \code{alpha.E} with
##'               regression coefficients for the spatio-temporal model and
##'               land-use covaraiates; variances are found in \code{gamma.V}
##'               and \code{alpha.V}; cross-covariance between gamma and alpha in
##'               \code{gamma.alpha.C}.}
##'   \item{beta}{A list with estimates of the beta-fields, including the
##'               regression mean \code{mu}, conditional expectations \code{EX},
##'               possibly variances \code{VX}, and the full covariance matrix
##'               \code{VX.full}.} 
##'   \item{EX.mu}{predictions based on the regression parameters, geographic
##'                covariates, and temporal trends. I.e. only the deterministic
##'                part of the spatio-temporal model.}
##'   \item{EX.mu.beta}{Predictions based on the latent-beta fields, but excluding
##'                     the residual nu field.}
##'   \item{EX}{Full predictions at the space-time locations in
##'             \code{object} and/or \code{STdata}.}
##'   \item{VX}{Pointwise variances for all locations in \code{EX}.}
##'   \item{VX.pred}{Pointwise prediction variances for all locations in
##'                  \code{EX}, i.e. including contribution from
##'                  \code{nugget.unobs}.}
##'   \item{VX.full}{A list with (number of locations) elements, each element is a
##'                  (number of timepoints) - by - (number of timepoints) temporal
##'                  covariance matrix for the timeseries at each location.}
##'   \item{I}{A vector with the locations of the observations in \code{object} or
##'            \code{STdata}. To extract predictions at the observations locations use
##'            \code{EX[I]}.}
##' 
##' @example Rd_examples/Ex_predict_STmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @family predictSTmodel methods
##' @importFrom stats predict
##' @method predict STmodel
##' @export
predict.STmodel <- function(object, x, STdata=NULL, Nmax=1000, only.pars=FALSE,
                            nugget.unobs=0, only.obs=FALSE, pred.var=TRUE,
                            pred.covar=FALSE, beta.covar=FALSE,
                            combine.data=FALSE, type="p", ...){
  ##check class belongings
  stCheckClass(object, "STmodel", name="object")
  if( !is.null(STdata) ){
    stCheckClass(STdata, c("STdata","STmodel"), name="STdata")
  }
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid
  stCheckType(type)

  if( inherits(x,"estimateSTmodel") ){
    x <- coef(x,"all")$par
  }##if( inherits(x,"estimateSTmodel") )
  
  ##figure out a bunch of dimensions
  ##dimensions$L!=0 is used to check for spatio-temporal covariates
  dimensions <- loglikeSTdim(object)
  
  ##check length of input parameters
  if( type=="f" && length(x)!=dimensions$nparam ){
    stop( paste("type=f, requires", dimensions$nparam,
                "parameters but length(x) =", length(x)) )
  }
  if( type!="f" && length(x)==dimensions$nparam ){
    ##drop the regression parameters
    x <- x[(dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam]
  }
  if( type!="f" && length(x)!=dimensions$nparam.cov ){
    stop( paste("type!=f, requires", dimensions$nparam.cov,
                "parameters but length(x) =", length(x)) )
  }
  ##only.pars is strange when type=="f"
  if( only.pars && type=="f"){
    warning("only.pars=TRUE and type=(f)ull only returns KNOWN parameters.",
            immediate. = TRUE)
  }
  ##only.pars, and we can ignore a bunch of things
  if( only.pars ){
    only.obs <- FALSE
    combine.data <- FALSE
    pred.covar <- FALSE
  }
  ##check if only.obs is valid
  if( only.obs && is.null(STdata) ){
    stop("only.obs=TRUE requires STdata.")
  }
  ##do we have an object to combine with
  if( combine.data && is.null(STdata) ){
    warning("No data to combine with; predicting for 'object'", immediate. = TRUE)
    combine.data <- FALSE
  }
  ##can't combine data if we're predicting at only obs.
  if(only.obs && combine.data){
    warning("only.obs=TRUE implies combine.data=FALSE.", immediate. = TRUE)
    combine.data <- FALSE
  }
  ##computing prediction covariates require !only.obs and pred.var
  if(pred.covar && !pred.var){
    warning("pred.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
    pred.var <- TRUE
  }
  if(pred.covar && only.obs){
    warning("only.obs=TRUE implies pred.covar=FALSE.", immediate. = TRUE)
    pred.covar <- FALSE
  }
  ##computing beta covariates require pred.var
  if(beta.covar && !pred.var){
    warning("beta.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
    pred.var <- TRUE
  }
  
  ##create the STdata used for predictions
  if( is.null(STdata) ){
    ##pure copy, predict in the data set
    STdata <- object
  }else{
    if(combine.data){
      ##combine the two datasets, and use for predictions
      STdata <- c(object, STdata)
    }else if( !inherits(STdata,"STmodel") ){
      ##predict only at STdata, and STdata not of class STmodel: need to cast.
      if( is.null(STdata$trend) ){
        warning("STdata lacking a trend element, using trend from object.",
                immediate.=TRUE)
        STdata$trend <- object$trend
      }
      ##Create an STmodel from STdata
      STdata <- createSTmodel(STdata, LUR=object$LUR.list, ST=object$ST.list,
                              cov.beta=object$cov.beta, cov.nu=object$cov.nu,
                              locations=object$locations.list,
                              scale=!is.null(object$scale.covars),
                              scale.covars=object$scale.covars)
    }else{
      ##STdata is an STmodel object ->
      ##test for consistent covariates and scaling (not allowed).
      areSTmodelsConsistent(object, STdata, "STdata")
    }
  }##if( is.null(STdata) ){...}else{...}
  
  ##extract parameters from x
  tmp <- loglikeSTgetPars(x, object)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  cov.pars.beta <- tmp$cov.beta
  cov.pars.nu <- tmp$cov.nu

  ##nugget for unobserved sites
  nugget.unobs <- internalFixNuggetUnobs(nugget.unobs, STdata,
                                         cov.pars.nu$nugget)
  
  ##Create the Xtilde = [M FX] matrix
  Xtilde <- calc.FX(object$F, object$LUR, object$obs$idx)
  ##Add the spatio-temporal covariate (if it exists)
  if( dimensions$L!=0 ){
    Xtilde <- cbind(object$ST, Xtilde)
  }
  ##create covariance matrices, beta-field
  i.sigma.B <- makeSigmaB(cov.pars.beta$pars, dist = object$D.beta,
                          type = object$cov.beta$covf,
                          nugget = cov.pars.beta$nugget)
  ##and nu-field
  i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu,
                          type = object$cov.nu$covf,
                          nugget = cov.pars.nu$nugget,
                          random.effect = cov.pars.nu$random.effect,
                          blocks1 = object$nt, ind1 = object$obs$idx)
  ##calculate block cholesky factor of the matrices, in-place to conserve memory
  i.sigma.B <- makeCholBlock(i.sigma.B, n.blocks=dimensions$m)
  i.sigma.nu <- makeCholBlock(i.sigma.nu, block.sizes=object$nt)
  ##invert the matrices, in-place to conserve memory
  i.sigma.B <- invCholBlock(i.sigma.B, n.blocks=dimensions$m)
  i.sigma.nu <- invCholBlock(i.sigma.nu, block.sizes=object$nt)
  ##F'*inv(sigma.nu)*F
  tF.iS.F <- calc.tFXF(object$F, i.sigma.nu, object$obs$idx,
                       block.sizes=object$nt, n.loc=dimensions$n.obs)
  ##inv(sigma.B|Y) = F'*inv(sigma.nu)*F + inv(sigma.B)
  R.i.sigma.B.Y <- tF.iS.F + i.sigma.B
  ##calculate cholesky factor of inv(sigma.B|Y)
  R.i.sigma.B.Y <- makeCholBlock(R.i.sigma.B.Y)

  ##First lets compute the regression parameters (alpha and gamma)
  if(type=="f"){
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
    iS.nu.X <- blockMult(i.sigma.nu, Xtilde, block.sizes=object$nt)
    ##compute F'*iSigma.nu*Xtilde
    tF.iS.nu.X <- calc.tFX(object$F, iS.nu.X, object$obs$idx,
                           n.loc=dimensions$n.obs)
    ##compute inv(sigma.B|Y)*F'*iSigma.nu*Xtilde
    iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, tF.iS.nu.X, transpose=TRUE)
    iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.X, transpose=FALSE)
    ##compute inv(Xtilde'*Sigma^-1*Xtilde)
    i.XSX <- t(Xtilde) %*% iS.nu.X
    i.XSX <- i.XSX - t(tF.iS.nu.X) %*% iSBY.tF.iS.X
  
    ##compute iSigma.nu*Y
    iS.nu.Y <- blockMult(i.sigma.nu, object$obs$obs, block.sizes=object$nt)
    ##compute F'*iSigma.nu*Y
    tF.iS.nu.Y <- calc.tFX(object$F, iS.nu.Y, object$obs$idx,
                           n.loc=dimensions$n.obs)
    ##compute t(Xtilde'*Sigma^-1*Y)
    tmp2 <- object$obs$obs %*% iS.nu.X
    tmp2 <- tmp2 - t(tF.iS.nu.Y) %*% iSBY.tF.iS.X
    ##compute alpha and gamma
    gamma.alpha <- solve(i.XSX, t(tmp2))

    ##remove variables not needed
    rm(iS.nu.X, tF.iS.nu.X, iSBY.tF.iS.X, iS.nu.Y, tF.iS.nu.Y, tmp2)

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
  }##if(type=="f"){...}else{...}
  
  ##force to matrices
  gamma.E <- as.matrix(gamma.E)
  alpha.E <- as.matrix(alpha.E)
  ##add names
  names.tmp <- loglikeSTnames(object, TRUE)
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
  
  ##construct return object
  out <- list()
  ##set class
  class(out) <- "predictSTmodel"
  ##save options used for predict.
  out$opts <- list(only.pars=only.pars, nugget.unobs=nugget.unobs,
                   only.obs=only.obs, pred.var=pred.var, pred.covar=pred.covar,
                   beta.covar=beta.covar, combine.data=combine.data, type=type)
  ##collect the regression results
  out$pars <- list(gamma.E=gamma.E, alpha.E=alpha.E, 
                   gamma.V=gamma.V, alpha.V=alpha.V,
                   gamma.alpha.C=gamma.alpha.C)

  if( out$opts$only.pars ){
    ##return the parameters
    return( out )
  }
  ##remove variables not needed, regression parameters...
  rm(gamma.E, alpha.E, gamma.V, alpha.V, gamma.alpha.C, names.tmp)
  ##and options stored elsewhere
  rm(only.pars, nugget.unobs, only.obs, pred.var, pred.covar,
     beta.covar, combine.data, type)
    
  ##observations minus the mean value of the field
  C.minus.mu <- object$obs$obs - (Xtilde %*% gamma.alpha)

  ##compute iSigma.nu*(C-mu)
  iS.nu.C <- blockMult(i.sigma.nu, C.minus.mu, block.sizes=object$nt)
  ##compute F'*iSigma.nu*(C-mu)
  tF.iS.nu.C <- calc.tFX(object$F, iS.nu.C, object$obs$idx,
                         n.loc=dimensions$n.obs)
  ##compute inv(sigma.B|Y)*F'*iSigma.nu*(C-mu)
  iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, tF.iS.nu.C, transpose=TRUE)
  iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.C, transpose=FALSE)
  ##compute F'*i.sigma.nu
  tF.iS <- calc.tFX(object$F, i.sigma.nu, object$obs$idx,
                    n.loc=dimensions$n.obs)
  ##compute iSigma.tilde*(C-mu)
  ##this is the same for all the conditionals...
  obs.mu <- iS.nu.C - t(tF.iS) %*% iSBY.tF.iS.C

  ##remove variables not needed
  rm(C.minus.mu, iS.nu.C, tF.iS.nu.C, iSBY.tF.iS.C)

  ##create matrices with site index and temporal trends for all
  ##OBS date has to be INCREASING, i.e. same sorting as for STdata$obs...
  if( out$opts$only.obs ){
    ##unobserved points
    idx.unobs <- STdata$obs$idx
    ##time of observation for each site
    T1 <- STdata$obs$date
    ##trend functions for all points matrix
    F <- STdata$F
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- STdata$ST
    }
    ##size of the unobserved locations
    N.unobs <- length( unique(idx.unobs) )
    T.unobs <- length( unique(T1) )
  }else{
    ##size of the unobserved locations
    N.unobs <- dim(STdata$locations)[1]
    T.unobs <- dim(STdata$trend)[1]
  
    ##unobserved points
    idx.unobs <- rep(1:N.unobs, each=T.unobs)
    ##time of observation for each site
    T1 <- rep(STdata$trend$date, N.unobs)
    ##trend functions for all points matrix
    F <- matrix(1, (T.unobs*N.unobs), dimensions$m)
    for(i in (1:dimensions$m)){
      ##if not constant, find relevant row in $trend
      if(colnames(STdata$F)[i] != "const"){
        F[,i] <- rep(STdata$trend[,colnames(STdata$F)[i]], N.unobs)
      }
    }##for(i in (1:dimensions$m))
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- matrix(STdata$ST.all, (T.unobs*N.unobs), dimensions$L)
    }##if( dimensions$L!=0 )
    if(out$opts$pred.covar){
      Nmax <- T.unobs
    }
  }##if( out$opts$only.obs ) ... else ...
  ##compute number of observations for each time-point,
  ##taking into account that there might be missing dates.
  date.all <- sort(unique( c(object$trend$date, STdata$trend$date) ))
  nt.unobs <- nt.obs <- double(length(date.all))
  for(i in c(1:length(date.all))){
    nt.obs[i] <- sum( object$obs$date==date.all[i] )
  }

  ##compute LUR times the temporal trends [M F*X] for unobserved
  Xtilde.unobs <- calc.FX(F, STdata$LUR.all, idx.unobs)
  if( dimensions$L!=0 ){
    Xtilde.unobs <- cbind(ST.unobs, Xtilde.unobs)
  }
  ##mean value of the field
  out$EX <- as.matrix( Xtilde.unobs %*% gamma.alpha )
  ##reshape out$EX to match T-by-N
  if( !out$opts$only.obs ){
    dim(out$EX) <- c(T.unobs, N.unobs)
  }
  ##lets save the contribution from the mean.
  out$EX.mu <- out$EX
  ##Create a matrix containing pointwise variance
  if(out$opts$pred.var){
    out$VX <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
    out$VX.pred <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
  }
  ##and a list of covariances.
  if(out$opts$pred.covar){
    out$VX.full <- list()
  }
  
  ##compute LUR*alpha (mean values for the beta fields)
  alpha <- vector("list", dimensions$m)
  offset <- 1
  for(i in 1:length(alpha)){
    alpha[[i]] <- out$pars$alpha.E[offset:sum(dimensions$p[1:i])]
    offset <- sum(dimensions$p[1:i])+1
  }
  ##we need the dimension of the prediction object.
  out$beta <- list()
  out$beta$mu <- calc.mu.B(STdata$LUR.all, alpha)
  ##remove variables not needed
  rm(alpha, gamma.alpha)
  
  ##create distance matrix between unobserved and observed elements
  ##unobserved locations
  loc.unobs.nu <- STdata$locations[, c("x.nu","y.nu"), drop=FALSE]
  loc.unobs.beta <- STdata$locations[, c("x.beta","y.beta"), drop=FALSE]

  ##observed locations
  I.obs <- match( colnames(object$D.nu), object$locations$ID)
  loc.obs.nu <- object$locations[I.obs, c("x.nu","y.nu"), drop=FALSE]
  loc.obs.beta <- object$locations[I.obs, c("x.beta","y.beta"), drop=FALSE]
  
  ##We also have that in the stripped data the indecies may not match so we need a
  ##second vector giving which idx in striped matches original idx.
  Ind.2.1 <- match(object$locations$ID, STdata$locations$ID, nomatch=0)
    
  ##precompute the cross-covariance for the beta-fields
  sigma.B.C <- makeSigmaB(cov.pars.beta$pars,
                          dist = crossDist(loc.unobs.beta, loc.obs.beta),
                          type = object$cov.beta$covf,
                          nugget = cov.pars.beta$nugget,
                          ind2.to.1=Ind.2.1)
  
  ##compute beta fields - right most part is F'*iSigma.tilde*(C-mu)
  out$beta$EX <- (c(out$beta$mu) + (sigma.B.C %*%
                                    calc.tFX(object$F, obs.mu,
                                             object$obs$idx, 
                                             n.loc=dimensions$n.obs) ))
  dim(out$beta$EX) <- dim(out$beta$mu)
  dimnames(out$beta$EX) <- dimnames(out$beta$mu)
  ##and variances
  if( out$opts$pred.var ){
    Sby.iSb.Sou <- i.sigma.B %*% t(sigma.B.C)
    Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, transpose=TRUE)
    Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, transpose=FALSE)

    if( out$opts$beta.covar ){
      ##sigma.B for unobserved locations
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars,
                               dist = crossDist(loc.unobs.beta),
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget)
      ##compute variance matrix
      tmp <- sigma.B.uu - sigma.B.C %*% tF.iS.F %*% Sby.iSb.Sou
      out$beta$VX.full <- list()
      for(i in 1:dim(out$beta$EX)[2]){
        Ind <- (1:dim(out$beta$EX)[1]) + (i-1)*dim(out$beta$EX)[1]
        out$beta$VX.full[[i]] <- tmp[Ind, Ind, drop=FALSE]
        rownames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
        colnames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
      }
      names(out$beta$VX.full) <- colnames(out$beta$EX)
      tmp <- diag(tmp)
    }else{
      ##compute only diagonal elements of the STATIONARY sigma.beta matrix
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = matrix(0,1,1),
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget)
      ##expand to full matrix
      sigma.B.uu <- matrix(diag(sigma.B.uu), ncol=dim(sigma.B.uu)[1],
                           nrow=dim(loc.unobs.beta)[1], byrow=TRUE)
      ##and compute the relevant covariance
      tmp <- c(sigma.B.uu) - rowSums(sigma.B.C * t(tF.iS.F %*% Sby.iSb.Sou))
    }
    
    out$beta$VX <- matrix(tmp, ncol=dim(out$beta$EX)[2])
    dimnames(out$beta$VX) <- dimnames(out$beta$EX)
    
    ##remove variables not needed
    rm(tmp, Sby.iSb.Sou, tF.iS.F, sigma.B.uu)
  }
  ##remove variables not needed
  rm(i.sigma.B)
  
  ##cross distance for the nu coordinates
  cross.D.nu <- crossDist(loc.unobs.nu, loc.obs.nu)
  
  ##now we need to split the conditional expectation into parts to
  ##reduce memory footprint
  for( i in 1:ceiling(length(out$EX)/Nmax) ){
    ##index of the points we want to get conditional expectations for
    Ind <- c((1+(i-1)*Nmax):min(i*Nmax,length(out$EX)))
    ##compute number of unobserved per block for this subset
    T1.Ind <- T1[Ind]
    for(j in c(1:length(date.all))){
      nt.unobs[j] <- sum( T1.Ind==date.all[j] )
    }
    ##order T1.Ind for computation of cross.covariance
    T1.order <- order(T1.Ind)
    
    ##create full cross-covariance matrix for all the observations
    sigma.B.full.C <- calc.FXtF2(F[Ind,,drop=FALSE], sigma.B.C,
                                 loc.ind=idx.unobs[Ind], F2=object$F,
                                 loc.ind2=object$obs$idx)
    ##parameter order is sill, nugget, range (always predict with nugget=0)
    sigma.nu.C <- makeSigmaNu(cov.pars.nu$pars, dist = cross.D.nu,
                              type = object$cov.nu$covf, nugget = 0,
                              random.effect = cov.pars.nu$random.effect,
                              ind1 = (idx.unobs[Ind])[T1.order],
                              ind2 = object$obs$idx,
                              blocks1 = nt.unobs, blocks2 = nt.obs,
                              ind2.to.1=Ind.2.1)
    sigma.nu.C[T1.order,] <- sigma.nu.C
    sigma.nu.C <- sigma.nu.C + sigma.B.full.C
    ##calculate conditional expectation
    out$EX[Ind] <- out$EX.mu[Ind] + sigma.nu.C %*% obs.mu
    
    ##calculate pointwise variance
    if(out$opts$pred.var){
###OLD VERSION
#      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars,
#                               dist = crossDist(loc.unobs.beta),
#                               type = object$cov.beta$covf,
#                               nugget = cov.pars.beta$nugget)
      ##first the unobserved covariance matrix
#      V.uu <- calc.FXtF2(F[Ind,,drop=FALSE], sigma.B.uu,
#                         loc.ind=idx.unobs[Ind])
#      unobs.D.nu <- crossDist(loc.unobs.nu)
#      V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu,
#                                 type = object$cov.nu$covf, nugget = 0,
#                                 random.effect = cov.pars.nu$random.effect,
#                                 ind1 = idx.unobs[Ind], blocks1 = nt.unobs)
###NEW VERSION
      ##pick out the observed locations in this iteration
      I.loc <- sort(unique( idx.unobs[Ind] ))
      ##compute relevant part of the covariance matrix
      unobs.D.beta <- crossDist(loc.unobs.beta[I.loc,,drop=FALSE])
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = unobs.D.beta,
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget)
      ##first the unobserved covariance matrix
      V.uu <- calc.FXtF2(F[Ind,,drop=FALSE], sigma.B.uu,
                         loc.ind=idx.unobs[Ind]-min(I.loc)+1)
      ##distance matrices for the unobserved nu locations
      unobs.D.nu <- crossDist(loc.unobs.nu[I.loc,,drop=FALSE])
      V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu,
                                 type = object$cov.nu$covf, nugget = 0,
                                 random.effect = cov.pars.nu$random.effect,
                                 ind1 = idx.unobs[Ind]-min(I.loc)+1,
                                 blocks1 = nt.unobs)
      
      ##transpose of corss-covariance matrix
      t.sigma.nu.C <- t(sigma.nu.C)
      ##compute iSigma.nu*Sigma.ou
      iS.Sou <- blockMult(i.sigma.nu, t.sigma.nu.C, block.sizes=object$nt)
      ##compute F'*iSigma.nu*Sigma.ou
      tF.iS.Sou <- tF.iS %*% t.sigma.nu.C
      ##compute chol(inv(sigma.B.Y))' * (F'*iSigma.nu*Sigma.ou)
      Sby.tF.iS.Sou <- solveTriBlock(R.i.sigma.B.Y, tF.iS.Sou, transpose=TRUE)
      if( out$opts$pred.covar ){
        ##full matrix
        tmp <- - sigma.nu.C %*% iS.Sou + t(Sby.tF.iS.Sou) %*% Sby.tF.iS.Sou
        V.cond <- V.cond.0 <- V.uu + tmp
        ##add nugget to the diagonal
        diag(V.cond) <- diag(V.cond) + out$opts$nugget.unobs[idx.unobs[Ind]]
      }else{
        ##only diagonal elements
        tmp <- (- colSums(t(sigma.nu.C) * iS.Sou) +
                colSums(Sby.tF.iS.Sou * Sby.tF.iS.Sou))
        V.cond <- V.cond.0 <- diag(V.uu) + tmp
        V.cond <- V.cond + out$opts$nugget.unobs[idx.unobs[Ind]]
      }
      if( out$opts$type=="r" ){
        stop("NOT YET IMPLEMENTED. Use type=\"p\" instead.")
#        tmp <- Xtilde.unobs[Ind,,drop=FALSE] - t(V.ou) %*% Xtilde
#        if(out$opts$pred.covar){
#          ##full matrix
#          V.cond <- V.cond + (tmp %*% i.XSX) %*% t(tmp)
#        }else{
#          ##only diagonal elements
#          V.cond <- V.cond + rowSums( (tmp %*% i.XSX) * tmp)
#        }
      }
      ##extract diagonal elements
      if(out$opts$pred.covar){
        out$VX[Ind] <- diag(V.cond.0)
        out$VX.pred[Ind] <- diag(V.cond)
        out$VX.full[[i]] <- V.cond.0
      }else{
        out$VX[Ind] <- V.cond.0
        out$VX.pred[Ind] <- V.cond
      }
    }##if(out$opts$pred.var)
  }##for(i in c(1:ceiling(length(out$EX)/Nmax)))
  ##ensure positive variances
  if(out$opts$pred.var){
    out$VX <- pmax(out$VX, 0)
    out$VX.pred <- pmax(out$VX.pred, 0)
  }

  ##compute predicitons based on the beta fields
  EX.beta.list <- list()
  for(i in 1:dim(out$beta$EX)[2]){
    EX.beta.list[[i]] <- out$beta$EX[,i,drop=FALSE]
  }
  ##compute beta-fields times the temporal trends
  out$EX.mu.beta <- calc.FX(F, EX.beta.list, idx.unobs)
  ##add the contributions from the different temporal trends
  out$EX.mu.beta <- matrix( rowSums(out$EX.mu.beta) )
  ##Add the spatio-temporal covariates, if any
  if( dimensions$L!=0 ){
    out$EX.mu.beta <- out$EX.mu.beta + (ST.unobs %*% out$pars$gamma.E)
  }
  ##reshape EX.mu.beta to match T-by-N
  if(!out$opts$only.obs){
    dim(out$EX.mu.beta) <- c(T.unobs, N.unobs)
  }
  
  if( !out$opts$only.obs ){
    ##add names to the EX and VX matrices
    colnames(out$EX) <- STdata$locations$ID
    rownames(out$EX) <- as.character(STdata$trend$date)
    dimnames(out$EX.mu.beta) <- dimnames(out$EX.mu) <- dimnames(out$EX)
    if(out$opts$pred.var){
      dimnames(out$VX) <- dimnames(out$VX.pred) <- dimnames(out$EX)
    }
    ##add names to full prediction variances
    if(out$opts$pred.covar){
      names(out$VX.full) <- colnames(out$EX)
      for(i in 1:length(out$VX.full))
        colnames(out$VX.full[[i]]) <- rownames(out$VX.full[[i]]) <- rownames(out$EX)
    }
  }##if( !out$opts$only.obs )
  ##index for the unobserved values into the predicted values
  if( length(STdata$obs$obs)!=0 ){
    if( out$opts$only.obs ){
      I <- 1:length(out$EX)
    }else{
      I <- (match(STdata$obs$ID,colnames(out$EX))-1)*dim(out$EX)[1] +
        match(STdata$obs$date, STdata$trend$date)
    }
    out$I <- data.frame(I=I, date=STdata$obs$date, ID=STdata$obs$ID,
                        stringsAsFactors=FALSE)
  }##if( length(STdata$obs$obs)!=0 )

  return( out )
}##function predict.STmodel

###################################
## S3 methods for predictSTmodel ##
###################################
##' \code{\link[base:print]{print}} method for class \code{predictSTmodel}.
##'
##' @title Print details for \code{predictSTmodel} object
##' @param x \code{predictSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##'
##' @examples
##'   ##load data
##'   data(pred.mesa.model)
##'   print(pred.mesa.model)
##'   
##' 
##' @family predictSTmodel methods
##' @method print predictSTmodel
##' @export
print.predictSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "predictSTmodel", name="x")

  cat("Prediction for STmodel.\n\n")
  if( x$opts$only.pars ){
    cat("Only computed regression parameters:\n")
  }else{
    cat("Regression parameters:\n")
  }
  if( length(x$pars$gamma.E)==0 ){
    cat("\t","No spatio-temporal covariate.\n")
  }else{
    cat("\t", length(x$pars$gamma), "Spatio-temporal covariate(s).\n")
  }
  cat("\t", length(x$pars$alpha.E),
      "beta-fields regression parameters in x$pars.\n\n")
  if( x$opts$only.pars ){
    return(invisible())
  }

  cat("Prediction of beta-fields, (x$beta):\n")
  str(x$beta,1)
  cat("\n")
  
  if( x$opts$only.obs ){
    cat("Predictions only for", length(x$EX), "observations.\n")
  }else{
    cat("Predictions for", dim(x$EX)[1], "times at",
        dim(x$EX)[2], "locations.\n")
  }
  str(x[c("EX.mu","EX.mu.beta","EX")],1)
  cat("\n")
  
  if( x$opts$pred.covar ){
    cat("Variances and temporal covariances for each location\n")
    cat("\thave been computed.\n")
    str(x[c("VX","VX.pred","VX.full")],1)
  }else if( x$opts$pred.var ){
    cat("Variances have been computed.\n")
    str(x[c("VX","VX.pred")],1)
  }else{
    cat("Variances have NOT been computed.\n")
  }
  cat("\n")
  
  return(invisible())
}##function print.predictSTmodel


##' \code{\link[graphics:plot]{plot}} method for classes \code{predictSTmodel}
##' and \code{predCVSTmodel}. Provides several different plots of the
##' data.  
##'
##' @title Plots for \code{predictSTmodel} and \code{predCVSTmodel} Objects
##' 
##' @param x \code{predictSTmodel} or \code{predCVSTmodel} object to plot.
##' @param y Plot predictions as a function of either \code{"time"} or
##'   \code{"obs"}ervations. 
##' @param STmodel \code{STdata}/\code{STmodel} object containing observations
##'   with which to compare the predictions (not used for
##'   \code{plot.predCVSTmodel}). 
##' @param ID The location for which we want to plot predictions. A
##'   string matching names in \code{colnames(x$EX)} (or \code{x$I$ID},
##'   number(s) which are used as \code{ID = colnames(x$EX)[ID]}, or
##'   \code{"all"} in which case all predictions are used.
##'   If several locations are given (or \code{"all"}) then
##'   \code{y} must be \code{"obs"}.
##' @param col A vector of three colours: The first is the colour of the
##'   predictions, second for the observations and third for the polygon
##'   illustrating the confidence bands.
##' @param pch,cex,lty,lwd Vectors with two elements giving the point type,
##'   size, line type and line width to use when plotting the predictions and
##'   observations respectively. Setting a value to \code{NA} will give no
##'   points/lines for the predictions/observations. \cr
##'   When plotting predictions
##'   as a function of observations \code{lty[2]} is used for the addition of
##'   \code{\link[graphics:abline]{abline}(0,1, lty=lty[2], col=col[2],
##'   lwd=lwd[2])}; \code{pch[2]} and \code{cex[2]} are ignored. 
##' @param p Width of the plotted confidence bands (as coverage percentage,
##'   used to find appropriate two-sided normal quantiles).
##' @param pred.type Which type of prediction to plot, one of
##'   \code{"EX"}, \code{"EX.mu"}, or \code{"EX.mu.beta"}, see the 
##'   output from \code{\link{predict.STmodel}}
##' @param pred.var Should we plot confidence bands based on prediction (TRUE)
##'   or confidence intrevalls (FALSE), see \code{\link{predict.STmodel}}.
##'   Only relevant if \code{pred.type="EX"}. \cr
##'   \strong{NOTE:} \emph{The default differs for \code{plot.predictSTmodel}
##'   and \code{plot.predCVSTmodel}!} 
##' @param add Add to existing plot?
##' @param ... Ignored additional arguments.
##' 
##' @return Nothing
##'
##' @example Rd_examples/Ex_plot_predictSTmodel.R
##'
##' @author Johan Lindström
##' 
##' @family predictSTmodel methods
##' @method plot predictSTmodel
##' @export
plot.predictSTmodel <- function(x, y="time", STmodel=NULL, ID=x$I$ID[1],
                                col=c("black","red","grey"), pch=c(NA,NA),
                                cex=c(1,1), lty=c(1,1), lwd=c(1,1), p=0.95,
                                pred.type="EX", pred.var=FALSE,
                                add=FALSE, ...){
  ##check class belonging
  stCheckClass(x, "predictSTmodel", name="x")
  if( !is.null(STmodel) ){
    stCheckClass(STmodel, c("STmodel","STdata"), name="STmodel")
  }
  ##check for observations
  if( x$opts$only.pars ){
    message("Prediction structure contains only parameters.")
    ##return
    return(invisible())
  }
  ##we have to use y, cast to resonable name
  plot.type <- y
  pred.var <- internalPlotPredictChecks(plot.type, pred.type, pred.var)
  ##check ID
  ID <- internalFindIDplot(ID, colnames(x$EX))
  ID.all <- internalCheckIDplot(ID, y)

  ##pick out data for plotting
  if( !ID.all ){
    ##pick out predictions, two cases
    if( x$opts$only.obs ){
      I1 <- x$I$ID==ID
      I2 <- 1
      date <- x$I$date[I1]
    }else{
      I1 <- 1:dim(x[[pred.type]])[1]
      I2 <- ID
      date <- rownames(x[[pred.type]])
    }
    sd <- x[[pred.var]][I1,I2]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x[[pred.type]][I1, I2], sd=sd, date=date,
                       stringsAsFactors=FALSE)
    ##pick out observations, if any
    obs <- STmodel$obs[STmodel$obs$ID==ID,,drop=FALSE]
    if( !is.null(obs) ){
      tmp <- obs
      obs <- matrix(NA, dim(pred)[1], 1)
      rownames(obs) <- as.character(pred$date)
      obs[as.character(tmp$date),] <- tmp$obs
    }
  }else{
    ##all only relevant when we have plot by observations
    ##pick out predictions, two cases
    if(length(ID)==1 && ID=="all"){
      ID <- unique(x$I$ID)
    }
    I.ID <- x$I$ID %in% ID
    I <- x$I$I[I.ID]
    sd <- x[[pred.var]][I]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x[[pred.type]][I], sd=sd, date=x$I$date[I.ID],
                       stringsAsFactors=FALSE)
    ##pick out observations, if any
    obs <- createDataMatrix(STmodel)
    I <- (match(x$I$ID[I.ID],colnames(obs))-1)*dim(obs)[1] +
        match(as.character(x$I$date[I.ID]), rownames(obs))
    obs <- obs[I]
    if( dim(pred)[1]!=length(obs) ){
      stop("Number of observed locations do not match no. predicted.")
    }
  }##if( !ID.all ){...}else{...}

  internalPlotPredictions(plot.type, ID, pred, obs, col, pch, cex, lty, lwd, p, add)
  
  ##return
  return(invisible())
}##function plot.predictSTmodel
