###################################
## Functions for crossvalidation ##
###################################
##Functions in this file:
## predictCV.STmodel            EX:with estimateCV.STmodel
## predictCV                    EX:S3 method
## print.predCVSTmodel          EX:ok
## summary.predCVSTmodel        EX:ok
## print.summary.predCVSTmodel  EX:with summary.predCVSTmodel
## plot.predCVSTmodel           EX:missing

##' @param silent Show status after each iteration?
##' @rdname estimateCV.STmodel
##' @family predCVSTmodel functions
##' @family predCVSTmodel methods
##' @method predictCV STmodel
##' @export
predictCV.STmodel <- function(object, x, Ind.cv, ..., silent=TRUE){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  ##check cross-validation groups
  Ind.cv <- stCheckInternalCV(Ind.cv)
  
  ##ensure that Ind.cv is a matrix
  Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    stop("Some observation(s) are left out in several Cv-groups.")
  }

  if( inherits(x, "estCVSTmodel") ){
    x <- coef(x, "all")
  }else if( inherits(x, "estimateSTmodel") ){
    x <- coef(x, "all")$par
  }##if( inherits(x,"estimateSTmodel") ){...}else if( inherits(x, "estimateSTmodel") ){...}
                         
  ##expand the parameter vector so that we have one
  ##vector for each of the cv-sets we want to try.
  x <- as.matrix(x)
  if(dim(x)[2]==1){
    x <- matrix(x, length(x), N.CV.sets)
  }else if(dim(x)[2] != N.CV.sets){
    stop("Number of parameters does not match the number of cv-sets.")
  }

  ##loop over all locations and predict
  pred <- list()
  for(i in 1:N.CV.sets){
    if(!silent)
      message( sprintf("Predicting cv-set %d/%d", i, N.CV.sets) )
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- as.logical( Ind.cv[,i] )
    }
    ##create data matrices that contains observed
    object.obs <- dropObservations(object, Ind.current)
    ##and locations to be predicted
    suppressWarnings( object.pred <- dropObservations(object, !Ind.current) )

    ##compute nugget for unobserved
    nugget.unobs <- loglikeSTgetPars(x[,i], object)$cov.nu$nugget
    nugget.unobs <- nugget.unobs[object.pred$locations$ID, , drop=FALSE]
    ##lets obtain the predection for this set
    pred[[i]] <- predict(object.obs, x[,i], STdata=object.pred,
                         nugget.unobs=nugget.unobs, only.pars=FALSE,
                         combine.data=FALSE, ...)
  }##for(i in 1:N.CV.sets){

  ##list containing results
  out <- list()
  #pred.obs=pred.obs, pred.by.cv=pred.by.cv, pred.all=pred.all.res)
  class(out) <- "predCVSTmodel"
  out$opts <- pred[[1]]$opts
  if( !is.null(out$opts$nugget.unobs) ){
    out$opts$nugget.unobs <- unlist(sapply(pred, function(x){
      x$opts$nugget.unobs}))
    names(out$opts$nugget.unobs) <- unlist(sapply(pred,function(x){
      rownames(x$opts$nugget.unobs)}))
  }
  out$Ind.cv <- Ind.cv

  ##collect results.
  ##construct a matrix matching mesa.data$obs$obs with the obs
  out$pred.obs <- object$obs[,c("obs","date","ID"), drop=FALSE]
  ##add some more columns
  out$pred.obs$EX.mu <- NA
  out$pred.obs$EX.mu.beta <- NA
  out$pred.obs$EX <- NA
  ##do we have prediction variances?
  if( out$opts$pred.var ){
    out$pred.obs$VX <- NA
    out$pred.obs$VX.pred <- NA
  }
  ##loop over the CV-sets
  for(i in 1:N.CV.sets){
    Ind.current <- Ind.cv==i
    I <- pred[[i]]$I$I
    out$pred.obs[Ind.current, c("EX.mu","EX.mu.beta","EX")] <-
      cbind( pred[[i]]$EX.mu[I], pred[[i]]$EX.mu.beta[I], pred[[i]]$EX[I])
    if( out$opts$pred.var ){
      out$pred.obs[Ind.current, c("VX","VX.pred")] <-
        cbind( pred[[i]]$VX[I], pred[[i]]$VX.pred[I])
    }
  }
  ##residuals
  out$pred.obs$res <- out$pred.obs$obs - out$pred.obs$EX
  ##normalised residuals
  if( out$opts$pred.var ){
    out$pred.obs$res.norm <- out$pred.obs$res / sqrt(out$pred.obs$VX.pred)
  }
  
  ##determine if the same location is in several cross-validation groups
  if( out$opts$only.obs ){
    any.duplicates <- any( duplicated( unlist( sapply(pred, function(x){
      unique(x$I$ID)}) ) ) )
  }else{
    any.duplicates <- any( duplicated( unlist( sapply(pred, function(x){
      colnames(x$EX)}) ) ) )
  }

  ##Extract by location predictions
  if( !any.duplicates ){
    out$pred.all <- list()

    ##matrices that contain predictions at all times and locations
    out$pred.all$EX.mu <- matrix(NA, dim(object$trend)[1],
                                 dim(object$locations)[1])
    colnames(out$pred.all$EX.mu) <- object$locations$ID
    rownames(out$pred.all$EX.mu) <- as.character(object$trend$date)
    out$pred.all$EX <- out$pred.all$EX.mu.beta <- out$pred.all$EX.mu
    if( out$opts$pred.var ){
      out$pred.all$VX <- out$pred.all$VX.pred <- out$pred.all$EX
    }

    ##matrices that contain beta predictions
    out$pred.all$beta <- list()
    out$pred.all$beta$mu <- matrix(NA, dim(object$locations)[1],
                                length(object$LUR))
    colnames(out$pred.all$beta$mu) <- names(object$LUR)
    rownames(out$pred.all$beta$mu) <- object$locations$ID
    out$pred.all$beta$EX <- out$pred.all$beta$mu
    if( out$opts$pred.var ){
      out$pred.all$beta$VX <- out$pred.all$beta$EX
    }

    for(i in 1:N.CV.sets){
      ##sites predicted in this CV-group
      if( out$opts$only.obs ){
        ID.names <- rownames(pred[[i]]$beta$EX)
        I <- ((match( pred[[i]]$I$ID, object$locations$ID)-1) *
              dim(out$pred.all$EX)[1] +
              match( as.character(pred[[i]]$I$date),
                    rownames(out$pred.all$EX)) )
      }else{
        ID.names <- colnames(pred[[i]]$EX)
        I <- (rep(match(ID.names, object$locations$ID)-1,
                  each=dim(out$pred.all$EX)[1]) * dim(out$pred.all$EX)[1] +
              rep(match(rownames(pred[[i]]$EX), rownames(out$pred.all$EX)),
                  length(ID.names)) )
      }
      ##extract predictions
      out$pred.all$EX.mu[I] <- pred[[i]]$EX.mu
      out$pred.all$EX.mu.beta[I] <- pred[[i]]$EX.mu.beta
      out$pred.all$EX[I] <- pred[[i]]$EX
      ##...variances
      if( out$opts$pred.var ){
        out$pred.all$VX[I] <- pred[[i]]$VX
        out$pred.all$VX.pred[I] <- pred[[i]]$VX.pred
      }
      ##...and beta-fields
      out$pred.all$beta$mu[ID.names,] <- pred[[i]]$beta$mu
      out$pred.all$beta$EX[ID.names,] <- pred[[i]]$beta$EX
      if( out$opts$pred.var ){
        out$pred.all$beta$VX[ID.names,] <- pred[[i]]$beta$VX
      }
    }##for(i in 1:N.CV.sets)
  }else{
    out$pred.all.by.cv <- vector("list", N.CV.sets)
    for(i in 1:N.CV.sets){
      if( out$opts$only.obs ){
        out$pred.all.by.cv[[i]] <- list()
        ##create data matrix for the different parts
        EX.mu <- createDataMatrix(obs=pred[[i]]$EX.mu, date=pred[[i]]$I$date,
                                ID=pred[[i]]$I$ID)
        EX.mu.beta <- createDataMatrix(obs=pred[[i]]$EX.mu.beta,
                                       date=pred[[i]]$I$date, 
                                       ID=pred[[i]]$I$ID)
        EX <- createDataMatrix(obs=pred[[i]]$EX, date=pred[[i]]$I$date,
                               ID=pred[[i]]$I$ID)
        if( out$opts$pred.var ){
          VX <- createDataMatrix(obs=pred[[i]]$VX, date=pred[[i]]$I$date,
                                 ID=pred[[i]]$I$ID)
          VX.pred <- createDataMatrix(obs=pred[[i]]$VX.pred,
                                      date=pred[[i]]$I$date,
                                      ID=pred[[i]]$I$ID)
        }
        
        ##create target matrices with ALL dates
        out$pred.all.by.cv[[i]]$EX.mu <- matrix(NA, dim(object$trend)[1],
                                                dim(EX.mu)[2])
        colnames(out$pred.all.by.cv[[i]]$EX.mu) <- colnames(EX.mu)
        rownames(out$pred.all.by.cv[[i]]$EX.mu) <- as.character(object$trend$date)
        out$pred.all.by.cv[[i]]$EX.mu.beta <- out$pred.all.by.cv[[i]]$EX.mu
        out$pred.all.by.cv[[i]]$EX <- out$pred.all.by.cv[[i]]$EX.mu
        if( out$opts$pred.var ){
          out$pred.all.by.cv[[i]]$VX <- out$pred.all.by.cv[[i]]$EX
          out$pred.all.by.cv[[i]]$VX.pred <- out$pred.all.by.cv[[i]]$EX
        }
        
        ##fit predictions into matrices
        out$pred.all.by.cv[[i]]$EX.mu[rownames(EX.mu),] <- EX.mu
        out$pred.all.by.cv[[i]]$EX.mu.beta[rownames(EX.mu.beta),] <- EX.mu.beta
        out$pred.all.by.cv[[i]]$EX[rownames(EX),] <- EX
        if( out$opts$pred.var ){
          out$pred.all.by.cv[[i]]$VX[rownames(VX),] <- VX
          out$pred.all.by.cv[[i]]$VX.pred[rownames(VX.pred),] <- VX.pred
        }
        
        ##just copy beta predictions
        out$pred.all.by.cv[[i]]$beta <- pred[[i]]$beta
      }else{
        ##just copy
        if( out$opts$pred.var ){
          pick.names <- c("EX.mu", "EX.mu.beta", "EX", "VX", "VX.pred",
                          "beta")
        }else{
          pick.names <- c("EX.mu", "EX.mu.beta", "EX", "beta")
        }
        out$pred.all.by.cv[[i]] <- pred[[i]][pick.names]
      }##if( out$opts$only.obs ){...}else{...}
    }##for(i in 1:N.CV.sets)
  }##if( !any(duplicated(...)) ){...}else{...}

  ##return results
  return( out )
}##function predictCV.STmodel

######################################
## General S3 methods for predictCV ##
######################################

##' @rdname estimateCV
##' @export
predictCV <- function(object, x, Ind.cv, ...){
  UseMethod("predictCV")
}

##################################
## S3-METHODS FOR predCVSTmodel ##
##################################

##' \code{\link[base:print]{print}} method for class \code{predCVSTmodel}.
##'
##' @title Print details for \code{predCVSTmodel} object
##' @param x \code{predCVSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @examples
##' ##load some data
##' data(CV.mesa.model)
##' ##print basic information for the CV-predictions
##' print(pred.cv.mesa)
##' 
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method print predCVSTmodel
##' @export
print.predCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "predCVSTmodel", name="x")

  N.cv <- max(x$Ind.cv)
  by.CV <- !is.null(x$pred.all.by.cv)

  N.loc <- dim(x$pred.all$EX)
  
  cat("Cross-validation prediction for STmodel with", N.cv, "CV-groups.\n")
  if( x$opts$only.obs ){
    cat("  Predictions for observations only.\n")
  }else{
    cat("  Predictions for", N.loc[2], "locations and",
        N.loc[1], "time points.\n")
  }
  if( x$opts$pred.var ){
    cat("  Variance have been computed.\n")
  }else{
    cat("  Variance have NOT been computed.\n")
  }
  if( by.CV ){
    cat("  Same location in several CV-groups.\n")
  }

  return(invisible())
}##function print.predCVSTmodel

##' \code{\link[base:summary]{summary}} method for class \code{predCVSTmodel}.
##'
##' Computes summary statistics for cross validation. Statistics that are
##' computed include RMSE, R2, and coverage of CI:s; both for all observations
##' and (possibly) stratified by date.
##' 
##' @title Computes summary details for \code{predCVSTmodel} object
##' 
##' @param object \code{predCVSTmodel} object to compute summary information
##'   for; the output from \code{\link{predictCV.STmodel}}.
##' @param pred.naive Result of naive prediction; used to compute
##'   modified R2 values. The output from \code{\link{predictNaive}}.
##' @param by.date Compute individual cross-validation statistics for each
##'   time-point. May lead to \emph{very many} statistics.
##' @param p Approximate coverage of the computed confidence bands; the
##'   confidence bands are used when computing coverage of the
##'   cross-validated predictions.
##' @param transform Transform observations and predictions \emph{before}
##'   computing statistics; see also \code{\link{computeLTA}}
##' @param LTA Compute cross-validation statistics for the long term averages at
##'   each site, uses \code{\link{computeLTA}} to compute the averages.
##'   \code{transform} is passed to \code{\link{computeLTA}}.
##' @param ... Ignored additional arguments.
##' 
##' @return A \code{summary.predCVSTmodel} object.
##'
##' @examples
##' ##load some data
##' data(CV.mesa.model)
##' 
##' ##basic summary statistics
##' summary(pred.cv.mesa)
##' ##computed on the natural scale, and with long term averages
##' summary(pred.cv.mesa, transform=exp, LTA=TRUE)
##'
##' ##store the results
##' cv.summary <- summary(pred.cv.mesa, transform=exp, LTA=TRUE)
##' ##and study the contents
##' str(cv.summary)
##' 
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method summary predCVSTmodel
##' @export
summary.predCVSTmodel <- function(object, pred.naive=NULL, by.date=FALSE,
                                  p=0.95, transform=function(x){return(x)},
                                  LTA=FALSE, ...){
  ##check class belonging
  stCheckClass(object, "predCVSTmodel", name="object")
  ##inputs
  if(p<=0 || p>=1){
    stop("'p' not a probability.")
  }
  if( !is.function(transform) ){
    stop("'transform' should be a function")
  }

  ##compute predictions by date?
  if( by.date ){
    pred.by.date <- split(object$pred.obs, object$pred.obs$date)
  }else{
    pred.by.date <- NULL
  }##if( by.date ){...}else{...}
  
  ##total number of summaryStats to compute
  Nstat.naive <- Nstat <- 1 + LTA + length(pred.by.date)
  if( !is.null(pred.naive) ){
     Nstat.naive <- Nstat.naive + (dim(pred.naive$pred)[2]-2)
   }

  ##allocate memory for each of these
  R2 <- matrix(NA, Nstat.naive, 3)
  RMSE <- matrix(NA, Nstat, 3)
  coverage <- matrix(NA, Nstat, 1)
  ##add names to matrix
  colnames(coverage) <- "EX"
  colnames(RMSE) <- colnames(R2) <- c("EX.mu", "EX.mu.beta", "EX")
  Rn <- "obs"
  if(LTA){
    Rn <- c(Rn,"average")
  }
  Rn <- c(Rn, names(pred.by.date))
  rownames(RMSE) <- rownames(coverage) <- Rn
  ##also add naive predictions for R2
  if( !is.null(pred.naive) ){
    Rn <- c(Rn, names(pred.naive$pred)[!(names(pred.naive$pred) %in%
                                         c("ID","date"))])
  }
  rownames(R2) <- Rn

  ##recompute p for a two-sided CI
  p.org <- p
  q <- qnorm( (1+p)/2 )

  ##extract observations and predictions
  obs <- transform( object$pred.obs$obs )
  EX.mu <- transform( object$pred.obs$EX.mu )
  EX.mu.beta <- transform( object$pred.obs$EX.mu.beta )
  EX <- transform( object$pred.obs$EX )
  if( object$opts$pred.var ){
    res <- object$pred.obs$res.norm
  }

  ##drop NA:s
  I <- is.na(EX.mu) | is.na(EX.mu.beta) | is.na(EX)

  ##compute stats for raw observations
  RMSE["obs", "EX.mu"] <- sqrt(mean( (obs[!I] - EX.mu[!I])^2 ))
  RMSE["obs", "EX.mu.beta"] <- sqrt(mean( (obs[!I] - EX.mu.beta[!I])^2 ))
  RMSE["obs", "EX"] <- sqrt(mean( (obs[!I] - EX[!I])^2 ))

  R2["obs", ] <- 1 - RMSE["obs", ]^2 / var( obs[!I] )
  
  if( object$opts$pred.var ){
    coverage["obs",1] <- mean( abs(res[!I]) < q )
  }

  ##compute stats for long term average
  if( LTA ){
    lta.tmp <- computeLTA(object, transform)
    ##drop NA:s
    lta.tmp <- lta.tmp[!apply(is.na(lta.tmp),1,any),,drop=FALSE]

    RMSE["average", "EX.mu"] <- sqrt(mean( (lta.tmp$obs - lta.tmp$EX.mu)^2 ))
    RMSE["average", "EX.mu.beta"] <- sqrt(mean( (lta.tmp$obs -
                                                 lta.tmp$EX.mu.beta)^2 ))
    RMSE["average", "EX"] <- sqrt(mean( (lta.tmp$obs - lta.tmp$EX)^2 ))

    R2["average", ] <- 1 - RMSE["average", ]^2 / var( lta.tmp$obs )
  }#if( LTA )

  ##compute stats for each date
  if( length(pred.by.date)>0 ){
    for( i in 1:length(pred.by.date) ){
      ##extract observations and predictions
      obs <- transform( pred.by.date[[i]]$obs )
      EX.mu <- transform( pred.by.date[[i]]$EX.mu )
      EX.mu.beta <- transform( pred.by.date[[i]]$EX.mu.beta )
      EX <- transform( pred.by.date[[i]]$EX )
      if( object$opts$pred.var ){
        res <- pred.by.date[[i]]$res.norm
      }

      ##drop NA:s
      I <- is.na(EX.mu) | is.na(EX.mu.beta) | is.na(EX)

      ##compute stats for raw observations
      I.n <- names(pred.by.date)[i]
      RMSE[I.n, "EX.mu"] <- sqrt(mean( (obs[!I] - EX.mu[!I])^2 ))
      RMSE[I.n, "EX.mu.beta"] <- sqrt(mean( (obs[!I] - EX.mu.beta[!I])^2 ))
      RMSE[I.n, "EX"] <- sqrt(mean( (obs[!I] - EX[!I])^2 ))

      R2[I.n, ] <- 1 - RMSE[I.n, ]^2 / var( obs[!I] )
      
      if( object$opts$pred.var ){
        coverage[I.n,1] <- mean( abs(res[!I]) < q )
      }
    }## for( i in 1:length(pred.by.date) )
  }## if( length(pred.by.date)>0 )

  ##compute modified R2 for the naive predictions
  if( !is.null(pred.naive) ){
    ##check that pred.naive matches my predictions
    if( any( object$pred.obs$date != pred.naive$pred$date ) ){
      stop("Missmatch between dates in 'object' and 'pred.naive'")
    }
    if( any( object$pred.obs$ID != pred.naive$pred$ID ) ){
      stop("Missmatch between ID:s in 'object' and 'pred.naive'")
    }

    ##extract naive predictions
    I.naive <- !(names(pred.naive$pred) %in% c("ID","date"))
    pred <- transform( pred.naive$pred[,I.naive,drop=FALSE] )
    ##extract observations
    obs <- transform( object$pred.obs$obs )
    
    I <- !apply(is.na(pred), 1, any)
    ##compute modified R2:s
    for(i in 1:dim(pred)[2]){
      V.naive <- mean( (pred[I,i] - obs[I])^2 )
      R2[names(pred)[i], ] <- 1 - RMSE["obs",]^2 / V.naive
    }
  }
  
  ##ensure that R2 is >= 0
  R2 <- pmax(R2, 0)

  ##some CV-statistics
  stats <- list(N.cv=max(object$Ind.cv), Npred=sum(!is.na(object$pred.obs$EX)),
                pred.var=object$opts$pred.var)
  ##return the computed stats
  out <- list(RMSE=RMSE, R2=R2, p=p.org, stats=stats)
  if( object$opts$pred.var ){
    out$coverage <- coverage
  }
  class(out) <- "summary.predCVSTmodel"
  return(out)
}##function summary.predCVSTmodel


##' \code{\link[base:print]{print}} method for class \code{summary.predCVSTmodel}.
##'
##' @title Print details for \code{summary.predCVSTmodel} object
##' @param x \code{summary.predCVSTmodel} object to print information for.
##' @param ... Additional arguments, passed to
##'   \code{\link[base:print]{print.table}}.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method print summary.predCVSTmodel
##' @export
print.summary.predCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.predCVSTmodel", name="x")

  cat("Cross-validation predictions for STmodel with", x$stats$N.cv,
      "CV-groups.\n")
  cat("  Predictions for", x$stats$Npred, "observations.\n\n")

  cat("RMSE:\n")
  print( x$RMSE )
  cat("\n")
  
  cat("R2:\n")
  print( x$R2 )
  cat("\n")

  if( x$stats$pred.var ){
    cat("Coverage of ", 100*x$p, "% prediction intervalls:\n", sep="")
    print( x$coverage )
    cat("\n")
  }else{
    cat("Coverage (i.e. prediciton variance) have NOT been computed.\n")
  }
  cat("\n")

  return(invisible())
}##function print.summary.predCVSTmodel

##' @rdname plot.predictSTmodel
##' @family predCVSTmodel methods
##' @importFrom graphics plot
##' @method plot predCVSTmodel
##' @export
plot.predCVSTmodel <- function(x, y="time", ID=colnames(x$pred.all$EX)[1],
                                col=c("black","red","grey"), pch=c(NA,NA),
                                cex=c(1,1), lty=c(1,1), lwd=c(1,1), p=0.95,
                                pred.type="EX", pred.var=TRUE,
                                add=FALSE, ...){
  ##check class belonging
  stCheckClass(x, "predCVSTmodel", name="x")
  ##we have to use y, cast to resonable name
  plot.type <- y
  pred.var <- internalPlotPredictChecks(plot.type, pred.type, pred.var)
  ##check ID
  ID <- internalFindIDplot(ID, colnames(x$pred.all$EX))
  ID.all <- internalCheckIDplot(ID, y)

  ##plotting options
  if( !ID.all ){
    ##pick out predictions, two cases
    date <- rownames(x$pred.all[[pred.type]])
    sd <- x$pred.all[[pred.var]][,ID]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x$pred.all[[pred.type]][,ID], sd=sd, date=date,
                       stringsAsFactors=FALSE)
    ##pick out observations, if any
    obs <- x$pred.obs[x$pred.obs$ID==ID,,drop=FALSE]
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
      ID <- unique(x$pred.obs$ID)
    }
    I.ID <- x$pred.obs$ID %in% ID
    sd <- x$pred.obs[I.ID,pred.var]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x$pred.obs[I.ID,pred.type], sd=sd,
                       date=x$pred.obs[I.ID,"date"], stringsAsFactors=FALSE)
    ##pick out observations, if any
    obs <- x$pred.obs[I.ID,"obs"]
    if( dim(pred)[1]!=length(obs) ){
      stop("Number of observed locations do not match no. predicted.")
    }
  }##if( !ID.all ){...}else{...}

  internalPlotPredictions(plot.type, ID, pred, obs, col, pch, cex, lty, lwd, p, add)
  
  return(invisible())
}##function plot.predCVSTmodel
