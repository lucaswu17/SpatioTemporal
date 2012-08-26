################################
## INTERNAL UTILITY FUNCTIONS ##
################################
##INTERNAL functions in this file:
## commonPrintST
## commonSummaryST
## internalPlotPredictions
## internalPlotPredictChecks
## internalFindIDplot
## internalCheckIDplot
## stCheckLoglikeIn
## stCheckType
## stCheckX
## internalFixNuggetUnobs
## checkDimInternal
## createSTmodelInternalDistance

########################################################################
## Parts of S3 print that are common for STdata and STmodel, INTERNAL ##
########################################################################
commonPrintST <- function(x, name, print.type, type=NULL){
  if( print.type==1 ){
    ##general information regarding number of observations
    cat( sprintf("%s-object with:\n", name) )
    cat(sprintf("\tNo. locations: %d (observed: %d)\n",
                dim(x$covars)[1], length(unique(x$obs$ID)) ))
    cat(sprintf("\tNo. time points: %d (observed: %d)\n",
                max(dim(x$trend)[1], length(unique(x$obs$date))),
                length(unique(x$obs$date)) ))
    cat(sprintf("\tNo. obs: %d\n\n",length(x$obs$obs) ))
    ##Trends
    if( is.null(x$trend) ){
      cat("No trend specified\n")
    }else{
      if( dim(x$trend)[2]==1 ){
        cat("Only constant temporal trend,")
      }else{
        cat( sprintf("Trend with %d basis function(s):\n",
                     dim(x$trend)[2]-1))
        print( names(x$trend)[ -which(names(x$trend)=="date") ] )
      }
      cat( sprintf("with dates:\n\t%s\n",
                   paste(range(x$trend$date), collapse=" to ") ))
    }##if( is.null(x$trend) ){...}else{...}
    if( !is.null(x$old.trend) ){
      cat( "Observations have been detrended.\n")
    }
    cat("\n")
  }##if( print.type==1 )
  if( print.type==2 && !is.null(type) ){
    if( length(type)!=dim(x$covars)[1] ){
      stop( paste("length(type) =",length(type),
                  "has to equal the number fo sites =",
                  dim(x$covars)[1]) )
    }
    if( !is.factor(type) ){
      warning("Attempting to coerc 'type' to factor")
      type <- as.factor(type)
    }
    cat( "All sites:")
    print( table(type,dnn="") )
    cat( "Observed:")
    print( table(type[x$covars$ID %in% unique(x$obs$ID)]) )
    cat("\n")
    for(i in  levels(type)){
      I <- (x$obs$ID %in% x$covars$ID[type==i])
      cat( sprintf("For %s:\n",i) )
      if( sum(I)!=0 ){
        cat( sprintf("  Number of obs: %d\n", sum(I)) )
        cat( sprintf("  Dates: %s\n", paste(range(x$obs$date[I]),
                                            collapse=" to ")) )
      }else{
        cat("  No observations\n")
      }
    }##for(i in  levels(type))
  }##if( print.type==2 && !is.null(type) )
}##function commonPrintST

##########################################################################
## Parts of S3 summary that are common for STdata and STmodel, INTERNAL ##
##########################################################################
commonSummaryST <- function(object, type=NULL){
  out <- list()
  if( dim(object$obs)[1]!=0 ){
    out$obs <- summary( object$obs[, c("obs","date"), drop=FALSE])
  }
  if( !is.null(object$trend) ){
    out$trend <- summary(object$trend)
  }
  ##and possibly observations by type
  if( !is.null(type) ){
    if( length(type)!=dim(object$covars)[1] ){
      stop( paste("length(type) =",length(type),
                  "has to equal the number fo sites =",
                  dim(object$covars)[1]) )
    }
    if( !is.factor(type) ){
      warning("Attempting to coerc 'type' to factor")
      type <- as.factor(type)
    }
    out$obs.by.type <- vector("list", length(levels(type)))
    names(out$obs.by.type) <- as.character(levels(type))
    for(i in levels(type)){
      I <- (object$obs$ID %in% object$covars$ID[type==i])
      if( sum(I)!=0 ){
        out$obs.by.type[[i]] <- summary(object$obs[I, c("obs","date"), drop=FALSE])
      }
    }##for(i in  levels(type))
  }##if( !is.null(type) )

  return(out)
}##function commonSummaryST

###########################################################################
## Parts of S3 plot that are common for predictSTmodel and predCVSTmodel ##
###########################################################################
internalPlotPredictions <- function(plot.type, ID, pred, obs, col, pch, cex,
                                    lty, lwd, p, add){
  ##compute the quantile
  q <- qnorm((1-p)/2, lower.tail=FALSE)
  ##ensure that lty, lwd, pch, and cex are of length==2
  if( length(lty)==1 ){
    lty = c(lty,lty)
  }
  if( length(lwd)==1 ){
    lwd = c(lwd,lwd)
  }
  if( length(pch)==1 ){
    pch = c(pch,pch)
  }
  if( length(cex)==1 ){
    cex = c(cex,cex)
  }

  if( plot.type=="obs" ){
    ##drop unobserved
    pred <- pred[!is.na(obs),,drop=FALSE]
    obs <- obs[!is.na(obs)]
    ##plot as a function of sorted observations
    I <- order(obs)
    ##reorder
    obs <- obs[I]
    pred <- pred[I,,drop=FALSE]
    ##use t for x-axis below
    t <- obs
    xlab <- "observations"
  }else{
    ##dates
    t <- convertCharToDate(pred$date)
    xlab <- "time"
  }
  if( length(pred$x)==0 ){
    stop( paste("No observations for ID =", paste(ID,collapse=", ")) )
  }
  
  if( length(ID)==1 ){
    ID.name <- ID
  }else{
    ID.name <- paste(range(ID),collapse=" to ")
  }
  ##plot the results
  if(!add){
    plot(t, pred$x, type="n", ylab="predictions", xlab=xlab, main=ID.name,
         ylim=range(c(pred$x+q*pred$sd, pred$x-q*pred$sd, obs, pred$x),
           na.rm=TRUE))
  }
  ##Plot the polygon, NA if missing -> no polygon
  polygon(c(t,rev(t)), c(pred$x+q*pred$sd, rev(pred$x-q*pred$sd)),
          border=col[3], col=col[3])
  ##plot the predictions
  if( !is.na(lty[1]) )
    lines(t, pred$x, col=col[1], lty=lty[1], lwd=lwd[1])
  if( !is.na(pch[1]) )
    points(t, pred$x, col=col[1], pch=pch[1], cex=cex[1])
  if( plot.type=="time" ){
    ##and the observations
    if( !is.na(lty[2]) )
      lines(t, obs, col=col[2], lty=lty[2], lwd=lwd[2])
    if( !is.na(pch[2]) )
      points(t, obs, col=col[2], pch=pch[2], cex=cex[2])
  }else{
    if( !is.na(lty[2]) )
      abline(0, 1, lty=lty[2], col=col[2], lwd=lwd[2])
  }
  
  ##return
  return(invisible())
}##function internalPlotPredictions


###########################################################################
## Common helper functions for S3 plot.predictSTmodel/plot.predCVSTmodel ##
###########################################################################
internalPlotPredictChecks <- function(plot.type, pred.type, pred.var){
  if( !(plot.type %in% c("time", "obs")) ){
    stop("Unknown option for 'y'; should be time, obs")
  }
  ##check pred.type to use
  if( !(pred.type %in% c("EX", "EX.mu", "EX.mu.beta")) ){
    stop("Unknown option for 'pred.type'; should be EX, EX.mu, EX.mu.beta")
  }
  if( pred.type!="EX" ){
    pred.var <- "DO NOT USE"
  }else{
    if( pred.var ){
      pred.var <- "VX.pred"
    }else{
      pred.var <- "VX"
    }
  }
  return(pred.var)
}##function internalPlotPredictChecks

internalFindIDplot <- function(ID, names){
  if( is.null(ID) ){
    ID <- 1
  }
  if( is.numeric(ID) ){
    ID <- names[ID]
  }
  if( !is.character(ID) ){
    stop("ID not character, or not inferable.")
  }
  return(ID)
}##function internalFindIDplot

internalCheckIDplot <- function(ID, y){
  if( ID=="all" || length(ID)!=1 ){
    ID.all <- TRUE
  }else{
    ID.all <- FALSE
  }
  if( ID.all && y=="time" ){
    stop("ID=all (or several locations) and y=time incompatable")
  }
  return(ID.all)
}##funciton internalCheckIDplot

#########################################################################
## Check parameter input to loglikeST, and related functions, INTERNAL ##
#########################################################################
stCheckLoglikeIn <- function(x, x.fixed, type){
  stCheckType(type)
  if( !is.null(x.fixed) ){
    if( length(x)!=sum(is.na(x.fixed)) ){
      stop("length(x) must match number of NA:s in x.fixed.")
    }
    x.fixed[ is.na(x.fixed) ] <- x
    x <- x.fixed
  }
  return(x)
}

stCheckType <- function(type){
  ##check if type is valid
  if( !(type %in% c("r","p","f")) ){
    stop("Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")
  }
}

############################################################
## Check parameter input to estimation and MCMC, INTERNAL ##
############################################################
stCheckX <- function(x, x.fixed, dimensions, type, object){
  ##Check if we need to truncate or expand
  if(dim(x)[1] != dimensions$nparam && dim(x)[1] != dimensions$nparam.cov){
    stop( paste("dim(x)[1] must be", dimensions$nparam,
                "or", dimensions$nparam.cov) )
  }
  ##check x.fixed
  if( !is.null(x.fixed) && (!is.vector(x.fixed) || !is.numeric(x.fixed)) ){
    stop("'x.fixed' must be a numeric vector.")
  }else if( is.null(x.fixed) ){
    ##expand
    x.fixed <- rep(NA, dim(x)[1])
  }else if(length(x.fixed) != dimensions$nparam &&
           length(x.fixed) != dimensions$nparam.cov){
    stop( paste("length(x.fixed) must be", dimensions$nparam,
                "or", dimensions$nparam.cov) )
  }
 
  if( dim(x)[1]==dimensions$nparam && type!="f"){
    ##requested REML or profile but provided full parameter
    ##starting points -> truncate
    I <- (dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam
    if( length(x.fixed)==dimensions$nparam ){
      x.fixed <- x.fixed[I]
    }
    x <- x[I,,drop=FALSE]
  }else if( dim(x)[1]==dimensions$nparam.cov && type=="f" ){
    ##requested full but provided parameters for REML or profile -> expand
    x.old <- x
    x <- matrix(NA, dimensions$nparam, dim(x)[2])
    for(i in 1:dim(x)[2]){
      ##compute alpha and gamma as cond. exp. given data
      tmp <- predict.STmodel(object, x.old[,i], only.pars=TRUE, type="p")$pars
      x[,i] <- c(tmp$gamma.E, tmp$alpha.E, x.old[,i])
    }
    if( length(x.fixed)==dimensions$nparam.cov ){
      x.fixed <- c(rep(NA,dimensions$nparam-dimensions$nparam.cov), x.fixed)
    }
  }
  ##add names to x.fixed.
  if( type!="f" ){
    names(x.fixed) <- loglikeSTnames(object, all=FALSE)
  }else{
    names(x.fixed) <- loglikeSTnames(object, all=TRUE)
  }

  ##reduce x
  x.all <- x
  x <- x[is.na(x.fixed),,drop=FALSE]
  x.all[!is.na(x.fixed),] <- matrix(x.fixed[!is.na(x.fixed)],
                                    ncol=dim(x.all)[2])
  return( list(x.all=x.all, x=x, x.fixed=x.fixed) )
}##fucntion stCheckX

###################################################
## Create nugget for all sites in STmodel adding ##
## nugget.unobs.in to unobserved locations.      ##
###################################################
internalFixNuggetUnobs <- function(nugget.unobs.in, STmodel, nugget){
  ##nugget fo unobserved sites
  nugget.unobs <- matrix(NA, length(STmodel$locations$ID), 1)
  I <- match(rownames(nugget), STmodel$locations$ID)
  if( !any(is.na(I)) ){
    nugget.unobs[ I ] <- nugget
  }
  if( length(nugget.unobs.in)!=1 &&
     length(nugget.unobs.in)!=sum(is.na(nugget.unobs)) ){
    stop( paste("Needs,", 1, "or", sum(is.na(nugget.unobs)),
                "elements in nugget.unobs") )
  }
  nugget.unobs[ is.na(nugget.unobs) ] <- nugget.unobs.in
  ##add names
  rownames(nugget.unobs) <- STmodel$locations$ID
  return( nugget.unobs )
}##function internalFixNuggetUnobs

###########################################################
## Check dimensions for block matrix functions, INTERNAL ##
###########################################################
checkDimInternal <- function(X){
  dim.X <- sapply(X, dim)
  if( any(dim.X[1,1]!=dim.X[1,]) ){
    stop("all elements in X must have same number of columns.")
  }
  dim <- list(n=dim.X[1,1], m=length(X), p=dim.X[2,])
  return(dim)
}##function checkDimInternal

###########################################################
## Compute distance-matrices for createSTmodel, INTERNAL ##
###########################################################
createSTmodelInternalDistance <- function(STmodel){
  if( dim(STmodel$obs)[1]!=0 ){
    I.idx <- unique(STmodel$obs$idx)
    ##calculate distance matrices (only for observed locations),
    ##different for nu and beta
    I.idx <- 1:max(I.idx)
    STmodel$D.nu <- crossDist(STmodel$locations[I.idx, c("x.nu","y.nu"),
                                                drop=FALSE])
    colnames(STmodel$D.nu) <- rownames(STmodel$D.nu) <- STmodel$locations$ID[I.idx]
    STmodel$D.beta <- crossDist(STmodel$locations[I.idx, c("x.beta","y.beta"),
                                                  drop=FALSE])
    colnames(STmodel$D.beta) <- rownames(STmodel$D.beta) <- STmodel$locations$ID[I.idx]
    
    ##count number of observations at each time point
    dates <- sort(unique(STmodel$obs$date))
    STmodel$nt <- double(length(dates))
    for(i in c(1:length(dates))){
      STmodel$nt[i] <- sum( STmodel$obs$date==dates[i] )
    }
  }
  return( STmodel )
}##function createSTmodelInternalDistance
