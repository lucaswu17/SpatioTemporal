#######################################################################
## FUNCTIONS THAT COMPUTE SVDs FOR TEMPORAL TRENDS WITH MISSING DATA ##
#######################################################################
##Functions in this file:
## SVDmiss           - EX:ok
## SVDsmooth         - EX:ok
## SVDsmoothCV       - EX:with SVDsmooth
## calcSmoothTrends  - EX:ok
## plot.SVDcv        - EX:with SVDsmooth
## print.SVDcv       - EX:with SVDsmooth

##' Function that completes a data matrix using iterative svd as described in
##' Fuentes et. al. (2006). The function iterates between computing the svd for
##' the matrix and replacing the missing values by linear regression of the
##' columns onto the first \code{ncomp} svd components. As initial replacement
##' for the missing values regression on the column averages are used. The
##' function \emph{will fail} if entire rows and/or columns are missing from the
##' data matrix.
##'
##' @title Missing Data SVD
##' @param X Data matrix, with missing values marked by \code{NA}.
##' @param niter Maximum number of iterations to run before exiting,
##'   \code{Inf} will run until the \code{conv.reldiff} criteria is met.
##' @param ncomp Number of SVD components to use in the reconstruction.
##' @param conv.reldiff Assume the iterative procedure has converged when
##'   the relative difference between two consecutive iterations is less than
##'   \code{conv.reldiff}.
##' @return A list with the following components:
##'   \item{Xfill}{The completed data matrix with missing values replaced by
##'                fitting the data to the \code{ncomp} most important svd
##'                components}
##'   \item{svd}{The result of svd on the completed data matrix, i.e.
##'              \code{svd(Xfill)}}
##'   \item{status}{A vector of status variables: \code{diff}, the
##'                 absolute difference between the two last iterations;
##'                 \code{rel.diff}, the relative difference;
##'                 \code{n.iter}, the number of iterations; and
##'                 \code{max.iter}, the requested maximum number of
##'                 iterations.}
##' 
##' @references
##' M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
##'  Analyze Space-Time Processes in Statistical methods for spatio-temporal
##'  systems (B. Finkenstädt, L. Held, V. Isham eds.) 77-150
##' 
##' @author Paul D. Sampson and Johan Lindström
##'
##' @example Rd_examples/Ex_SVDmiss.R
##' @family SVD for missing data
##' @family data matrix
##' @export
SVDmiss <- function(X, niter=25, ncomp=min(4,dim(X)[2]), conv.reldiff=0.001)
{
  ##First find missing values
  Ina <- is.na(X)
  if( all(!Ina) ){
    ##if X has no missing data this is simple
    svd0 <- svd(X)
    XF <- X
    i <- diff <- reldiff <- 0
  }else{
    ##X has missing data, use iterative method
    ##Iterative svd calculation with missing data.
    ##Initial first element of U matrix is average curve.
    u1 <- rowMeans(X, na.rm = TRUE)
    XM <- matrix(1, nrow(X), ncol(X))
    XM[Ina] <- 0
    XZ <- X
    # v1 is proportional to X'u1/(u1'u1), but with calculations
    # filling in zeros for missing values in the sums.
    XZ[Ina] <- 0.
    # Fill in missing values for initial complete SVD calculation.
    # Then iterate  using only complete data.
    v1 <- diag(t(XZ) %*% (XM * u1))/diag(t(XM * u1) %*% (XM * u1))
    XF <- X
    XF[Ina] <- (matrix(u1, ncol = 1) %*% matrix(v1, nrow = 1))[Ina]
    if( any(is.na(XF)) )
      stop("Unable to complete matrix, too much missing data")
    reldiff <- conv.reldiff+1
    i <- 0
    while(i<niter && reldiff>conv.reldiff){
      svd0 <- svd(XF)
      Xnew <- X
      Xnew[Ina] <- (svd0$u[, 1:ncomp] %*%
                    diag(svd0$d[1:ncomp],nrow=length(svd0$d[1:ncomp])) %*%
                    t(svd0$v[,1:ncomp]))[Ina]
      diff <- max(abs(Xnew - XF))
      reldiff <- diff/max(abs(XF[Ina]))
      XF <- Xnew
      i <- i+1
    }
  }#if( all(!is.na(X)) ) ... else ...
  final.diff <- c(diff,reldiff,i,niter)
  names(final.diff) <- c("diff","rel.diff","n.iter","max.iter")
  return(list(svd=svd0, Xfill=XF, status=final.diff))
}##function SVDmiss

##' Function that computes smooth functions for a data matrix with missing
##' values, as described in Fuentes et. al. (2006), or does cross validation to
##' determine a suitable number of basis functions. The function uses
##' \code{\link{SVDmiss}} to complete the matrix and then computes smooth basis
##' functions by applying \code{\link[stats:smooth.spline]{smooth.spline}} to the
##' SVD of the completed data matrix.
##' 
##' \code{SVDsmoothCV} uses leave-one-column-out cross-validation; holding one column
##' out from \code{X}, calling \code{SVDsmooth}, and then regressing
##' the held out column on the resulting smooth functions. Cross-validation
##' statistics computed include RMSE, R-squared and BIC.
##' 
##' @title Smooth Basis Functions for Data Matrix with Missing Values
##' 
##' @param X Data matrix, with missing values marked by \code{NA}. Rows and/or columns
##'   that are completely missing will be dropped (with a message), for the rows the
##'   smooths will be interpolated using
##'   \code{\link[stats:predict.smooth.spline]{predict.smooth.spline}}.
##' @param n.basis Number of smooth basis functions to compute, will be passed
##'   as \code{ncomp} to \code{\link{SVDmiss}}; for \code{SVDsmoothCV} a
##'   vector with the different number of basis functions to evaluate.
##' @param date.ind Vector giving the observation time of each row in
##'   \code{X}, used as \code{x} in \cr
##'   \code{\link[stats:smooth.spline]{smooth.spline}} when computing the smooth
##'   basis functions. If missing \code{\link{convertCharToDate}} is used to
##'   coerce the \code{rownames(X)}.
##' @param scale If \code{TRUE}, will use \code{\link[base:scale]{scale}} to scale
##'   \code{X} before calling \code{\link{SVDmiss}}.
##' @param niter,conv.reldiff Controls convergence, passed to \code{\link{SVDmiss}}.
##' @param df,spar The desired degrees of freedom/smoothing parameter for the
##'   spline, \cr see \code{\link[stats:smooth.spline]{smooth.spline}}
##' @param ... Additional parameters passed to \code{SVDsmooth}; i.e. \code{date.ind},
##'   \code{scale}, \code{niter}, \code{conv.reldiff}, \code{df}, and \code{spar}.
##' 
##' @return Depends on the function:
##'   \item{SVDsmooth}{A matrix where each column is a smooth basis function
##'                    based on the SVD of the completed data matrix. The left
##'                    most column contains the smooth of the most important SVD.}
##'   \item{SVDsmoothCV}{A list of class \code{SVDcv} with components:
##'     \describe{
##'       \item{CV.stat}{A \code{data.frame} with statistics for each of the
##'                      number of basis functions evaluated.}
##'       \item{BIC.all}{A \code{data.frame} with the individual BIC values for
##'                      each column in the data matrix and for each number of
##'                      basis functions evaluated.}
##'       \item{smoothSVD}{A list with \code{length(n.basis)} components. Each
##'                        component contains an array where \code{smoothSVD[[j]][,,i]}
##'                        is the result of \code{SVDsmooth} applied to
##'                        \code{X[,-i]} with \code{n.basis[j]} smooth functions.}
##'     }
##'   }
##' 
##' @references
##' M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
##'  Analyze Space-Time Processes in Statistical methods for spatio-temporal
##'  systems (B. Finkenstädt, L. Held, V. Isham eds.) 77-150
##' 
##' @author Paul D. Sampson and Johan Lindström
##'
##' @example Rd_examples/Ex_SVDsmooth.R
##' @family SVD for missing data
##' @family data matrix
##' @family SVDcv methods
##' @export
SVDsmooth <- function(X, n.basis=min(2,dim(X)[2]), date.ind=NULL, scale=TRUE,
                      niter=100, conv.reldiff=0.001, df=NULL, spar=NULL){

  ##date.ind contains NA/missing, use rownames of X
  if( missing(date.ind) || is.null(date.ind) || any(is.na(date.ind)) ){
    date.ind <- convertCharToDate( rownames(X) )
  }
  if( is.null(date.ind) ){
    date.ind <- 1:dim(X)[1]
  }
  ##check length od date.ind
  if( length(date.ind)!=dim(X)[1] ){
    stop("length(date.ind)!=dim(X)[1]")
  }
  ##check number of columns
  if( max(n.basis)>dim(X)[2] ){
    stop("Number of basis functions cannot exceed dim(X)[2]")
  }
  ##scale data to zero mean and unit variance
  if( scale ){
    X <- scale(X)
  }
  ##drop missing columns and rows
  Icol <- colSums(is.na(X))!=dim(X)[1]
  Irow <- rowSums(is.na(X))!=dim(X)[2]
  X.reduced <- X[Irow,Icol]
  if( any(!Icol) ){
    message( paste("Dropping column(s)", paste(which(!Icol), collapse=", ")) )
  }
  if( any(!Irow) ){
    message( paste("Smooth interpolated at time(s):",
                   paste(as.character(date.ind[!Irow]), collapse=", ")) )
  }
  ##missing data SVD
  X.svd <- SVDmiss(X.reduced, niter=niter, ncomp=n.basis,
                   conv.reldiff=conv.reldiff)
  ##calculate basis functions
  X.comps <- matrix(NA, length(date.ind), n.basis)
  for(j in 1:n.basis){
    if( is.null(df) && is.null(spar) ){
      spline <- smooth.spline(date.ind[Irow], X.svd$svd$u[,j])
    }else{
      spline <- smooth.spline(date.ind[Irow], X.svd$svd$u[,j], df=df, spar=spar)
    }
    X.comps[,j] <- predict(spline, as.double(date.ind))$y
    #scale components to unit variance and zero mean
    X.comps[,j] <- scale(X.comps[,j],  center=mean(X.comps[Irow,j]),
                         scale=sd(X.comps[Irow,j]))
    #Ensure that components have alternating sign
    X.comps[,j] <- (-1)^j*X.comps[,j]*sign(X.comps[1,j])
  }
  rownames(X.comps) <- as.character(date.ind)
  colnames(X.comps) <- paste("V",1:n.basis,sep="")
  return(X.comps)
}##function SVDsmooth

##' @rdname SVDsmooth
##' @export
SVDsmoothCV <- function(X, n.basis, ...){
  ##check number of columns
  if( max(n.basis)>(dim(X)[2]-1) )
    stop("Number of basis functions cannot exceed dim(X)[2]-1.")
  
  ##matrix holding CV statistics
  CV.stat <- matrix(NA,length(n.basis),3)
  colnames(CV.stat) <- c("RMSE", "R2", "BIC")
  rownames(CV.stat) <- paste("n.basis", as.character(n.basis), sep=".")
  
  ##matrix that holds the individual BIC:s
  BIC <- matrix(NA, length(n.basis), dim(X)[2])
  colnames(BIC) <- colnames(X)
  rownames(BIC) <- rownames(CV.stat)
  ##matrix with the residual sum of squares
  RSS <- matrix(NA, length(n.basis), dim(X)[2])

  ##loop over the number of basis functions requested
  trend <- list()
  for(i in 1:length(n.basis)){
    ##for each basis functions do leave one out CV
    trend[[i]] <- array(NA, c(dim(X)[1], n.basis[i], dim(X)[2]))
    err <- matrix(NA, dim(X)[1], dim(X)[2])
    for(j in 1:dim(X)[2]){
      if(i==1 && j==1){
        trend[[i]][,,j] <- SVDsmooth(X[,-j,drop=FALSE], n.basis[i], ...)
      }else{
        ##suppress message given above to avoid noisy CV output
        suppressMessages( trend[[i]][,,j] <- SVDsmooth(X[,-j,drop=FALSE],
                                                       n.basis[i], ...))
      }
      I <- !is.na(X[,j])
      if( any(I) ){
        ##compute CV-error
        err[I,j] <- lm(X[,j]~trend[[i]][,,j])$residuals
        ##compute BIC
        RSS[i,j] <- sum( err[I,j]^2 )
        BIC[i,j] <- sum(I)*log(RSS[i,j]/sum(I)) + n.basis[i]*log(sum(I))
      }else{
        RSS[i,j] <- BIC[i,j] <- NA
      }
    }
    #compute CV statistics
    CV.stat[i,"RMSE"] <- sqrt( mean(err^2,na.rm=TRUE) )
    CV.stat[i,"R2"] <- 1 - CV.stat[i,"RMSE"]^2 / var(c(X),na.rm=TRUE)
    dimnames(trend[[i]]) <- list(rownames(X), paste("V", 1:n.basis[i], sep=""),
                                 colnames(X))
  }##for(i in 1:length(n.basis))
  ##compute total BIC for each CV-set
  n.all <- sum(!is.na(X))
  RSS.all <- rowSums(RSS, na.rm=TRUE)
  p.all <- (1+n.basis)*dim(X)[2]
  CV.stat[,"BIC"] <- n.all*log(RSS.all/n.all) + p.all*log(n.all)
  ##names for the list containing all the trends
  names(trend) <- rownames(CV.stat)
  ##return a class object
  out <- list(CV.stat=as.data.frame(CV.stat), BIC.all=as.data.frame(t(BIC)),
              smoothSVD=trend)
  class(out) <- "SVDcv"
  return( out )
}##function SVDsmoothCV


##' A front end function for calling \code{\link{SVDsmooth}} (and
##' \code{\link{SVDsmoothCV}}), with either a \code{STdata} object
##' or vectors containing observations, dates and locations.
##' 
##' The function uses \code{\link{createDataMatrix}} to create a
##' data matrix which is passed to \code{\link{SVDsmooth}} (and
##' \code{\link{SVDsmoothCV}}). The output can be used as \cr
##' \code{STdata$trend = calcSmoothTrends(...)$trend}, or \cr
##' \code{STdata$trend = calcSmoothTrends(...)$trend.cv[[i]]}.
##' However, it is recommended to use \code{\link{updateSTdataTrend}}.
##'
##' @title Smooth Basis Functions for a STdata Object
##' 
##' @param STdata A \code{STdata}/\code{STmodel} data structure containing
##'   observations, see \code{\link{mesa.data}}. Use either this or the \code{obs},
##'   \code{date}, and \code{ID} inputs.
##' @param obs A vector of observations.
##' @param date A vector of observation times.
##' @param ID A vector of observation locations.
##' @param subset A subset of locations to extract the data matrix for. A warning
##'   is given for each name not found in \code{ID}.
##' @param extra.dates Additional dates for which smooth trends should be computed.
##' @param n.basis Number of basis functions to compute, see
##'   \code{\link{SVDsmooth}}.
##' @param cv Also compute smooth functions using leave one out
##'   cross-validation, \cr see \code{\link{SVDsmoothCV}}.
##' @param ... Additional parameters passed to \code{\link{SVDsmooth}}
##'   and \code{\link{SVDsmoothCV}}.
##' 
##' @return Returns a list with
##'   \item{trend}{A data.frame containing the smooth trends and the dates.
##'                This can be used as the \code{trend} in \code{STdata$trend}.}
##'   \item{trend.cv}{If \code{cv==TRUE} a list of data.frames; each one containing
##'                   the smooth trend obtained when leaving one site out.
##'                   Similar to \cr \code{SVDsmoothCV(data)$smoothSVD[[1]]}).}
##' 
##' @author Johan Lindström and Paul D. Sampson
##' 
##' @example Rd_examples/Ex_calcSmoothTrends.R
##' @family SVD for missing data
##' @family STdata
##' @export
calcSmoothTrends <- function(STdata=NULL, obs=STdata$obs$obs,
                             date=STdata$obs$date, ID=STdata$obs$ID, subset=NULL,
                             extra.dates=NULL, n.basis=2, cv=FALSE, ...){
  ##add extra dates
  date <- c(date, extra.dates)
  ##and expand obs and ID
  obs <- c(obs, rep(NA,length(extra.dates)) )
  ID <- c(ID, rep(ID[1],length(extra.dates)) )
  ##create data matrix
  data <- createDataMatrix(obs=obs, date=date, ID=ID, subset=subset)

  ##internal function
  extractTrend <- function(x){
    x <- as.data.frame(x)
    x$date <- convertCharToDate( rownames(x) )
    rownames(x) <- NULL
    return(x)
  }
  ##now let's do SVD
  data.comps <- SVDsmooth(data, n.basis, ...)
  data.comps <- extractTrend(data.comps)
  ##and cross-validation
  if(cv){
    ##Message of interpolation alreadt displayed above...
    suppressMessages( svd.cv <- SVDsmoothCV(data, n.basis, ...) )
    svd.cv <- svd.cv$smoothSVD[[1]]
    svd.tmp <- list()
    for(i in 1:dim(svd.cv)[3]){
      svd.tmp[[i]] <- extractTrend(svd.cv[,,i])
    }
    names(svd.tmp) <- dimnames(svd.cv)[[3]]
  }else{
    svd.tmp <- NULL
  }
  return( list(trend=data.comps, trend.cv=svd.tmp) )
}##function calcSmoothTrends


##########################
## S3-METHODS FOR SVDcv ##
##########################

##' \code{\link[graphics:plot]{plot}} method for class \code{SVDcv}.
##' Plots summary statistics for the cross-validation. Plots include
##' RMSE, R2, BIC, and scatter plots of BIC for each column.
##'
##' @title Plot cross-validation statistics for \code{SVDcv} object
##' @param x \code{SVDcv} object to plot.
##' @param y Not used
##' @param pairs \code{TRUE}/\code{FALSE} plot cross-validation statistics,
##'   or scatter plot of individual BIC:s.
##' @param ... Additional parameters passed to \code{\link[graphics:plot]{plot}} or
##'   \code{\link[graphics:plot]{pairs}}.
##' @return Nothing
##'
##' @examples
##'   ## end of SVDsmooth example
##' 
##' @author Johan Lindström
##' 
##' @family SVDcv methods
##' @family SVD for missing data
##' @method plot SVDcv
##' @export
plot.SVDcv <- function(x, y=NULL, pairs=FALSE, ...){
  stCheckClass(x, "SVDcv", "'x'")

  if( !pairs ){
    ##plot cross-validation statistics
    par(mfcol=c(2,2),mar=c(4,4,.5,.5))
    plot(x$CV.stat$RMSE, type="l", ylab="RMSE", ...)
    plot(x$CV.stat$R2, type="l", ylab="R2", ...)
    plot(x$CV.stat$BIC, type="l", ylab="BIC", ...)
  }else{
    ##plot the BIC for each column
    pairs(x$BIC.all, panel=function(x,y){points(x,y); abline(0,1)}, ...)
  }
  return(invisible())
}##function plot.SVDcv

##' \code{\link[base:print]{print}} method for class \code{SVDcv}, prints
##' cross-validation statistics.
##'
##' @title Print details for \code{SVDcv} object
##' @param x \code{SVDcv} object to print information for.
##' @param ... ignored additional arguments.
##' @return Nothing
##' 
##' @examples
##'   ## end of SVDsmooth example
##' 
##' @author Johan Lindström
##' 
##' @family SVDcv methods
##' @family SVD for missing data
##' @method print SVDcv
##' @export
print.SVDcv <- function(x, ...){
  stCheckClass(x, "SVDcv", "'x'")
  cat("Result of SVDsmoothCV, summary of cross-validation:\n")
  print(x$CV.stat)
  cat("\nIndividual BIC:s for each column:\n")
  print( apply(x$BIC.all, 2, summary) )
  return(invisible())
}##function print.SVDcv 
