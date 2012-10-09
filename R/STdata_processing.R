##############################################
## FUNCTIONS THAT MANIPULATE STdata OBJECTS ##
##############################################
##Functions in this file:
## createDataMatrix    EX:ok
## removeSTcovarMean   EX:ok
## updateSTdataTrend   EX:ok
## detrendSTdata       EX:ok

##' Creates a data matrix from a \code{STdata}/\code{STmodel} object, e.g.
##' \code{\link{mesa.data}} or \code{\link{mesa.model}}. Missing
##' observations are marked as \code{NA}.
##' 
##' @title Create a Data Matrix
##' @param STdata A \code{STdata}/\code{STmodel} object containing
##'   observations, see \code{\link{mesa.data}}. Use either this or the \code{obs},
##'   \code{date}, and \code{ID} inputs.
##' @param obs A vector of observations.
##' @param date A vector of observation times.
##' @param ID A vector of observation locations.
##' @param subset A subset of locations to extract the data matrix for. A warning
##'   is given for each name not found in \code{ID}.
##' @return Returns a matrix with dimensions (number of timepoints)-by-(number
##'   of locations). Row and column names of the matrix are taken as \code{ID}
##'   and \code{sort(unique(date))} respectively.
##' @author Johan Lindström
##' 
##' @example Rd_examples/Ex_createDataMatrix.R
##' 
##' @family data matrix
##' @family STdata functions
##' @family STmodel functions
##' @export
createDataMatrix <- function(STdata=NULL, obs=STdata$obs$obs,
                             date=STdata$obs$date, ID=STdata$obs$ID,
                             subset=NULL){
  if( !is.null(STdata) ){
    ##check class belonging
    stCheckClass(STdata, c("STdata","STmodel"), name="STdata")
  }
  if( any(length(obs)!=c(length(date),length(ID))) ){
    stop("'obs', 'date' and 'ID' are of unequal length")
  }
  ##cast ID to character (just to be safe)
  ID <- as.character(ID)
  ##find unique indecies
  date.ind <- sort(unique(date))
  ID.ind <- sort(unique(ID))
  ##construct the data matrix
  data <- matrix(NA, length(date.ind), length(ID.ind))
  ##drop NA values from obs, and the other vectors
  ID <- ID[!is.na(obs)]
  date <- date[!is.na(obs)]
  obs <- obs[!is.na(obs)]
  ##find indecies in the matrix
  I <- (match(ID,ID.ind)-1)*length(date.ind) + match(date,date.ind)
  data[I] <- obs
  ##add col- and rownames
  colnames(data) <- ID.ind
  rownames(data) <- as.character(date.ind)
  
  if( !is.null(subset) ){
    ##subset the data - first cast to character
    subset <- as.character(subset)
    ##check that all requested names exist in the colnames
    if( any(!(subset %in% colnames(data))) )
      warning( paste("Subset names not found in ID:",
                     paste(subset[!(subset %in% colnames(data))],collapse=", ")) )
    ##convert from a character vector to a index vector
    subset <- which(colnames(data) %in% subset)
    data <- data[,subset,drop=FALSE]
  }
  return(data)
}## createDataMatrix



##' Removes the temporal mean at each location for the spatio-temporal
##' covariares. The means are added to the \code{covar} field in the returned
##' object and can be used as geographic covariates.
##' 
##' @title Mean-Centre the Spatio-Temporal Covariate
##' @param STdata A \code{STdata} object, see \code{\link{mesa.data}}.
##' @return Returns a modfied version of the input, where the spatio-temporal
##'   covariates have been expanded to include covariates where the site by
##'   site temporal average has been removed. The averages are seen as geographic
##'   covariates and added to \code{STdata$covars}.
##' 
##' @example Rd_examples/Ex_removeSTmean.R
##' @author Johan Lindström
##' @family STdata functions
##' @export
removeSTcovarMean <- function(STdata){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")
    
  if( is.null(STdata$SpatioTemporal) ){
    ##no ST-covariates
    return(STdata)
  }
    
  ##compute mean at each locaiont for spatio-temporal covariate
  ST.mean <- apply(STdata$SpatioTemporal, 2:3, mean)
  ##add names
  colnames(ST.mean) <- paste("mean", colnames(ST.mean), sep=".")
  ##expand  the average
  tmp <- array(rep(ST.mean, each=dim(STdata$SpatioTemporal)[1]),
               dim(STdata$SpatioTemporal))
  
  ##subtract the means from the spatio-temporal data array
  tmp <- STdata$SpatioTemporal - tmp
  ##create new array that combines mean 0 and original ST-covar
  new <- array(NA, c(dim(tmp)[1:2],2*dim(tmp)[3]))
  new[,,1:dim(tmp)[3]] <- STdata$SpatioTemporal
  new[,,(dim(tmp)[3]+1):(2*dim(tmp)[3])] <- tmp
  ##names for the new array
  rownames(new) <- rownames(STdata$SpatioTemporal)
  colnames(new) <- colnames(STdata$SpatioTemporal)
  dimnames(new)[[3]] <- c(dimnames(tmp)[[3]],
                          paste("mean.0", dimnames(tmp)[[3]], sep="."))
  ##use this as a new spatio-temporal covariate
  STdata$SpatioTemporal <- new
  
  ##add the mean to the LUR data, matching locations
  IND <- match(STdata$covars$ID, rownames(ST.mean))
  names <- c(colnames(STdata$covars), colnames(ST.mean))
  STdata$covars <- cbind(STdata$covars, ST.mean[IND,])
  colnames(STdata$covars) <- names

  return(STdata)
}## removeSTcovarMean

##' Updates/sets the temporal trend for \code{STdata} objects. It also checks
##' that the spatio-temporal covariate exists for all dates in the trend (mainly an
##' issue if \code{extra.dates!=NULL} adding additional times at which to do
##' predictions.
##'
##' @title Update Trend in \code{STdata} Object
##' 
##' @param STdata A \code{STdata} object, see \code{\link{mesa.data}}.
##' @param n.basis number of basis functions for the temporal trend
##' @param extra.dates Additional dates for which smooth trends should be computed.
##' @param ... Additional parameters passed to \code{\link{calcSmoothTrends}}.
##' 
##' @return Returns a modfied version of the input, with and added/altered
##'   smooth trend.
##' 
##' @example Rd_examples/Ex_updateSTdataTrend.R
##' 
##' @author Johan Lindström
##' @family STdata functions
##' @family SVD for missing data
##' @export
updateSTdataTrend <- function(STdata, n.basis=0, extra.dates=NULL, ...){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")
  
  ##is there an existing trend
  if( !is.null(STdata$trend) ){
    message("Replacing existing trend.")
  }
  
  if(n.basis==0){
    ##no temporal trend, just add dates for prediction
    STdata$trend <- data.frame( date=sort(unique( c(STdata$obs$date,extra.dates) )) )
  }else{
    STdata$trend <- calcSmoothTrends(STdata=STdata, n.basis=n.basis,
                                  extra.dates=extra.dates, ...)$trend
  }
  ##ensure that we have all dates (mainly checks that added dates are also in ST-covariate.
  trend.dates <- STdata$trend$date
  ST.dates <- convertCharToDate( rownames(STdata$SpatioTemporal) )
  if( any(!(trend.dates %in% ST.dates)) && !is.null(ST.dates) ){
    stop( paste("trend dates not found in rownames(STdata$SpatioTemporal):",
                paste(trend.dates[!(trend.dates %in% ST.dates)],
                      collapse=", ")) )
  }
  
  ##return modified object
  return(STdata)
}

##' Removes an estimated time-trend from the observations in a
##' \code{STdata} object. Returns a modifed \code{STdata} object with no trend;
##' the new object can be used to fit a simpler model.
##' 
##' Sometimes there is no apparent spatial structure to the time-trend
##' amplitude, or there is not enough identifiability in the data to properly
##' model the structure. In that case, it is possible, at least as a sensitivity
##' analysis, to de-trend the observations and run a model with a spatial field
##' for the intercept only (apart from the spatio-temporal residual field).
##' 
##' \code{detrendSTdata} will remove the trends from the observations, using
##' \code{STdata$trend}. 'method' is applied as \cr
##' \code{metod(STdata$obs$obs ~ F,...)} \cr
##' where \code{F} is the temporal trend from \code{STdata$trend} for each
##' observation; or as \cr
##' \code{metod(STdata$obs$obs ~ F*obs.region,...)} \cr
##' where \cr
##' \code{obs.region = factor(region[ match(STdata$obs$ID, STdata$covars$ID)])}.
##' allowing for different trends in different region (i.e. interaction between
##' the time trend(s) and region identifiers).
##' \cr
##' \code{ predict( method(...) )} is then subtracted from
##' \code{STdata$obs$obs}, detrending the data.
##' 
##' @title Removes Temporal Trend from Observations in a \code{STdata} Object
##' 
##' @param STdata A \code{STdata} object, see \code{\link{mesa.data}}.
##' @param region Vector of the same length and order as \code{STdata$covars$ID}.
##'   Indicates region(s) in which different trends are to be fitted and removed.
##' @param method Method for fitting the trend (set to \code{method=lm} if
##'   \code{is.null(method)}); should produce output that allows
##'   the use of \code{\link[stats:predict]{predict}}. Possible options include
##'   \code{\link[stats:lm]{lm}}, \code{\link[MASS:rlm]{rlm}}, or
##'   \code{\link[MASS:lqs]{lqs}}.
##' @param ... Additional parameters passed to \code{method}.
##' 
##' @return Returns a modfied version of the input, with detrended observations
##'   and some changes:
##'     \item{STdata$obs}{Has an additional column, \code{removed.trend}, with
##'                       the amount subtracted per observation.}
##'     \item{STdata$trend}{Is reduced to only the \code{date} column, indicating
##'                         a constant trend.}
##'     \item{STdata$old.trend}{The previous \code{STdata$trend}, which was used
##'                             for detrending.}
##'     \item{STdata$fit.trend}{The result of \code{method}; the trend component
##'                             removed for each observations can be obtained as
##'                             \code{predict(STdata$fit.trend)}. NOTE: Aditional
##'                             functions, such as \code{\link{createSTmodel}}, might
##'                             reorder \code{STdata$obs} implying that\cr
##'            \code{STdata$obs$removed.trend != predict(STdata$fit.trend)}.}
##' 
##' @example Rd_examples/Ex_detrendSTdata.R
##' 
##' @author Assaf P. Oron and Johan Lindström
##' @family STdata functions
##' @export
detrendSTdata <- function(STdata, region=NULL, method=lm, ...){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")
  
  ##default method
  if( is.null(method) ){
    method <- lm
  }
  ##method must be a function
  if( !is.function(method) ){
    stop("'method' must be a function, e.g. method=lm")
  }
  ##if no trend, we can't detrend
  if( is.null(STdata$trend) || dim(STdata$trend)[2]==1 ){
    warning("No trend present in 'STdata', unable to detrend!.")
    return(STdata)
  }
  ##check if data already detrended
  if( !is.null(STdata$old.trend) ){
    warning("Data already detrended once. Second detrending might be a bad idea.")
  }

  ##Should the removed trend depend on the location
  if( !is.null(region) ){
    if( length(region) != length(STdata$covars$ID) ){
      stop("length(region) not equal to number of locations in STdata$covars$ID")
    }
    obs.region <- factor(region[ match(STdata$obs$ID, STdata$covars$ID)])
  }else{
    obs.region <- factor(rep(1, length(STdata$obs$obs)))
  }

  ##locate the date column in the trend object
  I.date <- which( names(STdata$trend)=="date" )
  ##create a matrix with the trend for the observation locations
  trend.obs <- as.matrix(STdata$trend[ match(STdata$obs$date, STdata$trend$date), -I.date])

  ##if only one level, then don't do factor regression
  if( nlevels(obs.region)==1 ){
    fit.trend <- method( STdata$obs$obs ~ trend.obs, ...)
  }else{
    fit.trend <- method( STdata$obs$obs ~ trend.obs*obs.region, ...)
  }
  ##save old trend, and information about the fit.
  STdata$old.trend <- STdata$trend
  STdata$fit.trend <- fit.trend
  ##remove fitted trends
  STdata$obs$removed.trend <- predict(STdata$fit.trend)
  STdata$obs$obs <- STdata$obs$obs - STdata$obs$removed.trend
  ##replace trend with a constant (just the dates)
  STdata$trend <- data.frame(date=STdata$trend$date)
  
  ##return modified object
  return(STdata)
}## detrendSTdata
