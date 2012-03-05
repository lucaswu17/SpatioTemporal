############################
# ALL AUXILIRARY FUNCTIONS #
############################

### Create the basic data structure, from which 'create.data.model' works

setupSTdataset=function(rawobs,covardat,covarnames,trendf,varnames=list(yraw="lac",date="intended_wednesday",idobs="site_id",idcov="site.id",xcoord="lambert.x",ycoord="lambert.y",long="longitude",lat="latitude"),x.to.km=1000,transform=log,scale=TRUE,mesa=TRUE)

# rawobs: data frame with the raw observations to be modeled. Should include columns with location ID and date
# covardat: data frame with model covariates. Should also include location ID, coordinates, and (preferably) lat/long. 
# covarnames: character vector with the names of the covariates in `covardat' that might potentially be used in the model.
# trendf: matrix with the ``skeleton'' time trends used in the model in columns. Row names should be dates.
# varnames: list with the variable names to identify locations, observations, date, coordinates, etc.
# x.to.km: conversion factor from the x coordinates to km. Set to 1 to skip any conversion.
# transform: link function to transform the raw observations to the modeling scale. Note that `trendf' is already assumed to be in the modeling scale.
# scale: logical, should the covariate matrix be scaled s.t. each covariate has mean 0 and variance 1? Recommended when building the model.
# mesa: logical, do site IDs follow the MESA-Air naming conventions? If true, a `type' column will be created in the location data frame. 
{

# Our location and LUR data frames only include locations with observations
covardat=covardat[covardat[,varnames$idcov] %in% rawobs[,varnames$idobs],]

lout=list(location=data.frame(ID=covardat[,varnames$idcov],x=covardat[,varnames$xcoord]/x.to.km,y=covardat[,varnames$ycoord]/x.to.km,long=covardat[,varnames$long],lat=covardat[,varnames$lat]))

if (mesa)
{
	lout$location$type=rep("HOME",dim(lout$location)[1])
	lout$location$type[grep("[-]",lout$location$ID)]="COMCO"
	lout$location$type[grep("^[0-9]",lout$location$ID)]="AQS"
	lout$location$type[grep("[A-Z]00",lout$location$ID)]="FIXED"
}
lout$LUR=covardat[,covarnames]
row.names(lout$LUR)=covardat[,varnames$idcov]
if(scale) lout$LUR=as.data.frame(scale(lout$LUR))

lout$obs=data.frame(ID=rawobs[,varnames$idobs],date=rawobs[,varnames$date],obs=transform(rawobs[,varnames$yraw]))

lout$trend=as.data.frame(trendf)
trendim=dim(lout$trend)

names(lout$trend)=paste("V",1:(trendim[2]),sep="")
lout$trend$date=as.Date(row.names(lout$trend))

### If there are dates in data that don't have SEOF values, we will interpolate them:
if(any(!(lout$obs$date %in% lout$trend$date)))
{
   addates=unique(lout$obs$date[!(lout$obs$date %in% lout$trend$date)])

   cat(length(addates)," data dates don't have trend curves! Interpolating...\n")
   if (min(addates)<min(lout$trend$date) || max(addates)>max(lout$trend$date)) cat("Unmatched dates start/end beyond Trend-Function range!\n")

   newtrend=as.data.frame(matrix(NA,nrow=length(addates),ncol=trendim[2]))
   for (a in 1:(trendim[2]))
   {
       newtrend[,a]=predict(smooth.spline(as.numeric(lout$trend$date),lout$trend[,a]),x=as.numeric(addates))$y
   }
   newtrend$date=addates
   names(newtrend)=names(lout$trend)
	row.names(newtrend)=addates
   lout$trend=rbind(lout$trend,newtrend)
}

lout
} #### \end{setupSTdataset}


#########################################

###### Convert a 'mesa.data' data structure to a similar structure with the observations detrended, using the ###### structure's "trend" component. The returned structure has no trend component.
detrend.data=function(mesa.data,subregion=NA,method="rlm")

{
require(MASS)


nt=dim(mesa.data$trend)[2]-1
obsdim=dim(mesa.data$obs)
dout=mesa.data

trendmatch=as.matrix(mesa.data$trend[match(mesa.data$obs$date,mesa.data$trend$date),1:nt])

if (method=="rlm")
{
	trendfit=rlm(mesa.data$obs$obs~trendmatch)

	if(!is.na(subregion[1]))
	{
		if (length(subregion)!=length(mesa.data$location$ID)) stop("Subregion IDs should be a vector of same length and order as # of locations.")
		obsregion=subregion[match(mesa.data$obs$ID,mesa.data$location$ID)]
		trendfit=rlm(mesa.data$obs$obs~trendmatch*factor(obsregion))
	}
	dout$obs$obs=mesa.data$obs$obs-predict(trendfit)
}

dout$trend=data.frame(date=dout$trend$date)
dout$oldtrend=mesa.data$trend
dout$obs$removedtrend=predict(trendfit)

return(dout)
}  ############## /detrend.data


## Set up data for optimisation, precomputes model data
create.data.model <- function(mesa.data, LUR=NA, ST.Ind=NA, strip=TRUE,
                              strip.loc=strip, strip.time=strip){
  LUR <- default.LUR.list(mesa.data, LUR)
  ST.Ind <- default.ST.list(mesa.data, ST.Ind)

  ##construct a matrix holding the temporal basis for each observation
  ##(including the constant)
  N <- length(mesa.data$obs$obs)
  if(N==0 && any(c(strip.loc,strip.time)))
    stop("In 'create.data.model': strip.loc=TRUE and/or strip.time=TRUE makes no sense when mesa.data contains no observations.")
  M <- dim(mesa.data$trend)[2]
  F <- matrix(0,N,M)
  F[,1] <- 1
  ##temporal trends
  trend.temp <- mesa.data$trend
  ##drop the date column
  trend.temp$date <- NULL
  if(M>1 && N>0){
    for(i in c(2:M) )
      F[,i] <- spline(x=mesa.data$trend$date, y=trend.temp[,i-1],
                      xout=mesa.data$obs$date)$y
  }
  ##check for extrapolation
  if( any(is.na(F)) )
    stop("In 'create.data.model': mesa.data$obs$date contains elements outside the range of mesa.data$trend$date.")
  
  F <- as.data.frame(F)
  names(F)[1] <- "const"
  if(M>1)
    names(F)[2:M] <- names(trend.temp)
  ##we need F to be a matrix, to support matrix calculations
  F <- data.matrix(F)
  rownames(F) <- as.character(mesa.data$obs$date)

  ##extract and interpolate the trend
  date <- sort(unique(c(mesa.data$trend$date, mesa.data$obs$date)))
  trend <- mesa.data$trend
  if( length(date) > dim(trend)[1]){
    ##warn for interpolation
    interp.date <- mesa.data$obs$date[!(mesa.data$obs$date %in%
                                        mesa.data$trend$date)]
    interp.date <- paste(as.character(unique(interp.date)), collapse=", ")
    warning(paste("In 'create.data.model': The temporal trend will been interpolated to the following observations dates:", interp.date),immediate. = TRUE)
    ##do the interpolation
    trend[(dim(trend)[1]+1):length(date),] <- NA
    trend$date <- date
    ##reset rownames
    rownames(trend) <- NULL
    ##interpolate to timepoints in trend+obs
    if(M>1){
      for(i in names(trend)[names(trend)!="date"] )
        trend[,i] <- spline(x=mesa.data$trend$date, y=mesa.data$trend[,i],
                            xout=trend$date)$y
    }
  }
  if(strip.time){
    ##drop time points that are unobserved
    IND <- trend$date %in% unique(mesa.data$obs$date)
    trend <- trend[IND,,drop=FALSE]
  }
  
  ##construct matrices of basis vectors for the landuse regression
  X <- construct.LUR.basis(mesa.data, LUR)
  
  ##extract the observations and locations
  location <- mesa.data$location
  IND <- (as.character(location$ID) %in%
          as.character(unique(mesa.data$obs$ID)))
  ##test for colocated monitoring sites
  test.coloc <- location[IND,c("x","y")]
  if( any(duplicated(test.coloc)) )
    stop("In 'create.data.model': Model does not allow for colocated monitoring sites. Merge collocated sites.")
  
  ##drop locations that don't have observations
  if(strip.loc){
    location <- location[IND,]
  }else{
    ##test for colocated sites, this time for non-monitoring sites
    test.coloc <- location[,c("x","y")]
    if( any(duplicated(test.coloc)) )
      stop("In 'create.data.model': Model does not allow for colocated sites. Merge collocated sites, or use strip.loc=TRUE.")
  }
  
  location$ID <- as.character(location$ID)
  ##calculate distance matrix
  D <- as.matrix(dist( cbind(location$x,
                             location$y) ))
  colnames(D) <- rownames(D) <- location$ID
  ##order the covariates so they match the ordering of the locations
  for(i in 1:length(X)){
    IND <- match(location$ID, rownames(X[[i]]))
    X[[i]] <- X[[i]][IND,,drop=FALSE]
  }
  if( any(is.na(IND)) )
    stop("In 'create.data.model': Unable to find LUR-covariates for some locations. Names in mesa.data$location$ID has to match those in rownames(mesa.data$LUR).")
  
  ##create a matrix with the spatio-temporal covariate
  ST <- construct.ST.basis(mesa.data, ST.Ind)
  ##and an array for all the space-time locations
  if( length(ST.Ind)==0 ){
    ST.all <- NULL
  }else{
    ST.all <- mesa.data$SpatioTemp[,,ST.Ind,drop=FALSE]
    ##and order the locations
    IND <- match(location$ID,colnames(ST.all))
    ST.all <- ST.all[,IND,,drop=FALSE]
    if( any(is.na(IND)) )
      stop("In 'create.data.model': Unable to find spatio-temporal covariates for some locations. Names in mesa.data$location$ID has to match those in colnames(mesa.data$SpatioTemp).")
    ##and extract the relevant time points (recall that trend$date
    ##allready contains the obs$date
    IND <- match(trend$date, as.Date(rownames(ST.all)))
    ST.all <- ST.all[IND,,,drop=FALSE]
    if( any(is.na(IND)) )
      stop("In 'create.data.model': Unable to find spatio-temporal covariates for some dates. Dates in unique(c(mesa.data$trend$date,mesa.data$obs$date)) has to match those in as.Date(rownames(mesa.data$SpatioTemp)).")
  }
  
  ##extract the observations
  obs <- mesa.data$obs[,c("obs","date","ID")]
  if( !is.null(obs) ){
    obs$ID <- as.character(obs$ID)
    obs$idx <- match(obs$ID, location$ID)
    if( any(is.na(obs$idx)) ){
      stop("In 'create.data.model': Unable to match some observation locations (mesa.data$obs$ID) to location names (mesa.data$location$ID).")
    }
    ##sort the observations
    IND <- order(obs$date,obs$idx)
    obs <- obs[IND,,drop=FALSE]
    ##make sure that F and ST matches
    F <- F[IND,,drop=FALSE]
    ST <- ST[IND,,drop=FALSE]
    ##test for dubplicate observations
    if( any(duplicated(obs[,c("date","idx")])) )
      warning("In 'create.data.model': Some space-time location(s) have multiple observations. This is not within the model framework; some functions may give bogus results.",immediate. = TRUE)
  }
  ##count number of observations at each time point
  if( N==0 ){
    dates <- NULL
    nt <- double(length(dates))
  }else{
    dates <- sort(unique(mesa.data$obs$date))
    nt <- double(length(dates))
    for(i in c(1:length(dates)))
      nt[i] <- length(which(mesa.data$obs$date==dates[i]))
  }
  
  ##check if
  return( list(obs=obs, location=location, trend=trend, F=F, X=X,
               dist=D, dates=dates, nt=nt, SpatioTemp=ST,
               SpatioTemp.all=ST.all, LUR.list=LUR, ST.Ind=ST.Ind) )
}

##create a data matrix (T-by-N) from the data structure
create.data.matrix <- function(mesa.data=NA, obs=mesa.data$obs$obs,
                               date=mesa.data$obs$date, ID=mesa.data$obs$ID,
                               subset=NA){
  if( any(length(obs)!=c(length(date),length(ID))) )
    stop("In 'create.data.matrix': 'obs', 'date' and 'ID' are of unequal length")
  ##cast ID to character
  ID <- as.character(ID)
  ##find unique indecies
  date.ind <- sort(unique(date))
  ID.ind <- sort(unique(ID))
  ##construct the data matrix
  data <- matrix(NA,length(date.ind),length(ID.ind))
  I <- (match(ID,ID.ind)-1)*length(date.ind) + match(date,date.ind)
  data[I] <- obs
  colnames(data) <- ID.ind
  rownames(data) <- as.character(date.ind)
  if( all(!is.na(subset)) ){
    ##subset the data - first cast to character
    subset <- as.character(subset)
    ##check that all requested names exist in the colnames
    if( any(!(subset %in% colnames(data))) )
      warning( sprintf("In 'create.data.matrix': Subset names not found in colnames(data): %s", paste(subset[!(subset %in% colnames(data))],collapse=", ")) )
    ##convert from a character vector to a index vector
    subset <- which(colnames(data) %in% subset)
    data <- data[,subset,drop=FALSE]
  }
  return(data)
}

## Remove mean from caline and add it to the LUR variables #
remove.ST.mean <- function(mesa.data){
  if( !is.null(mesa.data$SpatioTemp) ){
    ##compute mean at each locaiont for spatio-temporal covariate
    ST.mean <- apply(mesa.data$SpatioTemp, 2:3, mean)
    ##add names
    colnames(ST.mean) <- paste("mean",colnames(ST.mean),sep=".")
    ##expand  the average
    tmp <- array(rep(ST.mean,each=dim(mesa.data$SpatioTemp)[1]),
                 dim(mesa.data$SpatioTemp))
    ##subtract the means from the spatio-temporal data array
    mesa.data$SpatioTemp <- mesa.data$SpatioTemp - tmp
    ##add the mean to the LUR data, matching locations
    IND <- match(rownames(mesa.data$LUR),rownames(ST.mean))
    names <- c(colnames(mesa.data$LUR),colnames(ST.mean))
    mesa.data$LUR <- cbind(mesa.data$LUR, ST.mean[IND,])
    colnames(mesa.data$LUR) <- names
  }
  return(mesa.data)
}

##extract parameters from x
get.params <- function(x, dimensions){
  if(length(x)!=dimensions$nparam && length(x)!=dimensions$nparam.cov)
    stop("In 'get.params': Parameter vector, 'x', has inconsistent size")
  Ind <- 0
  if(length(x)==dimensions$nparam){
    ##regression coefficients for model output
    if(dimensions$L!=0){ 
      gamma <- as.matrix(x[c((Ind+1):(Ind+dimensions$L))])
      dim(gamma) <- c(dimensions$L,1)
      Ind <- Ind + dimensions$L
    }else{
      gamma <- double(0)
    }
    ##land use regression coefficients for each of the temporal basis functions
    alpha <- vector("list",dimensions$m)
    for( i in c(1:dimensions$m)){
      alpha[[i]] <- matrix(x[c((Ind+1):(Ind+dimensions$p[i]))],ncol=1)
      Ind <- Ind + dimensions$p[i]
    }
  }else{
    alpha <- double(0)
    gamma <- double(0)
  }
  ##extract log range and sill for each of the residuals
  trend.range <- double(dimensions$m)
  trend.sill <- double(dimensions$m)
  for( i in c(1:dimensions$m)){
    trend.range[i] <- x[Ind+1]
    trend.sill[i] <- x[Ind+2]
    Ind <- Ind + 2
  }
  ##transform from log
  trend.range <- exp(trend.range)
  trend.sill <- exp(trend.sill)
  ##range, sill, and nugget for the residuals
  phi.residuals <- exp( x[c((Ind+1):length(x))] )

  return(list(gamma=gamma,alpha=alpha,
              range=trend.range,
              sill=trend.sill,
              phi.nu=phi.residuals))
}


##INTERNALS:

##################################################
# Functions that create default LUR and ST lists #
##################################################
default.LUR.list <- function(mesa.data,LUR=NA){
  if( is.null(mesa.data$LUR) )
    stop("In 'default.LUR.list': LUR object lacking from mesa.data, e.g is.null(mesa.data$LUR)==TRUE.")
  m <- dim(mesa.data$trend)[2] ##number of LURS, including intercept
  p <- dim(mesa.data$LUR)[2]
  if(is.null(LUR) || p==0){
    LUR <- list()
    for(i in 1:m)
      LUR[[i]] <- character(0)
    return(LUR)
  }
  if(!is.list(LUR) && any(is.na(LUR))){
    ##set the LUR indicator to use all land use regressors provided.
    LUR <- list()
    for(i in 1:m)
      LUR[[i]] <- c(1:p)
  }else if(!is.list(LUR)){
    ##LUR list is just a single vector replicate for use with all the trends.
    LUR.in <- LUR
    LUR <- list()
    for(i in 1:m)
        LUR[[i]] <- LUR.in
  }else if(length(LUR)!=m){
    ##add or drop extra elements and throw a warning
    warning("In 'default.LUR.list': Length of LUR list does not match number of temporal trends")
    LUR.in <- LUR
    LUR <- LUR.in[1:min(m,length(LUR.in))]
    if(length(LUR)<m)
      for(i in (length(LUR)+1):m)
        LUR[[i]] <- c(1:p)
  }
  ##convert LUR list to list of names
  for(i in 1:m){
    if(is.null(LUR[[i]])){
      LUR[[i]] <- character(0)
    }else if(!is.character(LUR[[i]])){
      if( any(LUR[[i]] > dim(mesa.data$LUR)[2]) )
        warning( sprintf("In 'default.LUR.list': Unable to find LUR with indecies larger than %d",dim(mesa.data$LUR)[2]) )
      LUR[[i]] <- LUR[[i]][LUR[[i]] <= dim(mesa.data$LUR)[2]]
      LUR[[i]] <- names(mesa.data$LUR)[LUR[[i]]]
    }
  }
  ##retain only matched names
  for(i in 1:m){
    if( !all(LUR[[i]] %in% names(mesa.data$LUR)) )
      warning("In 'default.LUR.list': Some requested names unmatched.")
    LUR[[i]] <- LUR[[i]][LUR[[i]] %in% names(mesa.data$LUR)]
    ##ensure no duplicates
    LUR[[i]] <- unique(LUR[[i]])
  }
  return(LUR)
}

default.ST.list <- function(mesa.data,ST.Ind=NA){
  if(is.null(mesa.data$SpatioTemp) || is.null(ST.Ind) ||
     dim(mesa.data$SpatioTemp)[3]==0)
    ST.Ind <- NULL
  else{
    if(any(is.na(ST.Ind)))
      ST.Ind <- c(1:dim(mesa.data$SpatioTemp)[3])
    ##transform ST.Ind to a character vector
    if(!is.character(ST.Ind)){
      if( any(ST.Ind > dim(mesa.data$SpatioTemp)[3]) )
      warning( sprintf("In 'default.ST.list': Unable to find SpatioTemp with indecies larger than %d.",dim(mesa.data$SpatioTemp)[3]) )
      ST.Ind <- ST.Ind[ST.Ind <= dim(mesa.data$SpatioTemp)[3]]
      ST.Ind <- dimnames(mesa.data$SpatioTemp)[[3]][ST.Ind]
    }
  }
  ##retain only matched names
  if( !all(ST.Ind %in% dimnames(mesa.data$SpatioTemp)[[3]]) )
    warning("In 'default.ST.list': Some requested names unmatched.")
  ST.Ind <- ST.Ind[ST.Ind %in% dimnames(mesa.data$SpatioTemp)[[3]]]
  ##ensure no duplicates
  ST.Ind <- unique(ST.Ind)
  return(ST.Ind)
}

##########################################
# Functions that create LUR and ST bases #
##########################################
construct.LUR.basis <- function(mesa.data, LUR){
  ##ensure safe parameters
  LUR <- default.LUR.list(mesa.data,LUR)
  ##construct matrices of basis vectors for the landuse regression
  X <- list()
  for(i in 1:length(LUR))
    X[[i]] <- mesa.data$LUR
  
  if(length(X)>1){
    trend.tmp <- mesa.data$trend
    trend.tmp$date <- NULL
    names(X)[2:length(X)] <- names(trend.tmp)
  }
  names(X)[1] <- "const"
  
  ##drop components and add intercept 
  for(i in 1:length(X)){
    if( is.null(LUR[[i]]) || length(LUR[[i]])==0 ){
      tmp <- as.matrix( data.frame("const"=rep(1,dim(X[[i]])[1])) )
      X[[i]] <- tmp
      rownames(X[[i]]) <- rownames(mesa.data$LUR)
    }else{
      X[[i]] <- cbind(matrix(1,ncol=1,nrow=dim(X[[i]])[1]),
                      data.matrix( X[[i]][,LUR[[i]],drop=FALSE] ) )
      colnames(X[[i]])[1] <- "const"
    }
  }
  return(X)
}

construct.ST.basis <- function(mesa.data,ST.Ind){
  ##ensure safe parameters
  ST.Ind <- default.ST.list(mesa.data,ST.Ind)
  ##find indecies for the names in ST.Ind
  ST.Ind <- match(ST.Ind, dimnames(mesa.data$SpatioTemp)[[3]])
  ##pick out spatio-temporal covariates.
  if( length(ST.Ind)==0 ){
    ST <- NULL
  }else{
    ##Create spatio-temporal matrix
    ST <- matrix(0, length(mesa.data$obs$obs), length(ST.Ind))

    ##dates and ID:s for the spatio-temporal covaraites
    ST.dates <- as.Date(rownames(mesa.data$SpatioTemp))
    ST.ID <- colnames(mesa.data$SpatioTemp)
    ##correspinding index
    I <- (match(mesa.data$obs$ID,ST.ID)-1)*length(ST.dates) +
      match(mesa.data$obs$date,ST.dates)
    ##extract these points from the Spatio temporal covariate
    for(i in 1:length(ST.Ind)){
      tmp <- mesa.data$SpatioTemp[,,ST.Ind[i]]
      ST[,i] <- tmp[I]
    }
    colnames(ST) <- dimnames(mesa.data$SpatioTemp)[[3]][ST.Ind]
    ST <- data.matrix(ST)
  }
  return(ST)
}

combineMesaData <- function(mesa.data.model, mesa.data){
  ##first drop duplicate locations
  IND <- !(as.character(mesa.data$location$ID) %in%
           as.character(mesa.data.model$location$ID))
  if( sum(IND)==0 ){
    warning("In 'combineMesaData': All locations in mesa.data already exist in mesa.data.model, returning mesa.data.model")
    return(mesa.data.model)
  }
  mesa.data$location <- mesa.data$location[IND,]
  ##Having dropped locations will drop corresponding LUR and ST elements in
  ##the 'create.data.model' call.
  
  ##drop observations in the mesa.data object
  mesa.data$obs <- NULL
  ##need the trend in mesa.data to match that oc mesa.data.model since
  ##these are the date that'll be used. Matching the trend also implies
  ##that extraction of ST.all elements will occur at the correct times.
  mesa.data$trend <- mesa.data.model$trend
  
  model.tmp <- create.data.model(mesa.data, LUR=mesa.data.model$LUR.list,
                                 ST.Ind=mesa.data.model$ST.Ind, strip=FALSE)
  ##combine the two data-sets
  mesa.data.model$location <- rbind(mesa.data.model$location,
                                    model.tmp$location)
  ##recompute distances
  mesa.data.model$dist <- as.matrix(dist(cbind(mesa.data.model$location$x,
                                               mesa.data.model$location$y) ))
  ##combine LUR
  for(i in 1:length(mesa.data.model$X))
    mesa.data.model$X[[i]] <- rbind(mesa.data.model$X[[i]],
                                    model.tmp$X[[i]])
  ##combine SpatioTemp
  if( !is.null(mesa.data.model$SpatioTemp.all) ){
    tmp <- array(NA,c(dim(mesa.data.model$SpatioTemp.all)[1],
                      dim(mesa.data.model$location)[1],
                      dim(mesa.data.model$SpatioTemp.all)[3]))
    dimnames(tmp) <- list(rownames(mesa.data.model$SpatioTemp.all),
                          as.character(mesa.data.model$location$ID),
                          dimnames(mesa.data.model$SpatioTemp.all)[[3]])
    tmp[,colnames(mesa.data.model$SpatioTemp.all),] <-
      mesa.data.model$SpatioTemp.all
    tmp[,model.tmp$location$ID,] <-
      mesa.data$SpatioTemp[rownames(mesa.data.model$SpatioTemp.all),
                           model.tmp$location$ID, dimnames(tmp)[[3]]]
    if( any(is.na(tmp)) )
      stop("In 'combineMesaData': Some SpatioTemporal values missing, probably from mesa.data. Ensure that mesa.data$SpatioTemp has all dates and locations needed.")
    mesa.data.model$SpatioTemp.all <- tmp
  }
  return(mesa.data.model)
}
