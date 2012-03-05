printMesaDataNbrObs <- function(mesa.data){
  ##add the ID field if missing
  if( is.null(mesa.data$obs$ID) )
    mesa.data$obs$ID <- mesa.data$location$ID[mesa.data$obs$idx]
  if( is.null(mesa.data$trend) ){
    NT <- 0
  }else{
    NT <- dim(mesa.data$trend)[1]
  }
  ##plot stats to console
  cat(sprintf("Nbr locations: %d (observed: %d)\nNbr time points: %d (observed: %d)\nNbr obs: %d\n",
              dim(mesa.data$location)[1], length(unique(mesa.data$obs$ID)),
              NT, length(unique(mesa.data$obs$date)),
              length(mesa.data$obs$obs)))
  if(NT==0){
    cat( sprintf("No trend\n") )
  }else{
    cat( sprintf("Trend dates: %s\n", paste(range(mesa.data$trend$date),
                                            collapse=" to ")) )
  }
  if( !is.null(mesa.data$obs$date) )
    cat( sprintf("Observed dates: %s\n", paste(range(mesa.data$obs$date),
                                               collapse=" to ")) )
  if( !is.null(mesa.data$location$type) ){
    cat( "All sites:")
    print( table(mesa.data$location$type) )
    cat( "Observed:")
    print( table(mesa.data$location$type[mesa.data$location$ID
                                         %in% unique(mesa.data$obs$ID)]) )
    for(i in  levels(mesa.data$location$type)){
      I <- (mesa.data$obs$ID %in%
            mesa.data$location$ID[mesa.data$location$type==i])
      cat( sprintf("For %s:\n",i) )
      if( sum(I)!=0 ){
        cat( sprintf("  Number of obs: %d\n", sum(I)) )
        cat( sprintf("  Dates: %s\n", paste(range(mesa.data$obs$date[I]),
                                            collapse=" to ")) )
      }else{
        cat( "  No observations\n")
      }
    }
  }
}

plotMonitoringLoc <- function(mesa.data, main="", col=NULL,
                          legend.loc="topleft", legend.names=NULL,
                          add=FALSE, obsplot=FALSE, pch=19, cex=.1,...){
  if( is.null(mesa.data$location$type) )
    ##add a constant type argument to enable plotting
    mesa.data$location$type <-
      factor( rep("Observations",dim(mesa.data$location)[1]) )
  if( is.null(mesa.data$obs$idx) )
    mesa.data$obs$idx <- match(as.character(mesa.data$obs$ID),
                               mesa.data$location$ID)
    
  ##find the types that have observations
  TYPE <- sort(unique(mesa.data$location$type[
                unique(mesa.data$obs$idx)]))

  if( is.null(col) )
    col=1:length(TYPE)
  if(length(col) != length(TYPE))
    stop( sprintf("In 'data.overview': Needs one colour for every one of %d types",length(TYPE)) )
  if( !is.null(legend.loc) ){
    if( is.null(legend.names) )
      legend.names <- as.character(TYPE)
    if(length(legend.names) != length(TYPE))
      stop( sprintf("In 'data.overview': Needs one name for every one of %d types",length(TYPE)) )
  }
  ##start plotting
  yvals=mesa.data$obs$idx
  ytext="Location ID"

  if(obsplot) #### Plotting observations rather than locations
 {
	yvals=mesa.data$obs$obs
	ytext="Observations (modeling scale)"
 }
  
 if(!add)
    plot(mesa.data$obs$date ,yvals, type="n",
         main=main,xlab="Date",ylab=ytext,...)
  j=1
  for(i in as.character(TYPE) ){
    Ind <- (mesa.data$obs$idx %in%
            which(mesa.data$location$type==i))
    points(mesa.data$obs$date[Ind], yvals[Ind],
           col=col[j], cex=cex, pch=pch,...)
    j=j+1
  }
  ##possibly add a legend
  if( !is.null(legend.loc) )
    legend(legend.loc, legend.names, pch=19, col=col, pt.cex=.5)
}
    
plotMesaData <- function(mesa.data, ID, add=FALSE, type="obs", lag.max=NULL){
  if( length(ID)!=1 )
    stop("In 'plotMesaData': Can only plot observation timeseries at a time.")
  if( !(type %in% c("obs","res","acf","pacf")) )
    stop("In 'plotMesaData': Unknown option for 'type' should be obs, res, or acf.")
  ##add the ID field if missing, allows calling using both data and data.model
  if( is.null(mesa.data$obs$ID) )
    mesa.data$obs$ID <- mesa.data$location$ID[mesa.data$obs$idx]
  if( !is.character(ID) )
    ID <- as.character(mesa.data$location$ID)[ID]
  ##pick out the site to use
  IND <- mesa.data$obs$ID==ID
  if( sum(IND)==0 )
    stop("In 'plotMesaData': No observations at that location.")
    
  date <- mesa.data$obs$date[IND]
  y <- mesa.data$obs$obs[IND]
  if( length(y) > dim(mesa.data$trend)[2]){
    if( dim(mesa.data$trend)[2]==1 ){ ##just a constant temporal trend
      y.p <- rep(mean(y,na.rm=TRUE),length(mesa.data$trend$date))
    }else{
      trend.tmp <- mesa.data$trend
      trend.tmp$date <- NULL
      if( is.null(mesa.data$F) ){
        F <- matrix(NA,length(date),dim(trend.tmp)[2])
        for(i in c(1:dim(F)[2]))
          F[,i] <- spline(x=mesa.data$trend$date, y=trend.tmp[,i],
                          xout=date)$y
      }else{
        F <- mesa.data$F[IND,,drop=FALSE]
        F <- F[,colnames(F)!="const",drop=FALSE]
      }
      y.fit <- lm(y~F)
      y.p <- as.matrix(trend.tmp) %*% y.fit$coefficients[-1] +
        y.fit$coefficients[1]
      if( type!="obs" ){
        y <- y.fit$residuals
        y.p <- NULL
      }
    }
  }else{
    stop("In 'plotMesaData': Unable to fit smooth at that location, insufficient number of observations")
  }

  if(type %in% c("acf","pacf")){
    date.exp <- seq(min(date),max(date),min(diff(date)))
    y.exp <- rep(NA,length(date.exp))
    names(y.exp) <- as.character(date.exp)
    y.exp[as.character(date)] <- y
    if(type == "acf"){
      acf(y,main=paste("ACF for",ID),ylab="",xlab="",lag.max=lag.max)
    }else{
      pacf(y,main=paste("PACF for",ID),ylab="",xlab="",lag.max=lag.max)
    }
  }else{
    if(add){
      points(date, y)
    }else{
      plot(date, y, ylab="",xlab="",
           xlim = range(c(date,mesa.data$trend$date),na.rm=TRUE),
           ylim = range(c(y,y.p),na.rm=TRUE),
           main=switch(type,"obs"=ID,"res"=paste("Residuals",ID)))
    }
    if(type=="obs"){
      lines(mesa.data$trend$date, y.p)
    }else{
      abline(h=0,col="grey")
    }
  }
}

plotPrediction <- function(cond.exp, ID, mesa.data=NA, add=FALSE, p=0.95,
                           col=c("black","red","grey"), lty=c(1,1),
                           pch=c(NA,NA),cex=c(1,1)){
  if( length(ID)!=1 )
    stop("In 'plotMesaData': Can only plot one timeseries at a time.")
  ##ensure that lty, pch, and cex are of length==2
  if( length(lty)==1 )
    lty = c(lty,lty)
  if( length(pch)==1 )
    pch = c(pch,pch)
  if( length(cex)==1 )
    cex = c(cex,cex)
  
  if( any(is.na(mesa.data)) ){
    if( !is.character(ID) )
      stop("In 'plotMesaData': If mesa.data=NULL, ID cannot be a number")
    y <- rep(NA, length(cond.exp$EX[,ID]))
  }else{
    if( !is.character(ID) )
      ID <- as.character(mesa.data$location$ID)[ID]
    ##add an idx field
    if( is.null(mesa.data$obs$idx) )
      mesa.data$obs$idx <- match(as.character(mesa.data$obs$ID),
                                 mesa.data$location$ID)
  
    ##pick out the site to use
    idx <- which(as.character(mesa.data$location$ID)==ID)
    IND <- (as.Date(rownames(cond.exp$EX)) %in%
            mesa.data$obs$date[mesa.data$obs$idx == idx])
    ##pick out variance and expectation
    y <- rep(NA,length(cond.exp$EX[,ID]))
    if( sum(IND)!=0 )
      y[IND] <- mesa.data$obs$obs[mesa.data$obs$idx==idx]
  }
  E <- cond.exp$EX[,ID]
  sd <- sqrt(pmax(cond.exp$VX[,ID],0))
  t <- as.Date(rownames(cond.exp$EX))
  
  ##compute the quantile
  q <- qnorm((1-p)/2, lower.tail=FALSE)
  
  ##plot the results
  if(!add){
    plot(t,E+q*sd, type="n", ylim=range(c(E+q*sd, E-q*sd, y, E),
                               na.rm=TRUE),
         ylab="",xlab="",main=ID)
  }
  ##Plot the polygon
  polygon(c(t,rev(t)), c(E+q*sd, rev(E-q*sd)), border=col[3], col=col[3])
  ##plot the predictions
  if( !is.na(lty[1]) )
    lines(t, E, col=col[1], lty=lty[1])
  if( !is.na(pch[1]) )
    points(t, E, col=col[1], pch=pch[1], cex=cex[1])
  ##and the observations
  if( !is.na(lty[2]) )
    lines(t, y, col=col[2], lty=lty[2])
  if( !is.na(pch[2]) )
    points(t, y, col=col[2], pch=pch[2], cex=cex[2])
}

plotCV <- function(pred.cv, ID, mesa.data=NA, add=FALSE, p=0.95,
                   col=c("black","red","grey"), lty=c(1,1),
                   pch=c(NA,NA),cex=c(1,1)){
  if( length(ID)!=1 )
    stop("In 'plotMesaData': Can only plot one timeseries at a time.")
  if( !is.character(ID) ){
    if( any(is.na(mesa.data)) )
      stop("In 'plotMesaData': mesa.data needed if ID not given by name.")
    ID <- as.character(mesa.data$location$ID[ID])
  }
  ##ensure that lty, pch, and cex are of length==2
  if( length(lty)==1 )
    lty = c(lty,lty)
  if( length(pch)==1 )
    pch = c(pch,pch)
  if( length(cex)==1 )
    cex = c(cex,cex)
  
  ##pick out the site to use
  IND <- sapply(pred.cv$pred.all, function(x) colnames(x)==ID)
  IND.cv <- which(sapply(IND,function(x) any(x)))[1]
  ##extract the predictions
  pred <- pred.cv$pred.all[[IND.cv]][,IND[[IND.cv]],]
  ##pick out the conditional expectation and the variance
  y <- pred[,"obs"]
  E <- pred[,"pred"]
  sd <- sqrt(pmax(pred[,"pred.var"],0))
  t <- as.Date(rownames(pred))
  ##compute the quantile
  q <- qnorm((1-p)/2, lower.tail=FALSE)
  ##plot the results
  if(!add){
    plot(t,E+q*sd,type="n",ylim=range(c(E+q*sd, E-q*sd,y),
                                na.rm=TRUE),
         ylab="",xlab="",main=ID)
  }
  ##The polygon
  polygon(c(t,rev(t)), c(E+q*sd, rev(E-q*sd)), border=col[3], col=col[3])
  ##plot the predictions
  if( !is.na(lty[1]) )
    lines(t, E, col=col[1], lty=lty[1])
  if( !is.na(pch[1]) )
    points(t, E, col=col[1], pch=pch[1], cex=cex[1])
  ##and the observations
  if( !is.na(lty[2]) )
    lines(t, y, col=col[2], lty=lty[2])
  if( !is.na(pch[2]) )
    points(t, y, col=col[2], pch=pch[2], cex=cex[2])
}

CVresiduals.qqnorm <- function(res, I.season=as.factor(rep("obs",length(res))),
                               I.type=NA, main="All Data", norm=FALSE, legend=TRUE,...){
  TYPE <- sort(unique(I.type))
  TYPE <- TYPE[!is.na(TYPE)]
  I.s <- as.double(I.season)
  qqnorm(res,col=I.s+1,xlab="",ylab="",
         main=main,...)
  if(norm)
    abline(0,1,lty=1)
  qqline(res,lty=2)
  if(legend){
    if(norm){
      legend("bottomright",c("abline(0,1)","all data",levels(I.season)),
             pch=c(NA,NA,rep(1,length(unique(I.s)))),
             lty=c(1,2,rep(NA,length(unique(I.s)))),
             col=c(1,1,unique(I.s)+1))
    }else{
      legend("bottomright",c("all data",levels(I.season)),
             pch=c(NA,rep(1,length(unique(I.s)))),
             lty=c(2,rep(NA,length(unique(I.s)))),
             col=c(1,unique(I.s)+1))
    }
  }
  if(length(TYPE)>1){
    for(j in 1:length(TYPE)){
      qqnorm(res[I.type==TYPE[j]],col=I.s[I.type==TYPE[j]]+1,
             xlab="",ylab="", main=paste(main,TYPE[j]))
      if(norm)
        abline(0,1,lty=1)
      qqline(res[I.type==TYPE[j]],lty=2)
    }
  }
}##CVresiduals.qqnorm

CVresiduals.scatter <- function(res, LUR, I.season =
                                as.factor(rep("obs",length(res))),
                                I.type=NA, xlab="", main="All Data",
                                legend=TRUE, df=10,...){
  ##temporarlly disable warnings
  warn.old <- options("warn")
  options(warn=-1)
  ##extract TYPE and seasons
  TYPE <- sort(unique(I.type))
  TYPE <- TYPE[!is.na(TYPE)]
  I.s <- as.double(I.season)
  ##scatterplots
  plot(LUR, res, col=I.s+1, xlab=xlab, ylab="residuals",
       main=main,...)
  tmp <- try(smooth.spline(LUR,res,df=df),silent=TRUE)
  if( class(tmp)!="try-error" )
    lines( tmp )
  if( length(unique(I.s))>1 ){
    for(i in unique(I.s)){
      tmp <- try(smooth.spline(LUR[I.s==i],res[I.s==i],df=df),silent=TRUE)
      if( class(tmp)!="try-error" )
        lines( tmp, col=i+1)
    }
  }
  if(legend)
    legend("bottomright",levels(I.season), pch=1, col=(unique(I.s)+1),
           bg="white",cex=.75)
  if(length(TYPE)>1){
    for(j in 1:length(TYPE)){
      plot(LUR[I.type==TYPE[j]], res[I.type==TYPE[j]],
           xlim=range(LUR,na.rm=TRUE),ylim=range(res,na.rm=TRUE),
           col=I.s[I.type==TYPE[j]]+1, xlab=xlab, ylab="residuals",
           main=paste(main,TYPE[j],sep=""))
      tmp <- try(smooth.spline(LUR[I.type==TYPE[j]], res[I.type==TYPE[j]],
                               df=df),silent=TRUE)
      if( class(tmp)!="try-error" )
        lines(tmp)
      if( length(unique(I.s[I.type==TYPE[j]]))>1 ){
        for(i in unique(I.s[I.type==TYPE[j]])){
          tmp <- try(smooth.spline(LUR[I.type==TYPE[j]][I.s[I.type==TYPE[j]]==i],
                                   res[I.type==TYPE[j]][I.s[I.type==TYPE[j]]==i],
                                   df=df),silent=TRUE)
          if( class(tmp)!="try-error" )
            lines(tmp, col=i+1)
        }
      }
    }
  }
  ##restore warning status
  options(warn.old)
}##CVresiduals.scatter
