###################################################
# Functions for analysing crossvalidation results #
###################################################

##Compute Naive predictions based on a few fixed sites
predictNaive <- function(mesa.data.model, location=NULL, type=NULL){
  ##naive predictions that provide something simple to compare the data to
  ##Tries a number of different approaches

  ##matrix holding the predictions
  pred <- matrix(NA,length(mesa.data.model$obs$obs),4)
  ##select sites to base the predictions on
  if( !is.null(location) ){
    m <- length(location)
    if( is.character(location) ){
      IND.loc <- which(as.character(mesa.data.model$location$ID) %in% location)
    }else{
      IND.loc <- as.double(location)
    }
    IND.loc <- IND.loc[ IND.loc <= dim(mesa.data.model$location)[1] ]
    IND.loc <- IND.loc[ IND.loc > 0 ]
    if( length(IND.loc) != m)
      warning("In 'predictNaive': Some 'location':s not found in data.")
  }else if( !is.null(type) ){
    IND.loc <- which(as.character(mesa.data.model$location$type) %in% type)
  }else{
    stop("In 'predictNaive': Both location and type are NULL, you must specify one.")
  }
  if( length(IND.loc)==0 )
    stop("In 'predictNaive': No sites to use for naive prediction found")
  ##indicator picking out the sites to use for naive prediction.
  IND <- (mesa.data.model$obs$idx %in% IND.loc)
  ##indicator for which is the closest fixed site
  IND.closest <- apply(mesa.data.model$dist[IND.loc,],2,order)
  
  ##1) fit temporal smooths to all the fixed sites
  pred[,1] <- mesa.data.model$F %*% lm(mesa.data.model$obs$obs[IND] ~
                                       mesa.data.model$F[IND,]-1)$coefficients
  
  ##2) use average of the fixed sites to predict
  for(t in mesa.data.model$dates)
    pred[mesa.data.model$obs$date==t,2] <-
      mean(mesa.data.model$obs$obs[IND & mesa.data.model$obs$date==t])
  
  ##3) fit temporal smooth to the closes fixed site to predict the data
  ##start by fitting the temporal to each of the fixed sites
  pred.tmp <- matrix(NA,length(mesa.data.model$obs$obs),length(IND.loc))
  for(i in 1:length(IND.loc)){
    IND.tmp <- mesa.data.model$obs$idx==IND.loc[i]
    pred.tmp[,i] <- mesa.data.model$F %*%
      lm(mesa.data.model$obs$obs[IND.tmp] ~
         mesa.data.model$F[IND.tmp,]-1)$coefficients
  }
  for(i in 1:dim(mesa.data.model$location)[1])
    pred[mesa.data.model$obs$idx==i,3] <-
      pred.tmp[mesa.data.model$obs$idx==i,IND.closest[1,i]]

  ##4) use observations from the closest available site
  pred.tmp <- matrix(NA,length(mesa.data.model$obs$obs),length(IND.loc))
  for(i in 1:length(IND.loc)){
    IND.tmp <- mesa.data.model$obs$idx==IND.loc[i]
    obs.tmp <- mesa.data.model$obs$obs[IND.tmp]
    pred.tmp[,i] <- obs.tmp[match(mesa.data.model$obs$date,
                                  mesa.data.model$obs$date[IND.tmp])]
  }
  for(i in 1:dim(mesa.data.model$location)[1])
    pred[mesa.data.model$obs$idx==i,4] <-
      apply(pred.tmp[mesa.data.model$obs$idx == i,
                     IND.closest[,i],drop=FALSE],1,
            function(x) switch(all(is.na(x))+1,x[min(which(!is.na(x)))],NA))
  
  ##attach names to the predict matric indicating which
  ##simple approach was used.
  colnames(pred) <- c("smooth.fixed", "avg.fixed",
                      "smooth.closest.fixed","closest.fixed")
  return( list(pred = pred, locations =
               as.character(mesa.data.model$location$ID[IND.loc])) )
}##predictNaive

##Computes summary statistics for the cross validation
##computes stats are RMSE, R2 and coverate of 95% CI:s
summaryStatsCV <- function(predCV, pred.naive=NULL, lta=FALSE,
                           by.date=FALSE, p=0.95, trans=NULL){
  if(p<=0 || p>=1)
    stop("In 'summaryStatsCV': 'p' not a probability.")
  ##first determine what to compute
  DATE <- NULL
  if(by.date){
    ##find all dates that have observations
    for(i in 1:length(predCV$pred.all))
      DATE <-
        c(as.Date(names( which( apply( !is.na(predCV$pred.all[[i]][,,"obs"]),
                                      1,any)))),DATE)
    DATE <- sort(unique(DATE))
  }
  ##total number of summaryStats to compute
  Nstat <- 1 + lta + length(DATE)
  if( !is.null(pred.naive) )
     Nstat <- Nstat + dim(pred.naive$pred)[2]

  ##allocate memory for each of these
  Stats <- matrix(NA,Nstat,3)
  ##add names to matrix
  colnames(Stats) <- c("RMSE","R2","coverage")
  Rn <- "obs"
  if(lta)
    Rn <- c(Rn,"lta")
  Rn <- c(Rn, as.character(DATE) )
  if( !is.null(pred.naive) )
    Rn <- c(Rn, colnames(pred.naive$pred) )
  rownames(Stats) <- Rn

  ##recompute p for a two-sided CI
  p.org <- p
  p <- (1+p)/2
  ##compute stats for raw observations
  if( !is.null(trans) ){
    if(trans==0){
      obs <- exp(predCV$pred.obs[,"obs"])
      pred <- exp(predCV$pred.obs[,"pred"])
    }else{
      obs <- predCV$pred.obs[,"obs"]^trans
      pred <- predCV$pred.obs[,"pred"]^trans
    }
  }else{
    obs <- predCV$pred.obs[,"obs"]
    pred <- predCV$pred.obs[,"pred"]
  }
  I <- !is.na(pred)
  Stats["obs","RMSE"] <- sqrt( mean( (obs[I]-pred[I])^2 ))
  Stats["obs","R2"] <- 1-Stats["obs","RMSE"]^2 / var(obs[I])
  res <- (predCV$pred.obs[,"obs"] - predCV$pred.obs[,"pred"])
  res.norm <- res / sqrt(predCV$pred.obs[,"pred.var"])
  Stats["obs","coverage"] <- mean(abs(res.norm[I]) < qnorm(p))

  ##compute stats for long term average
  if(lta){
    lta.tmp <- compute.ltaCV(predCV, trans)
    I <- apply(!is.na(lta.tmp),1,all)
    Stats["lta","RMSE"] <- sqrt( mean( (lta.tmp[I,"obs"] -
                                        lta.tmp[I,"pred"])^2 ))
    Stats["lta","R2"] <- 1-Stats["lta","RMSE"]^2 / var(lta.tmp[I,"obs"])
  }else{
    lta.tmp <- NULL
  }#if(lta) ... else ...

  ##compute stats for each date
  if( length(DATE)>0 ){
    ##first find all sites
    loc <- sapply(predCV$pred.all,function(x) dimnames(x)[[2]])
    loc <- sort(unique(unlist(loc)))
    for(i in 1:length(DATE)){
      cDate <- as.character(DATE[i])
      ##create a matrix of results
      tmp <- matrix(NA,length(loc),3)
      rownames(tmp) <- loc
      colnames(tmp) <- c("obs","pred","pred.var")
      ##extract observations for each date
      for(k in 1:length(predCV$pred.all))
        for(l in 1:dim(predCV$pred.all[[k]])[2]){
          tmp[dimnames(predCV$pred.all[[k]])[[2]][l],
              c("obs","pred","pred.var")] <-
            predCV$pred.all[[k]][cDate,l,c("obs","pred","pred.var")]
      }
      I <- apply(!is.na(tmp),1,all)
      tmp <- tmp[I,,drop=FALSE]
      ##compute t-statistic before any possible transform
      T <- (tmp[,"obs"] - tmp[,"pred"]) / sqrt(tmp[,"pred.var"])
      if( !is.null(trans) ){
        if(trans==0){
          tmp <- exp(tmp)
        }else{
          tmp <- tmp^trans
        }
      }
        
      ##compute stats
      Stats[cDate,"RMSE"] <- sqrt( mean( (tmp[,"obs"] - tmp[,"pred"])^2 ))
      Stats[cDate,"R2"] <- 1-Stats[cDate,"RMSE"]^2 / var(tmp[,"obs"])
      Stats[cDate,"coverage"] <- mean(abs(T) < qnorm(p))
    }
  }#if( length(DATE)>0 )

  ##compute modified R2 for the naive predictions
  if( !is.null(pred.naive) ){
    ##extract observations, and naive predictions
    if( !is.null(trans) ){
      if(trans==0){
        pred.naive$pred <- exp(pred.naive$pred)
        obs <- exp(predCV$pred.obs[,"obs"])
      }else{
        pred.naive$pred <- pred.naive$pred^trans
        obs <- predCV$pred.obs[,"obs"]^trans
      }
    }else{
      obs <- predCV$pred.obs[,"obs"]
    }
    I <- !is.na(pred)
    ##compute modified R2:s
    for(i in 1:dim(pred.naive$pred)[2]){
      Stats[colnames(pred.naive$pred)[i],"R2"] <-
        1 - Stats["obs","RMSE"]^2 /
          mean( (pred.naive$pred[I,i] - obs[I])^2 )
    }
  }
  ##ensure that R2 is >= 0
  Stats[,"R2"] <- pmax(Stats[,"R2"],0)
  
  ##return the computed stats
  return( list(Stats=Stats, res=res, res.norm=res.norm,
               lta=lta.tmp, p=p.org) )
}#summaryStatsCV

##Computes the long term average for each of the sites in the
##cross-validated predictions
compute.ltaCV <- function(predCV, trans=NULL){
  ##first find all sites
  loc <- sapply(predCV$pred.all,function(x) dimnames(x)[[2]])
  loc <- sort(unique(unlist(loc)))
  ##create a matrix of results
  lta <- matrix(NA,length(loc),2)
  rownames(lta) <- loc
  colnames(lta) <- c("obs","pred")
  ##compute long term averages
  for(i in 1:length(predCV$pred.all)){
    for(j in 1:dim(predCV$pred.all[[i]])[2]){
      pred <- predCV$pred.all[[i]][,j,"pred"]
      obs <- predCV$pred.all[[i]][,j,"obs"]
      if( !is.null(trans) ){
        if(trans==0){
          pred <- exp(pred)
          obs <- exp(obs)
        }else{
          pred <- pred^trans
          obs <- obs^trans
        }
      }
      ##find values actually observed
      I <- !is.na(obs)
      ##find location for these obs and preds
      loc <- colnames(predCV$pred.all[[i]])[j]
      ##warn if location previously computed
      if( any(!is.na(lta[loc,])) )
        warning(sprintf("In 'compute.ltaCV': %s computed before, location present in more than one cv-set, latest result used.",loc))
      ##compute lta obs and pred
      lta[loc,"pred"] <- mean(pred[I])
      lta[loc,"obs"] <- mean(obs[I])
    }
  }
  return(lta)
}##compute.ltaCV
