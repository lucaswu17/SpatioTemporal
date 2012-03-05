###################################
## Functions for crossvalidation ##
###################################
##Function that creates three default schemes for cv,
createCV <- function(mesa.data.model, groups=10, min.dist=.1, random=FALSE,
                     subset=NA, option="all"){
  ##options - for internal MESA usage:
  ## all      (just divide over all sites)
  ## fixed    (only AQS and FIXED sites)
  ## comco    (only comco, keeping gradients together, ignores min.dist)
  ## snapshot (same as comco)
  ## home     (only home)
  option <- tolower(option) ##ensure lower case
  if( !(option %in% c("all","fixed","comco","snapshot","home")) )
    stop("Unknown option")
  
  if(option=="all"){
    Ind <- 1:dim(mesa.data.model$location)[1]
    if( all(!is.na(subset)) ){
      ##subset the data - first cast to character
      subset <- as.character(subset)
      ##check that all requested names exist in the colnames
      if( any(!(subset %in% mesa.data.model$location$ID)) )
        warning( sprintf("In 'create.data.matrix': Subset names not found: %s",
                         paste(subset[!(subset %in% mesa.data.model$location$ID)],
                               collapse=", ")) )
      Ind <- which(mesa.data.model$location$ID %in% subset)
    }
  }else if(option=="fixed"){
    Ind <- which(toupper(mesa.data.model$location$type) %in% c("AQS","FIXED"))
  }else if(option=="home"){
    Ind <- which(toupper(mesa.data.model$location$type) %in% c("HOME"))
  }else{ ##remaining options are comco/snapshot
    Ind <- which(toupper(mesa.data.model$location$type) %in% c("COMCO"))
  }
  if(length(Ind)==0)
    stop( paste("No sites found for option:",option) )
  
  if(!random){
    ##set the seed to ensure that repeated calls returns the same grouping
    if( exists(".Random.seed") )
      old.seed <- .Random.seed ##save current seed so it can be re-set
    else
      old.seed <- NULL
    set.seed(0,kind="Mersenne-Twister")
  }

  if( option %in% c("all","fixed","home") ){
    ##extract distance matrix and find stations closer than min.dist
    D <- mesa.data.model$dist[Ind,Ind]
    D[lower.tri(D,diag=TRUE)] <- NA
    Ind.coloc <- apply(D,1,function(x){which(x<min.dist)})
    Ind.G <- list()
    if( any(sapply(Ind.coloc,length)>0) ){
      ##find Ind.coloc's that contain colocated sites
      Ind.coloc.I <- unname(which(sapply(Ind.coloc,length)>0))
      i <- 1
      while( length(Ind.coloc.I)>0 ){
        ##extract the i:th colocated sites and add it to the grouping vector
        Ind.G[[i]] <- unique(c(Ind.coloc.I[1],Ind.coloc[[Ind.coloc.I[1]]]))
        ##drop the site from Ind.coloc.I
        Ind.coloc.I <- Ind.coloc.I[-1]
        ##see if this leads to additional matches
        ##find all Ind.coloc elemets that contain things from Ind.G[[i]]
        ##This is done both looking for matches
        Ind.tmp <- which(sapply(Ind.coloc,function(x)
                                any(x %in% Ind.G[[i]])))
        ##and adding the Ind.G[[i]] elements since the "diagonal"
        ##isn't found thtough matches
        Ind.tmp <- unique(c(Ind.tmp,Ind.G[[i]]))
        ##Finnaly let's only consider things not yet used.
        Ind.tmp <- Ind.tmp[Ind.tmp %in% Ind.coloc.I]
        ##Add and repeat until no more new elements match
        while( length(Ind.tmp)>0 ){
          Ind.G[[i]] <- unique(c(Ind.G[[i]], Ind.tmp,
                                 unlist(Ind.coloc[Ind.tmp])))
          Ind.coloc.I <- Ind.coloc.I[!(Ind.coloc.I %in% Ind.tmp)]
          Ind.tmp <- which(sapply(Ind.coloc,function(x)
                                  any(x %in% Ind.G[[i]])))
          Ind.tmp <- unique(c(Ind.tmp,Ind.G[[i]]))
          Ind.tmp <- Ind.tmp[Ind.tmp %in% Ind.coloc.I]
        }
        i <- i+1
      }
      for(i in 1:length(Ind.G))
        Ind.G[[i]] <- sort(Ind[Ind.G[[i]]])
    }
    ##extract points that are not colocated
    Ind.R <- Ind[!(Ind %in% unlist(Ind.G))]
  }else{ ##the snapshot sites
    ##drop starting non-numerical and then numerical characters
    ##from site ID, if the first character of whats left=="R" then
    ##this is a random site. -> seprate indecies into random and gradient
    ##random sites are on the format [A-Z]-[A-Z][0-9]R[0-9]
    ##so remove the leading characters and see if we find an R.
    Ind.R <- Ind[substr(sub("^[A-Z]+-[A-Z]+[0-9]+","",
                            mesa.data.model$location$ID[Ind]),1,1)=="R"]
    ##non random sites are the rest.
    Ind <- Ind[!(Ind %in% Ind.R)]
    ##now lets group the gradient sites
    grad.id <- unique(sub("[A-Z]+$","",
                          as.character(mesa.data.model$location$ID[Ind])))
    Ind.G <- list()
    for(i in 1:length(grad.id))
      Ind.G[[i]] <- Ind[grep(grad.id[i],mesa.data.model$location$ID[Ind])]
  }##if( option %in% c("all","fixed","home") ){...}else{...}

  ##we have no created an Ind.G and a Ind.R variable, lets permute and
  ##create the CV-indicator matrix
  if(groups==0){ ##leave one out CV, place each group seperately
    Ind.fin <- Ind.G
    if( length(Ind.R)!=0 )
      for(i in 1:length(Ind.R))
        Ind.fin[[i+length(Ind.G)]] <- Ind.R[i]
    groups <- length(Ind.fin)
  }else{    
    groups <- min(c(groups,length(Ind.R)+length(Ind.G)))
    ##randomly permute the sites
    Ind.R <- sample(Ind.R,length(Ind.R),replace=FALSE)
    if(length(Ind.G)!=0){
      Ind.G <- Ind.G[sample(length(Ind.G),length(Ind.G),replace=FALSE)]
    }
    ##collect sites to be in each group, first random sites
    Ind.R <- matrix(c(Ind.R,rep(NA,ceiling(length(Ind.R)/groups)*groups -
                                length(Ind.R))), nrow=groups)
    ##and then collect the grouped sites
    if(length(Ind.G)==0){
      Ind <- matrix(NA,groups,1)
    }else{
      Ind <- 1:length(Ind.G)
      Ind <- matrix(c(Ind,rep(NA,ceiling(length(Ind)/groups)*groups -
                              length(Ind))), nrow=groups)
      ##flips the matrix so that we have more groups at the bottom, since
      ##there are more random sites at the top.
      Ind <- Ind[seq(dim(Ind)[1],1,-1),,drop=FALSE]
    }
    ##compose a list, each element containing one cv-group
    Ind.fin <- list()
    for(i in 1:groups){
      Ind.fin[[i]] <- c(Ind.R[i,],unlist(Ind.G[Ind[i,]]))
      Ind.fin[[i]] <- Ind.fin[[i]][!is.na(Ind.fin[[i]])]
    }
  }##if(groups==0) else
  if(!random){
    ##re-set the random number generator
    if( !is.null(old.seed) )
      .Random.seed <- old.seed
  }
  if( any(sapply(Ind.fin,length)==0) )
    stop("Some CV-groups contain NO sites, try reducing the number of requested groups.")
  Ind.cv <- matrix(FALSE,dim(mesa.data.model$obs)[1],groups)
  for(i in 1:groups)
    Ind.cv[,i] <- mesa.data.model$obs$idx %in% Ind.fin[[i]]
  
  ##create a suitable cross-validation matrix, each colum
  ##marks the elements to drop for a specific CV.
  return(Ind.cv)
}##createCV

####################################

###Function that does the prediction part of the cross validation

####################################
predictCV <- function(par, mesa.data.model, Ind.cv, type="p",express=FALSE,Nmax=1000, silent=TRUE){
  ##ensure lower case
  type <- tolower(type) 
  ##first check that type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'predictCV': Unknown option for type. Should be: (f)ull, (r)eml, (p)rofile")
  ##ensure that Ind.cv is a matrix
  if(is.vector(Ind.cv))
    Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    N.CV.sets <- dim(Ind.cv)[2]
  }
  
  ##expand the parameter vector so that we have one
  ##vector for each of the cv-sets we want to try.
  if(is.vector(par))
    par <- as.matrix(par)
  if(dim(par)[2]==1)
    par <- par %*% matrix(1,1,N.CV.sets)
  if(dim(par)[2] != N.CV.sets)
    stop("In 'predictCV': Number of parameters does not match the number of cv-sets")

	
if (express) { # "express" version producing only a data frame with pred-obs pairs

	dout=NULL
	for (a in 1:N.CV.sets)
	{
	    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==a
    }else{
      Ind.current <- Ind.cv[,a]
    }
    ##create data matrices that contains observed
    mesa.data.obs <- drop.observations(mesa.data.model, Ind.current)
    ##and locations to be predicted (the "LUR" component needs to be modified)
    mesa.data.pred <- drop.observations(mesa.data.model, !Ind.current)	 
		LUR <- mesa.data.pred$X[[1]]
		if( length(mesa.data.pred$X) > 1){
		  for(j in 2:length(mesa.data.pred$X)) LUR <- cbind(LUR,mesa.data.pred$X[[j]])
		}
		mesa.data.pred$LUR=as.data.frame(LUR)

		tmpred=cond.expectation(par[,a],mesa.data.obs,mesa.data.pred,no.nugget=TRUE,only.obs=TRUE,pred.var=FALSE)

		dout=rbind(dout,data.frame(ID=mesa.data.pred$obs$ID,date=mesa.data.pred$obs$date,raw=mesa.data.pred$obs$obs,lur=tmpred$EX.mu,krig=tmpred$EX))
	}
	return(dout)  

} ############# /express	
	
  pred <- list()
  for(i in 1:N.CV.sets){
    if(!silent)
      print( sprintf("Predicting cv-set %d/%d",i,N.CV.sets) )
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- Ind.cv[,i]
    }
    ##create data matrices that contains observed
    mesa.data.obs <- drop.observations(mesa.data.model, Ind.current)
    ##and locations to be predicted
    mesa.data.pred <- drop.observations(mesa.data.model, !Ind.current)

    ##combine gegraphic covariates
    LUR <- mesa.data.pred$X[[1]]
    if( length(mesa.data.pred$X) > 1){
      for(j in 2:length(mesa.data.pred$X))
        LUR <- cbind(LUR,mesa.data.pred$X[[j]])
    }
    LUR <- LUR[,colnames(LUR)!="const",drop=FALSE]
    LUR <- LUR[,unique(colnames(LUR)),drop=FALSE]
    ##create a data structure
    mesa.data.pred <- list(location = mesa.data.pred$location,
                           LUR = as.data.frame(LUR),
                           trend = mesa.data.pred$trend,
                           obs = mesa.data.pred$obs,
                           SpatioTemp = mesa.data.pred$SpatioTemp.all)
    ##lets obtain the predection for this set
    pred[[i]] <- cond.expectation(par[,i], mesa.data.obs, mesa.data.pred,
                                  Nmax=Nmax, no.nugget=FALSE, 
                                  only.obs=FALSE, compute.beta=FALSE,
                                  pred.var=TRUE, pred.covar=FALSE,
                                  type=type)
  }##for(i in 1:N.CV.sets){
  if(dim(Ind.cv)[2]==1 || max(rowSums(Ind.cv))==1){
    ##check that each location is only predicted once
    ##construct a matrix matching mesa.data$obs$obs with the obs
    pred.obs <- matrix(NA,length(mesa.data.model$obs$obs),3)
    colnames(pred.obs) <- c("obs","pred","pred.var")
    pred.obs[,"obs"] <- mesa.data.model$obs$obs
    for(i in 1:N.CV.sets){
      if( dim(Ind.cv)[2]==1 ){
        Ind.current <- Ind.cv==i
      }else{
        Ind.current <- Ind.cv[,i]
      }
      pred.obs[Ind.current,"pred"] <- pred[[i]]$EX[pred[[i]]$I]
      pred.obs[Ind.current,"pred.var"] <- pred[[i]]$VX[pred[[i]]$I]
    }
    pred.by.cv <- NULL
  }else{
    ##construct a list where each component contains a n-by-3 matrix
    ##with observations, predictions and variance
    pred.by.cv <- list()
    for(i in 1:N.CV.sets){
      pred.by.cv[[i]] <- cbind(mesa.data.model$obs$obs[Ind.cv[,i]],
                               pred[[i]]$EX[pred[[i]]$I],
                               pred[[i]]$VX[pred[[i]]$I])
      colnames(pred.by.cv[[i]]) <- c("obs","pred","pred.var")
    }
    pred.obs <- NULL
  }##if(max(rowSums(Ind.cv))==1) ... else ...
  pred.all.res <- list()
  for(i in 1:N.CV.sets){
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- Ind.cv[,i]
    }
    pred.all.res[[i]] <- array(NA,c(dim(pred[[i]]$EX)[1],
                                    dim(pred[[i]]$EX)[2],3))
    pred.all.res[[i]][,,2] <- pred[[i]]$EX
    pred.all.res[[i]][,,3] <- pred[[i]]$VX
    dimnames(pred.all.res[[i]]) <- list(rownames(pred[[i]]$EX),
                                        colnames(pred[[i]]$EX),
                                        c("obs","pred","pred.var"))
    obs.tmp <- mesa.data.model$obs[Ind.current,]
    ind.tmp <- 
      (match(as.character(mesa.data.model$location$ID[obs.tmp$idx]),
             colnames(pred[[i]]$EX))-1)*dim(pred[[i]]$EX)[1] + 
               match(as.character(obs.tmp$date),rownames(pred[[i]]$EX))
    tmp <- pred.all.res[[i]][,,"obs"]
    tmp[ind.tmp] <- obs.tmp$obs
    pred.all.res[[i]][,,"obs"] <- tmp
  }
  return( list(pred.obs=pred.obs, pred.by.cv=pred.by.cv,
               pred.all=pred.all.res) )
}##predictCV


##Function that does the estimation part of the cross validation
##Function takes additional parameters that are used by the optimisation
estimateCV <- function(par.init=3, mesa.data.model, Ind.cv,
                       type="p", h=1e-3, diff.type=1, lower=-15,
                       upper=15, hessian.all=FALSE,
                       control=list(trace=3, maxit=1000)){
  ##cast the initial parameters into a matrix
  if(is.vector(par.init))
    par.init <- as.matrix(par.init)
  if(length(par.init)==1){
    ##since par.init has length one we need to construct a matrix
    ##of initialvalues
    nparam.cov <- loglike.dim(mesa.data.model)$nparam.cov
    par.init <- matrix(NA,nparam.cov,par.init)
    init <- seq(-5, 5, length.out=dim(par.init)[2])
    for(i in 1:dim(par.init)[2])
      par.init[,i] <- init[i]
  }
  ##ensure that Ind.cv is a matrix
  if(is.vector(Ind.cv))
    Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    N.CV.sets <- dim(Ind.cv)[2]
  }
  
  res <- list()
  for(i in 1:N.CV.sets){
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- Ind.cv[,i]
    }
    ##create data matrices that contains observations
    mesa.data.aux <- drop.observations(mesa.data.model,Ind.current)
    
    ##lets estimate parameters for this set
    if( !is.null(control$trace) && control$trace!=0){
      print( sprintf("In: 'estimateCV': Estimation of cv-set %d/%d",
                     i,dim(Ind.cv)[2]) )
    }
    res[[i]] <- fit.mesa.model(par.init, mesa.data.aux, type=type,
                               h=h, diff.type=diff.type, lower=lower,
                               upper=upper, hessian.all=hessian.all,
                               control=control)
  }##for(i in 1:dim(Ind.cv)[2])
  par <- matrix(NA, length(res[[1]]$res.best$par), N.CV.sets)
  res.all <- list()
  for(i in 1:N.CV.sets){
    par[,i] <- res[[i]]$res.best$par
    res.all[[i]] <- res[[i]]$res.best
  }
  rownames(par) <- names(res[[1]]$res.best$par)
  ##Return the estimated parameters from the cross-validation
  return( list(par=par, res.all=res.all) )
}##estimateCV

##function that drops observations from mesa.data.model
##and recomputes relevant parts, used in the cross-validation.
drop.observations <- function(mesa.data.model, Ind.cv){
  ##check that the size of the indicator is consistent.
  N <- dim(mesa.data.model$obs)[1]
  if((is.vector(Ind.cv) && length(Ind.cv)!=N) ||
     (!is.vector(Ind.cv) && any(dim(Ind.cv) != c(N,1))))
    stop("In 'drop.observations': Length of Ind.cv must match dim(mesa.data.model$obs)[1]; also Ind.cv must be a vector or a column matrix.")
  
  ##drop observations, temporal trends and spatio-temporal covars.
  mesa.data.model$obs <- mesa.data.model$obs[!Ind.cv,,drop=FALSE]
  mesa.data.model$F <- mesa.data.model$F[!Ind.cv,,drop=FALSE]
  mesa.data.model$SpatioTemp <-
    mesa.data.model$SpatioTemp[!Ind.cv,,drop=FALSE]
  
  ##find which locations that remain
  IND <- (1:dim(mesa.data.model$location)[1] %in%
          unique(mesa.data.model$obs$idx))
  ##store the old locations
  ID.old <- as.character(mesa.data.model$location$ID)
  ##start dropping stuff
  mesa.data.model$location <- mesa.data.model$location[IND,,drop=FALSE]
  for(i in 1:length(mesa.data.model$X))
    mesa.data.model$X[[i]] <- mesa.data.model$X[[i]][IND,,drop=FALSE]
  mesa.data.model$dist <- mesa.data.model$dist[IND,IND,drop=FALSE]
  mesa.data.model$SpatioTemp.all <-
    mesa.data.model$SpatioTemp.all[,IND,,drop=FALSE]
  ##recompute somethings
  mesa.data.model$dates <- sort(unique(mesa.data.model$obs$date))
  mesa.data.model$nt <- double(length(mesa.data.model$dates))
  for(i in c(1:length(mesa.data.model$dates)))
    mesa.data.model$nt[i] <- length(which(mesa.data.model$obs$date ==
                                          mesa.data.model$dates[i]))
  ##and recompute the indecies
  mesa.data.model$obs$idx <-
    match(ID.old[mesa.data.model$obs$idx],
          as.character(mesa.data.model$location$ID))
  return(mesa.data.model)
}
