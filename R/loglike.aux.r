###############################################
## FUNCTIONS THAT DETERMINES NAMES AND SIZES ##
###############################################
##figure out the dimensions of a number of interesting objects.
loglike.dim <- function(mesa.data.model){
  ##nbr of time points
  T <- as.integer( length(mesa.data.model$dates) )
  ##nbr of temporal basis fnct. incld. the constant
  m <- as.integer( dim(mesa.data.model$F)[2] )
  ##nbr of sites
  n <- as.integer( dim(mesa.data.model$X[[1]])[1] )
  ##nbr land use regressions for each temporal basis fnct.
  p <- integer(m) 
  for(i in c(1:m)){
    p[i] <- as.integer( dim(mesa.data.model$X[[i]])[2] )
  }
  ##number of spatio-temporal covariates
  if( any(is.null(mesa.data.model$SpatioTemp)) ){
    ##No model output
    L <- as.integer( 0 )
  }else{
    L <- as.integer( dim(mesa.data.model$SpatioTemp)[2] )
  }
  return( list(T=T,m=m,n=n,p=p,L=L,nparam=L+sum(p)+2*m+3,nparam.cov=2*m+3) )
}

##names of the variabels expected as input in x
loglike.var.names <- function(mesa.data.model, all=TRUE){
  dimensions <- loglike.dim(mesa.data.model)
  ##first figure out total number of variables
  if(all)
    N <- dimensions$nparam
  else
    N <-  dimensions$nparam.cov
  ##create a blank vector of the right length
  x <- double(N)
  ##add appropriate names to said vector
  Ind <- 1
  if(all){
    if(dimensions$L!=0){ ##regression coefficients for the model
      for(i in c(1:dimensions$L)){
        names(x)[Ind] <- paste("gamma",
                               colnames(mesa.data.model$SpatioTemp)[i],sep=".")
        Ind <- Ind + 1
      }
    }
    for(i in c(1:dimensions$m)){
      for(j in c(1:dimensions$p[i])){
        names(x)[Ind] <- paste("alpha",names(mesa.data.model$X)[i],
                               colnames(mesa.data.model$X[[i]])[j],sep=".")
        Ind <- Ind + 1
      }
    }
  }
  for(i in c(1:dimensions$m)){
    names(x)[Ind] <- paste("log.range",names(mesa.data.model$X)[i],sep=".")
    names(x)[Ind+1] <- paste("log.sill",names(mesa.data.model$X)[i],sep=".")
    Ind <- Ind + 2
  }
  names(x)[Ind] <- "log.range.nu"
  names(x)[Ind+1] <- "log.sill.nu"
  names(x)[Ind+2] <- "log.nugget.nu"
  return(names(x))
}

##INTERNALS:
##########################################################
## AUXILIRARY FUNCTIONS FOR THE LOGLIKELIHOOD FUNCTIONS ##
##########################################################
##compute the expected spatial values as mu_i(s) = X_i(s)*alpha_i
calc.mu.B <- function(dim,X,alpha){
  ##calculate the land use regression for the temporal trends
  mu.B <- double(dim$n * dim$m)
  dim(mu.B) <- c(dim$n, dim$m)
  for( i in c(1:dim$m) ) ##loop over diff. trends
    mu.B[,i] <- X[[i]] %*% alpha[[i]]
  return(mu.B)
}

##calculate block-matrix %*% X (LUR basis)
calc.iS.X <- function(X,iS, dim){
  iS.X <- matrix(0, dim(iS)[1], sum(dim$p))
  Ind <- c(0,0)
  ##populate the matrices, taking the block structure into account
  for(i in 1:dim$m){
    iS.X[(Ind[1]+1):(Ind[1]+dim$n),
           (Ind[2]+1):(Ind[2]+dim$p[i])] <-
      iS[(Ind[1]+1):(Ind[1]+dim$n),
             (Ind[1]+1):(Ind[1]+dim$n)] %*% X[[i]]
    ##increase the indecies
    Ind[1] <- Ind[1]+dim$n
    Ind[2] <- Ind[2]+dim$p[i]
  }
  return(iS.X)
}

##calculate X' %*% block-matrix %*% X (LUR basis)
calc.X.iS.X <- function(X, iS.X, dim){
  X.iS.X <- matrix(0, sum(dim$p), sum(dim$p))
  Ind <- c(0,0)
  ##populate the matrices, taking the block structure into account
  for(i in 1:dim$m){
    X.iS.X[(Ind[2]+1):(Ind[2]+dim$p[i]),
           (Ind[2]+1):(Ind[2]+dim$p[i])] <- t(X[[i]]) %*%
             iS.X[(Ind[1]+1):(Ind[1]+dim$n),
                  (Ind[2]+1):(Ind[2]+dim$p[i])]
    ##increase the indecies
    Ind[1] <- Ind[1]+dim$n
    Ind[2] <- Ind[2]+dim$p[i]
  }
  return(X.iS.X)
}
