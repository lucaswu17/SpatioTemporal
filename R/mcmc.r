#########################################################
# MCMC for the parameters either alpha+Psi, or just Psi #
#########################################################
run.MCMC <- function(x, mesa.data.model, type="f", N=1000,
                     Hessian.prop=NA, Sigma.prop=NA, silent=TRUE){
  #first ensure that type is lower case
  type <- tolower(type)
  #check if type is valid
  if( !(type %in% c("r","p","f")) )
    stop("In 'run.MCMC': Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")

  if( type=="p" ){
    warning("In 'run.MCMC': Profile likelihood not a proper likelihood, running MCMC on the full parameter set (setting type=\"f\")", immediate.=TRUE)
    type="f"
  }

  dimensions <- loglike.dim(mesa.data.model)
  if(length(x) != dimensions$nparam && length(x) != dimensions$nparam.cov)
    stop("In 'run.MCMC': Size missmatch for x, number of parameters is incorrect")
  if(length(x) == dimensions$nparam && type!="f"){
    #requested REML but provided full parameter
    #starting points -> truncate
    x <- x[(dimensions$nparam-dimensions$nparam.cov+1):dim(x)[1],]
  }else if(length(x) == dimensions$nparam.cov && type=="f"){
    #requested full but provided parameters for REML -> expand
    tmp <- cond.expectation(x, mesa.data.model, only.pars=TRUE,
                            type="p")$pars
    x <- c(tmp$gamma.E, tmp$alpha.E, x)
  }
  
  if( any(is.na(Sigma.prop)) ){ #no proposal given
    #see if we need to compute the Hessian, either due to missing proposal
    if( any(is.na(Hessian.prop)) ){
      if(!silent)
        print("In 'run.MCMC': Computing Hessian")
      Hessian.prop <- loglike.hessian(x, mesa.data.model, type=type)
    }else if( dim(Hessian.prop)[1]!=length(x) || dim(Hessian.prop)[2]!=length(x) ){
      #or due to size missmatch
      warning("In 'run.MCMC': Size missmatch, recomputing Hessian", immediate.=TRUE)
      Hessian.prop <- loglike.hessian(x, mesa.data.model, type=type)
    }
    #compute proposal matrix from the Hessian ( inv(H)*2.38/d )
    #first check that Hessian is positive definite
    if( !all(eigen(Hessian.prop)$value < -1e-10) )
      stop("In 'run.MCMC': Hessian not negative definite (eigs < -1e-10)")
    Sigma.prop <- -solve(Hessian.prop)*2.38*2.38/dim(Hessian.prop)[1]
  }else{
    #check size of given proposal
    if( dim(Sigma.prop)[1]!=length(x) || dim(Sigma.prop)[2]!=length(x) )
      stop("In 'run.MCMC': Size missmatch, 'Sigma.prop' does not match parameter vector")
  }
  if( !all(eigen(Sigma.prop)$value > 1e-10) )
    stop("In 'run.MCMC': Proposal matrix not positive definite (eigs > 1e-10)")
  S <- t(chol(Sigma.prop))

  ##decide on steppping for the display
  disp.step <- min(max(round(N/50),1),100)
  #allocate storage
  l <- double(N)
  alpha <- double(N)
  par <- matrix(NA,N,length(x))
  #initial values
  par[1,] <- x
  l[1] <- loglike(par[1,], mesa.data.model, type=type)
  #run MCMC
  for(i in 2:N){
    #propose new value
    par.new <- par[i-1,] + S %*% rnorm(dim(S)[1])
    #evaluate likelihood for the new value
    l[i] <- loglike(par.new, mesa.data.model, type=type)
    alpha[i] <- min(exp(l[i]-l[i-1]),1)
    if(runif(1)<alpha[i]){ #accept
      par[i,] <- par.new
    }else{#reject
      par[i,] <- par[i-1,]
      l[i] <- l[i-1]
    }
    if(!silent && (i %% disp.step)==0)
      print( sprintf("MCMC %d/%d",i,N) )
  }#for(i in 2:N)
  ##add names to the results
  colnames(par) <- loglike.var.names(mesa.data.model,
                                     all = (type=="f"))
  colnames(S) <- rownames(S) <- colnames(par)
  ##return the results
  return( list(par=par,log.like=l,acceptance=alpha,chol.prop=S) )
}#run.MCMC
