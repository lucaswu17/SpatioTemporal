##load the data
data(mesa.model)

##Compute dimensions for the data structure
dim <- loglikeSTdim(mesa.model)

##Find out in which order parameters should be given
loglikeST(NULL, mesa.model)

##Let's create random vectors of values
x <- runif( dim$nparam.cov )
x.all <- runif( dim$nparam )

##Evaluate the log-likelihood for these values
loglikeST(x.all, mesa.model, "f")
loglikeST(x, mesa.model, "p")
loglikeST(x, mesa.model, "r")

\dontshow{
  ##check that profile and full give the same results
  data(est.mesa.model)
  x.all <- coef(est.mesa.model)$par
  x <- coef(est.mesa.model, "cov")$par
  if( abs(loglikeST(x.all, mesa.model, "f") -
          loglikeST(x, mesa.model, "p")) > 1e-8 ){
    stop("loglike: full and profile not equal")
  }
}

