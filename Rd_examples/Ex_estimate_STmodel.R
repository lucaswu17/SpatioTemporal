##load a model object
data(mesa.model)

##examine the model
print(mesa.model)

##Important dimensions of the model
dim <- loglikeSTdim(mesa.model)
print(dim)

##Set up initial parameter values for optimization
x.init <- cbind(rep(2,dim$nparam.cov), c(rep(c(1,-3),dim$m+1),-3))
##and add names to the initial values
rownames(x.init) <- loglikeSTnames(mesa.model, all=FALSE)
print(x.init)
\dontrun{
  ##estimate parameters
  ##This may take a while...
  est.mesa.model <- estimate(mesa.model, x.init, type="p", hessian.all=TRUE)
}
##Let's load precomputed results instead.
data(est.mesa.model)

##Optimisation status message
print(est.mesa.model)

##compare the estimated parameters for the two starting points
est.mesa.model$summary$par.all
##and values of the likelihood (and convergence info)
est.mesa.model$summary$status

##extract the estimated parameters and approximate uncertainties
x <- coef(est.mesa.model)

##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x$par, ylim=range(c( x$par-1.96*x$sd, x$par+1.96*x$sd )),
     xlab="", xaxt="n")
points(x$par - 1.96*x$sd, pch=3)
points(x$par + 1.96*x$sd, pch=3)
axis(1, 1:length(x$par), rownames(x), las=2)

\dontrun{
  ##fixed pars
  x.fixed <- coef(est.mesa.model)$par
  x.fixed[c(1,2,5:9)] <- NA
  est.fix <- estimate(mesa.model, x.init, x.fixed, type="p")
}
