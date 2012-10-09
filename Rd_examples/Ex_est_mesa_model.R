##load data
data(mesa.data)
data(mesa.model)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))
\dontrun{
  ##estimate parameters
  est.mesa.model <- estimate(mesa.model, x.init, type="p", hessian.all=TRUE)
}

##create model structure with ST-covariates
mesa.model.ST <- createSTmodel(mesa.data, LUR=mesa.model$LUR.list,
                               ST="lax.conc.1500")
mesa.data.0 <- removeSTcovarMean(mesa.data)
##add intercept to the mean model
LUR.ST0 <- c( list( c(all.vars(mesa.model$LUR.list[[1]]),
                      "mean.lax.conc.1500")), mesa.model$LUR.list[-1] )
mesa.model.ST.0 <- createSTmodel(mesa.data.0, LUR=LUR.ST0,
                                 ST="mean.0.lax.conc.1500")

\dontrun{
  ##estimate parameters
  est.mesa.model.ST <- estimate(mesa.model.ST, x.init, hessian.all=TRUE)
  est.mesa.model.ST.0 <- estimate(mesa.model.ST.0, x.init, hessian.all=TRUE)
}

##time consuming estimation, load pre-computed results instead
data(est.mesa.model)

#estimation results
print(est.mesa.model)

##extract the estimated parameters and approximate uncertainties
x <- coef(est.mesa.model)

##compare estimated parameters
##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x$par, ylim=range(c( x$par-1.96*x$sd, x$par+1.96*x$sd )),
     xlab="", xaxt="n")
points(x$par - 1.96*x$sd, pch=3)
points(x$par + 1.96*x$sd, pch=3)
##add axis labels
axis(1, 1:length(x$par), rownames(x), las=2)

##estimation results for the Spatio-temporal covariates
print(est.mesa.model.ST)
print(est.mesa.model.ST.0)

##extract the estimated parameters and approximate uncertainties
x.ST <- coef(est.mesa.model.ST)
x.ST0 <- coef(est.mesa.model.ST.0)

##plot the estimated parameters
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(c(1:5,7:19), x.ST$par, ylab="", xlab="", xaxt="n")
points(1:19, x.ST0$par, pch=3, col=2)
points(c(2:5,7:19), x$par, pch=4, col=3)
abline(h=0, col="grey")
axis(1, 1:length(x.ST0$par), rownames(x.ST0), las=2)
legend("bottomleft", col = c(3,2,1), pch = c(4,3,1), 
       legend=c("mesa.model", "mesa.model.ST.mean0", "mesa.model.ST"))
