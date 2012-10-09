##load data
data(mesa.model)

\dontrun{
  ##estimate parameters
  dim <- loglikeSTdim(mesa.model)
  x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))
  est.mesa.model <- estimate(mesa.model, x.init, type="p", hessian.all=TRUE)
  
  ##estimate parameters with different nuggets
  x.init <- c(coef(est.mesa.model)$par,0)
  mesa.model2 <- updateCovf(mesa.model, cov.nu=list(covf="exp",
                                          random.effect=FALSE, nugget="type"))
  est2.mesa.model <- estimate(mesa.model2, x.init, type="p")
}

##time consuming estimation, load pre-computed results instead
data(est.mesa.model)

#estimation results
print(est.mesa.model)
print(est2.mesa.model)

##extract the estimated parameters and approximate uncertainties
x <- coef(est.mesa.model, "cov")
x2 <- coef(est2.mesa.model, "cov")

##compare estimated parameters
##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x$par, xlim=c(1, max(length(x$par), length(x2$par))),
     ylim=range(c( x$par-1.96*x$sd, x$par+1.96*x$sd,
       x2$par-1.96*x2$sd, x2$par+1.96*x2$sd)),
     xlab="", xaxt="n")
points(x$par - 1.96*x$sd, pch=3)
points(x$par + 1.96*x$sd, pch=3)
##parameters from second estimation
points(x2$par, pch=1, col=2)
points(x2$par - 1.96*x2$sd, pch=3, col=2)
points(x2$par + 1.96*x2$sd, pch=3, col=2)
##add axis labels
axis(1, 1:length(x2$par), rownames(x2), las=2)
