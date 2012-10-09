\docType{data}
\name{MCMC.mesa.model}
\alias{MCMC.mesa.model}
\title{Example of a \code{mcmcSTmodel} Structure}
\format{A list with elements, see the return description in
  \code{\link{MCMC.STmodel}}.}
\source{
  Contains parametere estimates for the Spatio-Temporal
  model applied to monitoring data from the \strong{MESA
  Air} project, see Cohen et.al. (2009) and
  \code{\link{mesa.data.raw}} for details.
}
\description{
  The output from a Metropolis-Hastings algorithm,
  implemented in \code{\link{MCMC.STmodel}}), run for the
  model in \code{\link{mesa.model}}
}
\examples{
##load data
data(mesa.model)
##and results of estimation
data(est.mesa.model)

##strating point
x <- coef(est.mesa.model)
##Hessian, for use as proposal matrix
H <- est.mesa.model$res.best$hessian.all
\dontrun{
  ##run MCMC
  MCMC.mesa.model <- MCMC(mesa.model, x$par, N = 5000, Hessian.prop = H)
}
##lets load precomputed results instead
data(MCMC.mesa.model)

####TODO: NEEDS S3-update
##Examine the results
print(MCMC.mesa.model)

##and contens of result vector
names(MCMC.mesa.model)
 
##Summary
summary(MCMC.mesa.model)

##MCMC tracks for four of the parameters
par(mfrow=c(5,1),mar=c(2,2,2.5,.5))
plot(MCMC.mesa.model$l, ylab="", xlab="", type="l",
     main="log-likelihood")
for(i in c(4,9,13,15)){
  plot(MCMC.mesa.model$par[,i], ylab="", xlab="", type="l",
       main=colnames(MCMC.mesa.model$par)[i])
}

##And estimated densities for the log-covariance parameters.
##The red line is the approximate normal distribution given by
##the maximum-likelihood estimates, e.g. ML-estimate and standard 
##deviation from the observed information matrix.
par(mfrow=c(3,3),mar=c(4,4,2.5,.5))
for(i in 9:17){
  xd <- sort(unique(MCMC.mesa.model$par[,i]))
  yd <- dnorm(xd, mean=x$par[i],sd=x$sd[i])
  dens <- density(MCMC.mesa.model$par[,i])
  plot(dens, ylim=c(0,max(c(dens$y,yd))), main =
       colnames(MCMC.mesa.model$par)[i])
  lines(xd,yd,col=2)
}
}
\references{
  M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L.
  Curl, T.  Gould, D. Hardie, A. Ho, P. Kinney, T. V.
  Larson, P. D. Sampson, L.  Sheppard, K. D. Stukovsky, S.
  S. Swan, L. S. Liu, J. D. Kaufman. (2009) Approach to
  Estimating Participant Pollutant Exposures in the
  Multi-Ethnic Study of Atherosclerosis and Air Pollution
  (MESA Air). Environmental Science & Technology: 43(13),
  4687-4693.
}
\seealso{
  \code{\link{createSTmodel}} for creation of the
  originating \code{STmodel} object.

  Other example data: \code{\link{est.cv.mesa}},
  \code{\link{est.mesa.model}}, \code{\link{mesa.data}},
  \code{\link{mesa.data.raw}}, \code{\link{mesa.model}},
  \code{\link{pred.mesa.model}}
}
\keyword{datasets}

