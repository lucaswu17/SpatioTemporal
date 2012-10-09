\docType{data}
\name{est.cv.mesa}
\alias{est.cv.mesa}
\alias{pred.cv.mesa}
\title{Example of \code{estCVSTmodel} and \code{predCVSTmodel} Structures}
\format{A list with elements, see the return description in
  \code{\link{estimateCV.STmodel}} and \code{\link{predictCV.STmodel}}.}
\source{
  Contains parametere estimates for the Spatio-Temporal
  model applied to monitoring data from the \strong{MESA
  Air} project, see Cohen et.al. (2009) and
  \code{\link{mesa.data.raw}} for details.
}
\description{
  Example of 10-fold cross-validated for the model in
  \code{\link{mesa.model}} using
  \code{\link{estimateCV.STmodel}} and
  \code{\link{predictCV.STmodel}}.
}
\examples{
##load data
data(mesa.model)

##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))

\dontrun{
  ##estimate different parameters for each CV-group
  est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)
  ##Do cross-validated predictions using the just estimated parameters
  pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa$par.cov, 
                            est.cv.mesa$Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)
##and results of estimation
data(est.mesa.model)

##examining the estimation results
print( est.cv.mesa )
names( est.cv.mesa )
 
##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(12,2.5,2,.5), las=2)
boxplot( est.cv.mesa, plot.type="cov")
points( coef(est.mesa.model,"cov")$par, pch=19, col=2)
 
##examining the prediction
print( pred.cv.mesa )
names( pred.cv.mesa )
 
##Plot observations with CV-predictions and prediction intervals
##TODO TODO TODO
#par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
#plotCV(mesa.data.res$pred.cv,  1, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 10, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 17, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 22, mesa.data.model)
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
  \code{\link{estimateCV.STmodel}} and
  \code{\link{predictCV.STmodel}} for cross-validation.
  \cr \code{\link{createSTmodel}} for creation of the
  originating \code{STmodel} object.

  Other example data: \code{\link{est.mesa.model}},
  \code{\link{MCMC.mesa.model}}, \code{\link{mesa.data}},
  \code{\link{mesa.data.raw}}, \code{\link{mesa.model}},
  \code{\link{pred.mesa.model}}
}
\keyword{datasets}

