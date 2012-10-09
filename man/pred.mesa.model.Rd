\docType{data}
\name{pred.mesa.model}
\alias{pred.mesa.model}
\alias{pred.mesa.model.obs}
\title{Example of a \code{predictSTmodel} Structure}
\format{A list with elements, see the return description in
  \code{\link{predict.STmodel}}.}
\source{
  Contains parametere estimates for the Spatio-Temporal
  model applied to monitoring data from the \strong{MESA
  Air} project, see Cohen et.al. (2009) and
  \code{\link{mesa.data.raw}} for details.
}
\description{
  Example of a predictions for the model in
  \code{\link{mesa.model}} using
  \code{\link{predict.STmodel}}. Two sets of predictions
  are presented, \code{pred.mesa.model} and
  \code{pred.mesa.model.obs}.
}
\examples{
##load data
data(mesa.model)
data(est.mesa.model)

\dontrun{
  ##compute predictions at all locations, including variances
  pred.mesa.model <- predict(mesa.model, est.mesa.model, pred.var=TRUE)
}

##separate data into AQS and FIXED sites
##first find indicators
I.aqs <- mesa.model$locations$ID[mesa.model$locations$type=="AQS"]
I.aqs <- mesa.model$obs$ID \%in\% I.aqs
##then drop relevant parts of STmodel
mesa.model.aqs <- dropObservations(mesa.model, !I.aqs)
mesa.model.fixed <- dropObservations(mesa.model, I.aqs)
                                     
\dontrun{
  ##compute predictions at only observed FIXED locations, based on observations
  ##at AQS locations.
  pred.mesa.model.obs <- predict(mesa.model.aqs, est.mesa.model,
                                 STdata=mesa.model.fixed, only.obs=TRUE)
}

##lets load precomputed results instead
data(pred.mesa.model)

##study results
print(pred.mesa.model)
print(pred.mesa.model.obs)

##simple plots
par(mfrow=c(2,1))
plot(pred.mesa.model)
plot(pred.mesa.model.obs, STmodel=mesa.model.fixed)

##scatter plot predicitons as function of observed
par(mfrow=c(1,1))
plot(pred.mesa.model.obs, y="obs", STmodel=mesa.model.fixed, ID="all",
     pred.var=TRUE, lty=c(NA,1),pch=c(19,NA))
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
  \code{\link{predict.STmodel}} for prediction. \cr
  \code{\link{createSTmodel}} for creation of the
  originating \code{STmodel} object.

  Other example data: \code{\link{est.cv.mesa}},
  \code{\link{est.mesa.model}},
  \code{\link{MCMC.mesa.model}}, \code{\link{mesa.data}},
  \code{\link{mesa.data.raw}}, \code{\link{mesa.model}}
}
\keyword{datasets}

