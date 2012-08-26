\name{predict.STmodel}
\alias{predict.STmodel}
\title{Computes Conditional Expectation at Unobserved Locations}
\usage{
  \method{predict}{STmodel} (object, x, STdata = NULL,
    Nmax = 1000, only.pars = FALSE, nugget.unobs = 0,
    only.obs = FALSE, pred.var = TRUE, pred.covar = FALSE,
    beta.covar = FALSE, combine.data = FALSE, type = "p",
    ...)
}
\arguments{
  \item{object}{\code{STmodel} object for which to compute
  predictions.}

  \item{x}{Model parameters for which to compute the
  conditional expectation. Either as a vector/matrix or an
  \code{estimateSTmodel} from
  \code{\link{estimate.STmodel}}.}

  \item{STdata}{\code{STdata}/\code{STmodel} object with
  locations/times for which to predict. If not given
  predictions are computed for locations/times in
  \code{object}}

  \item{Nmax}{Limits the size of matrices constructed when
  computing expectations. Use a smaller value if memory
  becomes a problem.}

  \item{only.pars}{Compute only the regression parameters
  (using GLS) along with the related variance.}

  \item{nugget.unobs}{Value of nugget at unonserved
  locations, either a scalar or a vector with one element
  per unobserved site. \strong{NOTE:} All sites in
  \code{STdata} are considered unobserved!}

  \item{only.obs}{Compute predictions at only locations
  specified by observations in \code{STdata}. Used to limit
  computations when doing cross-validation.
  \code{only.obs=TRUE} \emph{implies}
  \code{pred.covar=FALSE} and \code{combine.data=FALSE}.
  Further \code{\link{createSTmodel}} will be called on any
  \code{STdata} input, possibly \emph{reordering the
  observations.}}

  \item{pred.var,pred.covar}{Compute point-wise prediction
  variances; or compute covariance matrices for the
  predicted time series at each location.
  \code{pred.covar=TRUE} \emph{implies}
  \code{pred.var=TRUE} and sets \code{Nmax} equal to the
  number of timepoints.}

  \item{beta.covar}{Compute the full covariance matrix for
  the latent beta-fields, otherwise only the diagonal
  elements of V(beta|obs) are computed.}

  \item{combine.data}{Combine \code{object} and
  \code{STdata} and predict for the joint set of points,
  see \code{\link{c.STmodel}}.}

  \item{type}{A single character indicating the type of
  log-likelihood to compute. Valid options are "f", "p",
  and "r", for \emph{full}, \emph{profile} or
  \emph{restricted maximum likelihood} (REML).}

  \item{...}{Ignored additional arguments.}
}
\value{
  The function returns a list containing (objects not
  computed will be missing): \item{opts}{Copy of options
  used in the function call.} \item{pars}{A list with
  regression parameters and related variances.  \code{pars}
  contain \code{gamma.E} and \code{alpha.E} with regression
  coefficients for the spatio-temporal model and land-use
  covaraiates; variances are found in \code{gamma.V} and
  \code{alpha.V}; cross-covariance between gamma and alpha
  in \code{gamma.alpha.C}.} \item{beta}{A list with
  estimates of the beta-fields, including the regression
  mean \code{mu}, conditional expectations \code{EX},
  possibly variances \code{VX}, and the full covariance
  matrix \code{VX.full}.} \item{EX.mu}{predictions based on
  the regression parameters, geographic covariates, and
  temporal trends. I.e. only the deterministic part of the
  spatio-temporal model.} \item{EX.mu.beta}{Predictions
  based on the latent-beta fields, but excluding the
  residual nu field.} \item{EX}{Full predictions at the
  space-time locations in \code{object} and/or
  \code{STdata}.} \item{VX}{Pointwise variances for all
  locations in \code{EX}.} \item{VX.pred}{Pointwise
  prediction variances for all locations in \code{EX}, i.e.
  including contribution from \code{nugget.unobs}.}
  \item{VX.full}{A list with (number of locations)
  elements, each element is a (number of timepoints) - by -
  (number of timepoints) temporal covariance matrix for the
  timeseries at each location.} \item{I}{A vector with the
  locations of the observations in \code{object} or
  \code{STdata}. To extract predictions at the observations
  locations use \code{EX[I]}.}
}
\description{
  Compute the conditional expectations (i.e. predictions)
  at the unobserved space-time locations. Predictions are
  computed for the space-time locations in \code{object}
  and/or \code{STdata}, conditional on the observations in
  \code{object} and parameters given in \code{x}.
}
\details{
  In addition to computing the conditional expectation at a
  number of space-time locations the function also computes
  predictions based on only the regression part of the
  model as well as the latent beta-fields.

  Prediction are computed as the conditional expectation of
  a latent field given observations. This implies that
  \code{E(X_i| Y_i) != Y_i}, with the difference being due
  to smoothing over the nugget. Further two possible
  variance are given, \code{V(X_i|Y_i)} and
  \code{V(X_i|Y_i)+nugget_i}. Here the nugget for
  unobserved locations needs to be specified as an
  additional argument \code{nugget.nobs}. The two variances
  correspond, losely, to confidence and prediction
  intervals.

  Predictions variances can also be computed. If
  \code{pred.var=TRUE} point-wise variances for the
  predictions (and the latent beta-fields) are computed. If
  instead \code{pred.covar=TRUE} the full covariance
  matrices for each predicted time series is computed; this
  implies that the covariances between temporal predictions
  at the same location are calculated but \emph{not}, due
  to memory restrictions, any covariances between
  locations. \code{beta.covar=TRUE} gives the full
  covariance matrices for the latent beta-fields.
}
\examples{
##load data
data(mesa.model)
data(est.mesa.model)

##find regression parameters using GLS
x.reg <- predict(mesa.model, est.mesa.model, only.pars = TRUE)
str(x.reg$pars)

##Split data into FIXED and AQS
I.aqs <- mesa.model$locations$ID[mesa.model$locations$type=="AQS"]
I.aqs <- mesa.model$obs$ID \%in\% I.aqs
mesa.model.aqs <- dropObservations(mesa.model, !I.aqs)
mesa.model.fixed <- dropObservations(mesa.model, I.aqs)

\dontrun{
  ##compute predictions at all locations, including beta-fields
  pred.mesa.model <- predict(mesa.model, est.mesa.model,
                             pred.var=TRUE)
  ##predict at FIXED using AQS-sites
  pred.mesa.model.obs <- predict(mesa.model.aqs, est.mesa.model,
                                 STdata=mesa.model.fixed, only.obs=TRUE)
}
##Let's load precomputed results instead.
data(pred.mesa.model)

##Compare the predictions at all locations and only obs
print(pred.mesa.model)
print(pred.mesa.model.obs)

##estimate beta from the observations for reference
##create data matrix
D <- createDataMatrix(mesa.model)
beta <- matrix(NA,dim(D)[2], dim(mesa.model$trend)[2])
##extact the temporal trends
F <- mesa.model$trend
##drop the date column
F$date <- NULL
##estimate the beta-coeficients at each location
for(i in 1:dim(D)[2]){
  beta[i,] <- coefficients( lm(D[,i] ~ as.matrix(F)) )
}

##Study the results
##Start by comparing beta fields
par(mfcol=c(1,1), mar=c(4.5,4.5,2,.5), pty="s")
plot(x=beta[,1], y=pred.mesa.model$beta$EX[,1],
     main="Temporal Intercept",
     xlab="Empirical estimate", ylab="Spatio-Temporal Model")
abline(0,1,col="grey")
##or just the regression part of the beta fields
points(x=beta[,1], y=pred.mesa.model$beta$mu[,1], col=2)

##plot predictions and observations for 4 locations
par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
plot(pred.mesa.model, ID=1, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=10, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=17, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=22, STmodel=mesa.model, pred.var=TRUE)
}
\author{
  Johan Lindström
}
\seealso{
  Other predictSTmodel methods:
  \code{\link{plot.predCVSTmodel}},
  \code{\link{plot.predictSTmodel}},
  \code{\link{print.predictSTmodel}}

  Other STmodel methods: \code{\link{c.STmodel}},
  \code{\link{createSTmodel}},
  \code{\link{estimate.STmodel}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{MCMC.STmodel}}, \code{\link{plot.STdata}},
  \code{\link{plot.STmodel}},
  \code{\link{predictCV.STmodel}},
  \code{\link{print.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}},
  \code{\link{summary.STmodel}}
}

