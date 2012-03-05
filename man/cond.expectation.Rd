\name{cond.expectation}
\encoding{latin1}
\Rdversion{1.1}
\alias{cond.expectation}

\title{
  Computes Conditional Expectation for Unobserved Locations
}
\description{
  Compute the conditional expectations (i.e. predictions) at the
  unobserved space-time locations.
}
\usage{
cond.expectation(x, mesa.data.model, mesa.data = NA, Nmax = 1000,
    only.pars = FALSE, compute.beta = FALSE, no.nugget = FALSE,
    only.obs = FALSE, pred.var = TRUE, pred.covar = FALSE,
    combine.data = FALSE, type = "p") 
}
\arguments{
  \item{x}{
    Model parameters for which to compute the conditional
    expectation.
  }
  \item{mesa.data.model}{
    Data structure holding observations and information regarding where
    to compute predictions. Predictions are computed \emph{conditional}
    on the observations in \code{\link{mesa.data.model}}. See further
    \code{\link{mesa.data.model}}.
  }
  \item{mesa.data}{
    Data structure holding additional locations, geographic covariates,
    and spatio-temporal covariates, smooth temporal trends, and
    \emph{optionally} additional observations. See further
    \code{\link{mesa.data}}.
  }
  \item{Nmax}{
    Limits the size of matrices constructed when computing
    expectations. Use a smaller value if memory becomes a problem.
  }
  \item{only.pars}{
    Compute only the regression parameters (using GLS) along with the
    related variance.

    \code{only.pars=TRUE} is \emph{incompatible} with \code{type="f"}
  }
  \item{compute.beta}{
    Compute the latent beta-fields at all locations in
    \code{mesa.data.model} and/or \code{mesa.data}.
  }
  \item{no.nugget}{
    Set the nugget for the residual nu field to zero when computing the
    conditional expectations. Makes no difference at unobserved
    locations but implies smoothing at locations with observations.
  }
  \item{only.obs}{
    Compute predictions at only the observations specified in
    \code{\link{mesa.data}}. Used to limit computations when doing 
    cross-validation.

    This \emph{requires} that \code{mesa.data!=NA} and that
    \code{mesa.data$obs$obs} exist.

    \code{only.obs = TRUE} is \emph{incompatible} with \code{combine.data =
      TRUE} and \cr \code{pred.covar = TRUE}.
  }
  \item{pred.var}{
    If \code{TRUE}, compute point-wise variances for the predictions.
  }
  \item{pred.covar}{
    If \code{TRUE}, compute covariance matrices for the predicted time
    series at each location.
    
    \code{pred.covar = TRUE} \emph{implies} \code{pred.var = TRUE} and
    sets \code{Nmax} equal to the number of timepoints.

    \code{pred.covar = TRUE} is \emph{incompatible} with \code{only.obs
      = TRUE}.
  }
  \item{combine.data}{
    If \code{TRUE} will combine locations and covariates from
    \code{\link{mesa.data.model}} and \code{\link{mesa.data}}, but will
    use temporal trends from \code{\link{mesa.data.model}}.
    Predictions are still conditional on \emph{only} observations in
    \code{\link{mesa.data.model}}. See further
    \code{\link{combineMesaData}}.

    \code{combine.data = TRUE} is \emph{incompatible} with
    \code{only.obs = TRUE}
  }
  \item{type}{
    A single character indicating the type of log-likelihood to
    use. Valid options are "f", "p", and "r", for \emph{full},
    \emph{profile} or \emph{restricted maximum likelihood}
    (REML).
  }
}
\details{
  Predictions are computed for the space-time locations in
  \code{mesa.data.model} and/or \cr \code{mesa.data}, conditional on the
  observations in \code{mesa.data.model} and parameters given in
  \code{x}.

  Depending on the the model requested in \code{type} the parameters
  specified in \code{x} can be either the \emph{log}-covariances or
  regression parameters and \emph{log}-covariances parameters.

  If \code{type="f"} the user has the freedom to specify
  \emph{all} model parameters -- including the regression coefficients
  -- rather than just the \emph{log}-covariance parameters.
  However, if \code{type="f"} and only \emph{log}-covariance
  parameters are given then the regression parameters are inferred
  using generalised least squares (GLS). If \code{type!="f"} the
  regression parameters are \emph{always} computed using GLS (any input
  is \emph{ignored}). The GLS can be obtained through
  \code{cond.expectation(x, mesa.data.model, only.pars = TRUE,
    type="p")}.

  The space-time locations at which to compute predictions are specified
  through \cr \code{mesa.data.model} and \code{mesa.data}. If
  \code{mesa.data=NA} predictions will be computed at all space-time
  locations given in \code{mesa.data.model}. If \code{mesa.data!=NA}
  three options exist:
  \enumerate{
    \item If \code{combine.data=FALSE} predictions will be given at
    \emph{only} the space-time locations in \code{mesa.data}. If
    \code{mesa.data} does not have a smooth temporal trend, the trend in
    \cr \code{mesa.data.model} will be used (with a
    warning). Predictions are still computed conditional on the
    observations in \code{mesa.data.model} and any observations in
    \code{mesa.data} are \emph{ignored}.

    \item If \code{combine.data=FALSE} and \code{only.obs=TRUE}
    predictions are computed as in (1) \cr above, but only at space-time
    locations given in \code{mesa.data$obs}. This is primarily used to
    limit the number of computations needed when doing
    cross-validation.
  
    \item If \code{combine.data=TRUE} locations and covariates from both
    \code{mesa.data.model} and \code{mesa.data} are combined using
    \code{\link{combineMesaData}}. The smooth temporal trends are taken
    from \code{mesa.data.model} with trends in \code{mesa.data} being
    \emph{ignored}. I.e. predictions are computed conditional on
    \emph{only} the observations in \code{mesa.data.model}, at locations
    given in both \code{mesa.data.model} and \code{mesa.data} and at
    \emph{only} the times specified in \code{mesa.data.model$trend}.
  }

  In addition to computing the conditional expectation at a number of
  space-time locations the function also computes predictions based on
  only the regression part of the model as well as the latent
  beta-fields (the latter only if \code{compute.beta=TRUE}).

  \strong{NOTE:} When computing the latent beta-fields, some matrix
  multiplications occur in a different order. this can lead to
  differences on the order of \code{1e-13} in some conditional
  expectations and variances.

  Predictions variances can also be computed. If \code{pred.var=TRUE}
  point-wise variances for the predictions (and the latent beta-fields)
  are computed. If instead \code{pred.covar=TRUE} the full covariance
  matrices for each predicted time series is computed (as well as full
  covariance matrices for the latent beta-fields); this implies that the
  covariances between temporal predictions at the same location are
  calculated but \emph{not}, due to memory restrictions, any covariances
  between locations.
}
\value{
  The function returns a list containing:
  \item{pars}{A list with regression parameters and related variances.
    The variances are set to zero if \code{type="f"} and the input
    \code{x} contains the regression parameters.

    \code{pars} contain \code{gamma.E} and \code{alpha.E} which are
    regression coefficients for the spatio-temporal model and land-use
    covaraiates respectively. The associated variances are found in
    \code{gamma.V} and \code{alpha.V}, with the cross-covariance between
    gamma and alpha in \code{gamma.alpha.C}
  }
  \item{EX.beta.mu}{A (number of locations) - by - (number of smooth
    trends) matrix with mean value component of the latent beta-fields
    (i.e. covariates * alpha).
    
    \code{NULL} if \code{compute.beta=FALSE}.
  }
  \item{EX.beta}{A (number of locations) - by - (number of smooth
    trends) matrix with predictions of the latent beta-fields. 

    \code{NULL} if \code{compute.beta=FALSE}.
  }
  \item{VX.beta}{A (number of locations) - by - (number of smooth trends)
    matrix or a (number of locations) - by - (number of locations) - by
    - (number of smooth trends) 3D-array containing the prediction
    (co)variances for the beta-fields.

    \code{NULL} if \code{compute.beta=FALSE} or \code{pred.var=FALSE}, a
    matrix if \cr \code{pred.covar=FALSE}, and a 3D-array if
    \code{pred.covar=TRUE}.
  }
  \item{EX.mu}{A (number of timepoints) - by - (number of locations)
    matrix with predictions based on only the regression parameters,
    geographic covariates, and temporal trends. I.e. only the
    deterministic part of the spatio-temporal model. Predictions are
    computed at the space-time locations defined by
    \code{mesa.data.model} and/or \code{mesa.data}.
    
    If \code{only.obs=TRUE} this is a vector with predictions at the
    space-time locations of the observations in \code{mesa.data$obs}.
  }
  \item{EX.mu.beta}{Predictions based on the latent-beta fields, but not
    on the residual nu field; i.e. the predictions in \code{EX.mu} with
    kriging for the geographic covariates.
    
    Same size as \code{EX.mu} above; \code{NULL} if
    \code{compute.beta=FALSE}.
  }
  \item{EX}{Full predictions at the space-time locations defined by
    \code{mesa.data.model} and/or \code{mesa.data}.

    Same size as \code{EX.mu} above.
  }
  \item{VX}{Pointwise prediction variances for all locations in
    \code{EX} above. \cr \code{NULL} if \code{pred.var=FALSE}.
  }
  \item{VX.full}{A list with (number of locations) elements, each
    element is a (number of timepoints) - by - (number of timepoints)
    temporal covariance matrix for the predicted timeseries at that
    location. \code{NULL} if \code{pred.covar=FALSE}.
  }
  \item{I}{A vector with the locations of the observations in
    \code{mesa.data.model} or \code{mesa.data}. To extract predictions
    at these observations locations use \code{EX[I]}. \code{NULL} if
    \code{only.obs=TRUE}.
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  For optimization functions see \code{\link{loglike}}, 
  \code{\link{loglike.var.names}}, \code{\link{create.data.model}}, \cr
  \code{\link{fit.mesa.model}}, and \code{\link{run.MCMC}}.
  
  For prediction see also \code{\link{fit.mesa.model}}, and 
  \code{\link{plotCV}} for plotting prediction results.

  \code{\link{combineMesaData}} is used to combine
  \code{\link{mesa.data}} and \code{\link{mesa.data.model}} structures.
}
\examples{
##load data
data(mesa.data)
data(mesa.data.model)
data(mesa.data.res)

##extract estimated parameters
x <- mesa.data.res$par.est$res.best$par

##find regression parameters using GLS
x.reg <- cond.expectation(x, mesa.data.model, only.pars = TRUE)
str(x.reg)

\dontrun{
##compute predictions at all locations, including beta-fields
EX <- cond.expectation(x, mesa.data.model, compute.beta=TRUE)
##compute predictions at only observations locations
EX.obs <- cond.expectation(x, mesa.data.model, mesa.data=mesa.data,
                           only.obs=TRUE, compute.beta=FALSE)
}
##Let's load precomputed results instead.
EX <- mesa.data.res$EX
EX.obs <- mesa.data.res$EX.obs

##Compare the predictions at all locations and only obs
dim(EX$EX)
dim(EX.obs$EX)

##estimate beta from the observations for reference
##create data matrix
D <- create.data.matrix(mesa.data)
beta <- matrix(NA,dim(D)[2], dim(mesa.data$trend)[2])
##extact the temporal trends
F <- mesa.data$trend
##drop the date column
F$date <- NULL
##estimate the beta-coeficients at each location
for(i in 1:dim(D)[2])
  beta[i,] <- summary(lm(D[,i] ~ as.matrix(F)))$coefficients[,1]

##Study the results
##Start by comparing beta fields
par(mfcol=c(1,1), mar=c(4.5,4.5,2,.5), pty="s")
plot(x=beta[,1], y=EX$EX.beta[,1], main="Temporal Intercept",
     xlab="Empirical estimate", ylab="Spatio-Temporal Model")
abline(0,1,col="grey")

##plot predictions and observations for 4 locations
par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
plotPrediction(EX,  1, mesa.data.model)
plotPrediction(EX, 10, mesa.data.model)
plotPrediction(EX, 17, mesa.data.model)
plotPrediction(EX, 22, mesa.data.model)

##compare the only obs predictions with what we can extract from EX
par(mfcol=c(1,1), mar=c(4.5,4.5,2,.5), pty="s")
plot(EX$EX[EX$I], EX.obs$EX)

##minor (1e-14) numerical differences in the results
print(range(EX$EX[EX$I]-EX.obs$EX))
}%\examples
