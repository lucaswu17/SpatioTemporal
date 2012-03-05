\name{compute.ltaCV}
\encoding{latin1}
\Rdversion{1.1}
\alias{compute.ltaCV}

\title{
  Computes the long term average for each sites.
}
\description{
  Computes the long term average of observations and cross-validated
  predictions for each of the sites in \code{predCV}. The long term
  averages are computed using \emph{only} timepoints that have
  observations, this applies to both the observed and predicted. Also
  the function allows for a transformation: if requested the
  transformation is applied \emph{before} the averaging.
}
\usage{
compute.ltaCV(predCV, trans = NULL)
}
\arguments{
  \item{predCV}{
    Result of a cross-validation. Should be the output from
    \code{\link{predictCV}}.
  }
  \item{trans}{
    Transform observations and predictions \emph{before} computing
    averages. \code{NULL} gives no transformation, \code{trans=0} takes
    the exponent of data (assuming that data was originally
    log-transformed) a non-zero value will raise the data to that power
    (assuming an original root-transform), e.g. \code{trans=2} gives
    the long term averages as \code{mean(obs^2)} and \code{mean(pred^2)}.
  }
}
\value{
  Returns a (number of locations) - by - 2 matrix with the observed and
  predicted value for each location.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
See \code{\link{createCV}} and \code{\link{estimateCV}} for cross-validation
set-up and estimation.

For computing CV statistics, see also \code{\link{predictNaive}}, 
and for further illustration see \code{\link{plotCV}}, 
\code{\link{CVresiduals.qqnorm}}, and \code{\link{summaryStatsCV}}.
}
\examples{
##load data
data(mesa.data.res)

##Extract pre-computed cross-validated predictions
pred.cv <- mesa.data.res$pred.cv

##compute long term averages of predictions and observations
pred.lta <- compute.ltaCV(pred.cv)

##we can now compare observed and predicted averages at each site
plot(pred.lta[,"obs"],pred.lta[,"pred"])
abline(0,1)
}
