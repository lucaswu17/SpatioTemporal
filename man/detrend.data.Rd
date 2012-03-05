\name{detrend.data}
\encoding{latin1}
\Rdversion{1.1}
\alias{detrend.data}
\title{
  Removes the estimated time-trend from the observations in a \code{\link{mesa.data}} object.  }
\description{
 Removes the estimated time-trend from the observations in a \code{\link{mesa.data}} object, returning a similar object with a "vacant" trend, that can then be used to fit a simpler model.}

\usage{detrend.data(mesa.data,subregion=NA,method="rlm")}

\arguments{
  \item{mesa.data}{A \code{\link{mesa.data}} object.}
  \item{subregion}{\code{NA} (default), or a factor vector of the same
    length and order of the observations, with levels corresponding to
    sub-regions that should be treated differntly. }

  \item{method}{The method used for trend removal. At present only "rlm"
    is implemented.}
}

\details{ 
Sometimes there is no apparent spatial structure to the time-trend amplitude, or there is not enough identifiability in the data to properly model the structure. In that case, it is possible, at least as a sensitivity analysis, to de-trend the observations and run a model with a spatial field for the intercept only (apart from the spatio-temporal residual field).

When the returned object is used as input to \code{\link{create.data.model}}, the \code{LUR} list should contain only one vector of variable IDs (for the intercept). Similarly, vectors with initial parameter values for \code{\link{fit.mesa.model}} should be of length 5.

\code{detrend.data} will remove the trends from the observations, using the \code{trend} component of the \code{\link{mesa.data}} input object. The "rlm" method simply performs robust regression on all the observations as a bulk, with the time trend(s) as predictor(s). The "observations" in the new returned object are the residuals from this regression.

If a vector with factor levels is provided for \code{subregion}, then the model will be an interaction between the time trend(s) and region identifiers. }

\seealso{
  \code{\link{mesa.data}}, \code{\link{mesa.data.model}}, \code{\link{create.data.model}}.
}
 
\author{Assaf P. Oron}

\value{
  A \code{\link{mesa.data}} object, with the following differences:
  \item{obs}{This data frame has one additional column labeled
    "removedtrend", with the amount subtracted per observation. This
    enables back-calculation of the original values and of
    original-scale predictions.}
  \item{trend}{The \code{trend} data frame only has a "date" vector. It
    is still needed in this format, for compatibility as input to
    \code{\link{create.data.model}} and other functions.}
  \item{oldtrend}{An additional data frame with the old trend(s) used
    for the detrending.}
}
