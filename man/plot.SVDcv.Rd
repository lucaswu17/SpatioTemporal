\name{plot.SVDcv}
\alias{plot.SVDcv}
\title{Plot cross-validation statistics for \code{SVDcv} object}
\usage{
  \method{plot}{SVDcv} (x, y = NULL, pairs = FALSE, ...)
}
\arguments{
  \item{x}{\code{SVDcv} object to plot.}

  \item{y}{Not used}

  \item{pairs}{\code{TRUE}/\code{FALSE} plot
  cross-validation statistics, or scatter plot of
  individual BIC:s.}

  \item{...}{Additional parameters passed to
  \code{\link[graphics:plot]{plot}} or
  \code{\link[graphics:plot]{pairs}}.}
}
\value{
  Nothing
}
\description{
  \code{\link[graphics:plot]{plot}} method for class
  \code{SVDcv}. Plots summary statistics for the
  cross-validation. Plots include RMSE, R2, BIC, and
  scatter plots of BIC for each column.
}
\examples{
## end of SVDsmooth example
}
\author{
  Johan Lindström
}
\seealso{
  Other SVDcv methods: \code{\link{print.SVDcv}},
  \code{\link{SVDsmooth}}, \code{\link{SVDsmoothCV}}

  Other SVD for missing data:
  \code{\link{calcSmoothTrends}},
  \code{\link{print.SVDcv}}, \code{\link{SVDmiss}},
  \code{\link{SVDsmooth}}, \code{\link{SVDsmoothCV}},
  \code{\link{updateSTdataTrend}}
}

