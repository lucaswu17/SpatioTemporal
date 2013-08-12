\name{plot.SVDcv}
\alias{boxplot.SVDcv}
\alias{plot.SVDcv}
\title{Plot and Boxplot cross-validation statistics for \code{SVDcv} object}
\usage{
  \method{plot}{SVDcv} (x,
    y = c("all", "MSE", "R2", "AIC", "BIC"), pairs = FALSE,
    sd = FALSE, ...)

  \method{boxplot}{SVDcv} (x,
    y = c("all", "MSE", "R2", "AIC", "BIC"), ...)
}
\arguments{
  \item{x}{\code{SVDcv} object to plot.}

  \item{y}{Which CV-statistic to plot. For pairs
  \code{"all"} implies \code{"BIC"}.}

  \item{pairs}{\code{TRUE}/\code{FALSE} plot
  cross-validation statistics, or scatter plot of
  individual BIC:s.}

  \item{sd}{\code{TRUE}/\code{FALSE} add uncertainty to
  each CV-statistic.}

  \item{...}{Additional parameters passed to
  \code{\link[graphics:plot]{plot}} or
  \code{\link[graphics:plot]{pairs}}.}
}
\value{
  Nothing
}
\description{
  \code{\link[graphics:plot]{plot}} and
  \code{\link[graphics:boxplot]{boxplot}} methods for class
  \code{SVDcv}. Plots summary statistics for the
  cross-validation. Plots include RMSE, R2, BIC, and
  scatter plots of BIC for each column.
}
\examples{
##See SVDsmooth example
}
\author{
  Johan Lindström
}
\seealso{
  Other SVD for missing data:
  \code{\link{calcSmoothTrends}},
  \code{\link{print.SVDcv}}, \code{\link{summary.SVDcv}},
  \code{\link{SVDmiss}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}},
  \code{\link{updateSTdataTrend}},
  \code{\link{updateTrend}},
  \code{\link{updateTrend.STdata}},
  \code{\link{updateTrend.STmodel}}

  Other SVDcv methods: \code{\link{print.SVDcv}},
  \code{\link{summary.SVDcv}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}}
}

