\name{summary.STdata}
\alias{summary.STdata}
\title{Computes summary details for \code{STdata} object}
\usage{
  \method{summary}{STdata} (object,
    type = object$covars$type, ...)
}
\arguments{
  \item{object}{\code{STdata} object to compute summary
  information for.}

  \item{type}{Factorial of \code{length(x$covars$ID)}, if
  not \code{NULL} summaries for the observations are
  computed per type of site.}

  \item{...}{Ignored additional arguments.}
}
\value{
  A \code{summary.STdata} object.
}
\description{
  \code{\link[base:summary]{summary}} method for class
  \code{STdata}.
}
\author{
  Johan Lindström
}
\seealso{
  Other STdata methods: \code{\link{createSTdata}},
  \code{\link{plot.STdata}}, \code{\link{plot.STmodel}},
  \code{\link{print.STdata}},
  \code{\link{print.summary.STdata}},
  \code{\link{qqnorm.predCVSTmodel}},
  \code{\link{qqnorm.STdata}},
  \code{\link{qqnorm.STmodel}},
  \code{\link{scatterPlot.predCVSTmodel}},
  \code{\link{scatterPlot.STdata}},
  \code{\link{scatterPlot.STmodel}}
}

