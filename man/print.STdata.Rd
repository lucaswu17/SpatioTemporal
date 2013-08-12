\name{print.STdata}
\alias{print.STdata}
\title{Print details for \code{STdata} object}
\usage{
  \method{print}{STdata} (x, type = x$covars$type, ...)
}
\arguments{
  \item{x}{\code{STdata} object to print information for.}

  \item{type}{Factorial of \code{length(x$covars$ID)}, if
  not \code{NULL} the output also presents summaries of
  number of sites and observations as well as time periods
  per type of site.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{STdata}.
}
\author{
  Johan Lindström
}
\seealso{
  Other STdata methods: \code{\link{createSTdata}},
  \code{\link{plot.STdata}}, \code{\link{plot.STmodel}},
  \code{\link{print.summary.STdata}},
  \code{\link{qqnorm.predCVSTmodel}},
  \code{\link{qqnorm.STdata}},
  \code{\link{qqnorm.STmodel}},
  \code{\link{scatterPlot.predCVSTmodel}},
  \code{\link{scatterPlot.STdata}},
  \code{\link{scatterPlot.STmodel}},
  \code{\link{summary.STdata}}
}

