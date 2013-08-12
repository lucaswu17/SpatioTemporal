\name{print.predCVSTmodel}
\alias{print.predCVSTmodel}
\title{Print details for \code{predCVSTmodel} object}
\usage{
  \method{print}{predCVSTmodel} (x, ...)
}
\arguments{
  \item{x}{\code{predCVSTmodel} object to print information
  for.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{predCVSTmodel}.
}
\examples{
##load some data
data(pred.cv.mesa)
##print basic information for the CV-predictions
print(pred.cv.mesa)
}
\author{
  Johan Lindström
}
\seealso{
  Other predCVSTmodel methods: \code{\link{estimateCV}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{plot.predCVSTmodel}},
  \code{\link{plot.predictSTmodel}},
  \code{\link{predictCV}}, \code{\link{predictCV.STmodel}},
  \code{\link{print.summary.predCVSTmodel}},
  \code{\link{qqnorm.predCVSTmodel}},
  \code{\link{qqnorm.STdata}},
  \code{\link{qqnorm.STmodel}},
  \code{\link{scatterPlot.predCVSTmodel}},
  \code{\link{scatterPlot.STdata}},
  \code{\link{scatterPlot.STmodel}},
  \code{\link{summary.predCVSTmodel}}
}

