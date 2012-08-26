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
data(CV.mesa.model)
##print basic information for the CV-predictions
print(pred.cv.mesa)
}
\author{
  Johan Lindström
}
\seealso{
  Other predCVSTmodel methods:
  \code{\link{estimateCV.STmodel}},
  \code{\link{plot.predCVSTmodel}},
  \code{\link{plot.predictSTmodel}},
  \code{\link{predictCV.STmodel}},
  \code{\link{print.summary.predCVSTmodel}},
  \code{\link{summary.predCVSTmodel}}
}

