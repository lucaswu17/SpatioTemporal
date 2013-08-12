\name{print.estCVSTmodel}
\alias{print.estCVSTmodel}
\title{Print details for \code{estCVSTmodel} object}
\usage{
  \method{print}{estCVSTmodel} (x, ...)
}
\arguments{
  \item{x}{\code{estCVSTmodel} object to print information
  for.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{estCVSTmodel}.
}
\examples{
##load some data
data(est.cv.mesa)
##print basic information for the CV-predictions
print(est.cv.mesa)
}
\author{
  Johan Lindström
}
\seealso{
  Other estCVSTmodel methods:
  \code{\link{boxplot.estCVSTmodel}},
  \code{\link{coef.estCVSTmodel}},
  \code{\link{estimateCV}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{predictCV}}, \code{\link{predictCV.STmodel}},
  \code{\link{print.summary.estCVSTmodel}},
  \code{\link{summary.estCVSTmodel}}
}

