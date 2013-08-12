\name{print.estimateSTmodel}
\alias{print.estimateSTmodel}
\title{Print details for \code{estimateSTmodel} object}
\usage{
  \method{print}{estimateSTmodel} (x, ...)
}
\arguments{
  \item{x}{\code{estimateSTmodel} object to print
  information for.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{estimateSTmodel}.
}
\examples{
##load data
  data(est.mesa.model)
  print(est.mesa.model)
}
\author{
  Johan Lindström
}
\seealso{
  Other estimateSTmodel methods:
  \code{\link{coef.estimateSTmodel}},
  \code{\link{estimate}}, \code{\link{estimate.STmodel}}
}

