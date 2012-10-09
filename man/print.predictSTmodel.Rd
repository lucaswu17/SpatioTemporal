\name{print.predictSTmodel}
\alias{print.predictSTmodel}
\title{Print details for \code{predictSTmodel} object}
\usage{
  \method{print}{predictSTmodel} (x, ...)
}
\arguments{
  \item{x}{\code{predictSTmodel} object to print
  information for.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{predictSTmodel}.
}
\examples{
##load data
  data(pred.mesa.model)
  print(pred.mesa.model)
}
\author{
  Johan Lindström
}
\seealso{
  Other predictSTmodel methods:
  \code{\link{plot.predCVSTmodel}},
  \code{\link{plot.predictSTmodel}},
  \code{\link{predict.STmodel}}
}

