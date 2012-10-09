\name{print.STmodel}
\alias{print.STmodel}
\title{Print details for \code{STmodel} object}
\usage{
  \method{print}{STmodel} (x, type = x$locations$type, ...)
}
\arguments{
  \item{x}{\code{STmodel} object to print information for.}

  \item{type}{Factorial of \code{length(x$locations$ID)},
  if not \code{NULL} the output also presents summaries of
  number of sites and observations as well as time periods
  per type of site.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{STmodel}.
}
\examples{
##load some data
data(mesa.model)
##print basic information regarding obs, locations, dates, etc
print(mesa.model)
}
\author{
  Johan Lindström
}
\seealso{
  Other STmodel methods: \code{\link{createSTmodel}},
  \code{\link{c.STmodel}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{estimate.STmodel}},
  \code{\link{MCMC.STmodel}}, \code{\link{plot.STdata}},
  \code{\link{plot.STmodel}},
  \code{\link{predictCV.STmodel}},
  \code{\link{predict.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}},
  \code{\link{summary.STmodel}}
}

