\name{summary.STmodel}
\alias{summary.STmodel}
\title{Computes summary details for \code{STmodel} object}
\usage{
  \method{summary}{STmodel} (object,
    type = object$covars$type, ...)
}
\arguments{
  \item{object}{\code{STmodel} object to compute summary
  information for.}

  \item{type}{Factorial of \code{length(x$locations$ID)},
  if not \code{NULL} the output also presents summaries of
  number of sites and observations as well as time periods
  per type of site.}

  \item{...}{Ignored additional arguments.}
}
\value{
  A \code{summary.STmodel} object.
}
\description{
  \code{\link[base:summary]{summary}} method for class
  \code{STmodel}.
}
\examples{
##load some data
data(mesa.model)
##Summary of data fields.
summary(mesa.model)
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
  \code{\link{print.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}}
}

