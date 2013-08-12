\name{summary.mcmcSTmodel}
\alias{summary.mcmcSTmodel}
\title{Computes summary details for \code{mcmcSTmodel} object}
\usage{
  \method{summary}{mcmcSTmodel} (object, burnIn = 0, ...)
}
\arguments{
  \item{object}{\code{mcmcSTmodel} object to compute
  summary information for.}

  \item{burnIn}{Number of initial iterations to drop.}

  \item{...}{Ignored additional arguments.}
}
\value{
  A \code{summary.mcmcSTmodel} object.
}
\description{
  \code{\link[base:summary]{summary}} method for class
  \code{mcmcSTmodel}.
}
\examples{
##load data
  data(MCMC.mesa.model)
  summary(MCMC.mesa.model)
}
\author{
  Johan Lindström
}
\seealso{
  Other mcmcSTmodel methods:
  \code{\link{density.mcmcSTmodel}}, \code{\link{MCMC}},
  \code{\link{MCMC.STmodel}},
  \code{\link{plot.density.mcmcSTmodel}},
  \code{\link{plot.mcmcSTmodel}},
  \code{\link{print.mcmcSTmodel}},
  \code{\link{print.summary.mcmcSTmodel}}
}

