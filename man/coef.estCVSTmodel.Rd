\name{coef.estCVSTmodel}
\alias{coef.estCVSTmodel}
\title{Returns estimated parameters for each CV-group.}
\usage{
  \method{coef}{estCVSTmodel} (object,
    pars = c("all", "cov", "reg"), ...)
}
\arguments{
  \item{object}{\code{estCVSTmodel} object from which to
  extract estimated parameters.}

  \item{pars}{One of "cov", "reg", "all"; which parameters
  to extract.}

  \item{...}{Ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[stats:coef]{coef}} method for class
  \code{estCVSTmodel}.
}
\examples{
##load data
  data(est.cv.mesa)
  ##extract all parameters
  coef(est.cv.mesa)
  ##extract only covariance parameters
  coef(est.cv.mesa, pars="cov")
}
\author{
  Johan Lindström
}
\seealso{
  Other estCVSTmodel methods:
  \code{\link{boxplot.estCVSTmodel}},
  \code{\link{estimateCV}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{predictCV}}, \code{\link{predictCV.STmodel}},
  \code{\link{print.estCVSTmodel}},
  \code{\link{print.summary.estCVSTmodel}},
  \code{\link{summary.estCVSTmodel}}
}

