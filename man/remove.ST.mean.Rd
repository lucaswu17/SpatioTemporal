\name{remove.ST.mean}
\encoding{latin1}
\Rdversion{1.1}
\alias{remove.ST.mean}

\title{
  Mean-Centre the Spatio-Temporal Covariate
}
\description{
  Removes the temporal mean at each location for the spatio-temporal
  covariares. The means are added to the returned data structure and
  can be used as geographic covariates.
}
\usage{
remove.ST.mean(mesa.data)
}
\arguments{
  \item{mesa.data}{
    Data structure holding observations, observation locations, smooth
    temporal trends, geographic covariates and spatio-temporal
    covariates. See \code{\link{mesa.data}}.
  }
}
\value{
  Returns a modfied version of the input \code{\link{mesa.data}}, where
  the site by site temporal average has been removed from each of the
  spatio-temporal covariates. The averages are seen as geographic
  covariates and added to \code{mesa.data$LUR}.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\examples{
##load data
data(mesa.data)

mesa.data.mean0 <- remove.ST.mean(mesa.data)

##compare the data structures
##geographic covariates
summary(mesa.data$LUR)
summary(mesa.data.mean0$LUR)

##mean of the spatio-temporal covariate
cbind(colMeans(mesa.data$SpatioTemp),
      colMeans(mesa.data.mean0$SpatioTemp))

##mean of the spatio-temporal covariate
##compared to the added mean covariate
plot(mesa.data.mean0$LUR$mean.lax.conc.1500,
     colMeans(mesa.data$SpatioTemp))
\dontshow{
if( max(abs(colMeans(mesa.data.mean0$SpatioTemp))) > 1e-10 ){
  stop("remove.ST.mean 1: mean(centred ST) != 0")
}
if( max(abs(mesa.data.mean0$LUR$mean.lax.conc.1500 - 
            colMeans(mesa.data$SpatioTemp))) > 1e-10 ){
  stop("remove.ST.mean 2: extracted mean != mean(ST)")
}
}%\dontshow
}%\examples
