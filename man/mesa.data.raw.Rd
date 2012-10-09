\docType{data}
\name{mesa.data.raw}
\alias{mesa.data.raw}
\title{Example of raw data}
\format{The structure contains observations, temporal trends, locations,
  geographic covariates, and spatio-temporal covariates. The data is
  stored as a list with elements:
  \describe{
    \item{X}{A data.frame containing names, locations, and (geographic)
             covariates for all the (observation) locations.}
    \item{obs}{A time-by-location matrix for the observed data, missing data
               marked as \code{NA}}
    \item{lax.conc.1500}{A time-by-location matrix of a spatio-temporal
                         covariate based on output from Caline3QHC.}
  }}
\source{
  Contains monitoring data from the \strong{MESA Air}
  project, see Cohen et.al. (2009) for details.
}
\description{
  The raw data that was used to create the
  \code{\link{mesa.data}} and \code{\link{mesa.model}}
  structures. \cr The data structure contains raw data from
  the \strong{MESA Air} project. The example below
  describes how to create the \code{\link{mesa.data}}
  structure from raw data.
}
\examples{
##load the data
data(mesa.data.raw)

##create a matrix of time-points
T <- t(matrix(rownames(mesa.data.raw$obs),
              nrow=dim(mesa.data.raw$obs)[2],
              ncol=dim(mesa.data.raw$obs)[1], byrow=TRUE))
##...and locations
ID <- matrix(colnames(mesa.data.raw$obs), 
             nrow=dim(mesa.data.raw$obs)[1],
             ncol=dim(mesa.data.raw$obs)[2], byrow=TRUE)
##let's study these matrices
print(T[1:5,1:5])
print(ID[1:5,1:5])

##combine with the observations
obs <- data.frame(obs=c(mesa.data.raw$obs), date=as.Date(T),
                  ID=c(ID))
##drop unmonitored locations
obs <- obs[!is.na(obs$obs),,drop=FALSE]
##sort the locations (strictly not needed)
obs <- obs[order(obs$date,obs$ID),,drop=FALSE]
##drop names
rownames(obs) <- NULL

##create a 3D-array for the spatio-temporal covariate
ST <- array(mesa.data.raw$lax.conc.1500, dim =
            c(dim(mesa.data.raw$lax.conc.1500),1))
dimnames(ST) <- list(rownames(mesa.data.raw$lax.conc),
                     colnames(mesa.data.raw$lax.conc),
                     "lax.conc.1500")
##create STdata object
mesa.data <- createSTdata(obs, mesa.data.raw$X,
                          SpatioTemporal=ST)

##optionally just give a list
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)
##works just as fine with the list
mesa.data.list <- createSTdata(obs, mesa.data.raw$X,
                               SpatioTemporal=ST.list)

##add  the smooth trends
mesa.data <- updateSTdataTrend(mesa.data, n.basis=2)

##This should be the same as the data in data(mesa.data)
}
\references{
  M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L.
  Curl, T.  Gould, D. Hardie, A. Ho, P. Kinney, T. V.
  Larson, P. D. Sampson, L.  Sheppard, K. D. Stukovsky, S.
  S. Swan, L. S. Liu, J. D. Kaufman. (2009) Approach to
  Estimating Participant Pollutant Exposures in the
  Multi-Ethnic Study of Atherosclerosis and Air Pollution
  (MESA Air). Environmental Science & Technology: 43(13),
  4687-4693.
}
\seealso{
  \code{\link{createSTdata}} for creation of \code{STdata}
  objects.

  Other data matrix: \code{\link{createDataMatrix}},
  \code{\link{SVDmiss}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}}

  Other example data: \code{\link{est.cv.mesa}},
  \code{\link{est.mesa.model}},
  \code{\link{MCMC.mesa.model}}, \code{\link{mesa.data}},
  \code{\link{mesa.model}}, \code{\link{pred.mesa.model}}
}
\keyword{datasets}

