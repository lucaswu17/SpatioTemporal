\name{mesa.data.raw}
\encoding{latin1}
\Rdversion{1.1}
\alias{mesa.data.raw}
\docType{data}
\title{
  Example of raw data
}
\description{
  The raw data that was used to create the \code{\link{mesa.data}} and
  \code{\link{mesa.data.model}} structures.
}
\usage{data(mesa.data.raw)}
\format{
  The structure contains observations, temporal trends, locations,
  geographic covariates, and spatio-temporal covariates. The data is
  stored in a list with elements: 
  \describe{
    \item{X}{A data.frame containing names, locations, and (geographic)
      covariates regarding all the (observation) locations.
    }
    \item{obs}{A time-by-location matrix for the observed data.
    }
    \item{trend}{A data.frame containing the smooth temporal trends.
    }
    \item{lax.conc.1500}{A time-by-location matrix of a spatio-temporal
      covariate based on output from Caline3QHC.
    }
  }
}
\details{
  The data structure contains raw data from the \strong{MESA Air}
  project. The example below describes how to create the
  \code{\link{mesa.data}} structure from raw data.
}
\source{
  Contains monitoring data from the \strong{MESA Air} project, see
  Cohen et.al. (2009) for details.
}
\references{
  M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T. Gould,
  D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L. Sheppard,
  K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009) Approach
  to Estimating Participant Pollutant Exposures in the Multi-Ethnic
  Study of Atherosclerosis and Air Pollution (MESA Air). Environmental
  Science & Technology: 43(13), 4687-4693.
}
\seealso{
 See \code{\link{mesa.data}} and \code{\link{create.data.model}}.

For other data functions, see \code{\link{mesa.data.model}},
\code{\link{create.data.matrix}}, \cr \code{\link{construct.LUR.basis}}
and \code{\link{default.LUR.list}}.
}
\examples{
##load the data
data(mesa.data.raw)

##create the location data.frame
location <- mesa.data.raw$X[,c("ID","x","y","long","lat","type")]
##ensure that it's a data.frame
location <- as.data.frame(location)
##ensure that ID and ID are factors
location$ID <- factor(as.character(location$ID))
location$type <- factor(as.character(location$type))

##extract the covariate LUR matrix
LUR <- mesa.data.raw$X
##ensure that it's a data.frame
LUR <- as.data.frame(LUR)
##set rownames to the IDs and drop location information
rownames(LUR) <- as.character(LUR$ID)
LUR <- LUR[,!(names(LUR) \%in\% c("ID","x","y","long","lat","type"))]

##create matrices of observation times
T <- t(matrix(rownames(mesa.data.raw$obs),
              nrow=dim(mesa.data.raw$obs)[2],
              ncol=dim(mesa.data.raw$obs)[1],byrow=TRUE))
##...and locations
ID <- matrix(colnames(mesa.data.raw$obs), nrow=dim(mesa.data.raw$obs)[1],
             ncol=dim(mesa.data.raw$obs)[2],byrow=TRUE)
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

##compute the smooth trends
trend <- calc.smooth.trends(obs=obs$obs, date=obs$date,
                            ID=obs$ID, n.basis = 2)$svd

##create a 3D-array for the spatio-temporal covariate
ST <- array(mesa.data.raw$lax.conc.1500, dim =
            c(dim(mesa.data.raw$lax.conc.1500),1))
##add names for time, location and which spatio-temporal covariate
dimnames(ST) <- list(rownames(mesa.data.raw$lax.conc),
                     colnames(mesa.data.raw$lax.conc),
                     "lax.conc.1500")

##combine to a mesa.data-object.
mesa.data.alt <- list(location=location, LUR=LUR, trend = trend,
                      obs=obs, SpatioTemp=ST)

##This should be the same as the data in mesa.data
}
\keyword{datasets}
