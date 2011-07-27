\name{mesa.data}
\encoding{latin1}
\Rdversion{1.1}
\alias{mesa.data}
\docType{data}
\title{
  Example of a mesa.data Structure
}
\description{
  Example of a data structure holding observations, observation
  locations, smooth temporal trends, geographic covariates and
  spatio-temporal covariates.
}
\usage{data(mesa.data)}
\format{
  A list with elements, a detailed description of each elements is given
  below (see details):
  \describe{
    \item{location}{A data.frame containing names, locations and possibly
      other information about the observation locations.
    }
    \item{LUR}{A data.frame with geographic covariates (land use
      regressors) for each of the observations locations. The rownames of
      LUR should match the names in \code{location$ID} so that the
      covariates can be matched to locations.
    }
    \item{trend}{A data.frame containing the smooth temporal trends.
    }
    \item{obs}{A data.frame containing the observations, along with
      information regarding where and when the observation was taken.
    }
    \item{SpatioTemp}{An array containing the spatio-temporal covariates,
      set to \code{NULL} if no covariates exist.
    }
  }
}
\details{
  This data structure provides an example of how the input data should
  be organised. The \cr \code{mesa.data} object is used as input to
  \code{\link{create.data.model}} which creates a \cr
  \code{\link{mesa.data.model}} that can be used as input to most of the
  functions in this package, \cr
  e.g. \code{\link{loglike}}, \code{\link{fit.mesa.model}},
  \code{\link{cond.expectation}}, \code{\link{run.MCMC}}, \cr
  \code{\link{simulateMesaData}}, etc.
  
  Each of the components contains a number of elements. Additional
  elements can be added to the end of each component.  Some common
  options are described below.

  \describe{
    \item{\code{location}}{Should contain:
      \describe{
	\item{ID}{A factor (or character) containing unique names for
	  each of the observation locations.}
	\item{x, y}{X- and Y-coordinates for the the observation
	  locations.}
      }
      The \code{location} element \emph{may also} contain:
      \describe{
	\item{long, lat}{Longitude and latitude for the the
	  (observation) locations; often needed to plot the locations on
	  a map.} 
	\item{type}{A factor describing the type of location, often
	  helpful when trying to distinguish different types of 
	  observations/locations.}
      }
    }
    \item{\code{LUR}}{Should contain as many geographic
      covariates as makes sense. The \code{rownames} of \code{LUR}
      should match the names in \code{location$ID} so that geographic
      covariates can be related to the correct monitor location. Note
      that all of these covariates are not necessarily used when fitting
      the model, instead \code{\link{create.data.model}} is used to
      select covariates.
    }
    \item{\code{trend}}{Should contain the temporal trends with
      \emph{one of the} columns of the data.frame being named
      \code{date} of class \code{\link{Date}} providing the time
      alignment for the temporal trends. Thus this \cr data.frame
      contains one column for each temporal trend plus one additional
      \code{date} column.
    }
    \item{\code{obs}}{Should contain:
      \describe{
	\item{obs}{The value of each observation.}
	\item{date}{The time at which the corresponding observation was
	  taken (as class \code{\link{Date}}).
	}
	\item{ID}{A factor (or character) giving location at which each 
	  observation was taken. These should be names that match the
	  elements in \code{location$ID}.
	}
      }
      The \code{obs} element \emph{could also} contain additional
      information about each observation.
    }
    \item{\code{SpatioTemp}}{Should contain a 3D-array of spatio-temporal
      covariates, or equal \code{NULL} if no spatio-temporal covariates
      are given. The array should be (number of timepoints) - by  -
      (number of locations) - by - (number of covariates) and provide
      spatio-temporal covariates for \emph{all} space-time locations,
      even unobserved ones (this is needed for prediction).

      The \code{rownames} of the array should represent dates and
      \code{colnames} should match the observation location names in
      \code{location$ID}. Not all of the spatio-temporal covariates are
      necessarily used when fitting the model, instead
      \code{\link{create.data.model}} is used to select spatio-temporal
      covariates.
    }
  }
}
\source{
  Contains monitoring data from the \strong{MESA Air} project, see
  Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
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
 See \code{\link{mesa.data.model}} and \code{\link{create.data.model}}.

For other data functions, see \code{\link{mesa.data.raw}},
\code{\link{create.data.matrix}}, \cr
\code{\link{construct.LUR.basis}}, and \code{\link{default.LUR.list}}.
}
\examples{
##load the data
data(mesa.data)

##Lets look at the data
names(mesa.data)

##Study the structure of the location data
head(mesa.data$location)

##...the covariates
head(mesa.data$LUR)

##...the smooth temporal trends
head(mesa.data$trend)

##...observations
head(mesa.data$obs)

##...and Spatio-temporal covariate
head(mesa.data$SpatioTemp)

##Look at the number of observations/locations
printMesaDataNbrObs(mesa.data)

##Let's plot the space-time monitoring locations
plotMonitoringLoc(mesa.data)

##Let's plot the observations as a function of time
plotMonitoringLoc(mesa.data, obsplot=TRUE, legend.loc="bottomleft")

##plot observations and residuals from the temporal trends
par(mfcol=c(3,2),mar=c(2.5,2.5,2,.5))
plotMesaData(mesa.data,5,type="obs")
plotMesaData(mesa.data,5,type="res")
plotMesaData(mesa.data,5,type="acf")
plotMesaData(mesa.data,18,type="obs")
plotMesaData(mesa.data,18,type="res")
plotMesaData(mesa.data,18,type="acf")
}
\keyword{datasets}
