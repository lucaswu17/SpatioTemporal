\name{mesa.data.model}
\encoding{latin1}
\Rdversion{1.1}
\alias{mesa.data.model}
\docType{data}
\title{
  Example of a mesa.data.model Structure
}
\description{
  Example of an auxiliary data structure holding information regarding
  which geographic and spatio-temporal covariates to use when fitting
  the model. Should always be created by calling \cr
  \code{\link{create.data.model}} on a \code{\link{mesa.data}} structure.
}
\usage{data(mesa.data.model)}
\format{
  A list with the following elements:
  \describe{
    \item{obs}{A data.frame with the observations, containing:
      \describe{
	\item{\code{obs}}{The value of each observation.}
	\item{\code{date}}{The time at which each observation was taken
	  (as class \code{\link{Date}}).}
	\item{\code{ID}}{The name of the station at which each
	  observation was taken.}
	\item{\code{idx}}{The location at which the observation were
	  taken, this acts as an index into the \cr \code{location}
	  elements, e.g. \code{idx=3} indicates that the observation was
	  taken at the third site in \code{location}. Thus, \cr
	  \code{ID == mesa.data.model$location[mesa.data.model$obs$idx]}.
	}
      }
      The data.frame is sorted by \code{date}, with ties broken by
      \code{idx}.
    }
    \item{location}{A data.frame with the names of each site in \code{ID} and
      the location in \code{x} and \code{y}.  Could also contain longitudinal 
      coordinates for plotting the locations on a map.
    }
    \item{trend}{The smooth temporal trends from \code{mesa.data$trend},
      possibly interpolated using \cr \code{\link{spline}} to match all
      dates in \code{mesa.data$trend$date} and
      \code{mesa.data$obs$date}.
    }
    \item{F}{A matrix contaning the intercept (first column) and smooth
      temporal trends for each observation. The number of rows equals the
      total number of observations. This is an expansion of the data
      in \code{mesa.data$trend}.
    }
    \item{X}{The geographic covariates to be used for model
      fitting. Consists of a list that contains one covariate matrix for
      each termporal trend (including the intercept). Each matrix consists
      of a subset of covariates from \code{mesa.data$LUR}, with an
      added column of ones for the intercept. For details on the selection
      of covariates see \code{\link{create.data.model}}.
    }
    \item{dist}{
      Distance matrix giving the distance between each observed
      location. Computed from \code{x} and \code{y} in
      \code{mesa.data$location}.
    }
    \item{dates}{The unique dates/timepoints of the observations,
      i.e. \cr \code{sort(unique(mesa.data$obs$date))}.
    }
    \item{nt}{The number of observations at each date/timepoint,
      same length as \code{date} above.
    }
    \item{SpatioTemp}{A matrix with the spatio-temporal covariate(s)
      for each observation, compare to \code{F}.

      \code{NULL} if no spatio-temporal covariate is used.
    }
    \item{SpatioTemp.all}{The full 3D-array of spatio-temporal
      covariate(s) for all space-time locations given by \code{location}
      and \code{trend$date}.

      \code{NULL} if no spatio-temporal covariate is used.
    }
    \item{LUR.list}{The \code{LUR} parameter requested in the call to
      \code{\link{create.data.model}} that created this data structure.
    }
    \item{ST.Ind}{The \code{ST.Ind} parameter requested in the call to
      \code{\link{create.data.model}} that created this data structure.
    }
  }
}
\details{  
  This data structure is required as input to most of
  the functions in this package, see \code{\link{loglike}},
  \code{\link{fit.mesa.model}}, \code{\link{cond.expectation}},
  \code{\link{simulateMesaData}}, \cr \code{\link{run.MCMC}}, etc.
}
\section{Warning}{
  This structure should always be created by calling
  \code{\link{create.data.model}} on a \code{\link{mesa.data}}
  structure. The user should \strong{never} attempt to create the
  structure from scratch.
}
\source{
  Contains monitoring data from the \strong{MESA Air} project, see
  Cohen et.al. (2009) and \code{\link{mesa.data.raw}} for details.
  
  Obtained from \code{\link{mesa.data}} as:\cr
  \code{mesa.data.model <- create.data.model(mesa.data,}\cr
  \code{LUR = list(c("log10.m.to.a1", "s2000.pop.div.10000",}\cr
  \code{"km.to.coast"),"km.to.coast", "km.to.coast"), ST.Ind=NULL)}\cr
  See \code{\link{mesa.data}} and \code{\link{create.data.model}} for
  details.
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

For other data functions, see \code{\link{mesa.data.raw}},
\code{\link{create.data.matrix}}, \cr \code{\link{construct.LUR.basis}}
and \code{\link{default.LUR.list}}.
}
\examples{
##load the data
data(mesa.data.model)

##examine components
names(mesa.data.model)
printMesaDataNbrObs(mesa.data.model)

##requested geographic and spatio-temporal covariates
mesa.data.model$LUR.list
mesa.data.model$ST.Ind

##covariates for the temporal intercept
head(mesa.data.model$X$const)
##...and the two smooth temporal trends
head(mesa.data.model$X$V1)
head(mesa.data.model$X$V2)

##Some important dimensions of the model
loglike.dim(mesa.data.model)
}
\keyword{datasets}
