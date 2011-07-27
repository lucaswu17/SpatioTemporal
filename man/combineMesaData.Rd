\name{combineMesaData}
\encoding{latin1}
\Rdversion{1.1}
\alias{combineMesaData}

\title{
  Merges two Data Objects.
}
\description{
  Function that merges a \code{\link{mesa.data.model}} with a
  \code{\link{mesa.data}}, creating a new \cr \code{\link{mesa.data.model}}
  object that combines the locations from both input objects.

  This function is mainly for internal use.
}
\usage{
combineMesaData(mesa.data.model, mesa.data)
}
\arguments{
  \item{mesa.data.model}{
    Data structure holding locations, covariates, temporal trends, and
    observations. See further \code{\link{mesa.data.model}}.
  }
  \item{mesa.data}{
    Data structure holding additional locations, geographic covariates,
    and spatio-temporal covariates, smooth temporal trends, and
    \emph{optionally} additional observations. See further
    \code{\link{mesa.data}}.
  }
}
\details{
  Function that combines locations in \code{\link{mesa.data.model}} and
  \code{\link{mesa.data}} creating a new \cr
  \code{\link{mesa.data.model}} object. The temporal trend and
  observations are taken from \cr \code{\link{mesa.data.model}}. Any
  trend and/or observation data in \code{\link{mesa.data}}
  \emph{are ignored}.
  
  This function is used internally by \code{\link{cond.expectation}} and
  \code{\link{simulateMesaData}}.
}
\value{
  Returns a new \code{\link{mesa.data.model}} combining the two inputs.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Used by \code{\link{cond.expectation}} and
  \code{\link{simulateMesaData}}.
}
\examples{
##load the data
data(mesa.data)
data(mesa.data.model)

##keep only observations from the AQS sites
ID.AQS <- mesa.data$location$ID[mesa.data$location$type=="AQS"]
mesa.data$obs <- mesa.data$obs[mesa.data$obs$ID \%in\% ID.AQS,]

##keep the original mesa.data.model
mesa.data.model.org <- mesa.data.model

##and create a model object, dropping unobserved locations.
mesa.data.model <- create.data.model(mesa.data,
    LUR = mesa.data.model$LUR.list, ST.Ind = mesa.data.model$ST.Ind)

##lets study the number of sites and observations in each object
printMesaDataNbrObs(mesa.data)
printMesaDataNbrObs(mesa.data.model)
printMesaDataNbrObs(mesa.data.model.org)

##combine mesa.data and mesa.data.model
printMesaDataNbrObs(combineMesaData(mesa.data.model, mesa.data))
##adding the 5 FIXED sites dropped by create.data.model

##combining mesa.data.model.org and mesa.data produces a warning
##since mesa.data contains no new sites.
tmp <- combineMesaData(mesa.data.model.org, mesa.data)
}
