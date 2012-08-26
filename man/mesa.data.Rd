\docType{data}
\name{mesa.data}
\alias{mesa.data}
\title{Example of a \code{STdata} Structure}
\format{A list with elements, a detailed description of each elements is
given below (see details):
  \describe{
    \item{covars}{A data.frame with geographic covariates, including location
      should include a column \code{covars$ID} that gives the names of
      observations locations, matching those in \code{obs$ID}.}
    \item{trend}{A data.frame containing the smooth temporal trends.}
    \item{obs}{A data.frame containing the observations, along with
       information regarding where and when the observation was taken.}
    \item{SpatioTemporal}{A 3D array containing the spatio-temporal covariates;
       \code{NULL} if no covariates exist.}
  }}
\source{
  Contains monitoring data from the \strong{MESA Air}
  project, see Cohen et.al. (2009) and
  \code{\link{mesa.data.raw}} for details.
}
\description{
  Example of a data structure holding observations,
  geographic covariates (including observation locations),
  smooth temporal trends, and spatio-temporal covariates.
}
\details{
  A \code{STdata} object consists of a list with, some or
  all of, the following elements: \describe{
  \item{covars}{Geographic covariates, locations and names
  of the observation locations (the later in
  \code{covars$ID}), \code{\link{createSTmodel}} will
  extract covariates (land use regressors), observations
  locations, etc from this data.frame when constructing the
  model specification.} \item{trend}{The temporal trends
  with \emph{one of the} columns being named \code{date},
  preferably of class \code{\link[base:Date]{Date}}
  providing the time alignment for the temporal trends.}
  \item{obs}{A data.frame with columns: \describe{
  \item{obs}{The value of each observation.}
  \item{date}{The observations time, preferably of class
  \code{\link[base:Date]{Date}}.} \item{ID}{A
  \code{character}-class giving observation locations;
  should match elements in \code{locations$ID}.} } }
  \item{SpatioTemporal}{A 3D-array of spatio-temporal
  covariates, or \code{NULL} if no covariates exist. The
  array should be (number of timepoints) - by - (number of
  locations) - by - (number of covariates) and provide
  spatio-temporal covariates for \emph{all} space-time
  locations, even unobserved ones (needed for prediction).
  The \code{rownames} of the array should represent
  dates/times and \code{colnames} should match the
  observation location names in \code{covars$ID}.}
  \item{old.trend,fit.trend}{Additional components added if
  the observations have been detrended, see
  \code{\link{detrendSTdata}}.} }
}
\examples{
##load the data
data(mesa.data)

##Look at the number of observations/locations
print(mesa.data)

##Look at the summary of observations, covariates, etc.
summary(mesa.data)

##Lets look at the data
names(mesa.data)

##Study the structure of the covariates data
head(mesa.data$covars)

##...the smooth temporal trends
head(mesa.data$trend)

##...observations
head(mesa.data$obs)

##...and Spatio-temporal covariate
mesa.data$SpatioTemporal[1:5,1:5,,drop=FALSE]

##Let's plot the space-time monitoring locations
plot(mesa.data, "loc")

##Let's plot the observations as a function of time
plot(mesa.data, "loc.obs", legend.loc="bottomleft")

##plot observations and residuals from the temporal trends
par(mfcol=c(3,2),mar=c(2.5,2.5,2,.5))
plot(mesa.data, "obs", ID=5)
plot(mesa.data, "res", ID=5)
plot(mesa.data, "acf", ID=5)
plot(mesa.data, "obs", ID=18)
plot(mesa.data, "res", ID=18)
plot(mesa.data, "acf", ID=18)

##create STmodel object
##define land-use covariates, for intercept and trends
LUR = list(c("log10.m.to.a1", "s2000.pop.div.10000", "km.to.coast"),
  "km.to.coast", "km.to.coast")
##and covariance model
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=TRUE, random.effect=FALSE)
##which locations to use
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")
##create object
mesa.model <- createSTmodel(mesa.data, LUR=LUR, cov.beta=cov.beta,
                            cov.nu=cov.nu, locations=locations)

##This should be the same as the data in data(mesa.model)
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
  \code{\link{createSTmodel}} for creation of
  \code{STmodel} objects. \cr \code{\link{createSTdata}}
  for creation of \code{STdata} objects.

  Other example data: \code{\link{est.cv.mesa}},
  \code{\link{est.mesa.model}},
  \code{\link{MCMC.mesa.model}},
  \code{\link{mesa.data.raw}}, \code{\link{mesa.model}},
  \code{\link{pred.mesa.model}}
}
\keyword{datasets}

