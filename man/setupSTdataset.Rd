\name{setupSTdataset}
\encoding{latin1}
\Rdversion{1.1}
\alias{setupSTdataset}
\title{
  Creates the data structure that is then used to build specific models via \code{\link{create.data.model}}.
  }
\description{
 Creates the data structure that is then used to build specific models via \code{\link{create.data.model}}.
 This includes observations, observation
  locations, smooth temporal trends, geographic covariates and spatio-temporal covariates (in the future).}
\usage{
setupSTdataset(rawobs, covardat, covarnames, trendf,
               varnames=list(yraw="lac",date="intended_wednesday",idobs="site_id",idcov="site.id",xcoord="lambert.x",ycoord="lambert.y",long="longitude",lat="latitude"), x.to.km=1000,transform=log,scale=TRUE,mesa=TRUE)
}
\arguments{
  \item{rawobs}{
    Data frame with the raw observations. Should have vectors whose
    names are identified by the items ``yraw'', ``date'',``idobs'' in
    \code{varnames}.}
  \item{covardat}{
    Data frame with the regression predictors. Should have vectors whose
    names are identified by the items ``idcov'', ``xcoord'', ``ycoord'',
    ``long'', ``lat'' in \code{varnames}.}
  \item{covarnames}{Character vector with the names of variables in
    \code{covardat} that should be included in the modeling dataset. If
    you are not sure, then just specify all names (i.e.,
    ``names(covardat)'').}
  \item{trendf}{Matrix or data frame with the time trend, on the
    modeling scale. Its format should be compatible with the output of
    \code{SVD.smooth}: columns represent the smooth time trends, and the
    row names are in the R \code{\link{Date}} format.}
  \item{varnames}{List of character strings denoting the variable names
    needed to identify and match the various data components. See
    above.}
  \item{x.to.km}{Numeric, conversion factor from the coordinate scale to
    km. Defaults to 1000. Set to 1 to ignore.} 
  \item{transform}{Function, the transformation link from the scale in
    which observations appear in \code{rawobs}, to the modeling
    scale. Defaults to \code{\link{log}}.} 
  \item{scale}{Logical: should the covariates in the LUR component be
    scaled to each have mean 0 and variance 1? Defaults to \code{TRUE}.}
  \item{mesa}{Logical: should we assume that location ID names follow
    the MESA-Air convention? Defaults to \code{TRUE}.}
}
 
\details{ 

This function is a time-saving aid, for creating the basic data structure for fitting spatiotemporal models, with some built-in error checks.

The output of this function can be fed as the \code{\link{mesa.data}} argument into \code{\link{create.data.model}}, to build specific models, or to prediction functions.

The function separates out the geographic information in \code{covardat} into the ``location'' and ``LUR'' components, ensuring that all location IDs with observations have location and covariate data. It also ensures that all observation dates have a trend function value -- if some do not exist, it calculates them using \code{\link{smooth.spline}}, giving out a warning.

If you have a high-dimensional covariate set but plan to use only a small subset, then \code{covarnames} helps produce a more compact data structure. Specify only the variables you really think you might experiment with. The \code{scale} option is recommended when building and checking the model, because it improves numerical stability of \code{\link{fit.mesa.model}}, and makes the effect estimate magnitudes comparable. When producing a final model for scientific use, you might want to set \code{scale} to \code{FALSE}, and manually divide the covariates into units that are scientifically communicated, but not too far in magnitude from unity.

Be sure that your variable names are correctly specified in \code{varnames}. 

If you are struggling with this function, you can always construct a \code{\link{mesa.data}} list manually by yourself, as explained in the package tutorial.

}

\seealso{
 \code{\link{mesa.data}}, \code{\link{mesa.data.model}}, \code{\link{create.data.model}}.
 }
 
 \author{Assaf P. Oron}
\note{At present, this function does not handle spatio-temporal covariates. If you have them in your model, you can run this function and then add the spatio-temporal array manually.}

\value{A list with data frames described in \code{\link{mesa.data}}, except for the ``SpatioTemp'' array. }

