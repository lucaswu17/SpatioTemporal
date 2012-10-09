##' Package for spatio-temporal modelling. Contains functions that estimate,
##' simulate and predict from the model described in (Szpiro et.al., 2010;
##' Sampson et.al., 2011; Lindström et.al., 2010). The package also
##' contains functions that handle missing data SVD in accordance with
##' (Fuentes et.al. 2006).
##' \cr
##' \tabular{ll}{
##'   Package: \tab SpatioTemporal\cr
##'   Type: \tab Package\cr
##'   Version: \tab 1.1.2\cr
##'   Date: \tab 2012-10-09\cr
##'   License: \tab GPL version 2 or newer\cr
##'   LazyLoad: \tab yes\cr
##' }
##' Examples in the package uses data from the Multi-Ethnic Study of
##' Atherosclerosis and Air Pollution (MESA Air), (Cohen et.al.,2009).
##' 
##' @title Spatio-Temporal Modelling
##' 
##' @name SpatioTemporal-package
##' @aliases SpatioTemporal-package SpatioTemporal
##' @docType package
##' @note Data used in the examples has been provided by the Multi-Ethnic Study
##' of Atherosclerosis and Air Pollution (MESA Air). Details regarding the data
##' can be found in Cohen et.al. (2009).
##' \cr
##' Although the research described in this article has been funded wholly or in
##' part by the United States Environmental Protection Agency through assistance
##' agreement CR-834077101-0 and grant RD831697 to the University of Washington,
##' it has not been subjected to the Agency's required peer and policy review
##' and therefore does not necessarily reflect the views of the Agency and no
##' official endorsement should be inferred.
##' \cr
##' Travel for J. Lindström has been paid by STINT (The
##' Swedish Foundation for International Cooperation in Research and Higher
##' Education) Grant IG2005-2047.
##' \cr
##' Additional funding was provided by grants to the University of Washington
##' from the Health Effects Institute (4749-RFA05-1A/06-10) and the National
##' Institute of Environmental Health Sciences (P50 ES015915).
##'
##' @section Changelog:
##' \describe{
##'   \item{1.1.2}{Minor changes/Bug fixes:}
##'   \itemize{
##'     \item{Fixed a bug in \code{\link{SVDsmooth}}, that caused the values in
##'           the temporal smooths to depend on the number of \emph{unobserved
##'           time points.}. This \emph{also affects}
##'           \code{\link{calcSmoothTrends}} and \code{\link{updateSTdataTrend}}
##'           when the option \code{extra.dates} is in use.}
##'     \item{Fixed bug in \code{\link{simulate.STmodel}} that caused \code{NA}
##'           values when simulating at unobserved sites.}
##'     \item{Fixed bug in \code{\link{predict.STmodel}} that could cause
##'           errors when predicting at unobserved sites.}
##'     \item{Fixed bug in \code{\link{predictCV.STmodel}} and
##'           \code{\link{predict.STmodel}}; these will now handle predictions
##'           at locations with incomplete nugget covariates.}
##'     \item{Updated \code{\link{c.STmodel}} and \code{\link{predict.STmodel}}
##'           to avoid errors/warnings due to more complex nugget models.}
##'     \item{Replaced warning in \code{\link{createSTdata}} when
##'           \code{extra.dates!=NULL} and \code{n.basis=NULL} with a message.}
##'   }
##'   \item{1.1.1}{Bug fixes:}
##'   \itemize{
##'     \item{\code{\link{c.STmodel}} will now combine \code{STmodel}
##'           objects with identical covariate scaling.} 
##'   }
##'   \item{1.1.0}{Major Changes:}
##'   \itemize{
##'     \item{Changed the return of the variances for \code{beta} in
##'           \code{\link{predict.STmodel}}.}
##'     \item{Reduced the memory footprint of \code{\link{predict.STmodel}}.´}
##'     \item{Error checks in \code{\link{c.STmodel}} and
##'           \code{\link{predict.STmodel}}, combination of \code{STmodel}
##'           objects with different covariate scaling is \strong{NOT}
##'           possible.} 
##'   }
##'   \item{1.0.7}{Added:}
##'   \itemize{
##'     \item{New plot function: \code{\link{plot.predCVSTmodel}}.}
##'     \item{\code{\link{coef.estimateSTmodel}} and
##'           \code{\link{coef.estCVSTmodel}} functions that extract estimated
##'            parameters.}
##'     \item{Parameters for \code{\link{predict.STmodel}} and
##'           \code{\link{predictCV.STmodel}} can be specified using
##'           \code{estimateSTmodel} or \code{estCVSTmodel} objects.}
##'     \item{An \code{lwd} option to \code{\link{plot.predictSTmodel}}.}
##'     \item{A short introductory vignette as complement to the full tutorial.}
##'   }
##'   \item{1.0.6}{Bug fixes:}
##'   \itemize{
##'     \item{\code{\link{predictNaive}} now works for only one locations.}
##'     \item{\code{\link{detrendSTdata}} now works for different regions.}
##'   }
##'   \item{1.0.5}{Added packages \code{maps} and \code{plotrix} to suggested packages.}
##'   \item{1.0.4}{Bug fixes:}
##'   \itemize{
##'     \item{prediction for leave-one-out CV.}
##'     \item{stop \link{updateCovf} crashing in Rscript/R CMD BATCH.}
##'   }
##'   \item{1.0.3}{Minor bug fixes}
##'   \item{1.0.2}{Updated documentation and vignette}
##'   \item{1.0.0}{Major change, most old functions are now deprecated. New features:}
##'   \itemize{
##'     \item{Different covariance functions}
##'     \item{Nuggets in the beta-fields}
##'     \item{Different nuggets for different locations in the nu-field.}
##'     \item{Different coordinates for beta and nu-fields, allowing for
##'           precomputed deformations}
##'     \item{Covariates can be specifed using
##'           \link[stats:formula]{formula}-objects}
##'   }
##'   \item{0.9.2}{Minor updates - no user visible changes}
##'   \item{0.9.0}{First CRAN-release}
##'   \item{0.1.0}{First released version, short course at TIES-2010}
##' }
##' 
##' @author Johan Lindström, Adam Szpiro, Paul D. Sampson,
##' Silas Bergen, Assaf P. Oron
##' 
##' @references
##' M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T.
##'   Gould, D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L.
##'   Sheppard, K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009)
##'   Approach to Estimating Participant Pollutant Exposures in the Multi-Ethnic
##'   Study of Atherosclerosis and Air Pollution (MESA Air). Environmental Science
##'   & Technology: 43(13), 4687-4693.
##' 
##' M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
##'  Analyze Space-Time Processes in Statistical methods for spatio-temporal
##'  systems (B. Finkenstädt, L. Held, V. Isham eds.) 77-150
##' 
##' J. Lindström, A. Szpiro, P. D. Sampson, L. Sheppard, A. Oron,
##'   M. Richards, and T. Larson T. (2010) A flexible spatio-temmporal model for
##'   air pollution: allowing for spatio-temporal covariates. Berkeley Electronic
##'   Press, University of Washington Biostatistics Working Paper Series, No. 370.
##'   \url{http://www.bepress.com/uwbiostat/paper370}
##'
##' A. Szpiro, P. D. Sampson, L. Sheppard, T. Lumley, S. D. Adar, and J. D.
##'   Kaufman. (2010) Predicting intra-urban variation in air pollution
##'   concentrations with complex spatio-temporal dependencies. Environmetrics:
##'   21, 606-631.
##'
##' P. D. Sampson, A. Szpiro, L. Sheppard, J. Lindström, J. D.  Kaufman. (2011)
##'   Pragmatic Estimation of a Spatio-temporal Air Quality Model with Irregular
##'   Monitoring Data. Atmospheric Environment: 45(36), 6593-6606.
##' 
##' @keywords package
##' @example Rd_examples/Ex_SpatioTemporal_package.R
NULL
