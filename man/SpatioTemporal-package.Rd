\name{SpatioTemporal-package}
\encoding{latin1}
\Rdversion{1.1}
\alias{SpatioTemporal-package}
\alias{SpatioTemporal}
\docType{package}
\title{
  Spatio-Temporal Modelling
}
\description{
  Package for spatio-temporal modelling. Contains functions that
  estimate, simulate and predict from the model described in (Szpiro
  et.al., 2010; \enc{Lindström}{Lindstrom} et.al., 2010). The package
  also contains functions that handle missing data SVD in accordance
  with (Fuentes et.al. 2006).
}
\details{
  \tabular{ll}{
    Package: \tab SpatioTemporal\cr
    Type: \tab Package\cr
    Version: \tab 0.9.2\cr
    Date: \tab 2012-03-05\cr
    License: \tab GPL version 2 or newer\cr
    LazyLoad: \tab yes\cr
  }
  Estimation, prediciton and cross-validation routines for the
  spatio-temporal model described in (Szpiro et.al., 2010;
  \enc{Lindström}{Lindstrom} et.al., 2010)

  Examples in the package uses data from the Multi-Ethnic Study of
  Atherosclerosis and Air Pollution (MESA Air), see Cohen et.al. (2009)
  for details).
}
\author{
\enc{Johan Lindström}{Johan Lindstrom}, Adam Szpiro, Paul D. Sampson,
Silas Bergen, Assaf P. Oron

Maintainer: \enc{Johan Lindström}{Johan Lindstrom} <johanl@maths.lth.se>
}
\note{
  Data used in the examples has been provided by the Multi-Ethnic
  Study of Atherosclerosis and Air Pollution (MESA Air). Details
  regarding the data can be found in Cohen et.al. (2009).

  Although the research described in this article has been funded 
  wholly or in part by the United States Environmental Protection 
  Agency through assistance agreement CR-834077101-0 and grant 
  RD831697 to the University of Washington, it has not been 
  subjected to the Agency's required peer and policy review and 
  therefore does not necessarily reflect the views of the Agency 
  and no official endorsement should be inferred.

  Travel for \enc{J. Lindström}{J. Lindstrom} has been paid by STINT (The
  Swedish Foundation for International Cooperation in Research and Higher
  Education) Grant IG2005-2047.

  Additional funding was provided by grants to the University of 
  Washington from the Health Effects Institute 
  (4749-RFA05-1A/06-10) and the National Institute of 
  Environmental Health Sciences (P50 ES015915).
}
\references{
  Bergen, S. and \enc{Lindström, J.}{Lindstrom, J.}. (2011)
  Tutorial for Spatio-Temporal R Package
  
  M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T. Gould,
  D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L. Sheppard,
  K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009) Approach
  to Estimating Participant Pollutant Exposures in the Multi-Ethnic
  Study of Atherosclerosis and Air Pollution (MESA Air). Environmental
  Science & Technology: 43(13), 4687-4693.
  
  M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
  Analyze Space-Time Processes in Statistical methods for
  spatio-temporal systems (B. \enc{Finkenstädt}{Finkenstadt}, L. Held,
  V. Isham eds.) 77-150

  J. \enc{Lindström}{Lindstrom}, A. Szpiro, P. D. Sampson, L. Sheppard,
  A. Oron, M. Richards, and T. Larson T. (2010) A flexible
  spatio-temmporal model for air pollution: allowing for spatio-temporal
  covariates. Berkeley Electronic Press, University of Washington
  Biostatistics Working Paper Series, No. 370.
  
  A. Szpiro, P. D. Sampson, L. Sheppard, T. Lumley, S. D. Adar, and
  J. D. Kaufman. (2010) Predicting intra-urban variation in air
  pollution concentrations with complex spatio-temporal
  dependencies. Environmetrics: 21, 606-631.
}
\keyword{ package }
\examples{
##For a worked out data-analysis see the tutorial.
\dontrun{
vignette("Tutorial",package="SpatioTemporal")
}
}%\examples
