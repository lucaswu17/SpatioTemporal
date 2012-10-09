\name{estimateCV.STmodel}
\alias{estimateCV.STmodel}
\alias{predictCV.STmodel}
\title{Cross-Validated Estimation and Prediction}
\usage{
  \method{estimateCV}{STmodel} (object, x, Ind.cv,
    control = list(trace = 3), ...)

  \method{predictCV}{STmodel} (object, x, Ind.cv, ...,
    silent = TRUE)
}
\arguments{
  \item{object}{\code{STmodel} object for which to perform
  cross-validation.}

  \item{x}{Either a vector or matrix of starting point(s)
  for the optimisation, see \code{\link{estimate.STmodel}};
  or a matrix with parameters, the i:th row being used for
  prediction of the i:th cross-validation set. For
  prediction either a \code{estCVSTmodel} or
  \code{estimateSTmodel} object, results from
  \code{\link{estimateCV.STmodel}} or
  \code{\link{estimate.STmodel}}, can be used.}

  \item{Ind.cv}{\code{Ind.cv} defines the cross-validation
  scheme.  Either a (number or observations) - by -
  (groups) logical matrix or an \emph{integer valued}
  vector with length equal to (number or observations). See
  further \code{\link{createCV}}.}

  \item{control}{A list of control parameters for the
  optimisation.  See \code{\link[stats:optim]{optim}} for
  details; setting \code{trace}=0 eliminates all ouput.}

  \item{...}{All additional parameters for
  \code{\link{estimate.STmodel}} or
  \code{\link{predict.STmodel}}.  For
  \code{\link{predict.STmodel}} a number of parameters are
  set in \code{predictCV.STmodel} and can \strong{NOT} be
  overriden, these are \code{nugget.unobs},
  \code{only.pars=FALSE}, and \code{combine.data=FALSE}.}

  \item{silent}{Show status after each iteration?}
}
\value{
  Either a \code{estCVSTmodel} object with elements:
  \item{status}{Data.frame with convergence information and
  best function value for each cross-validation group.}
  \item{Ind.cv}{The cross-validation grouping.}
  \item{x.fixed}{Fixed parameters in the estimation, see
  \code{\link{estimate.STmodel}}.} \item{par.all,
  par.cov}{Matrices with estimated parameters for each
  cross-validation group.} \item{par.all.sd,
  par.cov.sd}{Standard deviations computed from the
  Hessian/information matrix for set of estimated
  parameters.} \item{res.all}{Estimation results for each
  cross-validation group, contains the \code{res.best} and
  \code{summary} fields from the
  \code{\link{estimate.STmodel}} results.} Or a
  \code{predCVSTmodel} object with elements:
  \item{opts}{Copy of the \code{opts} field in the output
  from \code{\link{predict.STmodel}}.} \item{Ind.cv}{The
  cross-validation grouping.} \item{pred.obs}{A data.frame
  with a copy of observations from \code{object$obs},
  predictions (for different model components), variances,
  and residuals. Variance field will be missing if
  \code{pred.var=FALSE}.} \item{pred.all}{A list with
  time-by-location data.frames containing predictions and
  variances for all space-time locations as well as
  predictions and variances for the beta-fields. Unobserved
  points are \code{NA} for the option
  \code{only.obs=TRUE}.}
}
\description{
  Functions that perform cross-validated parameter
  estimation and prediction for the spatio-temporal model.
}
\details{
  For \code{predictCV.STmodel} the parameters used to
  compute predictions for the left out observations can be
  either a single vector or a matrix. For a single vector
  the same parameter values will be used for all
  cross-validation predictions; for a matrix the parameters
  in \code{x[,i]} will be used for the predictions of the
  i:th cross-validation set (i.e. for \code{Ind.cv[,i]}).
  Suitable matrices are provided in the output from
  \code{estimateCV.STmodel}.

  The cross-validation groups are given by \code{Ind.cv}.
  \code{Ind.cv} should be either a (number of observations)
  - by - (groups) logical matrix or an \emph{integer
  valued} vector with length equal to (number of
  observations). If a matrix then each column defines a
  cross-validation set with the \code{TRUE} values marking
  the observations to be left out. If a vector then
  \code{1}:s denote observations to be dropped in the first
  cross-validation set, \code{2}:s observations to be
  dropped in the second set, etc. Observations marked by
  values \code{<=0} are never dropped. See
  \code{\link{createCV}} for details.
}
\examples{
##load data
data(mesa.model)

################
## estimateCV ##
################
##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))

\dontrun{
  ##estimate different parameters for each CV-group
  est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)

##examine the estimation results
print( est.cv.mesa )
##estimated parameters for each CV-group
est.cv.mesa$par.cov

##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(12,2.5,2,.5), las=2)
boxplot( est.cv.mesa, plot.type="cov")

###############
## predictCV ##
###############
\dontrun{
  ##Do cross-validated predictions using the just estimated parameters
  pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, est.cv.mesa$Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)

##examine the prediction results
names( pred.cv.mesa )
print( pred.cv.mesa )
##compute CV-statistics
pred.cv.stats <- summary( pred.cv.mesa, LTA=TRUE)
print( pred.cv.stats )
 
##Plot observations with CV-predictions and prediction intervals
##TODO TODO TODO
#par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
#plotCV(mesa.data.res$pred.cv,  1, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 10, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 17, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 22, mesa.data.model)

##Residual QQ-plots
#par(mfrow=c(1,2),mar=c(3,2,1,1),pty="s")
##Raw Residual QQ-plot
#CVresiduals.qqnorm(pred.cv.stats$res)
##Normalized Residual QQ-plot
#CVresiduals.qqnorm(pred.cv.stats$res.norm, norm=TRUE)

\dontrun{
  ##A faster option is to only consider the observations and not to compute
  ##variances
  pred.cv.fast <- predictCV(mesa.model, est.cv.mesa, est.cv.mesa$Ind.cv,
                            only.obs=TRUE, pred.var=FALSE)
  print( pred.cv.fast )
  summary( pred.cv.fast )
}
}
\author{
  Johan Lindström
}
\seealso{
  Other cross-validation functions:
  \code{\link{computeLTA}}, \code{\link{createCV}},
  \code{\link{dropObservations}},
  \code{\link{predictNaive}}

  Other estCVSTmodel methods:
  \code{\link{boxplot.estCVSTmodel}},
  \code{\link{coef.estCVSTmodel}},
  \code{\link{print.estCVSTmodel}},
  \code{\link{print.summary.estCVSTmodel}},
  \code{\link{summary.estCVSTmodel}}

  Other predCVSTmodel functions: \code{\link{computeLTA}}

  Other predCVSTmodel methods:
  \code{\link{plot.predCVSTmodel}},
  \code{\link{plot.predictSTmodel}},
  \code{\link{print.predCVSTmodel}},
  \code{\link{print.summary.predCVSTmodel}},
  \code{\link{summary.predCVSTmodel}}

  Other STmodel methods: \code{\link{createSTmodel}},
  \code{\link{c.STmodel}}, \code{\link{estimate.STmodel}},
  \code{\link{MCMC.STmodel}}, \code{\link{plot.STdata}},
  \code{\link{plot.STmodel}},
  \code{\link{predict.STmodel}},
  \code{\link{print.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}},
  \code{\link{summary.STmodel}}
}

