\name{estimateCV}
\encoding{latin1}
\Rdversion{1.1}
\alias{estimateCV}
\alias{predictCV}

\title{
  Cross-Validated Estimation and Prediction
}
\description{
  Functions that perform cross-validated parameter estimation and
  prediction for the spatio-temporal model.
}
\usage{
estimateCV(par.init = 3, mesa.data.model, Ind.cv, type = "p",
           h = 0.001, diff.type = 1, lower = -15, upper = 15,
           hessian.all = FALSE, control =
           list(trace = 3, maxit = 1000))

predictCV(par, mesa.data.model, Ind.cv, type = "p",express=FALSE, Nmax = 1000,
          silent = TRUE)
}
\arguments{
  \item{par.init}{
    Vector or matrix of starting point(s) for the parameter
    estimation. See further \code{\link{fit.mesa.model}}.
  }
  \item{par}{Vector or matrix of parameters. Cross-validated predictions
    are carried out using these parameter values. If \code{par} is a
    vector then the same parameter values will be used for all
    cross-validation sets. If \code{par} is a matrix the parameters in
    \code{par[,i]} will be used for the predictions of the
    cross-validation set given by \code{Ind.cv[,i]}.
  }
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding the
    observation locations. See \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
  \item{Ind.cv}{
    \code{Ind.cv} defines the cross-validation scheme.
    Either a (number or observations) - by - (groups) logical matrix
    or an \emph{integer valued} vector with length equal to (number or
    observations). See further \code{\link{createCV}}.
  }
  \item{type}{
    A single character indicating the type of log-likelihood to
    use. Valid options are "f", "p", and "r", for \emph{full},
    \emph{profile} or \emph{restricted maximum likelihood}
    (REML).
  }
  \item{express}{
	logical: should an ``express'' version of predictions be returned as a single data frame, without variances? 
	For more datails see under ``Value'' (defaults to \code{FALSE}). 
  }
  \item{h, diff.type}{
    Step length and type of finite difference to use when computing
    gradients. See \code{\link{loglike.grad}} and
    \code{\link{gen.gradient}} for details.
  }
  \item{lower, upper}{
    Bounds on the variables, see \code{\link{optim}}.
  }
  \item{hessian.all}{
    If \code{type!="f"} computes hessian (and uncertainties) for both
    regression and \emph{log}-covariance parameters, not only for
    \emph{log}-covariance parameters. \cr
    See \code{\link{fit.mesa.model}} for details.
  }
  \item{control}{
    A list of control parameters for the optimisation.
    See \code{\link{optim}} for details.
  }
  \item{Nmax}{
    Limits the size of matrices constructed when computing the
    expectations. See \code{\link{cond.expectation}}.
  }
  \item{silent}{
    If \code{FALSE} outputs brief progress information. If \code{TRUE}
    no output.
  }
}
\details{
  For \code{estimateCV} the inital parameters for the optimisation
  can be given as a vector of single starting point, as a matrix of
  multiple starting point, or as a single value. If \code{par.init} is a
  single integer then multiple starting points will be created as
  a set of constant vectors with the values of each starting point taken
  as \code{seq(-5, 5, length.out=par.init)}.

  For \code{predictCV} the parameters used to compute predictions for
  the left out observations can be either a single vector or a matrix
  with \code{dim(par)[2]=dim(Ind.cv)[2]} (or \cr
  \code{dim(par)[2]=max(Ind.cv)}). For a single vector the same
  parameter values will be 
  used for all cross-validation predictions, for a matrix the parameters
  in \code{par[,i]} will be used for the predictions of the
  cross-validation set given by \code{Ind.cv[,i]} (or \code{Ind.cv==i}).
  A suitable matrix is provided in the output from \code{estimateCV}.

  The cross-validation groups are given by \code{Ind.cv}. \code{Ind.cv}
  should be either a (number of observations) - by - (groups) logical
  matrix or an \emph{integer valued} vector with length equal to (number
  of observations). If a matrix then each column defines a
  cross-validation set with the \code{TRUE} values marking the
  observations to be left out. If a vector then \code{1}:s denote
  observations to be dropped in the first cross-validation set,
  \code{2}:s observations to be dropped in the second set,
  etc. Observations marked by values \code{<=0} are never dropped. See
  \code{\link{createCV}} for details.
}
\value{
  \code{estimateCV} returns a list containing:
  \item{par}{
    A (number of parameters) - by - (number of groups) matrix containing
    the estimated parameters for each cross-validation group. This can
    be used as \code{par}-input to \code{predictCV}.
  }
  \item{res.all}{
    A list with (number of groups) elements. Each element is the best
    estimation result for that cross-validation group, \cr
    e.g. \code{fit.mesa.model(...)$res.best}.
  }

  If ``express''=\code{FALSE} (default), \code{predictCV} returns a list containing:
  \item{pred.obs}{
    A data.frame with the observations, predictions and prediction
    variances. The elements in this data.frame are in the same order as
    those in \cr \code{mesa.data.model$obs}.

    \code{NULL} if any observation occurs in more than one
    cross-validation set, \cr
    e.g. \code{max(rowSums(Ind.cv))!=1}.
  }
  \item{pred.by.cv}{
    A list with (number of groups) elements. Each element contains a
    data.frame with the observations, predictions and prediction
    variances for that \cr cross-validation group. The elements in each
    data.frame are in the same order as those in
    \code{mesa.data.model$obs}.

    \code{NULL} if no observation occurs in more than one
    cross-validation set, \cr
    e.g. \code{max(rowSums(Ind.cv))==1}.
  }
  \item{pred.all}{
    A list with (number of groups) elements. Each element contains a
    (number of time points) - by - (number of locations in that group) -
    by - (3) 3D-array. The arrays contain observations, predictions and
    prediction variances for the locations in each each cross-validation
    group.
  }
  
  If ``express''=\code{TRUE}, \code{predictCV} returns a data frame with information pertaining only to 
  observations predicted during the CV cycle (i.e., observations left out and predicted upon at least once).
  That data frame has vectors:
  
  \item{ID}{The location ID, as it appears in the ``location'' component of ``mesa.data.model''}
  \item{date}{The observation's date}
  \item{raw}{The actual observation as it appears in the ``obs'' component of ``mesa.data.model''}
  \item{lur}{The prediction (on the model scale) without any spatial smoothing: see ``EX.mu'' under \code{cond.expectation}}
  \item{krig}{The prediction (on the model scale) of the full model including spatial smoothing: see ``EX'' under \code{cond.expectation}}
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} (\code{predictCV} ``express'' option by Assaf Oron)
}
\seealso{
  See \code{\link{createCV}} for creating CV-schemes.
  
  The functions use \code{\link{fit.mesa.model}} for estimation and
  \code{\link{cond.expectation}} for prediction.
  
For computing CV statistics, see also \code{\link{compute.ltaCV}}, 
\code{\link{predictNaive}}, and for further illustration see
\code{\link{plotCV}}, \code{\link{CVresiduals.qqnorm}}, 
and \code{\link{summaryStatsCV}}.
}
\examples{
##Load data
data(mesa.data.model)
data(mesa.data.res)

##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.data.model, groups=10, min.dist=.1)

\dontrun{
##estimate different parameters for each CV-group
dimm <- loglike.dim(mesa.data.model)
x.init <- as.matrix(cbind(rep(2,dimm$nparam.cov),c(rep(c(1,-3),dimm$m+1),-3)))
##This may take a while...
par.est.cv <- estimateCV(par.est$res.best$par, 
      mesa.data.model, Ind.cv)
}
##Get the pre-computed results instead.
par.est.cv <- mesa.data.res$par.est.cv

##let's examine estimation of the results
names(par.est.cv)
##estimated parameters for each CV-group
par.est.cv$par

##Convergence of the optimisation for each group
sapply(par.est.cv$res.all,function(x) x$conv)

##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(7,2.5,2,.5), las=2)
boxplot(t(par.est.cv$par))

\dontrun{
##Do cross-validated predictions using the newly 
##This may take a while...
pred.cv <- predictCV(par.est.cv$par, mesa.data.model,
                     Ind.cv, silent = FALSE)
}
##Get the pre-computed results instead.
pred.cv <- mesa.data.res$pred.cv
##examine the results
names(pred.cv)

##Alternatively, one can use the ``express'' option which should be faster and returns a compact dataset:
pred0.cv <- predictCV(par.est.cv$par, mesa.data.model,Ind.cv, silent = FALSE,express=TRUE)
#Instead of a list, a single data frame is returned:
head(pred0.cv)



##compute CV-statistics
pred.cv.stats <- summaryStatsCV(pred.cv, lta=TRUE, trans=0)
pred.cv.stats$Stats

##Plot observations with CV-predictions and prediction intervals
par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
plotCV(pred.cv,  1, mesa.data.model)
plotCV(pred.cv, 10, mesa.data.model)
plotCV(pred.cv, 17, mesa.data.model)
plotCV(pred.cv, 22, mesa.data.model)

##Residual QQ-plots
par(mfrow=c(1,2),mar=c(3,2,1,1),pty="s")
##Raw Residual QQ-plot
CVresiduals.qqnorm(pred.cv.stats$res)
##Normalized Residual QQ-plot
CVresiduals.qqnorm(pred.cv.stats$res.norm, norm=TRUE)
}