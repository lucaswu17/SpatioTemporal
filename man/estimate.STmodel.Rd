\name{estimate.STmodel}
\alias{estimate.STmodel}
\title{Estimation of the Spatio-Temporal Model}
\usage{
  \method{estimate}{STmodel} (object, x, x.fixed = NULL,
    type = "p", h = 0.001, diff.type = 1,
    hessian.all = FALSE, lower = -15, upper = 15,
    method = "L-BFGS-B",
    control = list(trace = 3, maxit = 1000), ...)
}
\arguments{
  \item{object}{\code{STmodel} object for which to estimate
  parameters.}

  \item{x}{Vector or matrix of starting point(s) for the
  optimisation. A vector will be treated as a single
  starting point. If \code{x} is a matrix the optimisation
  will be run using each column as a separate starting
  point.  If \code{x} is a single integer then multiple
  starting points will be created as a set of constant
  vectors with the values of each starting point taken as
  \code{seq(-5, 5, length.out=x)}. See \code{details}
  below.}

  \item{x.fixed}{Vector with parameter to be held fixed;
  parameters marked as \code{NA} will still be estimated.}

  \item{type}{A single character indicating the type of
  log-likelihood to use.  Valid options are "f", "p", and
  "r", for \emph{full}, \emph{profile} or \emph{restricted
  maximum likelihood} (REML).}

  \item{h,diff.type}{Step length and type of finite
  difference to use when computing gradients, see
  \code{\link{loglikeSTGrad}}.}

  \item{hessian.all}{If \code{type!="f"} computes hessian
  (and uncertainties) for both regression and
  \emph{log}-covariance parameters, not only for
  \emph{log}-covariance parameters. See \code{value}
  below.}

  \item{lower,upper,method}{Parameter bound and
  optimisation method, passed to
  \code{\link[stats:optim]{optim}}.}

  \item{control}{A list of control parameters for the
  optimisation.  See \code{\link[stats:optim]{optim}} for
  details; setting \code{trace}=0 eliminates all ouput.}

  \item{...}{Ignored additional arguments.}
}
\value{
  \code{estimateSTmodel} object containing:
  \item{res.best}{A list containing the best optimisation
  result; elements are described below. Selection of the
  best result is described in \code{details} above.}
  \item{res.all}{A list with all the optimisations results,
  each element contains (almost) the same information as
  \code{res.best}, e.g.  \code{res.all[[i]]} contains
  optimisation results for the i:th starting point.}
  \item{summary}{A list with parameter estimates and
  convergence information for all starting points.}
}
\description{
  Estimates parameters of the spatio-temporal model using
  maximum-likelihood, profile maximum likelihood or
  restricted maximum likelihood (REML). The function uses
  the \emph{L-BFGS-B} method in
  \code{\link[stats:optim]{optim}} to maximise
  \code{\link{loglikeST}}.
}
\details{
  The starting point(s) for the optimisation can either
  contain both regression parameters and log-covariances
  parameters for a total of
  \code{loglikeSTdim(object)$nparam} parameters or only
  contain the log-covariances covariances parameters \cr
  i.e. \code{loglikeSTdim(object)$nparam.cov} parameters.
  \cr If regression parameters are given but not needed
  (\code{type!="f"}) they are dropped; if they are needed
  but not given they are inferred through a generalised
  least squares (GLS) computation, obtained by calling
  \code{\link{predict.STmodel}}.

  If multiple starting points are used this function
  returns all optimisation results, along with an
  indication of the best result. The best result is
  determined by first evaluating which of the optimisations
  have converged. Convergence is determined by checking
  that the output from \code{\link[stats:optim]{optim}} has
  \code{convergence==0} and that the \code{hessian} is
  negative definite, \cr i.e.
  \code{all(eigen(hessian)$value < -1e-10)}. \cr Among the
  converged optimisations the one with the highest
  log-likelihood value is then selected as the best result.

  If none of the optimisations have converged the result
  with the highest log-likelihood value is selected as the
  best result.

  Most of the elements in \code{res.best} (and in
  \code{res.all[[i]]}) are obtained from
  \code{\link[stats:optim]{optim}}. The following is a
  brief description: \describe{ \item{par}{The best set of
  parameters found.} \item{value}{Log-likelihood value
  corresponding to \code{par}.} \item{counts}{The number of
  function/gradient calls.} \item{convergence}{\code{0}
  indicates successful convergence, see
  \code{\link[stats:optim]{optim}}.}
  \item{message}{Additional information returned by
  \code{\link[stats:optim]{optim}}.} \item{hessian}{A
  symmetric matrix giving the finite difference Hessian of
  the function \code{par}.} \item{conv}{A logical variable
  indicating convergence; \code{TRUE} if
  \code{convergence==0} and \code{hessian} is negative
  definite, see \code{details} above.} \item{par.init}{The
  initial parameters used for this optimisation.}
  \item{par.all}{All parameters (both regression and
  \emph{log}-covariance).  Identical to \code{par} if
  \code{type="f"}.} \item{hessian.all}{The hessian for all
  parameters (both regression and
  \emph{log}-covariance).\cr \strong{NOTE:} Due to
  computational considerations \code{hessian.all} is
  computed \emph{only} for \cr \code{res.best}.} }
}
\examples{
##load a model object
data(mesa.model)

##examine the model
print(mesa.model)

##Important dimensions of the model
dim <- loglikeSTdim(mesa.model)
print(dim)

##Set up initial parameter values for optimization
x.init <- cbind(rep(2,dim$nparam.cov), c(rep(c(1,-3),dim$m+1),-3))
##and add names to the initial values
rownames(x.init) <- loglikeSTnames(mesa.model, all=FALSE)
print(x.init)
\dontrun{
  ##estimate parameters
  ##This may take a while...
  est.mesa.model <- estimate(mesa.model, x.init, type="p", hessian.all=TRUE)
}
##Let's load precomputed results instead.
data(est.mesa.model)

##Optimisation status message
print(est.mesa.model)

##compare the estimated parameters for the two starting points
est.mesa.model$summary$par.all
##and values of the likelihood (and convergence info)
est.mesa.model$summary$status

##extract the estimated parameters and approximate uncertainties
x <- coef(est.mesa.model)

##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x$par, ylim=range(c( x$par-1.96*x$sd, x$par+1.96*x$sd )),
     xlab="", xaxt="n")
points(x$par - 1.96*x$sd, pch=3)
points(x$par + 1.96*x$sd, pch=3)
axis(1, 1:length(x$par), rownames(x), las=2)

\dontrun{
  ##fixed pars
  x.fixed <- coef(est.mesa.model)$par
  x.fixed[c(1,2,5:9)] <- NA
  est.fix <- estimate(mesa.model, x.init, x.fixed, type="p")
}
}
\author{
  Johan Lindström
}
\seealso{
  Other estimateSTmodel methods:
  \code{\link{coef.estimateSTmodel}},
  \code{\link{print.estimateSTmodel}}

  Other STmodel methods: \code{\link{createSTmodel}},
  \code{\link{c.STmodel}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{MCMC.STmodel}}, \code{\link{plot.STdata}},
  \code{\link{plot.STmodel}},
  \code{\link{predictCV.STmodel}},
  \code{\link{predict.STmodel}},
  \code{\link{print.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}},
  \code{\link{summary.STmodel}}
}

