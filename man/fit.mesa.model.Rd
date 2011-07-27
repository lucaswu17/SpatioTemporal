\name{fit.mesa.model}
\encoding{latin1}
\Rdversion{1.1}
\alias{fit.mesa.model}

\title{
  Estimation of the Spatio-Temporal Model
}
\description{
  Estimates parameters of the spatio-temporal model using
  maximum-likelihood, profile maximum likelihood or restricted maximum 
  likelihood (REML). The function uses the \emph{L-BFGS-B} method in 
  \code{\link{optim}} to maximise \code{\link{loglike}}.
}
\usage{
fit.mesa.model(x, mesa.data.model, type = "p", h = 0.001,
    diff.type = 1, lower = -15, upper = 15, hessian.all = FALSE,
    control = list(trace = 3, maxit = 1000))
}
\arguments{
  \item{x}{
    Vector or matrix of starting point(s) for the optimisation. A vector
    will be treated as a single starting point. If \code{x} is a matrix
    the optimisation will be run using each column as a separate
    starting point. If \code{x} is a single integer then multiple
    starting points will be created as a set of constant vectors with
    the values of each starting point taken as \code{seq(-5, 5,
      length.out=par.init)}. See \code{details} below.
  }
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding which
    \cr geographic and spatio-temporal covariates to use when fitting
    the model. \cr
    See \code{\link{create.data.model}} and \code{\link{mesa.data.model}}.
  }
  \item{type}{
    A single character indicating the type of log-likelihood to
    use. Valid options are "f", "p", and "r", for \emph{full},
    \emph{profile} or \emph{restricted maximum likelihood}
    (REML).
  }
  \item{h, diff.type}{
    Step length and type of finite difference to use when computing
    gradients. \cr See \code{\link{loglike.grad}} and
    \code{\link{gen.gradient}} for details.
  }
  \item{lower, upper}{
    Bounds on the variables, see \code{\link{optim}}.
  }
  \item{hessian.all}{
    If \code{type!="f"} computes hessian (and uncertainties) for both
    regression and \emph{log}-covariance parameters, not only for
    \emph{log}-covariance parameters. See \code{value} below.
  }
  \item{control}{
    A list of control parameters for the optimisation.
    See \code{\link{optim}} for details.
  }
}
\details{
  The function estimates parameters for the spatio-temporal model using
  the full likelihood formulation, profile likelihood, or restricted
  maximum likelihood (REML). In principal full likelihood and profile
  likelihood should give the same results, corresponding to the maximum
  likelihood estimate, with the full likelihood approach being
  \emph{slower}.

  The starting point(s) for the optimisation can either contain both
  regression parameters and log-covariances parameters for a
  total of \code{loglike.dim(mesa.data.model)$nparam} parameters or only
  contain the log-covariances covariances parameters \cr
  i.e. \code{loglike.dim(mesa.data.model)$nparam.cov} parameters. \cr If
  regression parameters are given but not needed (\code{type!="f"}) they
  are dropped; if they are needed but not given they are inferred
  through a generalised least squares (GLS) computation, obtained by
  calling \code{cond.expectation(x, mesa.data.model, only.pars=TRUE)}.

  If \code{x} is a single integer then that many starting points will be
  created as vectors with the constant values for the log-covariance
  parameters being in the sequence \code{seq(-5, 5, length.out=x)}; if
  needed the corresponding regression parameters are computed using GLS.

  The estimation uses the \emph{L-BFGS-B} method in \code{\link{optim}}
  to maximise \code{\link{loglike}}. Gradient calculations are done by
  \code{\link{loglike.grad}}; when specifying the type of finite
  differences to use note that central differences (\code{diff.type=0})
  will \strong{drastically} increase estimation time, often with
  \strong{little or no} benefit to the optimisation.

  If multiple starting points are used this function returns all
  optimisation results, along with an indication of the best result. The
  best result is determined by first evaluating which of the
  optimisations have converged. Convergence is determined by checking
  that the output from \code{\link{optim}} has \code{convergence==0} and
  that the \code{hessian} is negative definite, \cr
  i.e. \code{all(eigen(hessian)$value < -1e-10)}. \cr
  Among the converged optimisations the one with the highest
  log-likelihood value is then selected as the best result.

  If none of the optimisations have converged the result with the highest
  log-likelihood value is selected as the best result.
}
\value{
  Returns a list containing:
  \item{res.best}{A list containing the best optimisation result;
    elements are described below. Selection of the best result is
    described in \code{details} above.
  }
  \item{res.all}{A list with all the optimisations results, each element
    contains (almost) the same information as \code{res.best}, e.g.
    \code{res.all[[i]]} contains optimisation results for the i:th
    starting point.
  }
  \item{message}{A text string with information regarding best value and
    number of converged optimisations.
  }

  Most of the elements in \code{res.best} (and in \code{res.all[[i]]}) are
  obtained from \code{\link{optim}}. The following is a brief
  description:
  \describe{
    \item{par}{The best set of parameters found.}
    \item{value}{Log-likelihood value corresponding to \code{par}.}
    \item{counts}{The number of calls to \code{\link{loglike}} and
      \code{\link{loglike.grad}} respectively.}
    \item{convergence}{\code{0} indicates successful convergence, see
      \code{\link{optim}} for other possible values.}
    \item{message}{Additional information returned by \code{\link{optim}}.}
    \item{hessian}{A symmetric matrix giving the finite difference
      Hessian of the log-likelihood at \code{par}.}
    \item{conv}{A logical variable indicating convergence; \code{TRUE} if
      \code{convergence==0} and \code{hessian} is negative definite, see
      \code{details} above.
    }
    \item{par.init}{The initial parameters used for this optimisation.}
    \item{par.all}{All parameters (both regression and
      \emph{log}-covariance). Identical to \code{par} if \code{type="f"},
      otherwise computed as \cr
      \code{cond.expectation(par, mesa.data.model, only.pars=TRUE)}.
    }
    \item{hessian.all}{The hessian for all parameters (both regression and
      \emph{log}-covariance). Identical to \cr \code{hessian} if
      \code{type="f"}.
      
      \strong{NOTE:} Due to computational considerations
      \code{hessian.all} is computed \emph{only} for \cr \code{res.best}.
    }
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Uses the \emph{L-BFGS-B} method in \code{\link{optim}} to maximise the
  log-likelihood given by \code{\link{loglike}} (with gradients from
  \code{\link{loglike.grad}}.

  Expected names for \code{x} are given by
  \code{\link{loglike.var.names}}.
  
  For optimization functions see \code{\link{loglike}}, 
  \code{\link{loglike.var.names}}, \code{\link{create.data.model}}, \cr
  \code{\link{run.MCMC}}, and \code{\link{cond.expectation}}.
  
  For prediction see also \code{\link{cond.expectation}}, and 
  \code{\link{plotCV}} for plotting prediction results.
}
\examples{
##load a model object
data(mesa.data.model)

##examine the model
printMesaDataNbrObs(mesa.data.model)

##covariates
mesa.data.model$LUR.list
mesa.data.model$ST.Ind

##Important dimensions of the model
dim <- loglike.dim(mesa.data.model)
print(dim)

##Set up initial parameter values for optimization
x.init <- cbind(rep(2,dim$nparam.cov), 
                c(rep(c(1,-3),dim$m+1),-3))
##and add names to the initial values
rownames(x.init) <- loglike.var.names(mesa.data.model,
                                      all=FALSE)
print(x.init)

\dontrun{
##estimate parameters
##This may take a while...
par.est <- fit.mesa.model(x.init, mesa.data.model, type="p",
      hessian.all=TRUE, control=list(trace=3,maxit=1000))
}
##Let's load precomputed results instead.
data(mesa.data.res)
par.est <- mesa.data.res$par.est

##Optimisation status message
par.est$message

##compare the estimated parameters for the two starting points
cbind(par.est$res.all[[1]]$par.all,
      par.est$res.all[[2]]$par.all)
##and values of the likelihood
cbind(par.est$res.all[[1]]$value,
      par.est$res.all[[2]]$value)

##extract the estimated parameters
x <- par.est$res.best$par.all

##and approximate uncertainties from the hessian
x.sd <- sqrt(diag(-solve(par.est$res.best$hessian.all)))
names(x.sd) <- names(x)

##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x, ylim=range(c(x-1.96*x.sd,x+1.96*x.sd)),
     xlab="", xaxt="n")
points(x-1.96*x.sd, pch=3)
points(x+1.96*x.sd, pch=3)
axis(1,1:length(x),names(x),las=2)
}
