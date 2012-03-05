\name{run.MCMC}
\encoding{latin1}
\Rdversion{1.1}
\alias{run.MCMC}
\title{
  MCMC Inference of Parameters in the Spatio-Temporal Model.
}
\description{
  Estimates parameters and parameter uncertainties for the
  spatio-temporal model using a Metropolis-Hastings based Markov Chain
  Monte Carlo (MCMC) algorithm.
}
\usage{
run.MCMC(x, mesa.data.model, type = "f", N = 1000, 
         Hessian.prop = NA, Sigma.prop = NA, silent = TRUE)
}
\arguments{
  \item{x}{
    Point at which to start the MCMC. Could be either only
    \emph{log}-covariance parameters or regression and
    \emph{log}-covariance parameters. If regression parameters are given but
    not needed they are dropped, if they are needed but not given they
    are inferred by calling \cr \code{cond.expectation(x, mesa.data.model,
      only.pars=TRUE)}.
  }
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding which
    \cr geographic and spatio-temporal covariates to use when fitting
    the model. \cr See \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
  \item{type}{
    A single character indicating the type of log-likelihood to
    compute. Valid options are "f" or "r", for \emph{full},
    or \emph{restricted maximum likelihood} (REML).

    Since \emph{profile} is not a proper likelihood \code{type="p"}
    will revert (with a warning) to using the \emph{full}
    log-likelihood.

    Computation of the \emph{full} log-likelihood needs both
    regression and \emph{log}-covariance parameters while
    \emph{REML} only needs the \emph{log}-covariance parameters.
  }
  \item{N}{
    Number of MCMC iterations to run.
  }
  \item{Hessian.prop}{
    Hessian (information) matrix for the log-likelihood, can be used to
    create a proposal matrix for the MCMC.
  }
  \item{Sigma.prop}{
    Proposal matrix for the MCMC.
  }
  \item{silent}{
    If \code{FALSE} outputs status information and brief progress
    information every 100:th iteration. If \code{TRUE} no output.
  }
}
\details{
  The function runs uses a Metropolis-Hastings algorithm (Hastings,
  1970) to sample from the parameters of the spatio-temporal model,
  assuming flat priors for all the parameters (flat on the log-scale for
  the covariance parameters).

  At each iteration of the MCMC new parameters are proposed using a
  random-walk with a proposal covariance matrix. The proposal matrix is
  determined as:
  \describe{
    \item{1}{If \code{Sigma.prop} is given then this is used.}
    \item{2}{If \code{Sigma.prop=NA} then we follow Roberts et.al. (1997)
      and compute \cr
      \code{c <- 2.38*2.38/dim(Hessian.prop)[1]}\cr
      \code{Sigma.prop <- -c*solve(Hessian.prop)}.
    }
    \item{3}{If both \code{Sigma.prop=NA} and \code{Hessian.prop=NA} then
      the hessian is computed as \cr
      \code{Hessian.prop <- loglike.hessian(x, mesa.data.model, type)} \cr
      and \code{Sigma.prop} is computed according to point 2.
    }
  }
  The resulting proposal matrix is checked to ensure that it is positive
  definite before proceeding, \cr
  \code{all(eigen(Sigma.prop)$value > 1e-10)}.
}
\value{
  Returns a list containing:
  \item{par}{A \code{N} - by - (number of parameters) matrix with
    the MCMC paths of the estimated parameters.
  }
  \item{log.like}{A vector of length \code{N} with the log-likelihood
    values at each iteration.}
  \item{acceptance}{A vector of length \code{N} with the acceptance
    probability of the MCMC proposal at each iteration. Can be used to
    study the overall acceptance rate.}
  \item{chol.prop}{Choleskey factor of the proposal matrix.}
}
\references{
  W.K. Hastings. (1970) Monte Carlo sampling methods using Markov chains
  and their applications, Biometrika, 57, 97--109.

  G.O. Roberts, A. Gelman, W.R. Gilks. (1997) Weak convergence and
  optimal scaling of random walk Metropolis algorithms, Annals of
  Probability, 7, 110--120. 
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See \code{\link{fit.mesa.model}} for ML-parameter estimation, both
  functions use log-likelihood given by \code{\link{loglike}}.

  Expected names for \code{x} are given by
  \code{\link{loglike.var.names}}.
  
  For further optimization functions see \code{\link{loglike}}, 
  \code{\link{loglike.var.names}}, \cr
  \code{\link{create.data.model}}, and \code{\link{cond.expectation}}.
}
\examples{
##load data
data(mesa.data.model)
data(mesa.data.res)

##Extract parameters,
par <- mesa.data.res$par.est$res.best$par.all
##... and the Hessian.
H <- mesa.data.res$par.est$res.best$hessian.all
\dontrun{
##run the MCMC, this may take a while...
MCMC.res <- run.MCMC(par, mesa.data.model, N = 5000, 
                     Hessian.prop = H, silent = FALSE)
}
##Get the precomputed results instead.
MCMC.res <- mesa.data.res$MCMC.res

##components of the MCMC results
names(MCMC.res)

##The acceptance probability (alpha) for each step 
##in the Metropolis-Hastings algorithm.
summary(MCMC.res$acceptance)

##The MCMC-estimated parameters
summary(MCMC.res$par)

##MCMC tracks for four of the parameters
par(mfrow=c(4,1),mar=c(2,2,2.5,.5))
for(i in c(4,9,13,15)){
  plot(MCMC.res$par[,i], ylab="", xlab="", type="l",
       main=colnames(MCMC.res$par)[i])
}
}

