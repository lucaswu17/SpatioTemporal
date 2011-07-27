\name{calc.smooth.trends}
\encoding{latin1}
\Rdversion{1.1}
\alias{calc.smooth.trends}
\title{
  Smooth Basis Functions for a mesa.data Structure
}
\description{
  A front end function for calling \code{\link{SVD.smooth}} (and 
  \code{\link{SVD.smooth.cv}}), with either a \cr \code{\link{mesa.data}} 
  structure or vectors containing observations, observations times and
  locations. The function uses \code{\link{create.data.matrix}} to
  create a data matrix which is passed to \code{\link{SVD.smooth}}
  (and \code{\link{SVD.smooth.cv}}).
}
\usage{
calc.smooth.trends(mesa.data = NA, obs = mesa.data$obs$obs,
    date = mesa.data$obs$date, ID = mesa.data$obs$ID, 
    subset = NA, n.basis = 2, cv = FALSE, niter = 100, 
    conv.reldiff = 0.001, df = NULL, spar = NULL)
}
\arguments{
  \item{mesa.data}{
    A data structure containing the observations, see
    \code{\link{mesa.data}}. Use either this or the \code{obs},
    \code{date}, and \code{idx} inputs.
  }
  \item{obs}{
    A vector of observations, see \code{\link{mesa.data}}.
  }
  \item{date}{
    A vector of observation times, see \code{\link{mesa.data}}.
  }
  \item{ID}{
    A vector of observation locations, see \code{\link{mesa.data}}.
  }
  \item{subset}{
    A subset of locations to use when computing the smooth basis
    functions, see \code{\link{create.data.matrix}}.
  }
  \item{n.basis}{
    Number of basis functions to compute, see \code{\link{SVD.smooth}}.
  }
  \item{cv}{
    Also compute smooth functions using leave one out cross-validation, \cr
    see \code{\link{SVD.smooth.cv}}.
  }
  \item{niter, conv.reldiff}{
    Controls convergence for \code{\link{SVD.miss}}.
  }
  \item{df, spar}{
    The desired degrees of freedom/smoothing parameter for the spline, \cr
    see \code{\link{smooth.spline}}
  }
}
\value{
  Returns a list with the following components
  \item{svd}{A data.frame containing the smooth trends and the dates.
    This can be used as the \code{trend} in \code{mesa.data$trend}.
  }
  \item{svd.cv}{If \code{cv==TRUE} this is list of dataframes. Each
    dataframe contains the smooth trend obtained when leaving one
    site out. Similar to \cr 
    \code{SVD.smooth.cv(data)$smooth.SVD[[1]]}).

    If \code{cv==FALSE} this is set to \code{NA}.
  }
}
\references{
  M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
  Analyze Space-Time Processes in Statistical methods for
  spatio-temporal systems (B. Finkenst\"adt, L. Held, V. Isham eds.)
  77-150
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Paul D. Sampson
}
\seealso{
  Front end for calling \code{\link{SVD.smooth}} and
  \code{\link{SVD.smooth.cv}}, using \code{\link{create.data.matrix}} to
  create a data matrix from a \code{\link{mesa.data}} structure.

  See also \code{\link{SVD.miss}}, \code{\link{svd}}, and
  \code{\link{smooth.spline}}.
}
\examples{
##let's load some data
data(mesa.data)

##let's compute two smooth trend functions
trend <- calc.smooth.trends(mesa.data, n.basis = 2)
##and study the trends
par(mfrow=c(2,1), mar=c(2.5,2.5,.5,.5))
plot(trend$svd$date, trend$svd$V1, type="l", ylab="",xlab="",
     ylim=range(c(trend$svd$V1,trend$svd$V2)))
lines(trend$svd$date, trend$svd$V2, col=2)

##Let's exclude locations with fewer than 100 observations
IND <- names(which(table(mesa.data$obs$ID) >= 100))
##now we also compute the CV trends.
trend2 <- calc.smooth.trends(mesa.data, n.basis = 2, subset=IND,
                             cv = TRUE)

##Let's compare to the previous result
lines(trend2$svd$date, trend2$svd$V1, lty=2)
lines(trend2$svd$date, trend2$svd$V2, lty=2, col=2)

##we can also study the cross validated results to examine the
##possible variation in the estimated trends.
plot(trend$svd$date, trend2$svd$V1, type="n", ylab="",xlab="",
     ylim=range(c(trend2$svd$V1,trend2$svd$V2)))
for(i in 1:length(trend2$svd.cv)){
  lines(trend2$svd.cv[[i]]$date, trend2$svd.cv[[i]]$V1, col=1)
  lines(trend2$svd.cv[[i]]$date, trend2$svd.cv[[i]]$V2, col=2)
}

##Let's use our new trend functions in the mesa.data
mesa.data$trend <- trend2$svd
}
