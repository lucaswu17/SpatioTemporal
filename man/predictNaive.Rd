\name{predictNaive}
\encoding{latin1}
\Rdversion{1.1}
\alias{predictNaive}

\title{
  Naive Temporal Predictions
}
\description{
  Computes naive predictions that are based on a few sites. These
  predictions can then be used, e.g. in \code{\link{summaryStatsCV}}, to
  evaluate how much better the spatial-temporal model performs compared
  to simple (temporal) predictions.
}
\usage{
predictNaive(mesa.data.model, location = NULL, type = NULL)
}
\arguments{
  \item{mesa.data.model}{
  }
  \item{location}{
    List of locations on which to base the naive predictions.
  }
  \item{type}{
    The type of sites to  base the predictions on, uses the (optional)
    field \cr \code{mesa.data.model$location$type}.
  }
}
\details{
  The function requires one of \code{location} and \code{type} to be
  specified, if both are given \code{location} \emph{will be used over}
  \code{type}. If only \code{type} is given then all locations such that
  \cr \code{as.character(mesa.data.model$location$type) \%in\% type} \cr
  will be used.

  Given a set of locations the function computes 4 sets of naive
  prediction for all observations in \code{mesa.data.model$obs$obs}.
  The compute predictions are:
  \describe{
    \item{smooth.fixed}{
      The smooth trend in \code{mesa.data.model$trend} is fit to
      \emph{all} observations at the sites in \code{location} using a
      standard linear regression, \code{\link{lm}}.
      The fitted smooth function is then used as a naive prediction for
      all locations. The resulting predictions are \emph{smooth in time}
      and \emph{constant in space}.
    }
    \item{avg.fixed}{
      The temporal average over sites in \code{location} is computed
      and used as a naive prediction. This gives predictions that
      \emph{vary in time} but are \emph{constant in space}.
    }
    \item{smooth.closest.fixed}{
      This fits the smooth trend in \code{mesa.data.model$trend} to each
      site in \code{location} and then uses the smooth fit at the
      closest fixed site as a naive prediction at all other sites.
      The resulting predictions are \emph{smooth in time}
      and \emph{vary in space}.
    }
    \item{closest.fixed}{
      This uses the observations at the closest site in \code{location}
      to predict observations at all other sites. This gives predictions
      that \emph{vary both in time and space}.      
    }
  }
}
\value{
  \item{locations}{
    A character vector with the locations used for the naive
    predictions.
  }
  \item{pred}{
    A (number of observations) - by - (4) matrix containing the four
    naive predictions described under \code{details} above.
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See \code{\link{createCV}} and \code{\link{estimateCV}} for cross-validation
set-up and estimation.

For computing CV statistics, see also \code{\link{compute.ltaCV}}, 
and for further illustration see \code{\link{plotCV}}, 
\code{\link{CVresiduals.qqnorm}}, and \code{\link{summaryStatsCV}}.
}
\examples{
##load data
data(mesa.data.model)

##naive predictions based on either AQS
pred.aqs <- predictNaive(mesa.data.model, type="AQS")
##or FIXED sites
pred.fixed <- predictNaive(mesa.data.model, type="FIXED")

##plot the predictions - The two cases that are constant in space
par(mfcol=c(3,1), mar=c(4.5,4.5,1,.5))
##observations as a function of date
plot(mesa.data.model$obs$date, mesa.data.model$obs$obs,
     col=mesa.data.model$obs$idx, pch=19, cex=.25,
     xlab="Date", ylab="observations")
##Add the predictions based on the smooth fitted to all sites
lines(mesa.data.model$obs$date,
      pred.aqs$pred[,"smooth.fixed"], col=1)
lines(mesa.data.model$obs$date,
      pred.fixed$pred[,"smooth.fixed"], col=2)

plot(mesa.data.model$obs$date, mesa.data.model$obs$obs,
     col=mesa.data.model$obs$idx, pch=19, cex=.25,
     xlab="Date", ylab="observations")
##Add predictions based on the temporal average
lines(mesa.data.model$obs$date,
      pred.aqs$pred[,"avg.fixed"], col=1)
lines(mesa.data.model$obs$date,
      pred.fixed$pred[,"avg.fixed"], col=2)

##plot the predictions - One of the cases that vary in space
##first extract a data matrix
D <- create.data.matrix(date=mesa.data.model$obs$date,
    obs=pred.aqs$pred[,"smooth.closest.fixed"],
    ID=mesa.data.model$obs$ID)

##observations as a function of date
##(only five sites for clarity)
IND <- (mesa.data.model$obs$idx \%in\% c(1:5))
plot(mesa.data.model$obs$date[IND], mesa.data.model$obs$obs[IND],
     col=mesa.data.model$obs$idx[IND], pch=19, cex=.25,
     xlab="Date", ylab="observations")
##Add the predictions based on the smooth
##fitted to the closest site
for(i in 1:5)
  lines(as.Date(rownames(D)), D[,mesa.data.model$location$ID[i]],
        col=i)
}