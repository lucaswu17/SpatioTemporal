\name{summaryStatsCV}
\encoding{latin1}
\Rdversion{1.1}
\alias{summaryStatsCV}

\title{
  Computes Summary Statistics for Cross-validation
}
\description{
Computes summary statistics for cross validation. Statistics that are
computed include RMSE, R2, and coverage of CI:s; both for all
observations and stratified by date.
}
\usage{
summaryStatsCV(predCV, pred.naive = NULL, lta = FALSE, 
               by.date = FALSE, p = 0.95, trans = NULL)
}
\arguments{
  \item{predCV}{
    Result of a cross-validation. Should be the output from
    \code{\link{predictCV}}.
  }
  \item{pred.naive}{
    Result of naive prediction, this is used to compute modified
    R2 values. Should be the output from \code{\link{predictNaive}}.
  }
  \item{lta}{
    Compute cross-validation statistics for the long term averages at
    each site. If \code{trans!=NULL} the transformation will be applied
    \emph{before} computation of the averages, see
    \code{\link{compute.ltaCV}} for details.
  }
  \item{by.date}{
    Compute individual cross-validation statistics for each time-point.
    May lead to the computation of \emph{very many} statistics.
  }
  \item{p}{
    Approximate coverage of the computed confidence bands. The
    confidence bands are used when computing the coverage of the
    cross-validated predictions.
  }
  \item{trans}{
    Transform observations and predictions \emph{before} computing
    statistics. Different values for \code{trans} give different
    transforms:
    \describe{
      \item{\code{NULL}}{Gives no transformation}
      \item{\code{0}}{Takes the exponent of data, essentially
	assuming that data was originally log-transformed.}
      \item{non-zero values}{
	Raises the data to that power, assuming an original
	\cr root-transform, e.g. \code{trans=2} gives the long term
	averages as \cr
	\code{mean(obs^2)} and \code{mean(pred^2)}.
      }
    }
  }
}
\value{
  Returns a list containing:
  \item{Stats}{
    A data.frame where the columns contain RMSE, R2 and coverage of the
    width \code{p} confidence intervall(s). At a minimum this is
    computed for all observations.

    If \code{pred.naive!=NULL} four additional rows are added to
    \code{Stats}. These rows contain adjusted R2 that compare
    cross-validated predictions to predictions computed using
    \code{predictNaive}. The adjusted R2 are computed as
    (1 - MSE_cv/MSE_naive). For this to make sense the locations used
    for the naive predictions \emph{should not} be among the locations
    that cross-validated predictions are computed for.

    If \code{lta=TRUE} one additional rows containing RMSE and R2 for
    the long term average predictions given by
    \code{\link{compute.ltaCV}} is added to \code{Stats}.

    If \code{by.date=TRUE} one additional rows containing RMSE, R2 and
    coverage is added to \code{Stats} for \emph{each unique observation
      date}.
  }
  \item{res, res.norm}{
    Residuals and normalised residuals from the cross-validated
    predictions. Two (nbr of observations) - by (1) vectors with
    residuals for the observations in \code{mesa.data.model$obs}.

    The residuals are computed as: \cr
    \code{res <- (predCV$pred.obs[,"obs"] -}\cr
    \code{predCV$pred.obs[,"pred"])} \cr
    \code{res.norm <- res / sqrt(predCV$pred.obs[,"pred.var"])} \cr
    Here the normalised residuals are divided by the prediction standard
    deviation.
  }
  \item{lta}{A data.frame with predicted and observed long term averages
    at each site, or \code{NULL} if \code{lta=FALSE}. If given this is
    the output from:\cr 
    \code{compute.ltaCV(predCV, trans)} \cr
    See \code{\link{compute.ltaCV}}
  }
  \item{p}{Approximate coverage of the computed confidence bands, same
    as \code{p} in the input.
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See \code{\link{createCV}} and \code{\link{estimateCV}} for cross-validation
set-up and estimation.

For computing CV statistics, see also \code{\link{predictNaive}} and
\code{\link{compute.ltaCV}}; for further illustration see \code{\link{plotCV}} 
and \code{\link{CVresiduals.qqnorm}}.
}
\examples{
##load data
data(mesa.data.model)
data(mesa.data.res)

##Extract pre-computed cross-validated predictions
pred.cv <- mesa.data.res$pred.cv

##Naive predictions based on AQS sites only
pred.N <- predictNaive(mesa.data.model, type="AQS")

##compute summary statistics
stat.CV <- summaryStatsCV(pred.cv, pred.naive=pred.N,
                          lta=TRUE, by.date=TRUE)

##study the summary statistics (for observations and long term average)
stat.CV$Stats[1:2,]

##adjusted R2 values, these are slightly strange since we
##(in this case) are basing the naive predictions on 
##things left out of the cross-validation.
stat.CV$Stats[(dim(stat.CV$Stats)[1]-3):dim(stat.CV$Stats)[1],]

##plot the RMSE for each date as a function of date
plot(as.Date(rownames(stat.CV$Stats[3:(dim(stat.CV$Stats)[1]-4),])),
     stat.CV$Stats[3:(dim(stat.CV$Stats)[1]-4),"RMSE"],
     xlab="Date",ylab="RMSE")
##add over all RMSE as reference
abline(h=stat.CV$Stats["obs","RMSE"])

##Some plots for the residuals
par(mfrow=c(2,2), mar=c(4.5,4.5,3,.5))
## residuals against observations
plot(mesa.data.model$obs$obs, stat.CV$res,
     ylab="Residuals", xlab="Observations")
## Norm-plot for the residuals
CVresiduals.qqnorm(stat.CV$res)
## Norm-plot and normalised residuals, these should be N(0,1).
CVresiduals.qqnorm(stat.CV$res.norm, norm=TRUE)
## normalised residuals against the first temporal trend
CVresiduals.scatter(stat.CV$res.norm, mesa.data.model$F[,2],
                    xlab="First temporal trend")
}
