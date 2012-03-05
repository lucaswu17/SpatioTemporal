\name{plotCV}
\encoding{latin1}
\Rdversion{1.1}
\alias{plotCV}
\alias{plotPrediction}

\title{
  Illustrates Predictions and Cross-validated Predictions
}
\description{
  Plots the result of (cross-validated) predictions. Given output from
  \code{\link{cond.expectation}} or \code{\link{predictCV}}
  this will plot the predicted values for one site along with
  approximate confidence intervalls, and observed data.
}
\usage{
plotCV(pred.cv, ID, mesa.data = NA, add = FALSE, p = 0.95,
       col = c("black", "red", "grey"), lty=c(1,1), pch=c(NA,NA),
       cex=c(1,1))

plotPrediction(cond.exp, ID, mesa.data = NA, add = FALSE,
               p = 0.95, col = c("black", "red", "grey"),
               lty=c(1,1), pch=c(NA,NA), cex=c(1,1))
}
\arguments{
  \item{pred.cv}{
    Result of a cross-validation. Should be the output from
    \code{\link{predictCV}}.
  }
  \item{cond.exp}{
    Result of conditional expectations (predictions). Should be the
    output from \code{\link{cond.expectation}}, with
    \code{only.obs=FALSE}.
  }
  \item{ID}{
    The location for which we want plots.
    Either a string matching the names in \code{mesa.data$location$ID}
    or an integer; if an integer the functions will plot data from 
    \code{ID=mesa.data.model$location$ID[ID]}.
  }
  \item{mesa.data}{
    A data structure containing the observations, either
    \code{\link{mesa.data}} or \cr \code{\link{mesa.data.model}}. This
    is \emph{only needed} if the location in \code{ID} is specified as
    an integer.

    If \code{mesa.data=NA} then \code{plotPrediction} will \emph{not}
    plot any observations.
  }
  \item{add}{
    Add to an existing plot, or start a new plot.
  }
  \item{p}{
    Approximate coverage of the plotted confidence bands.
  }
  \item{col}{
    A vector of three colours to use when plotting the cross-validated
    predictions. The first is the colour of the predictions, second for
    the observations and third for the polygon illustrating the
    confidence bands.
  }
  \item{lty, pch}{
    Line and/or point type to use when plotting the predictions and
    observations. The first value is for the predictions, the second for
    the observations. If any value is given as \code{NA} this implies no
    lines and/or points. Thus the default gives \emph{solid lines} for
    both predictions and observations and \emph{no points}. See
    \code{\link{points}}, \code{\link{lines}}, and \code{\link{par}}.

    If only one value is given the function will use that value for both
    predictions and observations, e.g. \code{lty=c(lty,lty)}.
  }
  \item{cex}{
    Size of the points to use if \code{pch!=NA}, see
    \code{\link{points}} and \code{\link{par}}.
  }
}
\value{
  Does not return anything.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See \code{\link{createCV}} and \code{\link{estimateCV}} for cross-validation
set-up and estimation.

For prediction, see \code{\link{fit.mesa.model}} and 
\code{\link{cond.expectation}}.

For computing CV statistics, see also \code{\link{predictNaive}} and
\code{\link{compute.ltaCV}}; for further illustration see 
\code{\link{CVresiduals.qqnorm}} and \code{\link{summaryStatsCV}}.
}
\examples{
##load data
data(mesa.data.model)
data(mesa.data.res)

##############
##  plotCV  ##
##############
##Extract pre-computed cross-validated predictions
pred.cv <- mesa.data.res$pred.cv

##Plot observations with CV-predictions and 
##95\% prediction intervals
par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
plotCV(pred.cv,  1, mesa.data.model)
##different colours
plotCV(pred.cv, 10, mesa.data.model,
    col=c("blue","magenta","light blue"))
##points and lines for the observations
plotCV(pred.cv, 17, mesa.data.model, lty=1, pch=c(NA,19), cex=.5)
##location by name
plotCV(pred.cv, "L002", mesa.data.model)

######################
##  plotPrediction  ##
######################
##Extract pre-computed conditional expectations
EX <- mesa.data.res$EX

##plot predictions and observations for 4 locations
par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
plotPrediction(EX,  1, mesa.data.model)
##different colours
plotPrediction(EX, 10, mesa.data.model,
    col=c("blue","black","light blue"))
##points for the observations
plotPrediction(EX, 17, mesa.data.model, lty=c(1,NA),
               pch=c(NA,19), cex=.5)
##location by name
plotPrediction(EX, "L002", mesa.data.model)
}
