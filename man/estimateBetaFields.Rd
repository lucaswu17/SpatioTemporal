\name{estimateBetaFields}
\alias{estimateBetaFields}
\title{Regression Estimates of beta-Fields}
\usage{
  estimateBetaFields(STdata = NULL, subset = NULL)
}
\arguments{
  \item{STdata}{A \code{STdata}/\code{STmodel} object
  containing observations. Use either this or the
  \code{obs}, \code{date}, and \code{ID} inputs.}

  \item{subset}{A subset of locations for which to estimate
  the beta-fields. A warning is given for each name not
  found in \code{ID}.}
}
\value{
  A list with two matrices; the estimated beta-coefficients
  and standard deviations of the estimates.
}
\description{
  Estimates the latent-beta fields for a
  \code{STdata}/\code{STmodel} object by regressing the
  observations for each site on the temporal trends.
}
\examples{
require(plotrix)

##load data
data(mesa.model)

##Regression based estimate of the beta-fields
beta <- estimateBetaFields(mesa.model)

##check regression coefficients
summary(beta$beta)

##or plot as a function of distance to coast,
##with uncertainties
par(mfrow=c(2,2))
for(i in 1:3){
  plotCI(mesa.model$LUR[[1]][,"log10.m.to.a1"], beta$beta[,i],
         uiw=1.96*beta$beta.sd[,i],
         ylab=colnames(beta$beta)[i])
}

##or compare to the fields from predict.STmodel
data(pred.mesa.model)

##Study the results
##Start by comparing beta fields
par(mfcol=c(1,1), mar=c(4.5,4.5,2,.5), pty="s")
plotCI(x=beta$beta[,1], y=pred.mesa.model$beta$EX[,1],
       uiw=1.96*sqrt(pred.mesa.model$beta$VX[,1]),
       main="Temporal Intercept",
       xlab="Empirical estimate",
       ylab="Spatio-Temporal Model")
plotCI(x=beta$beta[,1], y=pred.mesa.model$beta$EX[,1],
       uiw=1.96*beta$beta.sd[,1], add=TRUE, err="x")
abline(0,1,col="grey")

##or just the regression part of the beta fields
points(x=beta$beta[,1], y=pred.mesa.model$beta$mu[,1],
       col=2, pch=19)
}
\author{
  Johan Lindström
}
\seealso{
  Other data matrix: \code{\link{createDataMatrix}},
  \code{\link{mesa.data.raw}}, \code{\link{SVDmiss}},
  \code{\link{SVDsmooth}}, \code{\link{SVDsmoothCV}}

  Other STdata functions: \code{\link{c.STmodel}},
  \code{\link{createDataMatrix}},
  \code{\link{createSTdata}}, \code{\link{createSTmodel}},
  \code{\link{detrendSTdata}},
  \code{\link{removeSTcovarMean}},
  \code{\link{updateSTdataTrend}},
  \code{\link{updateTrend}},
  \code{\link{updateTrend.STdata}},
  \code{\link{updateTrend.STmodel}}

  Other STmodel functions: \code{\link{createCV}},
  \code{\link{createDataMatrix}},
  \code{\link{createSTmodel}},
  \code{\link{dropObservations}}, \code{\link{loglikeST}},
  \code{\link{loglikeSTdim}}, \code{\link{loglikeSTnaive}},
  \code{\link{predictNaive}},
  \code{\link{processLocation}}, \code{\link{processLUR}},
  \code{\link{processST}}, \code{\link{updateCovf}},
  \code{\link{updateSTdataTrend}},
  \code{\link{updateTrend}},
  \code{\link{updateTrend.STdata}},
  \code{\link{updateTrend.STmodel}}
}

