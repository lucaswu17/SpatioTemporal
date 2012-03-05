\name{CVresiduals.qqnorm}
\encoding{latin1}
\Rdversion{1.1}
\alias{CVresiduals.qqnorm}
\alias{CVresiduals.scatter}

\title{
  QQ-norm and Scatterplots for Data and Residual Analysis
}
\description{
  Does QQ-norm and scatterplots of observations or residuals from a
  cross-validation. The scatterplots are of observations against any
  reasonable covariate; examples of how todo scatter plots of
  observations against temporal smooths as well as geographic
  covariates are given.
}
\usage{
CVresiduals.qqnorm(res, I.season =
    as.factor(rep("obs", length(res))), I.type = NA,
    main = "All Data", norm = FALSE, legend = TRUE, ...) 

CVresiduals.scatter(res, LUR, I.season =
    as.factor(rep("obs",length(res))), I.type = NA,
    xlab = "", main = "All Data", legend = TRUE,
    df = 10, ...)
}
\arguments{
  \item{res}{
    Vector of observations or residuals to plot. Could be one of \cr
    \code{mesa.data$obs$obs}, \code{mesa.data.model$obs$obs}, \cr
    \code{predCV$pred.obs[,"pred"]}, \code{statsCV$res}, or \cr
    \code{statsCV$res.norm}. \cr
    Here \code{predCV} and \code{statsCV} are outputs from
    \code{\link{predictCV}} and \cr
    \code{\link{summaryStatsCV}}, respectively.
  }
  \item{LUR}{
    A vector of land-use-regressors or other variable of interest
    against which to plot the residuals. Can be obtained from
    \code{\link{mesa.data.model}} as, e.g. \cr
    \code{mesa.data.model$X[[1]][mesa.data.model$obs$idx, 2]}. To plot
    the residuals against the temporal trends take \cr
    \code{LUR=mesa.data.model$F[,2]}.
  }
  \item{I.season}{
    A vector of factors with the same length as \code{res}. The factors
    are used to colour code the points in the plots. The default, with
    only one factor, gives the same colour to all points.
  }
  \item{I.type}{
    A vector of characters or factors with the same length as
    \code{res}. The factors are used to produce one plot with all the
    data and then subsequent plots including only data of one type.

    Each of these plots are still colour coded according to
    \code{I.season}. The default, \code{NA}, assumes no
    breakdown. Elements that are \code{NA} in \code{I.type} will be
    included in the initial all data plot, but not in any of the
    subsequent plots.
  }
  \item{xlab}{
    Name on the x-axis, typically the name of the \code{LUR} component
    against which the data is being plotted.
  }
  \item{main}{
    Title of the plots.
  }
  \item{norm}{
    If \code{TRUE} assumes standardised data/residuals (i.e. Gaussian
    with mean=0 and var=1) and adds \code{abline(0,1)} to the
    \code{\link{qqnorm}}-plot, in addition to \code{\link{qqline}}.
  }
  \item{legend}{
    If \code{TRUE} add a legend describing the colour coding due to
    \code{I.season} in the lower right hand corner.
  }
  \item{df}{
    Desired equivalent number of degrees of freedom for the fitted
    smooth spline, see \code{\link{smooth.spline}}.
  }
  \item{...}{
    Arguments passed on to the plotting function (\code{\link{qqnorm}}
    or \code{\link{plot}}).
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

For computing CV statistics, see also \code{\link{predictNaive}} and
\code{\link{compute.ltaCV}}; for further illustration see \code{\link{plotCV}}
and \code{\link{summaryStatsCV}}.
}
\examples{
##load data
data(mesa.data.model)

##create a vector dividing data into four seasons
I.season <- matrix(NA,length(mesa.data.model$obs$date),1)
I.season[months(mesa.data.model$obs$date) \%in\%
         c("December","January","February"),1] <- "DJF"
I.season[months(mesa.data.model$obs$date) \%in\%
         c("March","April","May"),1] <- "MAM"
I.season[months(mesa.data.model$obs$date) \%in\%
         c("June","July","August"),1] <- "JJA"
I.season[months(mesa.data.model$obs$date) \%in\%
         c("September","October","November"),1] <- "SON"
I.season <- factor(I.season,levels=c("DJF","MAM","JJA","SON"))

##create a vector dividing data into type (AQS or FIXED)
I.type <- as.character(mesa.data.model$location$type[
                       mesa.data.model$obs$idx])

##qq-plot of the observations
par(mfrow=c(2,2),mar=c(2,2,2,.5),mgp=c(1.7,.6,0),pty="s")
CVresiduals.qqnorm(mesa.data.model$obs$obs, I.season, I.type)

##scatter plot of the observations against some of the LUR-covariates
LUR <- mesa.data.model$X[[1]][mesa.data.model$obs$idx,]

par(mfrow=c(2,3),mar=c(3,2.5,2,.5),mgp=c(1.7,.6,0),pty="s")
##i=1 is the constant LUR, not very interesting
for(i in 2:4)
  CVresiduals.scatter(mesa.data.model$obs$obs, LUR[,i], I.season,
                      xlab = colnames(LUR)[i])
##and against the temporal trends
for(i in 2:3)
  CVresiduals.scatter(mesa.data.model$obs$obs, mesa.data.model$F[,i],
                      I.season, xlab = colnames(mesa.data.model$F)[i])

}