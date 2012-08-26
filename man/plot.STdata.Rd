\name{plot.STdata}
\alias{plot.STdata}
\alias{plot.STmodel}
\title{Different Plots for \code{STdata}/\code{STmodel} object}
\usage{
  \method{plot}{STdata} (x, y = "obs", ID = x$covars$ID[1],
    type = x$covars$type, col = NULL, pch = 19, cex = 0.1,
    legend.loc = "topleft", legend.names = NULL,
    add = FALSE, ...)

  \method{plot}{STmodel} (x, y = "obs",
    ID = x$locations$ID[1], type = x$locations$type, ...)
}
\arguments{
  \item{x}{\code{STdata}/\code{STmodel} object to plot.}

  \item{y}{Type of plot, options are \code{"obs"},
  \code{"res"}, \code{"acf"}, \code{"pacf"}, \code{"loc"},
  or \code{"loc.obs"}, see details below.}

  \item{ID}{The location for which we want to plot
  observations. Either a string matching the names in
  \code{x$covars$ID} or an integer; if an integer the
  functions will plot data from \code{ID=x$covars$ID[ID]}.}

  \item{type}{Factorial of \code{length(x$covars$ID)}, used
  by \code{"loc"} and for \code{"loc.obs"} to determine how
  many groups should be plotted and colour/type coded.}

  \item{col,pch}{Colour and size of points for \code{"loc"}
  and \code{"loc.obs"}, either one constant or a vector or
  \code{length(levels(type))}.}

  \item{cex}{Size of points for \code{"loc"} and
  \code{"loc.obs"}.}

  \item{legend.loc}{The location of the legend, for
  \code{"loc"} and \code{"loc.obs"}. See
  \code{\link[graphics:legend]{legend}}.}

  \item{legend.names}{A vector of character strings to be
  used in the legend, for \code{"loc"} and for
  \code{"loc.obs"}}

  \item{add}{Add to existing plot, only relevant if
  \code{y} is \code{"obs"}, \code{"res"}, \code{"loc"}, or
  \code{"loc.obs"}.}

  \item{...}{Additional parameters passed to
  \code{\link[graphics:plot]{plot}},
  \code{\link[stats:acf]{acf}}, or
  \code{\link[stats:pacf]{pacf}},}
}
\value{
  Nothing
}
\description{
  \code{\link[graphics:plot]{plot}} method for class
  \code{STdata} or \code{STmodel}. Provides several
  different plots of the data. When calles for
  \code{STmodel}, \code{STmodel$locations} acts as
  \code{STdata$covars}.
}
\details{
  Performs a variety of different plots determined by
  \code{y}: \describe{ \item{"obs"}{Plot observations for
  location \code{ID}, along with the fitted temporal
  trend.} \item{"res"}{Plot residuals for the fitted
  temporal trend at location \code{ID}; adds the \code{y=0}
  line for reference.} \item{"acf"}{Plot autocorrelation
  function for the residuals from the fitted temporal trend
  at location \code{ID}.} \item{"pacf"}{Plot partial
  autocorrelation function for the residuals from the
  fitted temporal trend at location \code{ID}.}
  \item{"loc"}{Plot the observation location index number
  as a function of the observation date, for all
  observations. Possibly coded by the \code{type} of
  observations locations.} \item{"loc.obs"}{Plot the
  observation value as a function of the observation date,
  for all observations. Possibly coded by the \code{type}
  of observations locations.} }

  For \code{y=c("loc","loc.obs")} a legend is added if
  \code{legend.loc!=NULL}. The vector \code{legend.names}
  should have length equal to the number of unique location
  types.  The default legend names are \code{levels(type)}.
  The function uses \code{col=1:length(levels(type))} if
  \code{col=NULL}.
}
\examples{
##load data
data(mesa.data)

##default plot
plot(mesa.data)

##plot monitor locations
plot(mesa.data, "loc")

##different names/colours/etc
plot(mesa.data, "loc", main="A nice plot", col=c("green","blue"),
    legend.names=c("Sites of one type", "..and of the other"),
    legend.loc="bottomleft")

##composite time-trend
plot(mesa.data, "loc.obs", legend.loc="bottomleft", cex=.5, pch=c(3,4))

##plot tim-series for the first site,
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID=1)
##residuals from the temporal trends,
plot(mesa.data, "res", ID=1)
##afc 
plot(mesa.data, "acf", ID=1)
##... and pafc for the residuals
plot(mesa.data, "pacf", ID=1)

##Same as above but calling the 2nd site by name
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID="60370016")
plot(mesa.data, "res", ID="60370016")
plot(mesa.data, "acf", ID="60370016")
plot(mesa.data, "pacf", ID="60370016")

##same, but with no temporal trend, first replace the trend with a constant
mesa.data <- updateSTdataTrend(mesa.data, n.basis=0)

par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID="60370016")
plot(mesa.data, "res", ID="60370016")
plot(mesa.data, "acf", ID="60370016")
plot(mesa.data, "pacf", ID="60370016")
}
\author{
  Johan Lindström and Assaf P. Oron
}
\seealso{
  Other STdata methods: \code{\link{createSTdata}},
  \code{\link{print.STdata}},
  \code{\link{print.summary.STdata}},
  \code{\link{summary.STdata}}

  Other STmodel methods: \code{\link{c.STmodel}},
  \code{\link{createSTmodel}},
  \code{\link{estimate.STmodel}},
  \code{\link{estimateCV.STmodel}},
  \code{\link{MCMC.STmodel}},
  \code{\link{predict.STmodel}},
  \code{\link{predictCV.STmodel}},
  \code{\link{print.STmodel}},
  \code{\link{print.summary.STmodel}},
  \code{\link{simulate.STmodel}},
  \code{\link{summary.STmodel}}
}

