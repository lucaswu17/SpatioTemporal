\name{simulateMesaData}
\encoding{latin1}
\Rdversion{1.1}
\alias{simulateMesaData}
\title{
  Simulate Data from the Spatio-Temporal Model
}
\description{
  Simulates data from the spatio-temporal model.
}
\usage{
simulateMesaData(x, mesa.data.model, mesa.data = NA, rep = 1,
                 combine.data = FALSE)
}
\arguments{
  \item{x}{
    Parameters for which to simulate data. If only \emph{log}-covariance
    parameters are given then the regression parameters are inferred by
    calling \cr
    \code{cond.expectation(x, mesa.data.model, only.pars=TRUE)}.
  }
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding where
    to simulate observations. See \code{\link{mesa.data.model}}.
  }
  \item{mesa.data}{
    Data structure holding observations, observation locations, smooth
    temporal trends, geographic covariates and spatio-temporal
    covariates. This can be used to obtain simulated values at different
    times and or locations than specified by
    \code{\link{mesa.data.model}}. See \code{\link{mesa.data}}. 
  }
  \item{rep}{
    Number of field replicates to simulate.
  }
  \item{combine.data}{
    If \code{TRUE} will combine locations and covariates from
    \code{\link{mesa.data.model}} and \code{\link{mesa.data}}, but will
    use temporal trends from \code{\link{mesa.data.model}}.
  }
}
\details{
  Data is simulated for the space-time locations in
  \code{mesa.data.model} and/or \code{mesa.data}, given the parameters
  in \code{x}. If only \emph{log}-covariance parameters are given then
  the regression parameters are inferred using generalised least squares
  (GLS), as \cr
  \code{cond.expectation(x, mesa.data.model, only.pars=TRUE, type="p")}.

  The space-time locations at which to simulate data are specified
  through \code{mesa.data.model} and/or \code{mesa.data}. If
  \code{mesa.data=NA} data will be simulated at all space-time
  locations given in \code{mesa.data.model}. If \code{mesa.data!=NA}
  two options exist:
  \enumerate{
    \item If \code{combine.data=FALSE} data will be simulated at
    \emph{only} the space-time locations in \cr \code{mesa.data}. If
    \code{mesa.data} does not have a smooth temporal trend, the trend in
    \cr \code{mesa.data.model} will be used (with a warning).

    \item If \code{combine.data=TRUE} locations and covariates from both
    \code{mesa.data.model} and \code{mesa.data} are combined using
    \code{\link{combineMesaData}}. The smooth temporal trends are taken
    from \code{mesa.data.model} with trends in \code{mesa.data} being
    \emph{ignored}. I.e. data is simulated at the combined set of
    locations in \code{mesa.data.model} and \code{mesa.data}; at
    \emph{only} the times specified in \code{mesa.data.model$trend}.
  }
  If the additional object \code{\link{mesa.data}} is given,
  \code{\link{mesa.data.model}} is used to specify the model (through
  the fields \code{LUR.list} and \code{ST.Ind}).
}
\value{
  Returns a list containing:
  \item{param}{Parameters used to simulate the data.}
  \item{B}{The simulated beta fields in a (number of locations) - by -
    (number of temporal trends) - by - (number of replicates) array.
  }
  \item{X}{The simulated spatio-temporal fields in a (number of
    timepoints) - by - (number of locations) - by - (number of
    replicates) array. Row and column names indicate the time and
    locations for each observation.
  }
  \item{obs}{A list with one element per replicate, containing the
    simulated observations extracted at space-time locations matching
    those in \code{mesa.data.model$obs} (or a sorted version of
    \code{mesa.data.model$obs}). Each element of the list is a 
    data.frame with \code{obs}, \code{date}, \code{idx}, and \code{ID}
    fields. To replace the observations with the i:th simulated values
    do:\cr \code{mesa.data.model$obs <- obs[[i]]}.
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  For estimation functions see \code{\link{fit.mesa.model}},
  \code{\link{run.MCMC}},  and \code{\link{estimateCV}}.


  For prediction see \code{\link{cond.expectation}} and
  \code{\link{predictCV}}
}
\examples{
##load the data
data(mesa.data)
data(mesa.data.model)
data(mesa.data.res)

##Get estimated parameters
x <- mesa.data.res$par.est$res.best$par.all

##Simulate 5 replicates from these parameters
sim.data <- simulateMesaData(x, mesa.data.model, rep = 5)

##plot the simulated observations as a function of time
par(mfrow=c(2,2),mar=c(4,4,.5,.5))
plot(sim.data$obs[[1]]$date, sim.data$obs[[1]]$obs,
     type="n", ylab="obs", xlab="Date")
for(i in 1:5)
  points(sim.data$obs[[i]]$date, sim.data$obs[[i]]$obs, col=i)

##and the latent beta-fields
for(i in 1:3){
  plot(sim.data$B[,i,1], ylim=range(sim.data$B[,i,]), type="n",
       xlab="loc", ylab=paste("beta",colnames(sim.data$B)[i]))
  for(j in 1:5)
    points(sim.data$B[,i,j], col=j)
}

###########################################
## A case with some unobserved locations ##
###########################################
##keep only observations from the AQS sites
ID.AQS <- mesa.data$location$ID[mesa.data$location$type=="AQS"]
mesa.data$obs <- mesa.data$obs[mesa.data$obs$ID \%in\% ID.AQS,]

##create a model object, dropping unobserved locations.
mesa.data.model <- create.data.model(mesa.data,
    LUR = mesa.data.model$LUR.list, ST.Ind = mesa.data.model$ST.Ind)

##simulate some replicates for this object
sim.data.obs <- simulateMesaData(x, mesa.data.model)
sim.data.all <- simulateMesaData(x, mesa.data.model,
    mesa.data=mesa.data)
sim.data.comb <- simulateMesaData(x, mesa.data.model,
    mesa.data=mesa.data, combine.data = TRUE)

##The first simulated object now only contains the AQS sites
mesa.data.model$location$ID
colnames(sim.data.obs$X)

##while the other two contain all locations
colnames(sim.data.all$X)
colnames(sim.data.comb$X)
\dontshow{
if( any(mesa.data.model$location$ID != colnames(sim.data.obs$X)) ){
  stop("simulateMesaData 1: Miss-match in simulated locations.")
}
if( any(colnames(sim.data.all$X) != colnames(sim.data.comb$X)) ){
  stop("simulateMesaData 2: Miss-match in simulated locations.")
}
}%\dontshow
}%\examples