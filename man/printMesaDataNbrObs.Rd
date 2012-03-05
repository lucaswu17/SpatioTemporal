\name{printMesaDataNbrObs}
\encoding{latin1}
\Rdversion{1.1}
\alias{printMesaDataNbrObs}

\title{
  Summary of Locations and Time Points
}
\description{
  Prints output summarizing locations and time points at both 
  observed and unobserved locations and dates in the input data structure.
}
\usage{
printMesaDataNbrObs(mesa.data)
}

\arguments{
  \item{mesa.data}{
    A data structure containing observations.  Could also contain
    locations and dates where data were not observed, i.e. if
    predictions were desired at those locations and dates.  Could be in
    the form of either \code{\link{mesa.data}} or \cr
    \code{\link{mesa.data.model}}.
  }
}
\details{
  The printed output gives the number of total locations 
  in \code{mesa.data} and specifies how many of those locations are for 
  observed data.  Similarly, the output gives the number of time points
  in \code{mesa.data} and the range of those dates, and specifies how many of 
  those time points are for observed data and the range of the observed time 
  points.  The number of total observations are given (over all locations
  and time points), as well as a summary of the range of dates and number of 
  observations for each monitoring type, if type is 
  specified in \code{mesa.data}.  If \emph{type is specified}, the
  output also summarizes the number of total monitoring types in
  \code{mesa.data} and the subset of which were for observed data.  
}
\value{
  Does not return anything.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  \code{\link{plotMesaData}}, \code{\link{plotMonitoringLoc}} for further
  data summaries.
}
\examples{
##load some data
data(mesa.data)
##print summary information regarding no. obs, locations, and dates
printMesaDataNbrObs(mesa.data)
}