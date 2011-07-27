\name{drop.observations}
\encoding{latin1}
\Rdversion{1.1}
\alias{drop.observations}

\title{
  Drop Observations from mesa.data.model
}
\description{
  Drops observations from \code{\link{mesa.data.model}}, removing
  marked observations along with the corresponding locations and
  recomputes a number of relevant elements. Used by the cross-validation
  functions \code{\link{estimateCV}} and \code{\link{predictCV}}.
  
}
\usage{
drop.observations(mesa.data.model, Ind.cv)
}
\arguments{
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding the
    observation locations. See \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
  \item{Ind.cv}{
    A logical vector with one element per observation in \cr
    \code{mesa.data.model$obs}. Observations marked with the \code{TRUE}
    will be \cr dropped from the data structure. Use
    \code{\link{createCV}} to create the logical vector.
  }
}
\value{
  Returns the \code{\link{mesa.data.model}} without the observations
  marked by \code{Ind.cv}. Any locations that lack observations as a
  result of this removal are also excluded from the returned structure.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See also \code{\link{createCV}}, \code{\link{estimateCV}}, and
  \code{\link{predictCV}}.
}
\examples{
##load data
data(mesa.data.model)

##Mark 30\% of observations
I <- runif(dim(mesa.data.model$obs)[1])<.3
##drop these observations
mesa.data.new <- drop.observations(mesa.data.model, I)

##This reduces the remaining number of observations
printMesaDataNbrObs(mesa.data.model)
printMesaDataNbrObs(mesa.data.new)
\dontshow{
if( (length(mesa.data.model$obs$obs)<=length(mesa.data.new$obs$obs)) ||
    (length(mesa.data.new$obs$obs)!=sum(!I)) ){
  stop("drop.observations 1: Observations not dropped")
}
}
##create cross validation structure
Icv <- createCV(mesa.data.model, groups = 10)

##drop observations from the second CV group
mesa.data.new <- drop.observations(mesa.data.model, Icv[,2])

##This reduces the remaining number of observations (and locations)
printMesaDataNbrObs(mesa.data.model)
printMesaDataNbrObs(mesa.data.new)
\dontshow{
if( (length(mesa.data.model$obs$obs)<=length(mesa.data.new$obs$obs)) ||
    (length(mesa.data.new$obs$obs)!=sum(!Icv[,2])) ||
    (dim(mesa.data.model$location)[1] <=
     dim(mesa.data.new$location)[1]) ){
  stop("drop.observations 2: Observations not dropped")
}
}%\dontshow
}%\examples
