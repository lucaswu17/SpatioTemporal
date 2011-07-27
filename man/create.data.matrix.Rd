\name{create.data.matrix}
\encoding{latin1}
\Rdversion{1.1}
\alias{create.data.matrix}
\title{
  Create a Data Matrix from a mesa.data Structure
}
\description{
  Creates a data matrix from a \code{\link{mesa.data}} structure, or
  from a set of vectors. Missing data are marked as \code{NA}.
}
\usage{
create.data.matrix(mesa.data = NA, obs = mesa.data$obs$obs,
    date = mesa.data$obs$date, ID = mesa.data$obs$ID,
    subset = NA) 
}
\arguments{
  \item{mesa.data}{
    A data structure containing the observations, see
    \code{\link{mesa.data}}. Use either this or the \code{obs},
    \code{date}, and \code{ID} inputs.
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
    A subset of locations to extract the data matrix for.  This would be a 
    vector of strings that matches names in \code{ID}. A warning will be given 
    for each name not found in \code{ID}.
  }
}
\value{
  Returns a matrix with dimensions (number of timepoints)-by-(number of locations.)
  Row and column names of the matrix are taken as \code{ID} and
  \code{sort(unique(date))} respectively.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Converts a \code{\link{mesa.data}} structure to a matrix, used
  internally by \code{\link{calc.smooth.trends}}. The data matrix can be
  used as input to \link{SVD.smooth}, \link{SVD.smooth.cv}, and
  \link{SVD.miss}.
  
  See also \code{\link{create.data.model}}.

For other data functions, see \code{\link{mesa.data.raw}},
\code{\link{mesa.data.model}}, \code{\link{default.LUR.list}}, \cr
and \code{\link{construct.LUR.basis}}.
}
\examples{
##load the data
data(mesa.data)

##create a data matrix
M1 <- create.data.matrix(mesa.data)
dim(M1)
head(M1)

##create data matrix for only a few locations
M2 <- create.data.matrix(mesa.data, subset =
     c("60370002","60370016","60370113","60371002",
       "60371103","60371201","L001","L002"))
dim(M2)
head(M2)
\dontshow{
if( (dim(M1)[1]!=dim(mesa.data$trend)[1]) ||
    (dim(M1)[2]!=dim(mesa.data$location)[1]) ){
  stop("create.data.matrix: dimension missmatch - M1")
}
if( (dim(M2)[1]!=dim(mesa.data$trend)[1]) || (dim(M2)[2]!=8) ){
  stop("create.data.matrix: dimension missmatch - M2")
}
if( max(abs(M1[,colnames(M2)]-M2),na.rm=TRUE) > 1e-10 ){
  stop("create.data.matrix: M1!=M2")
}
}
}
