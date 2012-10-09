\name{createDataMatrix}
\alias{createDataMatrix}
\title{Create a Data Matrix}
\usage{
  createDataMatrix(STdata = NULL, obs = STdata$obs$obs,
    date = STdata$obs$date, ID = STdata$obs$ID,
    subset = NULL)
}
\arguments{
  \item{STdata}{A \code{STdata}/\code{STmodel} object
  containing observations, see \code{\link{mesa.data}}. Use
  either this or the \code{obs}, \code{date}, and \code{ID}
  inputs.}

  \item{obs}{A vector of observations.}

  \item{date}{A vector of observation times.}

  \item{ID}{A vector of observation locations.}

  \item{subset}{A subset of locations to extract the data
  matrix for. A warning is given for each name not found in
  \code{ID}.}
}
\value{
  Returns a matrix with dimensions (number of
  timepoints)-by-(number of locations). Row and column
  names of the matrix are taken as \code{ID} and
  \code{sort(unique(date))} respectively.
}
\description{
  Creates a data matrix from a \code{STdata}/\code{STmodel}
  object, e.g. \code{\link{mesa.data}} or
  \code{\link{mesa.model}}. Missing observations are marked
  as \code{NA}.
}
\examples{
##load the data
data(mesa.data)

##create a data matrix
M1 <- createDataMatrix(mesa.data)
dim(M1)
head(M1)

##create data matrix for only a few locations
M2 <- createDataMatrix(mesa.data, subset =
                         c("60370002","60370016","60370113","60371002",
                           "60371103","60371201","L001","L002"))
dim(M2)
head(M2)
\dontshow{
  if( (dim(M1)[1]!=dim(mesa.data$trend)[1]) ||
     (dim(M1)[2]!=dim(mesa.data$covars)[1]) ){
    stop("createDataMatrix: dimension missmatch - M1")
  }
  if( (dim(M2)[1]!=dim(mesa.data$trend)[1]) || (dim(M2)[2]!=8) ){
    stop("createDataMatrix: dimension missmatch - M2")
  }
  if( max(abs(M1[,colnames(M2)]-M2),na.rm=TRUE) > 1e-13 ){
    stop("createDataMatrix: M1!=M2")
  }
}

}
\author{
  Johan Lindström
}
\seealso{
  Other data matrix: \code{\link{mesa.data.raw}},
  \code{\link{SVDmiss}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}}

  Other STdata functions: \code{\link{createSTdata}},
  \code{\link{createSTmodel}}, \code{\link{c.STmodel}},
  \code{\link{detrendSTdata}},
  \code{\link{removeSTcovarMean}},
  \code{\link{updateSTdataTrend}}

  Other STmodel functions: \code{\link{createCV}},
  \code{\link{createSTmodel}},
  \code{\link{dropObservations}}, \code{\link{loglikeST}},
  \code{\link{loglikeSTdim}}, \code{\link{loglikeSTnaive}},
  \code{\link{predictNaive}},
  \code{\link{processLocation}}, \code{\link{processLUR}},
  \code{\link{processST}}, \code{\link{updateCovf}}
}

