\name{create.data.model}
\encoding{latin1}
\Rdversion{1.1}
\alias{create.data.model}
\title{
  Creates the mesa.data.model structure; selection of model covariates
}
\description{
  Creates the auxiliary data structure \code{\link{mesa.data.model}};
  the structure holds information regarding which geographic and
  spatio-temporal covariates to use for model fitting.
}
\usage{
create.data.model(mesa.data, LUR = NA, ST.Ind = NA,
    strip = TRUE, strip.loc = strip, strip.time = strip)
}
\arguments{
  \item{mesa.data}{
    Data structure holding observations, observation locations, smooth
    temporal trends, geographic covariates and spatio-temporal
    covariates. See \code{\link{mesa.data}}.
  }
  \item{LUR}{
    A vector or list of vectors indicating which geographic covariates
    to use. \code{LUR = NA} uses all available covariates are used for
    all the temporal trends; a vector of integers or characters (names)
    specifies covariates to use for all temporal trends; and a list of
    vectors is used this allows for different covariates for each
    temporal trend.

    See further \code{\link{default.LUR.list}}.
  }
  \item{ST.Ind}{
    A vector indicating which spatio-temporal covariates to use.
    \code{ST.Ind = NA} uses all spatio-temporal covariates, \code{ST.Ind
      = NULL} uses \emph{no} spatio-temporal covariates.
    
    See further \code{\link{default.ST.list}}.
  }
  \item{strip, strip.loc, strip.time}{Drop times and/or locations
    without observations, \code{strip} is a shortcut for setting both
    \code{strip.loc} and \code{strip.time}.

    If \code{strip.loc=TRUE} all locations that do not have any
    observations are \cr dropped. If \code{strip.time=TRUE} all dates
    that do not have any observations are dropped from
    \code{mesa.data$trend}.
  }
}
\details{
  The function creates the auxiliary data structure
  \code{\link{mesa.data.model}} which is used in model fitting and
  prediction, see \code{\link{fit.mesa.model}},
  \code{\link{cond.expectation}}, and \code{\link{run.MCMC}}.

  The resulting structure holds information regarding which geographic
  and spatio-temporal covariates to use for model fitting, as well as a
  number of precomputed fields that speed up log-likelihood evaluations.

  For any observations that occur at times (dates) \emph{not in}
  \code{mesa.data$trend$date} the smooth temporal trends are
  interpolated to these times using \code{\link{spline}}.

  When selecting geographic covariates the code allows for different
  covariates for different temporal trends. \code{LUR = NA} gives all
  covariates are used for all the temporal trends. If a vector of
  integers or characters (names) is used to specify covariates then
  these covariates will be used for all the temporal trends
  (e.g. \code{LUR = c(1,2,3)}). If \code{LUR} instead is given as a list
  of vectors this allows for different covariates for each temporal
  trend. In this case the list needs to contain one vector for each of
  the temporal trends, starting with coefficients for the intercept,
  e.g.\cr \code{LUR = list(c(1,2,3),c(1,2),c(2))}. \cr
  \code{LUR = NULL} gives \emph{only an intercept}, no covariates.

  The option in \code{strip}, \code{strip.loc}, and \code{strip.time} to
  drop locations and times without observations can be used to reduce
  the dataset, thereby (drastically) speeding up the optimisation.
  To obtain predictions and simulations at the unobserved locations the
  original \code{mesa.data} structure can be passed to
  \code{\link{cond.expectation}} and \code{\link{simulateMesaData}}.
}
\value{
  Returns a \code{\link{mesa.data.model}} structure.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\section{Warning}{
  The structure \code{\link{mesa.data.model}} should \strong{always} be
  created by running \code{create.data.model} on a \code{\link{mesa.data}}
  structure.
}
\seealso{
  Used to create \code{\link{mesa.data.model}} from
  \code{\link{mesa.data}}. Uses \code{\link{default.LUR.list}}, \cr
  \code{\link{default.ST.list}}, \code{\link{construct.LUR.basis}}, and
  \code{\link{construct.ST.basis}} to extract covariates from
  \code{\link{mesa.data}}.
  
  For optimization functions see \code{\link{loglike}}, 
  \code{\link{loglike.var.names}}, \code{\link{run.MCMC}}, \cr
  \code{\link{fit.mesa.model}}, and \code{\link{cond.expectation}}.
  

  For other data functions, see \code{\link{mesa.data}},
  \code{\link{mesa.data.raw}}, \code{\link{create.data.matrix}}, \cr
  \code{\link{construct.LUR.basis}} and \code{\link{default.LUR.list}}.
}
\examples{
##load the data
data(mesa.data)

##create a mesa.data.model with all covariates
model.all <- create.data.model(mesa.data)

##no covariates
model.none <- create.data.model(mesa.data, LUR=NULL, ST.Ind=NULL)

##with different covariates for the different temporal trends
model.diff <- create.data.model(mesa.data, LUR=list(c(1,4,5,6),5:6,5:6),
                                ST.Ind=NULL)

##with the same covariates for all temporal trends, calling by name
model.same <- create.data.model(mesa.data, LUR=c("log10.m.to.a1",
    "log10.m.to.road","km.to.coast","s2000.pop.div.10000"),
    ST.Ind="lax.conc.1500") 
}
