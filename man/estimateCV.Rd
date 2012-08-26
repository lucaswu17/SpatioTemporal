\name{estimateCV}
\alias{estimateCV}
\alias{predictCV}
\title{Model cross-validation}
\usage{
  estimateCV(object, x, Ind.cv, ...)

  predictCV(object, x, Ind.cv, ...)
}
\arguments{
  \item{object}{model to do cross-validation for}

  \item{x}{parameter values for estimation (initial) /
  predictions (predict given x)}

  \item{Ind.cv}{Structure definig the cross-validation
  groups.}

  \item{...}{additional parameters}
}
\value{
  Estimation / prediction result
}
\description{
  General cross-validation, estimation and prediction, for
  a model defined by object.
}
\author{
  Johan Lindström
}

