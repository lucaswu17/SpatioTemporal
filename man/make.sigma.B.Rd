\name{make.sigma.B}
\alias{block.mult}
\alias{calc.F.times.X}
\alias{calc.smooth.trends}
\alias{calc.tF.mat.F}
\alias{calc.tF.times.mat}
\alias{combineMesaData}
\alias{compute.ltaCV}
\alias{cond.expectation}
\alias{construct.LUR.basis}
\alias{construct.ST.basis}
\alias{create.data.matrix}
\alias{create.data.model}
\alias{CVbasics}
\alias{CVresiduals.qqnorm}
\alias{CVresiduals.scatter}
\alias{default.LUR.list}
\alias{default.ST.list}
\alias{detrend.data}
\alias{dot.prod}
\alias{drop.observations}
\alias{fit.mesa.model}
\alias{gen.gradient}
\alias{gen.hessian}
\alias{get.params}
\alias{loglike}
\alias{loglike.dim}
\alias{loglike.grad}
\alias{loglike.hessian}
\alias{loglike.naive}
\alias{loglike.naive.grad}
\alias{loglike.naive.hessian}
\alias{loglike.var.names}
\alias{make.sigma.B}
\alias{make.sigma.B.full}
\alias{make.sigma.nu}
\alias{make.sigma.nu.cross.cov}
\alias{plotCV}
\alias{plotMesaData}
\alias{plotMonitoringLoc}
\alias{plotPrediction}
\alias{printMesaDataNbrObs}
\alias{remove.ST.mean}
\alias{run.MCMC}
\alias{setupSTdataset}
\alias{simulateMesaData}
\alias{summaryStatsCV}
\alias{SVD.miss}
\alias{SVD.smooth}
\alias{SVD.smooth.cv}
\alias{tstat}
\title{Deprecated functions, use replacements!}
\usage{
  make.sigma.B(...)

  make.sigma.B.full(...)

  make.sigma.nu(...)

  make.sigma.nu.cross.cov(...)

  calc.tF.times.mat(...)

  calc.F.times.X(...)

  calc.tF.mat.F(...)

  block.mult(...)

  dot.prod(...)

  SVD.miss(...)

  SVD.smooth(...)

  SVD.smooth.cv(...)

  calc.smooth.trends(...)

  setupSTdataset(...)

  printMesaDataNbrObs(...)

  plotMonitoringLoc(...)

  plotMesaData(...)

  create.data.matrix(...)

  remove.ST.mean(...)

  detrend.data(...)

  create.data.model(...)

  default.LUR.list(...)

  default.ST.list(...)

  construct.LUR.basis(...)

  construct.ST.basis(...)

  loglike.dim(...)

  loglike.var.names(...)

  get.params(...)

  loglike(...)

  loglike.naive(...)

  gen.gradient(...)

  gen.hessian(...)

  loglike.grad(...)

  loglike.hessian(...)

  loglike.naive.grad(...)

  loglike.naive.hessian(...)

  fit.mesa.model(...)

  cond.expectation(...)

  simulateMesaData(...)

  combineMesaData(...)

  drop.observations(...)

  plotPrediction(...)

  tstat(...)

  compute.ltaCV(...)

  CVbasics(...)

  summaryStatsCV(...)

  run.MCMC(...)

  plotCV(...)

  CVresiduals.qqnorm(...)

  CVresiduals.scatter(...)
}
\arguments{
  \item{...}{Unused, for compability.}
}
\value{
  Does not return.
}
\description{
  Deprecated functions, use replacements!
}
\details{
  Functions have been rename/replaced as: \describe{
  \item{block.mult}{\code{\link{blockMult}}}

  \item{calc.F.times.X}{\code{\link{calc.FX}}}
  \item{calc.smooth.trends}{\code{\link{calcSmoothTrends}}
  and \code{\link{updateTrend}}}
  \item{calc.tF.mat.F}{\code{\link{calc.tFXF}}}
  \item{calc.tF.times.mat}{\code{\link{calc.tFX}}, see also
  \code{\link{expandF}}.}
  \item{combineMesaData}{\code{\link{c.STmodel}}}
  \item{compute.ltaCV}{\code{\link{computeLTA}}}
  \item{cond.expectation}{\code{\link{predict.STmodel}}}
  \item{construct.LUR.basis}{\code{\link{createLUR}}}
  \item{construct.ST.basis}{\code{\link{createST}}}
  \item{create.data.matrix}{\code{\link{createDataMatrix}}}
  \item{create.data.model}{\code{\link{createSTmodel}}, see
  also \code{\link{updateCovf}} and
  \code{\link{processLocation}}.} \item{CVbasics}{Included
  in \code{\link{estimateCV.STmodel}}.}
  \item{CVresiduals.qqnorm}{\code{\link{qqnorm.STdata}},
  \code{\link{qqnorm.STmodel}}, or
  \code{\link{qqnorm.predCVSTmodel}}}
  \item{CVresiduals.scatter}{\code{\link{scatterPlot.STdata}},
  \code{\link{scatterPlot.STmodel}}, or
  \code{\link{scatterPlot.predCVSTmodel}}}

  \item{default.LUR.list}{\code{\link{processLUR}}}
  \item{default.ST.list}{\code{\link{processST}}}
  \item{detrend.data}{\code{\link{detrendSTdata}}}
  \item{dot.prod}{\code{\link{dotProd}}}
  \item{drop.observations}{\code{\link{dropObservations}}}

  \item{fit.mesa.model}{\code{\link{estimate.STmodel}}}

  \item{gen.gradient}{\code{\link{genGradient}}}
  \item{gen.hessian}{\code{\link{genHessian}}}
  \item{get.params}{\code{\link{loglikeSTgetPars}}}

  \item{loglike}{\code{\link{loglikeST}}}
  \item{loglike.dim}{\code{\link{loglikeSTdim}}}
  \item{loglike.grad}{\code{\link{loglikeSTGrad}}}
  \item{loglike.hessian}{\code{\link{loglikeSTHessian}}}
  \item{loglike.naive}{\code{\link{loglikeSTnaive}}}
  \item{loglike.naive.grad}{\code{\link{loglikeSTnaiveGrad}}}
  \item{loglike.naive.hessian}{\code{\link{loglikeSTnaiveHessian}}}
  \item{loglike.var.names}{\code{\link{loglikeSTnames}}}

  \item{make.sigma.B}{\code{\link{makeSigmaB}}, see
  \code{\link{parsCovFuns}}, \code{\link{namesCovFuns}},
  and \code{\link{evalCovFuns}} for new covariance
  specifications.}
  \item{make.sigma.B.full}{\code{\link{makeSigmaB}} and
  \code{\link{calc.FXtF2}}}
  \item{make.sigma.nu}{\code{\link{makeSigmaNu}}}
  \item{make.sigma.nu.cross.cov}{\code{\link{makeSigmaNu}}}

  \item{printMesaDataNbrObs}{\code{\link{print.STdata}},
  \code{\link{summary.STdata}},
  \code{\link{print.STmodel}}, and
  \code{\link{summary.STmodel}}}
  \item{plotCV}{\code{\link{plot.predCVSTmodel}}}
  \item{plotMesaData}{\code{\link{plot.STdata}} and
  \code{\link{plot.STmodel}}}
  \item{plotMonitoringLoc}{\code{\link{plot.STdata}} and
  \code{\link{plot.STmodel}}}
  \item{plotPrediction}{\code{\link{plot.predictSTmodel}}}

  \item{remove.ST.mean}{\code{\link{removeSTcovarMean}}}
  \item{run.MCMC}{\code{\link{MCMC.STmodel}}}

  \item{setupSTdataset}{\code{\link{createSTdata}}}
  \item{simulateMesaData}{\code{\link{simulate.STmodel}}}
  \item{summaryStatsCV}{\code{\link{summary.predCVSTmodel}}}
  \item{SVD.miss}{\code{\link{SVDmiss}}}
  \item{SVD.smooth}{\code{\link{SVDsmooth}}}
  \item{SVD.smooth.cv}{\code{\link{SVDsmoothCV}}, see also
  \code{\link{plot.SVDcv}} and \code{\link{print.SVDcv}}.}

  \item{tstat}{Included in \code{\link{predict.STmodel}}} }
}
\author{
  Johan Lindström
}

