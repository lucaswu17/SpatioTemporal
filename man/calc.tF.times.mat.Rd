\name{calc.tF.times.mat}
\encoding{latin1}
\Rdversion{1.1}
\alias{calc.tF.times.mat}
\alias{calc.tF.mat.F}
\alias{calc.F.times.X}

\title{
  Compute Matrix Productes Involving the Smooth Trends
}
\description{
  Computes the matrix products between a sparse matrix \code{F}
  containing the temporal basis functions and other matrices.
}
\usage{
calc.tF.times.mat(X, F, loc.ind, n.x = dim(X)[2],
     n.loc = max(loc.ind))

calc.tF.mat.F(mat, F, loc.ind, n.blocks = 1, block.sizes =
     rep(dim(mat)[1]/n.blocks, n.blocks), n.loc = max(loc.ind))

calc.F.times.X(LUR, F, loc.ind)
}
\arguments{
  \item{X}{
    A vector or matrix.
  }
  \item{F}{
    A (number of obs.) - by - (number of temporal trends) matrix
    containing the temporal trends. Usually \code{mesa.data.model$F},
    where \code{\link{mesa.data.model}} is obtained from
    \code{\link{create.data.model}}.
  }
  \item{loc.ind}{
    A vector indicating which location each row in \code{F} corresponds
    to, usually \cr \code{mesa.data.model$obs$idx}.
  }
  \item{n.loc}{
    Number of locations.
  }
  \item{n.x}{
    Number of columns in \code{X}, defaults correctly if \code{X} is a
    matrix.
  }
  \item{mat}{
    A block diagonal, square matrix.
  }
  \item{n.blocks}{
    Number of diagonal blocks in \code{mat}. Defaults to 1
    (i.e. a full matrix) if not given.
  }
  \item{block.sizes}{
    A vector of length \code{n.blocks} with the size of each of the
    diagonal blocks. If not given it will assume equal size blocks.
  }
  \item{LUR}{
    A list of matrices, usually \code{mesa.data.model$X}. Each matrix in
    the list should have the same number of rows, but the number of
    columns may vary.
  }
}
\details{
  \code{calc.tF.times.mat} computes \code{F} times a matrix.
  \code{calc.tF.mat.F} computes the quadratic form \code{F}' * 
  \code{mat} * \code{F} where \code{mat} is a block diagonal, square
  matrix. And \code{calc.F.times.X} computes \code{F} times a block
  matrix represented as a list of matrices.

  See the examples for details.
}
\value{
  Returns a matrix with the result of the relevant matrix computation(s).
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\seealso{
  Block matrices can be created by \code{\link{make.sigma.B}}, 
  \code{\link{make.sigma.B.full}}, and \cr \code{\link{make.sigma.nu}}.
  
  Other block matrix functions \code{\link{makeCholBlock}}, 
  \code{\link{block.mult}}, \code{\link{calc.iS.X}}, 
  \code{\link{dot.prod}}, and \cr \code{\link{sumLogDiag}}.
  
  This function is called by \code{\link{loglike}}.
}
\examples{
##This starts with a couple of simple examples, more elaborate examples
##with real data can be found further down.

##create a trend
trend <- cbind(1:5,sin(1:5))
##an index of locations
idx <- c(rep(1:3,3),1:2,2:3)
##a list of time points for each location/observation
T <- c(rep(1:3,each=3),4,4,5,5)
##create a random observations matrix
obs <- rnorm(length(T))

##expand the F matrix to match the locations/times in idx/T.
F <- trend[T,]

##compute tF \%*\% obs
tFobs <- calc.tF.times.mat(obs, F, idx, n.x=1)

##alternatievly this can be computed as observtions for each location
##multiplied by the trend function at the corresponding time points.
tFobs.alt <- double(length(tFobs))
for(j in 1:dim(trend)[2]){
  for(i in 1:max(idx)){
    tFobs.alt[i+(j-1)*max(idx)] <- sum(obs[idx==i]*trend[T[idx==i],j])
  }
}
##compare results
print( cbind(tFobs,tFobs.alt) )
\dontshow{
if( max(abs(tFobs-tFobs.alt)) > 1e-13 ){
  stop("calc.tF.times.mat 1: Results not equal")
}
}
##create a covariance matrix for locations and each of 
C <- make.sigma.nu(1, 0, 1, c(3,3,3,2,2), idx,
                   as.matrix(dist(1:max(idx))) )
##compute F' \%*\% C \%*\% F
tFmatF <- calc.tF.mat.F(C, F, idx, block.sizes = c(3,3,3,2,2),
                        n.blocks = 5)
##which is equivalent of
tFmatF.alt <- calc.tF.times.mat(t(calc.tF.times.mat(C, F, idx)), F, idx)

print(tFmatF-tFmatF.alt)
\dontshow{
if( max(abs(tFmatF-tFmatF.alt)) > 1e-13 ){
  stop("calc.tF.times.mat 2: Results not equal")
}
}

##some examples using the mesa.data
data(mesa.data.model)

##Some information about the size(s) of the model.
dim <- loglike.dim(mesa.data.model)

##compute F' \%*\% obs
tFobs <- calc.tF.times.mat(mesa.data.model$obs$obs, mesa.data.model$F,
                           mesa.data.model$obs$idx,n.x=1)

##The resulting matrix contains 75 elements (3 temporal trend at 25
##sites). The first element are the observations at the first site
##multiplied by the constant temporal trend, e.g.
print(tFobs[1])
print( sum(mesa.data.model$obs$obs[mesa.data.model$obs$idx==1]) )
\dontshow{
if( abs(sum(mesa.data.model$obs$obs[mesa.data.model$obs$idx==1]) -
        tFobs[1]) > 1e-10 ){
  stop("calc.tF.times.mat - real data 1: Results not equal")
}
}

##The 27:th element are the observations at the second site (27-25)
##multiplied by the first temporal trend (second element in F)
print(tFobs[dim$n+2])
print(sum(mesa.data.model$obs$obs[mesa.data.model$obs$idx==2] *
          mesa.data.model$F[mesa.data.model$obs$idx==2,2]))
\dontshow{
if( abs(sum(mesa.data.model$obs$obs[mesa.data.model$obs$idx==2] *
            mesa.data.model$F[mesa.data.model$obs$idx==2,2]) -
        tFobs[dim$n+2]) > 1e-10 ){
  stop("calc.tF.times.mat - real data 2: Results not equal")
}
}

##compute F \%*\% X
FX <- calc.F.times.X(mesa.data.model$X, mesa.data.model$F,
                     mesa.data.model$obs$idx)
##The resulting matrix is
##(number of time points) - by - (number of land use covariates)
##where the number of land use covariates are computed over all the
##two + intercept temporal trends.
##Each column contains the temporal trend for the observations
##multiplied by the corresponding LUR-covariate
par(mfrow=c(3,1))
plot(FX[,2])
points(mesa.data.model$X[[1]][mesa.data.model$obs$idx,2] *
       mesa.data.model$F[,1],col=2,pch=3)

plot(FX[,dim$p[1]+1])
points(mesa.data.model$X[[2]][mesa.data.model$obs$idx,1] *
       mesa.data.model$F[,2],col=2,pch=3)

plot(FX[,dim$p[1]+dim$p[2]+2])
points(mesa.data.model$X[[3]][mesa.data.model$obs$idx,2] *
       mesa.data.model$F[,3],col=2,pch=3)

##If the regression parameters, alpha, are known (or estimated)
##The intercept part of the model is given by FX \%*\% alpha
\dontshow{
if( max(abs(FX[,2] - mesa.data.model$X[[1]][mesa.data.model$obs$idx,2] *
        mesa.data.model$F[,1])) > 1e-10){
  stop("calc.tF.times.mat - real data 3a: Results not equal")
}
if( max(abs(FX[,dim$p[1]+1] - mesa.data.model$F[,2] *
            mesa.data.model$X[[2]][mesa.data.model$obs$idx,1])) > 1e-10){
  stop("calc.tF.times.mat - real data 3b: Results not equal")
}
if( max(abs(FX[,dim$p[1]+dim$p[2]+2] - mesa.data.model$F[,3] *
            mesa.data.model$X[[3]][mesa.data.model$obs$idx,2])) > 1e-10){
  stop("calc.tF.times.mat - real data 3c: Results not equal")
}
}%\dontshow
}%\examples
