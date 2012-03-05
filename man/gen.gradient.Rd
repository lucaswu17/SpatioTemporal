\name{gen.gradient}
\encoding{latin1}
\Rdversion{1.1}
\alias{gen.gradient}
\alias{gen.hessian}
\title{
  Compute Finite Difference Gradient and Hessians.
}
\description{
  Computes finite difference gradient and/or
  hessian. \code{gen.gradient} function does forward, backward or
  central differences, the \code{gen.hessian} function uses only central
  differences.
}
\usage{
gen.gradient(x, func, h = 0.001, diff.type = 0)

gen.hessian(x, func, h = 0.001)
}
\arguments{
  \item{x}{
    Point at which to compute the gradient or hessian.
  }
  \item{func}{
    A function that takes only \code{x} as an input argument. Use \cr
    \code{function(x)\{my.func(x,other.input)\}} to create a temporary
    function, see the example.
  }
  \item{h}{
    Step length for the finite difference.
  }
  \item{diff.type}{
    Type of finite difference, positive values gives forward
    differences, \code{0} gives central differences, and negative values
    gives backward differences.
  }
}
\value{
  \code{gen.gradient} returns the gradient as a vector.
  \code{gen.hessian} returns the hessian matrix.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Used by \code{\link{loglike.grad}}, \code{\link{loglike.hessian}},
  \code{\link{loglike.naive.grad}}, and \cr
  \code{\link{loglike.naive.hessian}}.
}
\examples{
#create a two variable function
f.test <- function(x)\{sin(x[1])*cos(x[2])\}

#compute the gradient using forward difference
gen.gradient(c(.5,.5), f.test, diff.type=1)
#and central difference
gen.gradient(c(.5,.5), f.test, diff.type=0)
#compared to the true value
c(cos(.5)*cos(.5),-sin(.5)*sin(.5))

#Compute the Hessian
gen.hessian(c(.5,.5), f.test, h=1e-4)
#and compare to the true value
matrix(c(-sin(.5)*cos(.5),-cos(.5)*sin(.5),
         -cos(.5)*sin(.5),-sin(.5)*cos(.5)),2,2)
\dontshow{
if( max(abs(gen.gradient(c(.5,.5), f.test, h=1e-5, diff.type=0) -
            c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-10 ){
  stop("gen.gradient, diff.type=0: Results not equal")
}
if( max(abs(gen.gradient(c(.5,.5), f.test, h=1e-5, diff.type=1) -
            c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-5 ){
  stop("gen.gradient, diff.type=1: Results not equal")
}
if( max(abs(gen.gradient(c(.5,.5), f.test, h=1e-5, diff.type=-1) -
            c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-5 ){
  stop("gen.gradient, diff.type=-1: Results not equal")
}
if( max(abs(gen.hessian(c(.5,.5), f.test, h=1e-4) -
        matrix(c(-sin(.5)*cos(.5),-cos(.5)*sin(.5),
               -cos(.5)*sin(.5),-sin(.5)*cos(.5)),2,2))) > 1e-5 ){
  stop("gen.hessian: Results not equal")
}
}%\dontshow
}%\examples
