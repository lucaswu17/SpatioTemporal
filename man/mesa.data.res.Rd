\name{mesa.data.res}
\encoding{latin1}
\Rdversion{1.1}
\alias{mesa.data.res}
\docType{data}
\title{
  Results of some time consuming code.
}
\description{
  Data structure holding results of some time consuming code. These
  precomputed results is used in internal examples and in the
  package-vignette (Bergen and \enc{Lindström}{Lindstrom}, 2011).
}
\usage{data(mesa.data)}
\format{
  A list with elements. The code for recomputing each element is given
  under source below.

  The list contains
  \describe{
    \item{par.est}{The result of fitting the model to
      \code{\link{mesa.data.model}} using \code{\link{fit.mesa.model}}.
    }
    \item{EX}{Predictions at all locations computed using
      \code{\link{cond.expectation}}; this includes the latent
      beta-fields.
    }
    \item{EX.obs}{Predictions at only observations locations computed
      using \code{\link{cond.expectation}}.
    }
    \item{Ind.cv}{Grouping for the 10-fold cross-validation}
    \item{par.est.cv}{The result of a 10-fold cross-validated
      parameter estimation using \code{\link{estimateCV}}.}
    \item{pred.cv}{The cross-validation predictions, obtainable through
      \code{\link{predictCV}}. 
    }
    \item{par.est.ST}{The result of fitting a model with a
      spatio-temporal covariate using \cr
      \code{\link{fit.mesa.model}}.
    }
    \item{par.est.ST.mean0}{The result of fitting a model with a
      mean sperated spatio-temporal covariate using
      \code{\link{fit.mesa.model}}.
    }
    \item{MCMC.res}{The output from a parameter estimation using a
      Metropolis-Hastings algorithm \cr
      (implemented in \code{\link{run.MCMC}}) on the data in
      \code{\link{mesa.data.model}}.
    }
  }
}
\source{
  Contains model fitting results for the data in
  \code{\link{mesa.data.raw}}. This uses monitoring data from the
  \strong{MESA Air} project, see Cohen et.al. (2009) and
  \code{\link{mesa.data.raw}} for details. The elements have been
  produced as
  
  \describe{
    \item{par.est}{
    Estimated parameters, obtained as:\cr
    \code{data(mesa.data.model)}\cr
    \code{##initial values}\cr
    \code{dim <- loglike.dim(mesa.data.model)}\cr
    \code{x.init<-cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))}\cr
    \code{par.est <- fit.mesa.model(x.init, mesa.data.model,}\cr
    \code{hessian.all=TRUE,control=list(trace=3,maxit=1000))}
  }
  \item{EX}{
    \code{EX <- cond.expectation(par.est$res.best$par,}\cr
    \code{mesa.data.model, compute.beta=TRUE)}
  }
  \item{EX.obs}{
    \code{EX.obs <- cond.expectation(par.est$res.best$par,}\cr
    \code{mesa.data.model, mesa.data=mesa.data,}\cr
    \code{only.obs=TRUE, compute.beta=FALSE)}
  }    
  \item{Ind.cv}{
    Grouping for the 10-fold cross-validation, obtained as:\cr
    \code{data(mesa.data.model)}\cr
    \code{##create the CV structure defining 10 different CV-groups}\cr
    \code{Ind.cv <- createCV(mesa.data.model, groups=10, min.dist=.1)}
  }
  \item{par.est.cv}{
    Estimated parameters from 10-fold cross-validation, obtained as:\cr
    \code{data(mesa.data.model)}\cr
    \code{##create the CV structure defining 10 different CV-groups}\cr
    \code{Ind.cv <- createCV(mesa.data.model, groups=10, min.dist=.1)}\cr
    \code{##initial values}\cr
    \code{dim <- loglike.dim(mesa.data.model)}\cr
    \code{x.init<-cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))}\cr
    \code{##estimate different parameters for each CV-group}\cr
    \code{par.est.cv <- estimateCV(x.init, mesa.data.model, Ind.cv)}
  }
  \item{pred.cv}{
    The cross-validation predictions, obtainable through:\cr
    \code{pred.cv <- predictCV(par.est.cv$par, mesa.data.model, Ind.cv)}
  }
  \item{par.est.ST}{
    Estimated parameters with a spatio-temporal covariate, obtained
    as:\cr
    \code{data(mesa.data.model)}\cr
    \code{##model structure with ST-covariate}\cr
    \code{mesa.data.model.ST <- create.data.model(mesa.data,}\cr
    \code{LUR = list(c("log10.m.to.a1", "s2000.pop.div.10000",}\cr
    \code{"km.to.coast"), "km.to.coast", "km.to.coast"),}\cr
    \code{ST.Ind="lax.conc.1500") }\cr
    \code{##initial values}\cr
    \code{dim <- loglike.dim(mesa.data.model.ST)}\cr
    \code{x.init <- c(rep(c(1,-3),dim$m+1),-3)}\cr
    \code{##estimate parameters}\cr
    \code{par.est.ST <- fit.mesa.model(x.init, mesa.data.model.ST,}\cr
    \code{hessian.all=TRUE, control=list(trace=3,maxit=1000))}
  }
  \item{par.est.ST.mean0}{
    Estimated parameters with a mean sperated spatio-temporal covariate,
    obtained as:\cr 
    \code{data(mesa.data.model)}\cr
    \code{##model structure with ST-covariate}\cr
    \code{mesa.data.mean0 <- remove.ST.mean(mesa.data) }\cr
    \code{mesa.data.model.ST.mean0 <- create.data.model(mesa.data.mean0, }\cr
    \code{LUR = list(c("log10.m.to.a1", "s2000.pop.div.10000", }\cr
    \code{"km.to.coast", "mean.lax.conc.1500"),  "km.to.coast",}\cr
    \code{"km.to.coast"), ST.Ind="lax.conc.1500")}\cr
    \code{##initial values}\cr
    \code{dim <- loglike.dim(mesa.data.model.ST.mean0)}\cr
    \code{x.init <- c(rep(c(1,-3),dim$m+1),-3)}\cr
    \code{par.est.ST.mean0 <- fit.mesa.model(x.init,}\cr
    \code{mesa.data.model.ST.mean0, hessian.all=TRUE,}\cr
    \code{control=list(trace=3,maxit=1000))}
  }
  \item{MCMC.res}{
    Output from a MCMC algorithm started at the estimate in \cr
    \code{mesa.data.res$par.est}. The MCMC results were obtained as:\cr
    \code{data(mesa.data.model)}\cr
    \code{data(mesa.data.res)}\cr
    \code{##initial parameters and Hessian}\cr
    \code{x <- mesa.data.res$par.est$res.best$par.all}\cr
    \code{H <- mesa.data.res$par.est$res.best$hessian.all}\cr
    \code{##run MCMC}\cr
    \code{MCMC.res <- run.MCMC(x, mesa.data.model, N = 5000,}\cr
    \code{Hessian.prop = H)}
  }
  }
}
\references{
  M. A. Cohen, S. D. Adar, R. W. Allen, E. Avol, C. L. Curl, T. Gould,
  D. Hardie, A. Ho, P. Kinney, T. V. Larson, P. D. Sampson, L. Sheppard,
  K. D. Stukovsky, S. S. Swan, L. S. Liu, J. D. Kaufman. (2009) Approach
  to Estimating Participant Pollutant Exposures in the Multi-Ethnic
  Study of Atherosclerosis and Air Pollution (MESA Air). Environmental
  Science & Technology: 43(13), 4687-4693.

  Bergen, S. and \enc{Lindström, J.}{Lindstrom, J.}. (2011)
  Tutorial for Spatio-Temporal R Package
}
\seealso{
  Results from the estimation functions \code{\link{fit.mesa.model}},
  \code{\link{createCV}}, \code{\link{estimateCV}},
  \code{\link{predictCV}}, and \code{\link{run.MCMC}} using data in
  \code{\link{mesa.data.raw}}.
}
\examples{
##load data
data(mesa.data.model)
data(mesa.data.res)

##examining which components we have
names(mesa.data.res)

##examining the results for the different components

#################
## For par.est ##
#################
names(mesa.data.res$par.est)

##Optimisation status message
mesa.data.res$par.est$message

##extract the estimated parameters
x <- mesa.data.res$par.est$res.best$par.all
##and approximate uncertainties from the hessian
x.sd <- sqrt(diag(-solve(mesa.data.res$par.est$res.best$hessian.all)))

##plot the estimated parameters
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(x,ylab="",xlab="",xaxt="n")
abline(h=0, col="grey")
axis(1,1:length(x),names(x),las=2)

########################
## For EX and EX.obs  ##
########################
names(mesa.data.res$EX)
names(mesa.data.res$EX.obs)


####################
## For par.est.cv ##
####################
names(mesa.data.res$par.est.cv)

##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(7,2.5,2,.5), las=2)
boxplot(t(mesa.data.res$par.est.cv$par))
points(mesa.data.res$par.est$res.best$par, pch=4, col=2)

##################
## For pred.cv  ##
##################
names(mesa.data.res$pred.cv)

##Plot observations with CV-predictions and prediction intervals
par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
plotCV(mesa.data.res$pred.cv,  1, mesa.data.model)
plotCV(mesa.data.res$pred.cv, 10, mesa.data.model)
plotCV(mesa.data.res$pred.cv, 17, mesa.data.model)
plotCV(mesa.data.res$pred.cv, 22, mesa.data.model)

#########################################
## For par.est.ST and par.est.ST.mean0 ##
#########################################
names(mesa.data.res$par.est.ST)
names(mesa.data.res$par.est.ST.mean0)

##Optimisation status message
mesa.data.res$par.est.ST$message
mesa.data.res$par.est.ST.mean0$message

##extract the estimated parameters
x.ST <- mesa.data.res$par.est.ST$res.best$par.all
x.ST0 <- mesa.data.res$par.est.ST.mean0$res.best$par.all

##plot the estimated parameters
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
plot(c(1:5,7:19), x.ST, ylab="",xlab="",xaxt="n")
points(1:19, x.ST0, pch=3, col=2)
points(c(2:5,7:19), x, pch=4, col=3)
abline(h=0, col="grey")
axis(1,1:length(x.ST0),names(x.ST0),las=2)
legend("bottomleft", col = c(1,2,3), pch = 1:3, 
       legend=c("par.est.ST","par.est.ST.mean0","par.est.ST"))

##################
## For MCMC.res ##
##################
names(mesa.data.res$MCMC.res)

##The MCMC-estimated parameters
summary(mesa.data.res$MCMC.res$par)

##MCMC tracks for four of the parameters
par(mfrow=c(4,1),mar=c(2,2,2.5,.5))
for(i in c(4,9,13,15)){
  plot(mesa.data.res$MCMC.res$par[,i], ylab="", xlab="", type="l",
       main=colnames(mesa.data.res$MCMC.res$par)[i])
}

##And estimated densities for the log-covariance parameters.
##The red line is the approximate normal distribution given by
##the maximum-likelihood estimates, e.g. ML-estimate and standard 
##deviation from the observed information matrix.
par(mfrow=c(3,3),mar=c(4,4,2.5,.5))
for(i in 9:17){
  xd <- sort(unique(mesa.data.res$MCMC.res$par[,i]))
  yd <- dnorm(xd, mean=x[i],sd=x.sd[i])
  dens <- density(mesa.data.res$MCMC.res$par[,i])
  plot(dens, ylim=c(0,max(c(dens$y,yd))), main =
       colnames(mesa.data.res$MCMC.res$par)[i])
  lines(xd,yd,col=2)
}
}
\keyword{datasets}
