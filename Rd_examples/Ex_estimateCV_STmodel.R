##load data
data(mesa.model)

################
## estimateCV ##
################
##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))

\dontrun{
  ##estimate different parameters for each CV-group
  est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)

##examine the estimation results
print( est.cv.mesa )
##estimated parameters for each CV-group
est.cv.mesa$par.cov

##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(12,2.5,2,.5), las=2)
boxplot( est.cv.mesa, plot.type="cov")

###############
## predictCV ##
###############
\dontrun{
  ##Do cross-validated predictions using the just estimated parameters
  pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, est.cv.mesa$Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)

##examine the prediction results
names( pred.cv.mesa )
print( pred.cv.mesa )
##compute CV-statistics
pred.cv.stats <- summary( pred.cv.mesa, LTA=TRUE)
print( pred.cv.stats )
 
##Plot observations with CV-predictions and prediction intervals
##TODO TODO TODO
#par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
#plotCV(mesa.data.res$pred.cv,  1, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 10, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 17, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 22, mesa.data.model)

##Residual QQ-plots
#par(mfrow=c(1,2),mar=c(3,2,1,1),pty="s")
##Raw Residual QQ-plot
#CVresiduals.qqnorm(pred.cv.stats$res)
##Normalized Residual QQ-plot
#CVresiduals.qqnorm(pred.cv.stats$res.norm, norm=TRUE)

\dontrun{
  ##A faster option is to only consider the observations and not to compute
  ##variances
  pred.cv.fast <- predictCV(mesa.model, est.cv.mesa, est.cv.mesa$Ind.cv,
                            only.obs=TRUE, pred.var=FALSE)
  print( pred.cv.fast )
  summary( pred.cv.fast )
}
