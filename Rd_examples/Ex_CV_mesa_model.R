##load data
data(mesa.model)

##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(rep(2,dim$nparam.cov),c(rep(c(1,-3),dim$m+1),-3))

\dontrun{
  ##estimate different parameters for each CV-group
  est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)
  ##Do cross-validated predictions using the just estimated parameters
  pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa$par.cov, 
                            est.cv.mesa$Ind.cv)
}
##lets load precomputed results instead
data(CV.mesa.model)
##and results of estimation
data(est.mesa.model)

##examining the estimation results
print( est.cv.mesa )
names( est.cv.mesa )
 
##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(12,2.5,2,.5), las=2)
boxplot( est.cv.mesa, plot.type="cov")
points( coef(est.mesa.model,"cov")$par, pch=19, col=2)
 
##examining the prediction
print( pred.cv.mesa )
names( pred.cv.mesa )
 
##Plot observations with CV-predictions and prediction intervals
##TODO TODO TODO
#par(mfcol=c(4,1),mar=c(2.5,2.5,2,.5))
#plotCV(mesa.data.res$pred.cv,  1, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 10, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 17, mesa.data.model)
#plotCV(mesa.data.res$pred.cv, 22, mesa.data.model)
