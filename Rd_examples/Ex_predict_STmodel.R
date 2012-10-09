##load data
data(mesa.model)
data(est.mesa.model)

##find regression parameters using GLS
x.reg <- predict(mesa.model, est.mesa.model, only.pars = TRUE)
str(x.reg$pars)

##Split data into FIXED and AQS
I.aqs <- mesa.model$locations$ID[mesa.model$locations$type=="AQS"]
I.aqs <- mesa.model$obs$ID %in% I.aqs
mesa.model.aqs <- dropObservations(mesa.model, !I.aqs)
mesa.model.fixed <- dropObservations(mesa.model, I.aqs)

\dontrun{
  ##compute predictions at all locations, including beta-fields
  pred.mesa.model <- predict(mesa.model, est.mesa.model,
                             pred.var=TRUE)
  ##predict at FIXED using AQS-sites
  pred.mesa.model.obs <- predict(mesa.model.aqs, est.mesa.model,
                                 STdata=mesa.model.fixed, only.obs=TRUE)
}
##Let's load precomputed results instead.
data(pred.mesa.model)

##Compare the predictions at all locations and only obs
print(pred.mesa.model)
print(pred.mesa.model.obs)

##estimate beta from the observations for reference
##create data matrix
D <- createDataMatrix(mesa.model)
beta <- matrix(NA,dim(D)[2], dim(mesa.model$trend)[2])
##extact the temporal trends
F <- mesa.model$trend
##drop the date column
F$date <- NULL
##estimate the beta-coeficients at each location
for(i in 1:dim(D)[2]){
  beta[i,] <- coefficients( lm(D[,i] ~ as.matrix(F)) )
}

##Study the results
##Start by comparing beta fields
par(mfcol=c(1,1), mar=c(4.5,4.5,2,.5), pty="s")
plot(x=beta[,1], y=pred.mesa.model$beta$EX[,1],
     main="Temporal Intercept",
     xlab="Empirical estimate", ylab="Spatio-Temporal Model")
abline(0,1,col="grey")
##or just the regression part of the beta fields
points(x=beta[,1], y=pred.mesa.model$beta$mu[,1], col=2)

##plot predictions and observations for 4 locations
par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
plot(pred.mesa.model, ID=1, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=10, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=17, STmodel=mesa.model, pred.var=TRUE)
plot(pred.mesa.model, ID=22, STmodel=mesa.model, pred.var=TRUE)
