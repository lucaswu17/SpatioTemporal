##load data
data(mesa.model)
data(est.mesa.model)

\dontrun{
  ##compute predictions at all locations, including variances
  pred.mesa.model <- predict(mesa.model, est.mesa.model, pred.var=TRUE)
}

##separate data into AQS and FIXED sites
##first find indicators
I.aqs <- mesa.model$locations$ID[mesa.model$locations$type=="AQS"]
I.aqs <- mesa.model$obs$ID %in% I.aqs
##then drop relevant parts of STmodel
mesa.model.aqs <- dropObservations(mesa.model, !I.aqs)
mesa.model.fixed <- dropObservations(mesa.model, I.aqs)
                                     
\dontrun{
  ##compute predictions at only observed FIXED locations, based on observations
  ##at AQS locations.
  pred.mesa.model.obs <- predict(mesa.model.aqs, est.mesa.model,
                                 STdata=mesa.model.fixed, only.obs=TRUE)
}

##lets load precomputed results instead
data(pred.mesa.model)

##study results
print(pred.mesa.model)
print(pred.mesa.model.obs)

##simple plots
par(mfrow=c(2,1))
plot(pred.mesa.model)
plot(pred.mesa.model.obs, STmodel=mesa.model.fixed)

##scatter plot predicitons as function of observed
par(mfrow=c(1,1))
plot(pred.mesa.model.obs, y="obs", STmodel=mesa.model.fixed, ID="all",
     pred.var=TRUE, lty=c(NA,1),pch=c(19,NA))
