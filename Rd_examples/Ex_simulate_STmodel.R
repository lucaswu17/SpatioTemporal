##load the data
data(mesa.data)
data(mesa.model)
data(est.mesa.model)

##Get estimated parameters
x <- coef(est.mesa.model)$par

##Simulate 5 replicates from these parameters
sim.data <- simulate(mesa.model, nsim=5, x=x)

##compute average beta fields
beta <- calc.mu.B(mesa.model$LUR, loglikeSTgetPars(x, mesa.model)$alpha)

##plot the simulated observations as a function of time
par(mfrow=c(2,2), mar=c(4,4,.5,.5))
plot(sim.data$obs[[1]]$date, sim.data$obs[[1]]$obs,
     type="n", ylab="obs", xlab="Date")
for(i in 1:5){
  points(sim.data$obs[[i]]$date, sim.data$obs[[i]]$obs, col=i)
}
##and the latent beta-fields
for(i in 1:3){
  plot(sim.data$B[,i,1], ylim=range(sim.data$B[,i,]), type="n",
       xlab="loc", ylab=paste("beta",colnames(sim.data$B)[i]))
  for(j in 1:5){
    points(sim.data$B[,i,j], col=j)
  }
  lines( beta[,i], col="grey")
}

###########################################
## A case with some unobserved locations ##
###########################################
##keep only observations from the AQS sites
I.aqs <- mesa.data$covars$ID[ mesa.data$covars$type=="AQS"]
mesa.data$obs <- mesa.data$obs[ mesa.data$obs$ID %in% I.aqs,]

##create a new model object
mesa.model <- createSTmodel(mesa.data, LUR = mesa.model$LUR.list)

##simulate some replicates for this object
sim.data2 <- simulate(mesa.model, nsim=1, x=x)

##The new object now only contains simulated observations for the AQS sites
table( sim.data2$obs[[1]]$ID )
##while the other two contain all locations
table( sim.data$obs[[1]]$ID )
