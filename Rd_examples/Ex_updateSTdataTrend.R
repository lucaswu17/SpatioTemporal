##load data
data(mesa.data)

##default data and time trend for one location
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data)

##let's try with no trend
mesa.data <- updateSTdataTrend(mesa.data, n.basis=0)
plot(mesa.data)

##...and just one basis function, based on only AQS sites
subset <- mesa.data$covars$ID[mesa.data$covars$type=="AQS"]
mesa.data <- updateSTdataTrend(mesa.data, n.basis=1, subset=subset)
plot(mesa.data)

##Five basis functions, based on only AQS sites and much less smooth
mesa.data <- updateSTdataTrend(mesa.data, n.basis=5, subset=subset, df=100)
plot(mesa.data)
