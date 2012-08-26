##load data
data(mesa.data)

##default plot
plot(mesa.data)

##plot monitor locations
plot(mesa.data, "loc")

##different names/colours/etc
plot(mesa.data, "loc", main="A nice plot", col=c("green","blue"),
    legend.names=c("Sites of one type", "..and of the other"),
    legend.loc="bottomleft")

##composite time-trend
plot(mesa.data, "loc.obs", legend.loc="bottomleft", cex=.5, pch=c(3,4))

##plot tim-series for the first site,
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID=1)
##residuals from the temporal trends,
plot(mesa.data, "res", ID=1)
##afc 
plot(mesa.data, "acf", ID=1)
##... and pafc for the residuals
plot(mesa.data, "pacf", ID=1)

##Same as above but calling the 2nd site by name
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID="60370016")
plot(mesa.data, "res", ID="60370016")
plot(mesa.data, "acf", ID="60370016")
plot(mesa.data, "pacf", ID="60370016")

##same, but with no temporal trend, first replace the trend with a constant
mesa.data <- updateSTdataTrend(mesa.data, n.basis=0)

par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID="60370016")
plot(mesa.data, "res", ID="60370016")
plot(mesa.data, "acf", ID="60370016")
plot(mesa.data, "pacf", ID="60370016")
