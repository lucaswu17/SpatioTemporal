##load the data
data(mesa.data)

##Look at the number of observations/locations
print(mesa.data)

##Look at the summary of observations, covariates, etc.
summary(mesa.data)

##Lets look at the data
names(mesa.data)

##Study the structure of the covariates data
head(mesa.data$covars)

##...the smooth temporal trends
head(mesa.data$trend)

##...observations
head(mesa.data$obs)

##...and Spatio-temporal covariate
mesa.data$SpatioTemporal[1:5,1:5,,drop=FALSE]

##Let's plot the space-time monitoring locations
plot(mesa.data, "loc")

##Let's plot the observations as a function of time
plot(mesa.data, "loc.obs", legend.loc="bottomleft")

##plot observations and residuals from the temporal trends
par(mfcol=c(3,2),mar=c(2.5,2.5,2,.5))
plot(mesa.data, "obs", ID=5)
plot(mesa.data, "res", ID=5)
plot(mesa.data, "acf", ID=5)
plot(mesa.data, "obs", ID=18)
plot(mesa.data, "res", ID=18)
plot(mesa.data, "acf", ID=18)

##create STmodel object
##define land-use covariates, for intercept and trends
LUR = list(c("log10.m.to.a1", "s2000.pop.div.10000", "km.to.coast"),
  "km.to.coast", "km.to.coast")
##and covariance model
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=TRUE, random.effect=FALSE)
##which locations to use
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")
##create object
mesa.model <- createSTmodel(mesa.data, LUR=LUR, cov.beta=cov.beta,
                            cov.nu=cov.nu, locations=locations)

##This should be the same as the data in data(mesa.model)
