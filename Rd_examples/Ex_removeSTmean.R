##load data
data(mesa.data)

mesa.data.mean0 <- removeSTcovarMean(mesa.data)

##compare the data structures
##geographic covariates
summary(mesa.data$covars)
summary(mesa.data.mean0$covars)

##mean of the spatio-temporal covariate, note that the new
##contains both mean-zero and original
cbind(colMeans(mesa.data$SpatioTemporal),
      colMeans(mesa.data.mean0$SpatioTemporal))

##mean of the spatio-temporal covariate
##compared to the added mean covariate
plot(mesa.data.mean0$covars$mean.lax.conc.1500,
     colMeans(mesa.data$SpatioTemporal))
\dontshow{
  if( max(abs(colMeans(mesa.data.mean0$SpatioTemporal[,,2]))) > 1e-10 ){
    stop("remove.ST.mean 1: mean(centred ST) != 0")
  }
  if( max(abs(mesa.data.mean0$covars$mean.lax.conc.1500 - 
              colMeans(mesa.data$SpatioTemporal))) > 1e-10 ){
    stop("remove.ST.mean 2: extracted mean != mean(ST)")
  }
}
