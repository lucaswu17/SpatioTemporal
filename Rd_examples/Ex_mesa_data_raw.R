##load the data
data(mesa.data.raw)

##create a matrix of time-points
T <- t(matrix(rownames(mesa.data.raw$obs),
              nrow=dim(mesa.data.raw$obs)[2],
              ncol=dim(mesa.data.raw$obs)[1], byrow=TRUE))
##...and locations
ID <- matrix(colnames(mesa.data.raw$obs), 
             nrow=dim(mesa.data.raw$obs)[1],
             ncol=dim(mesa.data.raw$obs)[2], byrow=TRUE)
##let's study these matrices
print(T[1:5,1:5])
print(ID[1:5,1:5])

##combine with the observations
obs <- data.frame(obs=c(mesa.data.raw$obs), date=as.Date(T),
                  ID=c(ID))
##drop unmonitored locations
obs <- obs[!is.na(obs$obs),,drop=FALSE]
##sort the locations (strictly not needed)
obs <- obs[order(obs$date,obs$ID),,drop=FALSE]
##drop names
rownames(obs) <- NULL

##create a 3D-array for the spatio-temporal covariate
ST <- array(mesa.data.raw$lax.conc.1500, dim =
            c(dim(mesa.data.raw$lax.conc.1500),1))
dimnames(ST) <- list(rownames(mesa.data.raw$lax.conc),
                     colnames(mesa.data.raw$lax.conc),
                     "lax.conc.1500")
##create STdata object
mesa.data <- createSTdata(obs, mesa.data.raw$X,
                          SpatioTemporal=ST)

##optionally just give a list
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)
##works just as fine with the list
mesa.data.list <- createSTdata(obs, mesa.data.raw$X,
                               SpatioTemporal=ST.list)

##add  the smooth trends
mesa.data <- updateSTdataTrend(mesa.data, n.basis=2)

##This should be the same as the data in data(mesa.data)
