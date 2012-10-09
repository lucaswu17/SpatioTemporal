##load the data
data(mesa.data)

##create a data matrix
M1 <- createDataMatrix(mesa.data)
dim(M1)
head(M1)

##create data matrix for only a few locations
M2 <- createDataMatrix(mesa.data, subset =
                         c("60370002","60370016","60370113","60371002",
                           "60371103","60371201","L001","L002"))
dim(M2)
head(M2)
\dontshow{
  if( (dim(M1)[1]!=dim(mesa.data$trend)[1]) ||
     (dim(M1)[2]!=dim(mesa.data$covars)[1]) ){
    stop("createDataMatrix: dimension missmatch - M1")
  }
  if( (dim(M2)[1]!=dim(mesa.data$trend)[1]) || (dim(M2)[2]!=8) ){
    stop("createDataMatrix: dimension missmatch - M2")
  }
  if( max(abs(M1[,colnames(M2)]-M2),na.rm=TRUE) > 1e-13 ){
    stop("createDataMatrix: M1!=M2")
  }
}

