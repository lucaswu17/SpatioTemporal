####################################################################
## FILE CONTAINING DERIVATIVE COMPUTATIONS FOR THE LOGLIKELIHOODS ##
####################################################################
##general function that computes gradients
gen.gradient <- function(x, func, h=1e-3, diff.type=0){
  diff.type <- sign(diff.type) ##ensure diff.type is [-1,0,1]
  if(diff.type==0){ ##central
    f.p <- double(length(x))
    f.m <- double(length(x))
    for(i in 1:length(x)){
      dx <- x
      dx[i] <- x[i]+h/2
      f.p[i] <- func(dx)
      dx[i] <- x[i]-h/2
      f.m[i] <- func(dx)
    }
    df <- (f.p-f.m)/h
  }else{ ##forward or backward
    f <- func(x)
    df <- double(length(x))
    for(i in 1:length(x)){
      dx <- x
      dx[i] <- dx[i]+diff.type*h
      df[i] <- func(dx)
    }
    df <- diff.type*(df-f)/h
  }
  return(df)
}

##functions that compute gradients for the loglikelihoods
loglike.grad <- function(x, mesa.data.model, type="p",
                         h=1e-3, diff.type=0){
  func <- function(x0){ loglike(x0, mesa.data.model, type) }
  df <- gen.gradient(x, func, h=h, diff.type=diff.type)
  return(df)
}

##general function that computes hessians
gen.hessian <- function(x, func, h=1e-3){
  ##uses only central differences
  N <- length(x)
  ##allocate matrix for the hessian
  hessian <- matrix(NA,N,N)
  ##start with the diagonal elements
  f <- double(5)
  for(i in 1:N){
    tmp <- x
    for(j in 1:5){
      tmp[i] <- x[i]+(j-3)*h
      f[j] <- func(tmp)
    }
    hessian[i,i] <- (-f[1]+16*f[2]-30*f[3]+16*f[4]-f[5])/(12*h^2)
  }
  ##off diagonal elements
  f <- double(4)
  Ind <- cbind(c(1,1,-1,-1),c(1,-1,1,-1))
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      tmp <- x
      for(k in 1:4){
        tmp[i] <- x[i] + Ind[k,1]*h
        tmp[j] <- x[j] + Ind[k,2]*h
        f[k] <- func(tmp)
      }
      hessian[i,j] <- (f[1]-f[2]-f[3]+f[4])/(4*h^2)
      ##use symmetri for the other half of the matrix
      hessian[j,i] <- hessian[i,j]
    }
  }
  return(hessian)
}

##functions that compute hessians for the loglikelihoods
loglike.hessian <- function(x, mesa.data.model, type="p", h=1e-3){
  func <- function(x0){ loglike(x0, mesa.data.model, type) }
  H <- gen.hessian(x, func, h=h)
  return(H)
}

##Naive functions, as timing alternatives.
loglike.naive.grad <- function(x, mesa.data.model, type="p",
                              h=1e-3, diff.type=0){
  func <- function(x0){ loglike.naive(x0, mesa.data.model, type) }
  df <- gen.gradient(x, func, h=h, diff.type=diff.type)
  return(df)
}

loglike.naive.hessian <- function(x, mesa.data.model, type="p", h=1e-3){
  func <- function(x0){ loglike.naive(x0, mesa.data.model, type) }
  H <- gen.hessian(x, func, h=h)
  return(H)
}
