###########################
## INTERFACES FOR C CODE ##
###########################

##############################################
## CODE THAT CONSTRUCTS COVARIANCE MATRICES ##
##############################################
##constructs the covariance matrix for the land use regression
##inputs are:
##number of temporal trends (incl. constant)
##number of locations
##one sill for each trend
##one range for each trend
##distance matrix with the distances between each location
make.sigma.B <- function(sill, range, dists){
  ##check dimensions
  if( length(range) != length(sill) )
    stop("In 'make.sigma.B': number of elements in 'range' and 'sill' differ")
  if(dim(dists)[1] != dim(dists)[2])
    stop("In 'make.sigma.B': 'dists' not a square matrix")
  ##call the c-function
  .Call("make_sigma_b", as.integer(length(sill)),
        as.integer(dim(dists)[1]), sill, range, dists)
}

##constructs the full cross-covariance matrix for the land use regression
##inputs are:
##number of observations for the first locations
##number of observations for the second locations
##number of temporal trends (incl. constant)
##number of locations for the first set of locations
##number of locations for the second set of locations
##location of each observation #1 (used to index the distance matrix)
##  (n_obs2-by-1)
##location of each observation #2 (used to index the distance matrix)
##  (n_obs2-by-1)
##one sill for each trend
##  (n_trends-by-1)
##one range for each trend
##  (n_trends-by-1)
##matrix with the temporal trends for the first set of locations 
##  (n_obs1-by-n_trends)
##matrix with the temporal trends for the second set of locations
##  (n_obs2-by-n_trends)
##distance matrix with the distances between each location
##  (n_loc1-by-n_loc2)
make.sigma.B.full <- function(sill, range, loc.ind1,
                              loc.ind2=loc.ind1, F1, F2=F1, dists){
  ##check dimensions
  if(max(loc.ind1) > dim(dists)[1])
    stop("In 'make.sigma.B.full': max(loc.ind1) > dim(dists)[1]")
  if(max(loc.ind2) > dim(dists)[2])
    stop("In 'make.sigma.B.full': max(loc.ind2) > dim(dists)[2]")
  if(dim(F1)[2] != dim(F2)[2] || length(sill) != length(range) ||
     dim(F1)[2] != length(sill))
    stop("In 'make.sigma.B.full': inconsistent number of components in 'F1', 'F2', 'sill', and/or 'range'")
  
  ##call the c-function
  .Call("make_sigma_b_full", as.integer(dim(F1)[1]),
        as.integer(dim(F2)[1]), as.integer(dim(F1)[2]),
        as.integer(dim(dists)[1]), as.integer(dim(dists)[2]),
        as.integer(loc.ind1), as.integer(loc.ind2), sill, range,
        F1, F2, dists)
}

##construct the covariance matrix for the spatio temporal residuals
##assuming the same spatial covariance for all time-points and no
##temporal dependence
##inputs are:
##sill, nugget and range
##total number of points
##number of blocks (time points)
##size of each block (no. obs per time point)
##list of the location index of each obs
##number of locations (size of distance matrix)
##distance matrix with distance between each location
make.sigma.nu <- function(sill, nugget, range, block.sizes, loc.index, dists){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##extract sizes
  n.tot <- length(loc.index)
  n.dist <- dim(dists)[1]
  n.blocks <- length(block.sizes)
  ##check dimensions
  if(dim(dists)[1] != dim(dists)[2])
    stop("In 'make.sigma.nu': 'dists' not a square matrix")
  if(max(loc.index) > n.dist)
    stop("In 'make.sigma.nu': max(loc.index) > n.dist")
  if(sum(block.sizes) != n.tot)
    stop("In 'make.sigma.nu': sum(block.sizes) != length(loc.index)")
  ##call the c-function
  .Call("make_sigma_nu", sill, nugget, range, as.integer(n.tot),
        as.integer(n.blocks), as.integer(block.sizes),
        as.integer(loc.index), as.integer(n.dist), dists)
}

##constructs the full cross-covariance matrix for the spatio temporal 
##residuals, assuming the same spatial covariance for all time-points 
##and no temporal dependence
##inputs are:
##sill, nugget and range
##number of observations for the first locations
##number of observations for the second locations
##number of locations for the first set of obs.
##number of locations for the second set of obs.
##location of each observation #1 (used to index the distance matrix)
##  (n_obs1-by-1)
##location of each observation #2 (used to index the distance matrix)
##  (n_obs2-by-1)
##observation time for observations #1
##  (n_obs1-by-1)
##observation time for observations #2
##  (n_obs2-by-1)
##distance matrix with the distances between each location
##  (n_loc1-by-n_loc2)
make.sigma.nu.cross.cov <- function(sill, nugget, range,
                                     loc.ind1, loc.ind2=loc.ind1,
                                     loc.ind2.to.1=1:max(loc.ind2),
                                     T1, T2=T1, dists){
  ##check dimensions - TODO
  if(length(loc.ind1) != length(T1))
    stop("In 'make.sigma.nu.cross.cov': length(loc.ind1) != length(T1)")
  if(length(loc.ind2) != length(T2))
    stop("In 'make.sigma.nu.cross.cov': length(loc.ind2) != length(T2)")
  if(max(loc.ind1) > dim(dists)[1])
    stop("In 'make.sigma.nu.cross.cov': max(loc.ind1) > dim(dists)[1]")
  if(max(loc.ind2) > dim(dists)[2])
    stop("In 'make.sigma.nu.cross.cov': max(loc.ind2) > dim(dists)[2]")
  if(max(loc.ind2) > length(loc.ind2.to.1))
    stop("In 'make.sigma.nu.cross.cov': max(loc.ind2) > length(loc.ind2.to.1)")
  ##call the c-function
  .Call("make_sigma_nu_cross_cov", sill, nugget, range,
        as.integer(length(loc.ind1)), as.integer(length(loc.ind2)),
        as.integer(dim(dists)[1]),as.integer(dim(dists)[2]),
        as.integer(loc.ind1), as.integer(loc.ind2),
        as.integer(loc.ind2.to.1), as.integer(T1), as.integer(T2),
        dists)
}

#############################################
## CODE THAT DOES SPARSE F TIMES SOMETHING ##
#############################################
##calculates F'*X, inputs:
##A matrix X
##The matrix F (for a given column each row corresponds to the value of the
##  temporal trend at the time of that observations, trends vary between columns)
##number of different locations
##Vector of locations, giving the location of each observation. Needed to know
##  which element of sigma.nu that we should multiply elements in F with.
##number of observations
##number of temporal trends
##returns a (n_loc*m)-by-dim(X)[2] matrix
calc.tF.times.mat <- function(X, F, loc.ind, n.x=dim(X)[2],
                              n.loc=max(loc.ind)){
  ##check dimensions
  if(length(X) != dim(F)[1]*n.x)
    stop("In 'calc.tF.times.mat': length(X) != dim(F)[1]*n.x")
  if(max(loc.ind) > n.loc)
    stop("In 'calc.tF.mat.F': max(loc.ind) > n.loc")
  ##call the c-function
  .Call("calc_tF_times_mat", X, F, as.integer(n.loc), as.integer(n.x),
        as.integer(loc.ind), as.integer(dim(F)[1]), as.integer(dim(F)[2]))
}

##calculates F*X, inputs:
##A list containing the LUR:s
##The matrix F (for a given column each row corresponds to the
##  value of the temporal trend at the time of that observations,
##  trends vary between columns)
##Vector of locations, giving the location of each observation. Needed to know
##  which element of sigma.nu that we should multiply elements in F with.
##number of different locations
##returns a column vector with (n_obs*m) elements
calc.F.times.X <- function(LUR, F, loc.ind){
  m <- length(LUR)
  p <- double(m)
  n.loc <- dim(LUR[[1]])[1]
  ##check dimensions
  if(max(loc.ind)>n.loc)
    stop("In 'calc.F.times.X': max(loc.ind) > number of locations")
  for(i in 1:m){
    p[i] <- dim(LUR[[i]])[2]
    if(n.loc != dim(LUR[[i]])[1])
      stop(sprintf("In 'calc.F.times.X': All components of 'LUR' does not have the same no. of locations (1 and %i not matching).",i))
  }
  n.obs <- dim(F)[1]
  FX <- matrix(0, n.obs, sum(p))
  Ind <- 0
  for(i in 1:m){
    tmp <- .Call("calc_F_part_X", LUR[[i]], F[,i],
                 as.integer(n.loc), as.integer(loc.ind),
                 as.integer(n.obs), as.integer(p[i]))
    FX[,(Ind+1):(Ind+p[i])] <- tmp
    Ind <- Ind + p[i]
  }
  return(FX)
}

##calculates F'*mat*F, inputs:
##The matrix inv(sigma.nu), ordinarily calculated by
##  make_chol_block and inv_chol_block
##The matrix F (each column corresponds to the temporal trend at the
##  time of that observations, trends vary between columns)
##Vector of locations, giving the location of each observation.
##number of different locations
##number of observations
##number of temporal trends
##returned matrix will be square with side n_loc*m (same as sigmaB)
calc.tF.mat.F <- function(mat, F, loc.ind, n.blocks=1,
                          block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks),
                          n.loc=max(loc.ind)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##number of obs
  n.obs <- dim(mat)[1]
  ##check dimensions
  if(dim(mat)[2] != n.obs)
    stop("In 'calc.tF.mat.F': 'mat' not a square matrix")
  if(sum(block.sizes) != n.obs)
    stop("In 'calc.tF.mat.F': sum(block.sizes) != dim(mat)[1]")
  if(dim(F)[1] != n.obs)
    stop("In 'calc.tF.mat.F': dim(F)[1] != dim(mat)[1]")
  if(max(loc.ind) > n.loc)
    stop("In 'calc.tF.mat.F': max(loc.ind) > n.loc")
  ##call the c-function
  .Call("calc_tF_mat_F", mat, F, as.integer(n.loc),
        as.integer(loc.ind), as.integer(n.blocks),
        as.integer(block.sizes), as.integer(n.obs),
        as.integer(dim(F)[2]))
}

##########################
## BASIC LINEAR ALGEBRA ##
##########################
##cholesky inverse of a block diagonal matrix using the diagonal structure.
##inputs is total matrix side, number of blocks, vector with the size of each
##block, maximum size of any block, the matrix it self.
makeCholBlock <- function(mat, n.blocks=1,
                            block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##find the largest block (needed for memory allocation in C-code)
  max.size <- max(block.sizes)
  ##check dimensions
  if(dim(mat)[1] != dim(mat)[2])
    stop("In 'makeCholBlock': 'mat' not a square matrix")
  if(sum(block.sizes) != dim(mat)[1])
    stop("In 'makeCholBlock': sum(block.sizes) != dim(mat)[1]")
  if(n.blocks==1){
    tmp <- try(chol(mat),silent=TRUE)
    if(class(tmp)=="try-error")
      tmp <- -1
  }else{
    ##call the c-function
    tmp <- .Call("make_chol_block", as.integer(dim(mat)[1]), as.integer(n.blocks),
                 as.integer(block.sizes), as.integer(max.size), mat)
  }
  ##check if matrix is pos.def
  if( tmp[1]==-1 )
    stop("In 'makeCholBlock': Matrix not positive definite.")
  return(tmp)
}

##calculate a mtrix inverse from a cholesky factor (requires a previous call
##to make_chol_block)
##inputs is total matrix side, number of blocks, vector with the size of each
##block, maximum size of any block, the matrix it self.
invCholBlock <- function(R, n.blocks=1,
                           block.sizes=rep(dim(R)[1]/n.blocks,n.blocks)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##find the largest block (needed for memory allocation in C-code)
  max.size <- max(block.sizes)
  ##check dimensions
  if(dim(R)[1] != dim(R)[2])
    stop("In 'invCholBlock': 'R' not a square matrix")
  if(sum(block.sizes) != dim(R)[1])
    stop("In 'invCholBlock': sum(block.sizes) != dim(R)[1]")

  if(n.blocks==1){
    return( chol2inv(R) )
  }else{
    ##call the c-function
    .Call("inv_chol_block", as.integer(dim(R)[1]), as.integer(n.blocks),
          as.integer(block.sizes), as.integer(max.size), R)
  }
}

##solves the system T*X=B where T is block upper triagular (output from
##make_chol_block)
##inputs is total matrix side, number of blocks, vector with the size of each
##block, maximum size of any block, number of system to solve (dim(B)[2]),
##indicator wether or not to transpose T (T/F), the triangular matrix, and the
##lhs (B) of the equation system.
solveTriBlock <- function(R, B, n.x=dim(B)[2], n.blocks=1,
                            block.sizes=rep(dim(R)[1]/n.blocks,n.blocks),
                            transpose=FALSE){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##find the largest block (needed for memory allocation in C-code)
  max.size <- max(block.sizes)
  ##check dimensions
  if(dim(R)[1] != dim(R)[2])
    stop("In 'solveTriBlock': 'R' not a square matrix")
  if(sum(block.sizes) != dim(R)[1])
    stop("In 'solveTriBlock': sum(block.sizes) != dim(R)[1]")
  B <- as.matrix(B)
  if(length(B) != (n.x*dim(R)[1]))
    stop("In 'solveTriBlock': length(B) != (n.x*dim(R)[1])")
  
  if(n.blocks==1){
    return( backsolve(R,as.matrix(B,block.sizes,n.x),transpose=transpose) )
  }else{
    ##call the c-function
    .Call("solve_tri_block", as.integer(dim(R)[1]), as.integer(n.blocks),
          as.integer(block.sizes), as.integer(max.size), as.integer(n.x),
          as.integer(transpose), R, B)
  }
}

##multiplication of block diagonal matrix with a matrix
##inputs are:
##block matrix size
##X matrix size
##number of blocks
##vector with the size of each block
##the block matrix it self
##the target matrix
block.mult <- function(mat, X, n.x=dim(X)[2], n.blocks=1,
                       block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##check dimensions
  if(dim(mat)[1] != dim(mat)[2])
    stop("In 'block.mult': 'mat' not a square matrix")
  if(sum(block.sizes) != dim(mat)[1])
    stop("In 'block.mult': sum(block.sizes) != dim(mat)[1]")
  if(length(X) != (n.x*dim(mat)[1]))
    stop("In 'block.mult': length(X) != (n.x*dim(mat)[1])")
  ##call the c-function
  .Call("block_mult", as.integer(dim(mat)[1]), as.integer(n.x),
        as.integer(n.blocks), as.integer(block.sizes), mat, X)
}

##calculate sum of the logarithm of diagonal elements, corresponds to the
##log determinant of a cholesky factor
##input matrix and matrix size
sumLogDiag <- function(mat){
  ##check dimensions
  if(dim(mat)[1] != dim(mat)[2])
    stop("In 'sumLogDiag': 'mat' not a square matrix")
  ##call the c-function
  .Call("sum_log_diag",as.integer(dim(mat)[1]),mat)
}
##calculate sum of the squared elements of a vector (or matrix) (|x|^2)
##input vector and vector size
norm2 <- function(v1){
  ##no need to check dimensions
  .Call("norm2_c",as.integer(length(v1)),v1)
}
##calculate dot product of two vectors
##input vectors and vector length
dot.prod <- function(v1,v2){
  if(length(v1)!=length(v2))
    warning("In 'dot.prod': vectors of unequal length, truncating the longer vector")
  ##no need to check dimensions
  .Call("dot_prod",as.integer(min(length(v1),length(v2))),v1,v2)
}
