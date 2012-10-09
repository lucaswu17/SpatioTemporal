#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate sum of the squared elements of a vector (|x|^2)
//input vector and vector size
SEXP C_norm2(SEXP matR){

  //extract input from R
  double *mat = REAL(matR);
  int n = length(matR);
  //loop variables
  int j;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++){
    result[0] += mat[j]*mat[j];
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
