#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate sum of the squared elements of a vector (|x|^2)
//input vector and vector size
SEXP norm2_c(SEXP nR, SEXP matR){

  int j;
  int n;
  double *mat, *result;
  SEXP resultR;

  n = *INTEGER(nR);
  mat = REAL(matR);
  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++){
    result[0] += mat[j]*mat[j];
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
