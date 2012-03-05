#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate sum of the logarithm of diagonal elements, corresponds to the
//log determinant of a cholesky factor
//input matrix and matrix size
SEXP sum_log_diag(SEXP nR, SEXP matR){

  int j;
  int n;
  double *mat, *result;
  SEXP resultR;

  n = *INTEGER(nR);
  mat = REAL(matR);
  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++)
    result[0] += log(mat[j+n*j]);
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
