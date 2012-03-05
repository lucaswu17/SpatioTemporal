#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate dot product of two vectors
//input vectors and vector length
SEXP dot_prod(SEXP nR, SEXP v1R, SEXP v2R){

  int j;
  int n;
  double *v1, *v2, *result;
  SEXP resultR;

  n = *INTEGER(nR);
  v1 = REAL(v1R);
  v2 = REAL(v2R);
  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++){
    result[0] += v1[j]*v2[j];
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
