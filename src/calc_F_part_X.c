#include "R.h"
#include "Rinternals.h"

//calculates F_i*X_i, inputs:
//A matrix X_i
//The vector F (one "column" of F, corresponds to the value of the
//  temporal trend at the time of that observations)
//number of different locations
//Vector of locations, giving the location of each observation. Needed to know
//  which element of X that we should multiply elements in F with.
//number of observations
//number of columns in X (number of LUR components)
//returns a column vector with (n_obs*m) elements
SEXP calc_F_part_X(SEXP XR, SEXP FR, SEXP n_locR, SEXP locIndR,
		       SEXP n_obsR, SEXP PiR){
  int *locInd, n_loc, n_obs, pi, i, i_pi;
  double *X, *F, *result;
  SEXP resultR;

  X = REAL(XR);
  F = REAL(FR);
  n_loc = *INTEGER(n_locR);
  locInd = INTEGER(locIndR);
  n_obs = *INTEGER(n_obsR);
  pi = *INTEGER(PiR);

  PROTECT(resultR = allocMatrix(REALSXP, n_obs, pi));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_obs*pi*sizeof(double));

  for(i_pi=0; i_pi<pi; ++i_pi)
    for(i=0; i<n_obs; ++i)
      result[i + i_pi*n_obs] += F[i] * X[locInd[i]-1 + n_loc*i_pi];
  UNPROTECT(1); //resultR
  return resultR;
}
