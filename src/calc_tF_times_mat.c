#include "R.h"
#include "Rinternals.h"

//calculates F'*v, inputs:
//A vector v
//The matrix F (for a given column each row corresponds to the value of the
//  temporal trend at the time of that observations, trends vary between columns)
//number of different locations
//Vector of locations, giving the location of each observation. Needed to know
//  which element of sigma.res that we should multiply elements in F with.
//number of observations
//number of temporal trends
//returns a column vector with (n_obs*m) elements
SEXP calc_tF_times_mat(SEXP XR, SEXP FR, SEXP n_locR, SEXP n_xR,
		      SEXP locIndR, SEXP n_obsR, SEXP mR){
  int *locInd, n_loc, n_obs, n_x, m, n_tot, j, i, i_m;
  double *X, *F, *result;
  SEXP resultR;

  X = REAL(XR);
  F = REAL(FR);
  n_loc = *INTEGER(n_locR);
  n_x = *INTEGER(n_xR);
  locInd = INTEGER(locIndR);
  n_obs = *INTEGER(n_obsR);
  m = *INTEGER(mR);
  n_tot = n_loc*m;

  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_tot*n_x*sizeof(double));

  for(j=0; j<n_x; ++j)
    for(i_m=0; i_m<m; ++i_m)
      for(i=0; i<n_obs; ++i)
	result[locInd[i]-1 + i_m*n_loc + j*n_tot] += 
	  F[i+i_m*n_obs] * X[i+j*n_obs];
  
  UNPROTECT(1); //resultR
  return resultR;
}
