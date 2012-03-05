#include "R.h"
#include "Rinternals.h"
#include <string.h>

//constructs the covariance matrix for the land use regression
//inputs are:
//number of temporal trends (incl. constant)
//number of locations
//one sill for each trend
//one range for each trend
//distance matrix with the distances between each location
SEXP make_sigma_b(SEXP n_trendR, SEXP n_locR,
		  SEXP trend_sillR, SEXP trend_rangeR,
		  SEXP distsR){

  int i, j, k;
  int n_trend, n_loc;
  double *trend_sill, *trend_range, *dist;
  double *result;
  SEXP resultR;

  n_trend = *INTEGER(n_trendR);
  n_loc = *INTEGER(n_locR);

  trend_sill = REAL(trend_sillR);
  trend_range = REAL(trend_rangeR);
  dist = REAL(distsR);

  PROTECT(resultR = allocMatrix(REALSXP, n_trend*n_loc, n_trend*n_loc));
  result = REAL(resultR);


  /* INITIALIZE TO 0 */
  memset(result,0,n_trend*n_loc * n_trend*n_loc * sizeof(double));
  /* FILL IN MATRIX (use symmetri)*/
  for(k=0; k < n_trend; k++){
    for (i=0; i < n_loc; i++){
      result[(k*n_loc+i) + (k*n_loc+i)*(n_loc*n_trend)] = trend_sill[k];
      for (j=i+1; j < n_loc; j++){
        result[(k*n_loc+j) + (k*n_loc+i)*(n_loc*n_trend)] =
	  trend_sill[k]*exp(-dist[j+n_loc*i]/trend_range[k]);
	result[(k*n_loc+i) + (k*n_loc+j)*(n_loc*n_trend)] =
	  result[(k*n_loc+j) + (k*n_loc+i)*(n_loc*n_trend)];
      }
    }
  } 

  UNPROTECT(1); /* resultR */
  return resultR;
}
