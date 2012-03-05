#include "R.h"
#include "Rinternals.h"
#include <string.h>
//constructs the full cross-covariance matrix for the land use regression
//inputs are:
//number of observations for the first locations
//number of observations for the second locations
//number of temporal trends (incl. constant)
//number of locations for the first set of locations
//number of locations for the second set of locations
//location of each observation #1 (used to index the distance matrix)
//  (n_obs1-by-1)
//location of each observation #2 (used to index the distance matrix)
//  (n_obs2-by-1)
//one sill for each trend
//  (n_trends-by-1)
//one range for each trend
//  (n_trends-by-1)
//matrix with the temporal trends for the first set of locations 
//  (n_obs1-by-n_trends)
//matrix with the temporal trends for the second set of locations
//  (n_obs2-by-n_trends)
//distance matrix with the distances between each location
//  (n_loc1-by-n_loc2)
SEXP make_sigma_b_full(SEXP n_obs1R, SEXP n_obs2R, SEXP n_trendsR,
		       SEXP n_loc1R, SEXP n_loc2R,
                       SEXP loc_ind1R, SEXP loc_ind2R,
		       SEXP trend_sillR, SEXP trend_rangeR,
		       SEXP F1R, SEXP F2R, SEXP distsR){
  int i, j, k;
  int n_obs1, n_obs2, n_trend, n_loc1, n_loc2;
  int *loc_ind1, *loc_ind2;
  double tmp;
  double *trend_sill, *trend_range, *F1, *F2, *dist;
  double *result;
  SEXP resultR;

  n_obs1 = *INTEGER(n_obs1R);
  n_obs2 = *INTEGER(n_obs2R);
  n_trend = *INTEGER(n_trendsR);
  n_loc1 = *INTEGER(n_loc1R);
  n_loc2 = *INTEGER(n_loc2R);
  loc_ind1 = INTEGER(loc_ind1R);
  loc_ind2 = INTEGER(loc_ind2R);

  trend_sill = REAL(trend_sillR);
  trend_range = REAL(trend_rangeR);
  F1 = REAL(F1R);
  F2 = REAL(F2R);
  dist = REAL(distsR);

  PROTECT(resultR = allocMatrix(REALSXP, n_obs1, n_obs2));
  result = REAL(resultR);

  /* INITIALIZE TO 0 */
  memset(result,0,n_obs1*n_obs2* sizeof(double));
  /* FILL IN MATRIX*/
  for(k=0; k < n_trend; k++){
    for (i=0; i < n_obs1; i++){
      for (j=0; j < n_obs2; j++){
	tmp = F1[i+k*n_obs1] * F2[j+k*n_obs2] * trend_sill[k] *
	  exp(-dist[(loc_ind1[i]-1)+n_loc1*(loc_ind2[j]-1)]/trend_range[k]);
        result[i + j*n_obs1] += tmp;
      }
    }
  } 

  UNPROTECT(1); /* resultR */
  return resultR;
}
