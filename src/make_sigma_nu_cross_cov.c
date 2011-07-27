#include "R.h"
#include "Rinternals.h"
//constructs the full cross-covariance matrix for the spatio temporal 
//residuals, assuming the same spatial covariance for all time-points 
//and no temporal dependence
//inputs are:
//sill, nugget and range
//number of observations for the first locations
//number of observations for the second locations
//number of locations for the first set of obs.
//number of locations for the second set of obs.
//location of each observation #1 (used to index the distance matrix)
//  (n_obs1-by-1)
//location of each observation #2 (used to index the distance matrix)
//  (n_obs2-by-1)
//observation time for observations #1
//  (n_obs1-by-1)
//observation time for observations #2
//  (n_obs2-by-1)
//distance matrix with the distances between each location
//  (n_loc1-by-n_loc2)
SEXP make_sigma_nu_cross_cov(SEXP res_sillR, SEXP res_nuggetR, 
			      SEXP res_rangeR,
			      SEXP n_obs1R, SEXP n_obs2R, 
			      SEXP n_loc1R, SEXP n_loc2R,
			      SEXP loc_ind1R, SEXP loc_ind2R,
			      SEXP loc_ind2to1R,
			      SEXP T1R, SEXP T2R, SEXP distsR){
  int i, j;
  int n_obs1, n_obs2, n_loc1, n_loc2;
  int *loc_ind1, *loc_ind2, *loc_ind2to1, *T1, *T2;
  double res_sill, res_nugget, res_range;
  double *dists, *result;
  SEXP resultR;

  res_sill = *REAL(res_sillR);
  res_nugget = *REAL(res_nuggetR);
  res_range = *REAL(res_rangeR);

  n_obs1 = *INTEGER(n_obs1R);
  n_obs2 = *INTEGER(n_obs2R);
  n_loc1 = *INTEGER(n_loc1R);
  n_loc2 = *INTEGER(n_loc2R);
  loc_ind1 = INTEGER(loc_ind1R);
  loc_ind2 = INTEGER(loc_ind2R);
  loc_ind2to1 = INTEGER(loc_ind2to1R);
  T1 = INTEGER(T1R);
  T2 = INTEGER(T2R);

  dists = REAL(distsR);

  PROTECT(resultR = allocMatrix(REALSXP, n_obs1, n_obs2));
  result = REAL(resultR);

  /* INITIALIZE TO 0 */
  memset(result, 0, n_obs1*n_obs2*sizeof(double));

  /* FILL IN MATRIX */
  for(i=0; i<n_obs1; ++i){
    for(j=0; j<n_obs2; ++j){
      //check if both observations are from the same time point 
      //(o.w we assume uncorrelated)
      if(T1[i]==T2[j]){
	if(loc_ind1[i]==loc_ind2to1[loc_ind2[j]-1]) //co-located use nugget
	  result[i+ j*n_obs1] = res_sill+res_nugget;
	else //not co-located calculate covariance function
	  result[i+ j*n_obs1] = res_sill * 
	    exp(-dists[loc_ind1[i]-1 + (loc_ind2[j]-1)*n_loc1]/res_range);
      }
    }//for(j=0; j<n_obs2; ++j)
  }//for(i=0; i<n_obs1; ++i)

  UNPROTECT(1); /* resultR */
  return resultR;
}
