#include "R.h"
#include "Rinternals.h"

//construct the covariance matrix for the spatio temporal residuals
//assuming the same spatial covariance for all time-points and no
//temporal dependence
//inputs are:
//sill, nugget and range
//total number of points
//number of blocks (time points)
//size of each block (nbr obs per time point)
//list of the location index of each obs
//number of locations (size of distance matrix)
//distance matrix with distance between each location
SEXP make_sigma_nu(SEXP res_sillR, SEXP res_nuggetR, SEXP res_rangeR,
		    SEXP n_totR, SEXP n_blocksR, SEXP block_sizesR,
		    SEXP loc_indexR, SEXP n_distR, SEXP distsR){
  int n_tot, n_blocks, n_dist, *block_sizes, *loc_index;
  int n_cum, n_this, i, j, k;
  double res_sill, res_nugget, res_range;
  double *dists, *result, *tmp_CF;
  SEXP resultR;

  res_sill = *REAL(res_sillR);
  res_nugget = *REAL(res_nuggetR);
  res_range = *REAL(res_rangeR);
  n_tot = *INTEGER(n_totR);
  n_blocks = *INTEGER(n_blocksR);
  block_sizes = INTEGER(block_sizesR);
  loc_index = INTEGER(loc_indexR);
  n_dist = *INTEGER(n_distR);
  dists = REAL(distsR);

  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  /* INITIALIZE TO 0 */
  memset(result, 0, n_tot*n_tot*sizeof(double));
  //temporary block that stores calculations so far, avoids double 
  //calculation of the covariance-function
  tmp_CF = malloc(n_dist*n_dist * sizeof(double));
  memset(tmp_CF, 0, n_dist*n_dist * sizeof(double));

  /* FILL IN MATRIX (use symmetri)*/
  n_cum = 0;
  for(k=0; k < n_blocks; ++k){
    n_this = block_sizes[k];
    for(i=0; i<n_this; ++i){
      //diagonal elements = nugget+sill
      result[i+n_cum + (i+n_cum)*n_tot] = res_sill+res_nugget;
      for(j=i+1; j<n_this; ++j){
	//check if we need to calculate for this distance
	if( tmp_CF[loc_index[j+n_cum]-1 + (loc_index[i+n_cum]-1)*n_dist]==0 )
	  tmp_CF[loc_index[j+n_cum]-1 + (loc_index[i+n_cum]-1)*n_dist] = 
	    res_sill * exp(-dists[loc_index[j+n_cum]-1 + 
				  (loc_index[i+n_cum]-1)*n_dist]/res_range);
	//pre-calculated by above if-statement, just use the value
	result[j+n_cum + (i+n_cum)*n_tot] =
	  tmp_CF[loc_index[j+n_cum]-1 + (loc_index[i+n_cum]-1)*n_dist];
	result[i+n_cum + (j+n_cum)*n_tot] = result[j+n_cum + (i+n_cum)*n_tot];
      }
    }
    n_cum += n_this;
  }//for(k=0; k < n_blocks; ++k)
  free(tmp_CF);

  UNPROTECT(1); /* resultR */
  return resultR;
}
