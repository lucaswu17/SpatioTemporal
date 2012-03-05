#include "R.h"
#include "Rinternals.h"

//calculates F'*inv(sigma.res)*F, inputs:
//The matrix inv(sigma.res), ordinarily calculated by calls to
//  make_sigma_res, make_chol_block, and inv_chol_block
//The matrix F (for a given column each row corresponds to the value of the
//  temporal trend at the time of that observations, trends vary between columns)
//number of different locations
//Vector of locations, giving the location of each observation. Needed to know
//  which element of sigma.res that we should multiply elements in F with.
//number of blocks
//vector with the size of each block
//number of observations
//number of temporal trends
//returned matrix will be square with side n_loc*m (same as sigmaB)
//  the returned matrix is band diagonal with bandwidth n_loc.
SEXP calc_tF_mat_F(SEXP iSigmaR, SEXP FR, SEXP n_locR, SEXP locIndR,
		       SEXP n_blocksR, SEXP block_sizesR,
		       SEXP n_obsR, SEXP mR){

  int *locInd, n_loc, n_obs, m, n_tot, n_blocks, *block_sizes;
  int i, j, k, l, i_m, loc1, loc2;
  double *iSigma, *F, *result, *iS_F;
  SEXP resultR;

  iSigma = REAL(iSigmaR);
  F = REAL(FR);
  n_loc = *INTEGER(n_locR);
  locInd = INTEGER(locIndR);
  n_blocks = *INTEGER(n_blocksR);
  block_sizes = INTEGER(block_sizesR);
  n_obs = *INTEGER(n_obsR);
  m = *INTEGER(mR);
  n_tot = n_loc*m;

  
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_tot*n_tot*sizeof(double));

  //temporary storage for inv(Sigma)*F
  iS_F = malloc( n_obs*n_tot * sizeof (double));
  memset(iS_F, 0, n_obs*n_tot * sizeof(double));

  //calculate inv(Sigma)*F
  for(i_m=0; i_m<m; ++i_m){
    loc1=0;
    for(i=0; i<n_blocks; ++i){
      for(k=0; k<block_sizes[i]; ++k){
	for(l=0; l<block_sizes[i]; ++l){
	  loc2 = locInd[loc1+l]-1;
	  iS_F[loc1+k + loc2*n_obs + i_m*n_obs*n_loc] +=
	    F[loc1+l+i_m*n_obs]*iSigma[loc1+k + (loc1+l)*n_obs];
	}
      }
      loc1 += block_sizes[i];
    }
  }
  //calculate F' * (inv(Sigma)*F)
  for(l=0; l<n_tot; ++l)
    for(k=0; k<m; ++k)
      for(i=0; i<n_obs; ++i){
	loc1 = locInd[i]-1;
	result[loc1+k*n_loc + l*n_tot] += F[i + k*n_obs] * iS_F[i + l*n_obs];
      }
  
  free(iS_F);

  UNPROTECT(1); //resultR
  return resultR;
}
